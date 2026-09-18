# Node keys glue species and gene together. The published workflow used
# paste0(Species, "-", Gene), which silently merges two distinct genes
# whenever a gene id contains a hyphen; \x01 cannot occur in an
# identifier, so the join is injective.
.gcg_sep <- "\x01"


#' Assert that identifiers are free of the node-key separator
#'
#' @param edges Edge table with species/gene columns.
#' @return Invisibly `TRUE`; errors when a separator is embedded.
#' @noRd
.gcg_check_ids <- function(edges) {
  cols <- c(edges$species1, edges$species2, edges$gene1, edges$gene2)
  if (any(grepl(.gcg_sep, cols, fixed = TRUE))) {
    stop(
      "species and gene identifiers must not contain \\x01: it ",
      "separates the two halves of a gene-graph node key"
    )
  }
  invisible(TRUE)
}


#' Refuse a clique table that row-binds two runs under one clique_id
#'
#' `gene_clique_graph()` numbers cliques positionally within a HOG, so
#' two runs at different `alpha_graph` produce the same ids unless the
#' caller gave them distinct `id_prefix` values. Row-binding them merges
#' unrelated genes into one oversized "clique" whose `n_members`,
#' `n_species`, `n_pairs` and `missing_species` are then arithmetic over
#' genes that were never adjacent -- silently, since nothing downstream
#' can tell a merged block from a real clique.
#'
#' Three independent signatures of such a merge are refused, since each
#' misses cases the others catch: the repeated (id, species, gene)
#' triple needs the colliding cliques to share a member, the shared id
#' needs them to come from different HOGs, and the row-count check needs
#' the input to carry `n_members` at all. Their union is not complete --
#' a table without `n_members` whose colliding cliques come from one HOG
#' and share no member passes all three. Give each
#' [gene_clique_graph()] run a distinct `id_prefix` rather than relying
#' on detection.
#'
#' @param cl_id Clique id of each row.
#' @param mk_all Node key (species + gene) of each row.
#' @param cl_hog HOG of each row.
#' @param n_members Declared clique size of each row, or `NULL` when the
#'   input does not carry the column.
#' @return Invisibly `TRUE`; errors on any merge signature.
#' @noRd
.gcg_check_clique_ids <- function(cl_id, mk_all, cl_hog, n_members) {
  hint <- ": give each gene_clique_graph() run a distinct id_prefix"
  if (anyDuplicated(paste(cl_id, mk_all, sep = .gcg_sep)) > 0L) {
    stop(
      "cliques contain repeated (clique_id, species, gene) rows", hint
    )
  }
  n_id <- length(unique(cl_id))
  if (length(unique(paste(cl_id, cl_hog, sep = .gcg_sep))) != n_id) {
    stop("cliques contain a clique_id spanning several hogs", hint)
  }
  if (!is.null(n_members)) {
    by_id <- split(as.integer(n_members), cl_id)
    # Only MORE rows than declared can be a merge -- stacking two blocks
    # adds rows, it never removes them. Fewer rows is a row-filtered
    # table, which is a legitimate thing to hand this function, and
    # refusing it was a regression that misdiagnosed the cause. A single
    # clique_id declaring two different sizes is still two blocks.
    bad <- vapply(by_id, function(v) {
      u <- unique(v[!is.na(v)])
      length(u) > 1L || (length(u) == 1L && length(v) > u)
    }, logical(1))
    if (any(bad)) {
      stop(
        "cliques have more member rows than n_members declares", hint
      )
    }
  }
  invisible(TRUE)
}


#' Empty result template for [gene_clique_graph()]
#'
#' @param has_effect Whether an `effect_size` column is carried.
#' @return A zero-row data frame with the full column set.
#' @noRd
.gcg_empty <- function(has_effect) {
  out <- data.frame(
    clique_id = character(0), hog = character(0),
    species = character(0), gene = character(0),
    n_members = integer(0), n_species = integer(0),
    n_edges = integer(0), mean_q = numeric(0), max_q = numeric(0),
    stringsAsFactors = FALSE
  )
  if (has_effect) out$mean_effect_size <- numeric(0)
  out
}


#' Attach the floor diagnostics to a [gene_clique_graph()] result
#'
#' Carried as columns as well as attributes: `subset()` drops attributes,
#' so a caller who filters the result would otherwise lose the record of
#' how little resolution `mean_q` has.
#'
#' @param out Result frame (possibly zero-row).
#' @param alpha_graph,min_size Call parameters.
#' @param q_floor,mean_q_floor,n_tied Floor diagnostics.
#' @return `out` with floor columns and attributes.
#' @noRd
.gcg_graph_attrs <- function(out, alpha_graph, min_size, q_floor,
                             mean_q_floor, n_tied) {
  out$mean_q_floor <- rep(mean_q_floor, nrow(out))
  out$n_cliques_at_q_floor <- rep(n_tied, nrow(out))
  attr(out, "alpha_graph") <- alpha_graph
  attr(out, "min_size") <- min_size
  attr(out, "q_floor") <- q_floor
  attr(out, "mean_q_floor") <- mean_q_floor
  attr(out, "n_cliques_at_q_floor") <- n_tied
  out
}


#' Maximal cliques of the per-orthogroup gene graph
#'
#' Builds, for each ortholog group, an undirected graph whose nodes are
#' (species, gene) pairs and whose edges are the co-expressolog calls in
#' `edges`, then enumerates every maximal clique of at least `min_size`
#' nodes. This is the computation used by Netotea-style clique workflows
#' and it differs from \code{\link{find_cliques}}, which cliques the
#' \emph{species} graph and returns a single best gene assignment per
#' species clique: here every paralog combination that forms a clique is
#' reported separately.
#'
#' Because co-expressolog edges are always cross-species, no two genes of
#' the same species can be adjacent, so a clique carries at most one gene
#' per species. The `n_species` column records that explicitly rather
#' than assuming it, so a caller that supplies within-species edges can
#' still detect the violation instead of getting silently wrong
#' species-pair arithmetic downstream. Only a row whose two endpoints are
#' the *same* (species, gene) node is dropped as a self-loop; a
#' within-species edge between two distinct paralogs is kept, and shows
#' up as `n_species < n_members`.
#'
#' @section Duplicate rows:
#' Two rows describing the same undirected pair within one HOG collapse
#' to the more significant of the two, broken on `effect_size`
#' descending when the q-values tie -- q-values saturate at the
#' permutation floor, so a tie there must not be settled by input row
#' order. This is a duplicate-row collapse, not a direction combine:
#' \code{\link{find_coexpressologs}} already merges the two comparison
#' directions with `pval_combine` (default `"max"`, the reciprocal
#' criterion) and emits one row per pair. A caller who row-binds the two
#' directional tables instead gets `"min"` semantics here, so combine
#' upstream if the reciprocal criterion is wanted.
#'
#' @param edges Data frame of co-expressolog calls, as returned by
#'   \code{\link{find_coexpressologs}}: columns `gene1`, `gene2`,
#'   `species1`, `species2`, `hog` and `q.value`. An `effect_size`
#'   column is used when present.
#' @param min_size Minimum number of nodes in a reported clique
#'   (default 3, matching the published workflow).
#' @param alpha_graph Edges with `q.value < alpha_graph` build the
#'   graph. Use the calling threshold (e.g. 0.1) for complete cliques,
#'   a permissive value (0.9) for partially-significant cliques, and a
#'   value **above** 1 (`Inf`) for an unfiltered graph. The comparison
#'   is strict, so `alpha_graph = 1` drops every edge whose `q.value` is
#'   exactly 1 -- not a corner case, since BH q-values cap at 1 and 20 of
#'   511 module q-values sit there on the package's own worked
#'   example. Passing several thresholds and combining
#'   the results is the intended way to feed
#'   \code{\link{classify_gene_cliques}}, since a clique that is maximal
#'   at one threshold need not be maximal at another.
#' @param id_prefix String prepended to every `clique_id`. Clique ids
#'   are `<prefix><hog>_<k>`, so runs at different `alpha_graph` values
#'   need distinct prefixes before they can be row-bound.
#'
#' @section rcomplex container:
#' The `.rcomplex` method builds the graph from `x$edges`, which in the
#' container workflow is usually the permutation table
#' ([find_coexpressologs()] with `method = "permutation"`). Its
#' Besag-Clifford q-values top out well below 1 (0.4 on the package's own
#' worked example), so a permissive `alpha_graph` set to relax against
#' that ceiling can end up with nothing to relax against. The `.rcomplex`
#' method warns when `alpha_graph` exceeds the largest `q.value` in
#' `x$edges`.
#'
#' @return A data frame with one row per clique member:
#'   \describe{
#'     \item{clique_id}{Clique identifier, unique within the run}
#'     \item{hog}{Ortholog group}
#'     \item{species, gene}{The member}
#'     \item{n_members}{Nodes in the clique}
#'     \item{n_species}{Distinct species in the clique (equal to
#'       `n_members` for cross-species edge tables)}
#'     \item{n_edges}{Graph edges inside the clique}
#'     \item{mean_q, max_q}{q-value summaries over those edges}
#'     \item{mean_effect_size}{Present only when `edges` carries
#'       `effect_size`. Prefer it over `mean_q` for ranking: q-values
#'       saturate at the permutation floor and cannot separate cliques
#'       once they get there}
#'     \item{mean_q_floor, n_cliques_at_q_floor}{Smallest clique
#'       `mean_q` in the run and how many cliques are tied there. A
#'       large tie count means `mean_q` cannot rank those cliques at
#'       all}
#'   }
#'   The same two floor diagnostics, plus `alpha_graph`, `min_size` and
#'   `q_floor` (smallest edge q-value in the graph), are also recorded as
#'   attributes. Both are always set, empty results included.
#'
#' @examples
#' edges <- data.frame(
#'   gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
#'   species1 = c("SP_A", "SP_A", "SP_B"),
#'   species2 = c("SP_B", "SP_C", "SP_C"),
#'   hog = "HOG1", q.value = c(0.01, 0.02, 0.03)
#' )
#' gene_clique_graph(edges)
#'
#' @seealso \code{\link{classify_gene_cliques}},
#'   \code{\link{find_cliques}}
#' @references
#' Rodriguez E, Birkeland S, Chapple ED, et al. (2026).
#' Comparative regulomics of wood formation across dicot and
#' conifer trees. \emph{Nature Communications} 17(1).
#' \doi{10.1038/s41467-026-75624-2}
#' @param ... Additional arguments passed to the default method.
#' @export
gene_clique_graph <- function(edges, ...) UseMethod("gene_clique_graph")

#' @rdname gene_clique_graph
#' @export
gene_clique_graph.default <- function(edges, min_size = 3L,
                                      alpha_graph = 0.1,
                                      id_prefix = "", ...) {
  rlang::check_dots_empty()
  required <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value"
  )
  absent <- setdiff(required, names(edges))
  if (length(absent) > 0L) {
    stop(
      "edges missing required columns: ",
      paste(absent, collapse = ", ")
    )
  }
  min_size <- as.integer(min_size)
  if (length(min_size) != 1L || is.na(min_size) || min_size < 2L) {
    stop("min_size must be a single integer >= 2")
  }
  ok_alpha <- is.numeric(alpha_graph) && length(alpha_graph) == 1L &&
    !is.na(alpha_graph)
  if (!ok_alpha) {
    stop("alpha_graph must be a single non-missing number")
  }
  if (!is.character(id_prefix) || length(id_prefix) != 1L) {
    stop("id_prefix must be a single string")
  }
  .gcg_check_ids(edges)

  has_effect <- "effect_size" %in% names(edges)
  keep <- !is.na(edges$q.value) & edges$q.value < alpha_graph
  edges <- edges[keep, , drop = FALSE]
  if (nrow(edges) == 0L) {
    return(.gcg_graph_attrs(
      .gcg_empty(has_effect), alpha_graph, min_size,
      NA_real_, NA_real_, 0L
    ))
  }

  key1 <- paste(edges$species1, edges$gene1, sep = .gcg_sep)
  key2 <- paste(edges$species2, edges$gene2, sep = .gcg_sep)
  qv <- as.numeric(edges$q.value)
  ev <- if (has_effect) {
    as.numeric(edges$effect_size)
  } else {
    rep(NA_real_, nrow(edges))
  }
  hog_chr <- as.character(edges$hog)

  by_hog <- split(seq_len(nrow(edges)), hog_chr)
  cid_l <- list()
  sp_l <- list()
  gn_l <- list()
  hog_v <- character(0)
  nm_v <- integer(0)
  ns_v <- integer(0)
  ne_v <- integer(0)
  mq_v <- numeric(0)
  xq_v <- numeric(0)
  me_v <- numeric(0)
  j <- 0L

  for (h in names(by_hog)) {
    idx <- by_hog[[h]]
    k1 <- key1[idx]
    k2 <- key2[idx]
    # Both endpoints on the same node is a self-loop, which igraph would
    # happily accept and which would inflate n_edges. A within-species
    # row between two *distinct* paralogs is a real edge and is kept;
    # it surfaces downstream as n_species < n_members.
    loop <- k1 == k2
    idx <- idx[!loop]
    if (length(idx) == 0L) next
    k1 <- key1[idx]
    k2 <- key2[idx]

    nodes <- c(k1, k2)
    node_sp <- c(edges$species1[idx], edges$species2[idx])
    node_gn <- c(edges$gene1[idx], edges$gene2[idx])
    first <- !duplicated(nodes)
    nodes <- nodes[first]
    node_sp <- node_sp[first]
    node_gn <- node_gn[first]
    if (length(nodes) < min_size) next

    i1 <- match(k1, nodes)
    i2 <- match(k2, nodes)
    # Duplicate rows for one undirected pair collapse to the most
    # significant copy. The q tie-break is effect_size descending, never
    # input row order: q saturates at the permutation floor, and the
    # surviving row is what sets mean_effect_size, the ranking key.
    pkey <- paste(pmin(i1, i2), pmax(i1, i2), sep = "-")
    ord <- order(pkey, qv[idx], -ev[idx], na.last = TRUE)
    sel <- ord[!duplicated(pkey[ord])]
    i1 <- i1[sel]
    i2 <- i2[sel]
    eq <- qv[idx][sel]
    ee <- ev[idx][sel]

    g <- igraph::make_graph(as.vector(rbind(i1, i2)),
      n = length(nodes), directed = FALSE
    )
    cliques <- igraph::max_cliques(g, min = min_size)
    if (length(cliques) == 0L) next

    for (k in seq_along(cliques)) {
      cc <- as.integer(cliques[[k]])
      inc <- logical(length(nodes))
      inc[cc] <- TRUE
      hit <- inc[i1] & inc[i2]
      j <- j + 1L
      cid_l[[j]] <- rep(paste0(id_prefix, h, "_", k), length(cc))
      sp_l[[j]] <- node_sp[cc]
      gn_l[[j]] <- node_gn[cc]
      hog_v[j] <- h
      nm_v[j] <- length(cc)
      ns_v[j] <- length(unique(node_sp[cc]))
      ne_v[j] <- sum(hit)
      mq_v[j] <- mean(eq[hit])
      xq_v[j] <- max(eq[hit])
      me_v[j] <- if (has_effect) mean(ee[hit]) else NA_real_
    }
  }

  if (j == 0L) {
    return(.gcg_graph_attrs(
      .gcg_empty(has_effect), alpha_graph, min_size,
      min(qv), NA_real_, 0L
    ))
  }

  out <- data.frame(
    clique_id = unlist(cid_l, use.names = FALSE),
    hog = rep(hog_v, times = nm_v),
    species = unlist(sp_l, use.names = FALSE),
    gene = unlist(gn_l, use.names = FALSE),
    n_members = rep(nm_v, times = nm_v),
    n_species = rep(ns_v, times = nm_v),
    n_edges = rep(ne_v, times = nm_v),
    mean_q = rep(mq_v, times = nm_v),
    max_q = rep(xq_v, times = nm_v),
    stringsAsFactors = FALSE
  )
  if (has_effect) {
    out$mean_effect_size <- rep(me_v, times = nm_v)
  }
  # Tie counts go through the tolerant comparison for the same reason
  # pvalue_resolution() does: mean_q is a mean over a different edge
  # subset per clique, so two mathematically equal values need not be
  # bitwise equal, and a tie count that exists to say "these cannot be
  # ranked" must not under-report the tie.
  mq_ties <- .tol_min_ties(mq_v)
  .gcg_graph_attrs(
    out, alpha_graph, min_size, min(qv), mq_ties$min, mq_ties$n_at_min
  )
}


#' State of one absent species relative to a clique
#'
#' Separates the three states the tier logic depends on: the species is
#' absent from the orthogroup, it is present but was never tested
#' against these clique members, or it was tested and failed. Only the
#' first two are annotation gaps; folding the third into them would let
#' `partial_present` claim cliques that were actually rejected. A species
#' whose every failed test had power below `min_power` is
#' `"underpowered"`: its failure is no rejection, so `partial_present`
#' counts it as a gap, but `lineage_specific` does not, because that call
#' would have to survive reading the species as conserved.
#'
#' @param s Species to classify.
#' @param rows Row indices of `edges` belonging to the clique's HOG.
#' @param edges Full edge table.
#' @param key1,key2 Node keys of `edges`, parallel to its rows.
#' @param mk Node keys of the clique members.
#' @param m Number of clique members.
#' @param alpha_call Significance threshold.
#' @param power Detection power of each row of `edges` (`NA` when the
#'   table carries none).
#' @param min_power Power below which a failed test is uninformative.
#' @return One of "absent", "untested", "tested_ns", "underpowered",
#'   "extendable".
#' @noRd
.gcg_missing_reason <- function(s, rows, edges, key1, key2, mk, m,
                                alpha_call, power, min_power) {
  if (length(rows) == 0L) {
    return("absent")
  }
  sp1 <- edges$species1[rows]
  sp2 <- edges$species2[rows]
  if (!any(sp1 == s) && !any(sp2 == s)) {
    return("absent")
  }
  k1 <- key1[rows]
  k2 <- key2[rows]
  # A row tests the species against the clique only when the species
  # sits on one side and a clique member on the other.
  hit1 <- sp1 == s & k2 %in% mk
  hit2 <- sp2 == s & k1 %in% mk
  if (!any(hit1) && !any(hit2)) {
    return("untested")
  }
  cand <- c(k1[hit1], k2[hit2])
  partner <- c(k2[hit1], k1[hit2])
  qs <- c(edges$q.value[rows][hit1], edges$q.value[rows][hit2])
  sig <- !is.na(qs) & qs < alpha_call
  if (any(sig)) {
    # A candidate significant against every member would enlarge the
    # clique: the input clique set was built on a looser graph, so this
    # is a threshold mismatch rather than an annotation gap. Count
    # distinct members, not rows -- duplicate rows for one pair would
    # otherwise fake a candidate up to full membership.
    seen <- paste(cand[sig], partner[sig], sep = .gcg_sep)
    uniq_cand <- cand[sig][!duplicated(seen)]
    if (max(table(uniq_cand)) >= m) {
      return("extendable")
    }
  }
  # Whatever significant edges the species has, they were not enough to
  # join the clique, so what keeps it out are the tests that failed --
  # whether or not some other test succeeded. If none of those failures
  # could have succeeded, that is absent evidence rather than a
  # rejection; NA power keeps the old reading, and a species whose every
  # test was significant (but which still cannot extend the clique) has
  # no failure to excuse and stays `tested_ns`.
  pw <- c(power[rows][hit1], power[rows][hit2])[!sig]
  if (length(pw) > 0L && all(!is.na(pw) & pw < min_power)) {
    return("underpowered")
  }
  "tested_ns"
}


#' Tier order of the classification waterfall
#' @noRd
.gcg_tiers <- c(
  "complete_conserved", "lineage_specific",
  "partial_significant", "partial_present",
  "differentiated", "underpowered", "unclassified"
)


#' Classify gene-graph cliques into conservation tiers
#'
#' Applies the published five-tier taxonomy to the cliques returned by
#' \code{\link{gene_clique_graph}}. Every threshold is derived from
#' the number of species actually supplied, so nothing is tied to the
#' six species and fifteen species-pairs of the original workflow.
#'
#' With `S = length(species)` species, `P = choose(S, 2)` species-pairs,
#' lineage sizes `n_l`, `W = sum(choose(n_l, 2))` within-lineage pairs
#' and `X = P - W` cross-lineage pairs, the tiers are:
#' \describe{
#'   \item{complete_conserved}{All `S` species present and all
#'     `choose(S, 2)` pairs significant at `alpha_call`.}
#'   \item{lineage_specific}{A complete clique over one entire lineage,
#'     with no species outside it testable against the clique. A species
#'     that *was* compared against every member and came back
#'     non-significant is evidence of a boundary rather than a gap, so
#'     it blocks this tier -- that clique is a candidate for
#'     `differentiated`, on the cliques of the unfiltered graph.}
#'   \item{partial_significant}{All `S` species present, every clique
#'     edge below `alpha_graph`, and at least `choose(S - 1, 2) + 1`
#'     pairs significant at `alpha_call`.}
#'   \item{partial_present}{`S - g` species present for
#'     `1 <= g <= max_gap`, all `choose(S - g, 2)` pairs significant,
#'     and every absent species an annotation gap or `underpowered`
#'     rather than a rejected test.}
#'   \item{differentiated}{All `S` species present, every pair tested,
#'     at least one lineage of two or more species, every such lineage
#'     fully significant within itself, and at most `cross_max`
#'     cross-lineage pairs significant.}
#'   \item{underpowered}{A clique that would be `lineage_specific` or
#'     `differentiated` but for tests that could not have succeeded,
#'     read from a `power` column in `edges` (see
#'     \code{\link{comparison_to_edges}}). It takes the place of
#'     `lineage_specific` when every outside species is a gap or
#'     `underpowered` and at least one is `underpowered`, and of
#'     `differentiated` when `n_sig_cross + n_underpowered_cross`
#'     exceeds `cross_max`: a specificity or divergence call must survive
#'     treating every underpowered pair as possibly significant. A
#'     low-degree gene cannot reach the call whatever its conservation,
#'     so without this its missing edges read as a lineage boundary.}
#' }
#' The waterfall is evaluated in that order, `underpowered` at the
#' position of the tier it replaces, and the first match wins.
#'
#' `choose(S - 1, 2) + 1` equals `choose(S, 2) - (S - 2)`, so the
#' `partial_significant` tolerance is `S - 2` non-significant edges.
#' That is the largest tolerance under which no member can be isolated:
#' cutting a member loose needs all `S - 1` of its edges removed, so
#' with `S - 2` removals every member still keeps a significant edge.
#'
#' `differentiated` demands that every member pair actually has a row in
#' `edges`. Without that, a clique whose cross-lineage pairs were never
#' tested would be scored as diverged on absent evidence -- the same
#' conflation `partial_present` refuses through `missing_reason`.
#'
#' @param cliques Data frame from \code{\link{gene_clique_graph}}, or
#'   any table with `clique_id`, `hog`, `species` and `gene` columns.
#'   Combine runs at several `alpha_graph` values (with distinct
#'   `id_prefix`) to expose every tier: a clique complete at
#'   `alpha_call` need not be maximal on a looser graph. Two runs that
#'   collide on a `clique_id` are refused where the collision is
#'   detectable; detection is not complete, so give each
#'   [gene_clique_graph()] run a distinct `id_prefix` rather than
#'   relying on it.
#' @param edges The full, unfiltered co-expressolog table. It must not
#'   be pre-filtered on `q.value`: the gap tier needs to see rows that
#'   were tested and failed in order to refuse them. An optional `power`
#'   column (from \code{\link{comparison_to_edges}}) enables the
#'   `underpowered` tier; without it, or where it is `NA`, the
#'   classification is unchanged.
#' @param species Character vector of every species in the analysis.
#'   Every species appearing in `cliques` must be listed; a stranger
#'   would be counted into the clique's species total while also being
#'   reported as missing.
#' @param lineage Optional named vector mapping each species to a
#'   lineage. Required for `lineage_specific` and `differentiated`;
#'   when `NULL` those two tiers are skipped.
#' @param alpha_call Significance threshold for calling a species pair
#'   conserved (default 0.1).
#' @param alpha_graph Loose threshold defining the `partial_significant`
#'   graph (default 0.9).
#' @param max_gap Largest number of absent species tolerated by
#'   `partial_present` (default 1).
#' @param cross_max Maximum number of significant cross-lineage pairs
#'   allowed by `differentiated`. Defaults to `choose(S - 1, 2) - W`,
#'   the generalisation of the published cut; the original six-species
#'   script used a looser hard-coded 6.
#' @param min_power Detection power below which a non-significant pair
#'   is read as uninformative rather than as evidence against
#'   conservation (default 0.8). Only used when `edges` has `power`.
#'
#' @section rcomplex container:
#' The `.rcomplex` method calls with `edges = x$edges`, which in the
#' container workflow is usually the permutation table
#' ([find_coexpressologs()] with `method = "permutation"`). Its
#' Besag-Clifford q-values top out well below 1 (0.4 on the package's own
#' worked example), so the default `alpha_graph = 0.9` has nothing to
#' relax against there. The `.rcomplex` method warns when `alpha_graph`
#' exceeds the largest `q.value` in `x$edges`.
#'
#' @return A data frame with one row per clique:
#'   \describe{
#'     \item{clique_id, hog}{Clique identity}
#'     \item{classification}{Tier, or `"unclassified"`}
#'     \item{hog_class}{Earliest tier reached by any clique of that HOG,
#'       reproducing the HOG-level precedence of the published scripts
#'       without discarding the per-clique detail}
#'     \item{n_members, n_species}{Clique size}
#'     \item{n_pairs, n_present, n_sig}{Member pairs, pairs with a row
#'       in `edges`, and pairs significant at `alpha_call`}
#'     \item{n_sig_within, n_sig_cross}{Significant pairs within and
#'       across lineages (`NA` without `lineage`)}
#'     \item{n_underpowered_cross}{Cross-lineage pairs present in
#'       `edges`, not significant, and with `power` below `min_power`
#'       (`NA` without `lineage`)}
#'     \item{n_missing, missing_species, missing_reason}{Species not in
#'       the clique, and why: `absent` (no row in the HOG),
#'       `untested` (in the HOG but never compared to a member),
#'       `tested_ns` (compared and not significant), `underpowered`
#'       (compared and not significant, every such test with `power`
#'       below `min_power`), `extendable`
#'       (significant against every member, so the clique came from a
#'       looser graph). Comma-separated and positionally aligned}
#'     \item{mean_q, max_q}{Recomputed from `edges` over present pairs}
#'     \item{mean_effect_size}{Present when `edges` carries
#'       `effect_size`; the ranking key to prefer, since `mean_q`
#'       saturates}
#'     \item{mean_q_floor, n_cliques_at_q_floor}{Smallest `mean_q` over
#'       the classified cliques and how many are tied there -- the
#'       resolution `mean_q` actually has as a ranking key}
#'   }
#'   The two floor diagnostics are also attributes, alongside
#'   `n_species_total`, `n_pairs_total`, `n_within_pairs`,
#'   `n_cross_pairs`, `cross_max`, `alpha_call`, `alpha_graph` and
#'   `max_gap`. All are set even when the result is empty.
#'
#' @examples
#' edges <- data.frame(
#'   gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
#'   species1 = c("SP_A", "SP_A", "SP_B"),
#'   species2 = c("SP_B", "SP_C", "SP_C"),
#'   hog = "HOG1", q.value = c(0.01, 0.02, 0.03)
#' )
#' cl <- gene_clique_graph(edges)
#' classify_gene_cliques(cl, edges, c("SP_A", "SP_B", "SP_C"))
#'
#' @seealso \code{\link{gene_clique_graph}}
#' @references
#' Rodriguez E, Birkeland S, Chapple ED, et al. (2026).
#' Comparative regulomics of wood formation across dicot and
#' conifer trees. \emph{Nature Communications} 17(1).
#' \doi{10.1038/s41467-026-75624-2}
#' @param ... Additional arguments passed to the default method.
#' @export
classify_gene_cliques <- function(cliques, ...) {
  UseMethod("classify_gene_cliques")
}

#' @rdname classify_gene_cliques
#' @export
classify_gene_cliques.default <- function(cliques, edges, species,
                                          lineage = NULL, alpha_call = 0.1,
                                          alpha_graph = 0.9, max_gap = 1L,
                                          cross_max = NULL,
                                          min_power = 0.8, ...) {
  rlang::check_dots_empty()
  need_cl <- c("clique_id", "hog", "species", "gene")
  absent <- setdiff(need_cl, names(cliques))
  if (length(absent) > 0L) {
    stop(
      "cliques missing required columns: ",
      paste(absent, collapse = ", ")
    )
  }
  need_ed <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value"
  )
  absent <- setdiff(need_ed, names(edges))
  if (length(absent) > 0L) {
    stop(
      "edges missing required columns: ",
      paste(absent, collapse = ", ")
    )
  }
  species <- as.character(species)
  if (length(species) < 2L || anyDuplicated(species) > 0L) {
    stop("species must be at least 2 unique species names")
  }
  # The only numeric parameters that were unchecked. An NA alpha makes
  # every `q < alpha` comparison NA, so n_sig is NA, the tier predicates
  # evaluate to NA, and the call dies inside a helper with base R's
  # "missing value where TRUE/FALSE needed" -- naming neither argument.
  for (nm in c("alpha_call", "alpha_graph")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(nm, " must be a single non-missing number")
    }
  }
  ok_power <- is.numeric(min_power) && length(min_power) == 1L &&
    !is.na(min_power) && min_power >= 0 && min_power <= 1
  if (!ok_power) {
    stop("min_power must be a single number in [0, 1]")
  }
  cl_sp <- as.character(cliques$species)
  # A clique species outside `species` is counted into the clique's own
  # species total but never into choose(S, 2), so the clique would be
  # scored complete while the same row reports a real species missing.
  stray <- setdiff(unique(cl_sp), species)
  if (length(stray) > 0L) {
    stop(
      "cliques contain species absent from `species`: ",
      paste(stray, collapse = ", ")
    )
  }
  max_gap <- as.integer(max_gap)
  if (length(max_gap) != 1L || is.na(max_gap) || max_gap < 0L) {
    stop("max_gap must be a single non-negative integer")
  }
  .gcg_check_ids(edges)

  n_sp <- length(species)
  n_pair <- choose(n_sp, 2)

  lin <- NULL
  if (!is.null(lineage)) {
    if (is.null(names(lineage))) {
      stop("lineage must be a named vector (names are species)")
    }
    miss <- setdiff(species, names(lineage))
    if (length(miss) > 0L) {
      stop(
        "lineage missing entries for: ",
        paste(miss, collapse = ", ")
      )
    }
    lin <- stats::setNames(as.character(lineage[species]), species)
    # An NA lineage passes the presence check above (the species is a
    # name in `lineage`, just with a missing value) but table(lin) drops
    # NA entries silently, so a clique containing that species would be
    # scored against undercounted lin_sizes / w_pairs / x_pairs instead of
    # failing loudly.
    if (anyNA(lin)) {
      stop(
        "lineage has missing values for: ",
        paste(species[is.na(lin)], collapse = ", ")
      )
    }
  }
  lin_sizes <- if (is.null(lin)) integer(0) else table(lin)
  w_pairs <- if (is.null(lin)) NA_real_ else sum(choose(lin_sizes, 2))
  x_pairs <- if (is.null(lin)) NA_real_ else n_pair - w_pairs
  if (is.null(cross_max)) {
    cross_max <- if (is.null(lin)) {
      NA_real_
    } else {
      max(0, choose(n_sp - 1L, 2) - w_pairs)
    }
  } else {
    ok_cross <- is.numeric(cross_max) && length(cross_max) == 1L &&
      !is.na(cross_max) && cross_max >= 0
    if (!ok_cross) {
      stop("cross_max must be a single non-negative number")
    }
  }

  has_effect <- "effect_size" %in% names(edges)
  ids <- unique(as.character(cliques$clique_id))
  const <- list(
    n_sp = n_sp, n_pair = n_pair, w_pairs = w_pairs,
    x_pairs = x_pairs, cross_max = cross_max,
    alpha_call = alpha_call, alpha_graph = alpha_graph,
    max_gap = max_gap
  )
  if (length(ids) == 0L) {
    return(.gcg_class_attrs(.gcg_empty_class(has_effect), const))
  }

  ekey1 <- paste(edges$species1, edges$gene1, sep = .gcg_sep)
  ekey2 <- paste(edges$species2, edges$gene2, sep = .gcg_sep)
  ehog <- as.character(edges$hog)
  # Duplicated undirected pairs collapse to their most significant row,
  # matching the graph the cliques were built on. Ties break on
  # effect_size descending, never on row order: q saturates at the
  # permutation floor and the survivor sets mean_effect_size.
  pkey <- paste(ehog, pmin(ekey1, ekey2), pmax(ekey1, ekey2),
    sep = .gcg_sep
  )
  ev_all <- if (has_effect) {
    as.numeric(edges$effect_size)
  } else {
    rep(NA_real_, nrow(edges))
  }
  ord <- order(pkey, edges$q.value, -ev_all, na.last = TRUE)
  uniq <- ord[!duplicated(pkey[ord])]
  lut_key <- pkey[uniq]
  lut_q <- as.numeric(edges$q.value[uniq])
  lut_e <- ev_all[uniq]
  pw_all <- if ("power" %in% names(edges)) {
    as.numeric(edges$power)
  } else {
    rep(NA_real_, nrow(edges))
  }
  lut_p <- pw_all[uniq]
  rows_by_hog <- split(seq_len(nrow(edges)), ehog)

  cl_hog <- as.character(cliques$hog)
  mk_all <- paste(cl_sp, as.character(cliques$gene), sep = .gcg_sep)
  cl_id <- as.character(cliques$clique_id)
  .gcg_check_clique_ids(cl_id, mk_all, cl_hog, cliques$n_members)

  cl_by_id <- split(seq_len(nrow(cliques)), factor(cl_id, levels = ids))
  cmb_l <- lapply(cl_by_id, function(rr) {
    if (length(rr) < 2L) {
      matrix(integer(0), nrow = 2L)
    } else {
      utils::combn(length(rr), 2L)
    }
  })
  # One match() over every clique's pair keys. Matching per clique
  # rehashes the whole edge table once per clique, which is O(cliques x
  # edges) and turns a real run into hours.
  qk_l <- lapply(seq_along(ids), function(i) {
    rr <- cl_by_id[[i]]
    cmb <- cmb_l[[i]]
    if (ncol(cmb) == 0L) {
      return(character(0))
    }
    a <- mk_all[rr[cmb[1L, ]]]
    b <- mk_all[rr[cmb[2L, ]]]
    paste(cl_hog[rr[1L]], pmin(a, b), pmax(a, b), sep = .gcg_sep)
  })
  hit_all <- match(unlist(qk_l, use.names = FALSE), lut_key)
  np <- vapply(qk_l, length, integer(1))
  ends <- cumsum(np)
  starts <- ends - np + 1L

  res <- lapply(seq_along(ids), function(i) {
    hit <- if (np[i] == 0L) integer(0) else hit_all[starts[i]:ends[i]]
    .gcg_classify_one(
      rr = cl_by_id[[i]], cmb = cmb_l[[i]], hit = hit,
      cl_hog = cl_hog, cl_sp = cl_sp, mk_all = mk_all, edges = edges,
      species = species, lin = lin, lin_sizes = lin_sizes,
      alpha_call = alpha_call, alpha_graph = alpha_graph,
      max_gap = max_gap, cross_max = cross_max, n_sp = n_sp,
      lut_q = lut_q, lut_e = lut_e, ekey1 = ekey1, ekey2 = ekey2,
      rows_by_hog = rows_by_hog, lut_p = lut_p, power = pw_all,
      min_power = min_power
    )
  })
  pick <- function(field, what) {
    vapply(res, function(z) z[[field]], what)
  }

  out <- data.frame(
    clique_id = ids,
    hog = pick("hog", character(1)),
    classification = pick("cls", character(1)),
    n_members = pick("m", integer(1)),
    n_species = pick("m_sp", integer(1)),
    n_pairs = pick("n_pairs", integer(1)),
    n_present = pick("n_present", integer(1)),
    n_sig = pick("n_sig", integer(1)),
    n_sig_within = pick("n_sig_w", integer(1)),
    n_sig_cross = pick("n_sig_x", integer(1)),
    n_underpowered_cross = pick("n_up_x", integer(1)),
    n_missing = pick("n_missing", integer(1)),
    missing_species = pick("missing_species", character(1)),
    missing_reason = pick("missing_reason", character(1)),
    mean_q = pick("mean_q", numeric(1)),
    max_q = pick("max_q", numeric(1)),
    stringsAsFactors = FALSE
  )
  if (has_effect) out$mean_effect_size <- pick("mean_e", numeric(1))

  # HOG-level precedence: the published scripts removed a whole
  # orthogroup from later tiers once any of its cliques matched.
  rank <- match(out$classification, .gcg_tiers)
  best <- vapply(split(rank, out$hog), min, numeric(1))
  out$hog_class <- .gcg_tiers[best[out$hog]]
  .gcg_class_attrs(out, const)
}


#' Attach the derived constants and floor diagnostics
#'
#' The floor diagnostics are columns as well as attributes: `subset()`
#' drops attributes, so filtering the result would otherwise discard the
#' record of how little resolution `mean_q` has as a ranking key.
#'
#' @param out Result frame (possibly zero-row).
#' @param const List of derived constants from the caller.
#' @return `out` with floor columns and attributes.
#' @noRd
.gcg_class_attrs <- function(out, const) {
  finite_q <- out$mean_q[is.finite(out$mean_q)]
  floor_q <- if (length(finite_q) == 0L) NA_real_ else min(finite_q)
  n_tied <- if (is.na(floor_q)) {
    0L
  } else {
    .tol_min_ties(out$mean_q)$n_at_min
  }
  out$mean_q_floor <- rep(floor_q, nrow(out))
  out$n_cliques_at_q_floor <- rep(n_tied, nrow(out))
  attr(out, "n_species_total") <- const$n_sp
  attr(out, "n_pairs_total") <- const$n_pair
  attr(out, "n_within_pairs") <- const$w_pairs
  attr(out, "n_cross_pairs") <- const$x_pairs
  attr(out, "cross_max") <- const$cross_max
  attr(out, "alpha_call") <- const$alpha_call
  attr(out, "alpha_graph") <- const$alpha_graph
  attr(out, "max_gap") <- const$max_gap
  attr(out, "mean_q_floor") <- floor_q
  attr(out, "n_cliques_at_q_floor") <- n_tied
  rownames(out) <- NULL
  out
}


#' Empty result template for [classify_gene_cliques()]
#' @noRd
.gcg_empty_class <- function(has_effect) {
  out <- data.frame(
    clique_id = character(0), hog = character(0),
    classification = character(0), n_members = integer(0),
    n_species = integer(0), n_pairs = integer(0),
    n_present = integer(0), n_sig = integer(0),
    n_sig_within = integer(0), n_sig_cross = integer(0),
    n_underpowered_cross = integer(0),
    n_missing = integer(0), missing_species = character(0),
    missing_reason = character(0), mean_q = numeric(0),
    max_q = numeric(0), stringsAsFactors = FALSE
  )
  if (has_effect) out$mean_effect_size <- numeric(0)
  out$hog_class <- character(0)
  out
}


#' Classify a single clique
#'
#' Split out of [classify_gene_cliques()] so the per-clique tier logic
#' stays readable; every argument is precomputed by the caller,
#' including `hit`, this clique's slice of the one global `match()`.
#'
#' @return A named list of scalars, one per output column.
#' @noRd
.gcg_classify_one <- function(rr, cmb, hit, cl_hog, cl_sp, mk_all,
                              edges, species, lin, lin_sizes,
                              alpha_call, alpha_graph, max_gap,
                              cross_max, n_sp, lut_q, lut_e, ekey1,
                              ekey2, rows_by_hog, lut_p, power,
                              min_power) {
  hog <- cl_hog[rr[1L]]
  msp <- cl_sp[rr]
  mk <- mk_all[rr]
  m <- length(rr)
  m_sp <- length(unique(msp))

  qp <- lut_q[hit]
  ep <- lut_e[hit]
  pp <- lut_p[hit]
  if (m >= 2L && !is.null(lin)) {
    l1 <- lin[msp[cmb[1L, ]]]
    l2 <- lin[msp[cmb[2L, ]]]
    within <- unname(l1 == l2)
    # A within-lineage pair is labelled by its (single) lineage;
    # cross-lineage pairs carry no label.
    pair_lin <- ifelse(within, unname(l1), NA_character_)
  } else {
    within <- rep(NA, length(qp))
    pair_lin <- rep(NA_character_, length(qp))
  }
  sig <- !is.na(qp) & qp < alpha_call
  n_pairs <- length(qp)
  n_present <- sum(!is.na(qp))
  n_sig <- sum(sig)
  n_sig_w <- if (is.null(lin)) NA_integer_ else sum(sig & within)
  n_sig_x <- if (is.null(lin)) NA_integer_ else sum(sig & !within)
  # Present, non-significant cross pairs that could not have been called.
  n_up_x <- if (is.null(lin)) {
    NA_integer_
  } else {
    sum(!is.na(qp) & !sig & !within & !is.na(pp) & pp < min_power)
  }
  max_q <- if (n_present == 0L) NA_real_ else max(qp, na.rm = TRUE)
  mean_q <- if (n_present == 0L) NA_real_ else mean(qp, na.rm = TRUE)
  mean_e <- if (all(is.na(ep))) NA_real_ else mean(ep, na.rm = TRUE)

  gone <- setdiff(species, msp)
  rows <- rows_by_hog[[hog]]
  if (is.null(rows)) rows <- integer(0)
  reason <- vapply(gone, .gcg_missing_reason, character(1),
    rows = rows, edges = edges, key1 = ekey1,
    key2 = ekey2, mk = mk, m = m,
    alpha_call = alpha_call, power = power,
    min_power = min_power, USE.NAMES = FALSE
  )
  # Only "absent" and "untested" are annotation gaps. Admitting
  # "tested_ns" would let partial_present claim a clique whose missing
  # species was in fact rejected, and would let lineage_specific claim
  # one whose outside species were compared against every member and
  # diverged -- evidence of a boundary, not a gap, and the case
  # `differentiated` exists to score. Folding it in here would make the
  # two tiers indistinguishable. The published workflow draws the same
  # line: its dicot- and conifer-specific sets require Cross == 0 in a
  # matrix whose 1s mark a *tested* species pair, not a significant one,
  # while its differentiated set requires both lineages present.
  gap_only <- length(gone) == 0L ||
    all(reason %in% c("absent", "untested"))
  # An underpowered outside species was tested, so it is no gap, but its
  # failure is no boundary either: a lineage-specific call has to survive
  # reading it as conserved, which it cannot.
  up_only <- length(gone) > 0L && any(reason == "underpowered") &&
    all(reason %in% c("absent", "untested", "underpowered"))
  # partial_present only asks that no missing species was *rejected*. An
  # underpowered one is unknown, not rejected, so it counts as a gap here;
  # lineage_specific keeps gap_only, because its call would have to
  # survive reading that species as conserved (up_only above).
  gap_pp <- length(gone) == 0L ||
    all(reason %in% c("absent", "untested", "underpowered"))

  cls <- "unclassified"
  # One gene per species is the invariant the pair arithmetic rests on;
  # a within-species edge in the input breaks it, so refuse to score.
  # A one-member "clique" has choose(1, 2) == 0 pairs, which every tier
  # count would vacuously satisfy; it carries no conserved edge at all.
  if (m_sp == m && m_sp >= 2L) {
    lin_m <- if (is.null(lin)) NULL else unique(lin[msp])
    n_l <- if (is.null(lin) || length(lin_m) != 1L) {
      NA_integer_
    } else {
      as.integer(lin_sizes[[lin_m]])
    }
    diff_ok <- FALSE
    # n_present == n_pairs: an untested cross-lineage pair is absent
    # evidence, not evidence of divergence.
    if (!is.null(lin) && m_sp == n_sp && n_present == n_pairs) {
      big <- names(lin_sizes)[lin_sizes >= 2L]
      full <- vapply(big, function(g) {
        sum(sig & !is.na(pair_lin) & pair_lin == g) ==
          choose(as.integer(lin_sizes[[g]]), 2L)
      }, logical(1))
      # all(logical(0)) is TRUE, so with every lineage a singleton the
      # tier would collapse to "few cross pairs significant" and pass a
      # clique with no significant pair at all.
      diff_ok <- length(big) >= 1L && all(full) && n_sig_x <= cross_max
    }
    gap <- n_sp - m_sp
    is_complete <- m_sp == n_sp && n_sig == choose(n_sp, 2)
    lin_core <- !is.na(n_l) && m_sp == n_l && n_sig == choose(n_l, 2)
    is_lineage <- lin_core && gap_only
    # choose(S - 1, 2) + 1 == choose(S, 2) - (S - 2): the tolerance is
    # S - 2 non-significant edges, the largest that cannot isolate a
    # member, since cutting one loose needs all S - 1 of its edges.
    is_part_sig <- m_sp == n_sp && n_present == choose(n_sp, 2) &&
      !is.na(max_q) && max_q < alpha_graph &&
      n_sig >= choose(n_sp - 1L, 2) + 1
    is_part_pres <- gap >= 1L && gap <= max_gap &&
      n_sig == choose(m_sp, 2) && gap_pp
    if (is_complete) {
      cls <- "complete_conserved"
    } else if (is_lineage) {
      cls <- "lineage_specific"
    } else if (lin_core && up_only) {
      cls <- "underpowered"
    } else if (is_part_sig) {
      cls <- "partial_significant"
    } else if (is_part_pres) {
      cls <- "partial_present"
    } else if (diff_ok) {
      # The call must survive treating every underpowered cross pair as
      # possibly significant.
      cls <- if (n_sig_x + n_up_x > cross_max) {
        "underpowered"
      } else {
        "differentiated"
      }
    }
  }

  list(
    hog = hog, cls = cls, m = m, m_sp = m_sp, n_pairs = n_pairs,
    n_present = n_present, n_sig = n_sig, n_sig_w = n_sig_w,
    n_sig_x = n_sig_x, n_up_x = n_up_x, n_missing = length(gone),
    missing_species = paste(gone, collapse = ","),
    missing_reason = paste(reason, collapse = ","),
    mean_q = mean_q, max_q = max_q, mean_e = mean_e
  )
}
