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


#' Empty result template for [gene_clique_graph()]
#'
#' @param has_effect Whether an `effect_size` column is carried.
#' @param has_score Whether a `score` column is carried.
#' @param sp_lv Species levels of the `members` table.
#' @return The two zero-row tables with their full column sets.
#' @noRd
.gcg_empty <- function(has_effect, has_score, sp_lv) {
  cliques <- data.frame(
    clique_id = integer(0), hog = character(0),
    n_members = integer(0), n_species = integer(0),
    n_edges = integer(0), mean_q = numeric(0), max_q = numeric(0)
  )
  if (has_effect) cliques$mean_effect_size <- numeric(0)
  if (has_score) cliques$score <- numeric(0)
  members <- data.frame(
    clique_id = integer(0), species = factor(character(0), sp_lv),
    gene = character(0)
  )
  list(cliques = cliques, members = members)
}


#' Attach the floor diagnostics to a [gene_clique_graph()] result
#'
#' Carried as columns as well as attributes: `subset()` drops attributes,
#' so a caller who filters the result would otherwise lose the record of
#' how little resolution `mean_q` has.
#'
#' @param res List of the `cliques` and `members` tables.
#' @param alpha_graph,min_size Call parameters.
#' @param q_floor,mean_q_floor,n_tied Floor diagnostics.
#' @return `res` of class `gene_cliques`, with floor columns and
#'   attributes on `cliques`.
#' @noRd
.gcg_graph_attrs <- function(res, alpha_graph, min_size, q_floor,
                             mean_q_floor, n_tied) {
  out <- res$cliques
  out$mean_q_floor <- rep(mean_q_floor, nrow(out))
  out$n_cliques_at_q_floor <- rep(n_tied, nrow(out))
  attr(out, "alpha_graph") <- alpha_graph
  attr(out, "min_size") <- min_size
  attr(out, "q_floor") <- q_floor
  attr(out, "mean_q_floor") <- mean_q_floor
  attr(out, "n_cliques_at_q_floor") <- n_tied
  res$cliques <- out
  structure(res, class = "gene_cliques")
}


#' Statistics of one block of cliques of one HOG
#'
#' The edges of each clique are visited in edge order, as the edge scan
#' of the reference did, so the sums come out in the same order. The
#' mean is the two-pass mean of `mean.default()`.
#'
#' @param cl Cliques as integer vectors of node indices.
#' @param em Node x node matrix of edge indices.
#' @param eq,ee,es q-value, effect size and score of each edge.
#' @param sp_h Species index of each node.
#' @param nodes Global node id of each node.
#' @return A list of per-clique statistics and the members' node ids.
#' @noRd
.gcg_clique_stats <- function(cl, em, eq, ee, es, sp_h, nodes) {
  n_k <- length(cl)
  sz <- lengths(cl)
  nd <- unlist(cl, use.names = FALSE)
  cq <- rep.int(seq_len(n_k), sz)
  n_right <- sz[cq] - sequence(sz)
  left <- rep.int(seq_along(nd), n_right)
  e_k <- em[cbind(nd[left], nd[left + sequence(n_right)])]
  o <- order(cq[left], e_k)
  e_h <- e_k[o]
  g <- cq[left][o]
  grp <- structure(g, N.groups = n_k, class = "qG")
  first_sp <- !duplicated((cq - 1L) * max(sp_h) + sp_h[nd])
  list(
    sz = sz, ns = tabulate(cq[first_sp], nbins = n_k),
    ne = tabulate(g, nbins = n_k), mq = .gcg_gmean(eq[e_h], g, n_k),
    xq = collapse::fmax(eq[e_h], grp, na.rm = FALSE, use.g.names = FALSE),
    me = .gcg_gmean(ee[e_h], g, n_k),
    sc = collapse::fsum(es[e_h], grp, na.rm = FALSE, use.g.names = FALSE),
    gi = nodes[nd]
  )
}


#' Grouped two-pass mean, as `mean.default()` computes it
#'
#' Values are summed in input order within each group, as `mean()` sums
#' them. An empty group gives `NA`. R's `mean()` sums in long double where
#' the platform has one, so the last bit can differ from this helper on
#' x86_64.
#'
#' @param x Values.
#' @param g Group of each value, in `seq_len(n_k)`.
#' @param n_k Number of groups.
#' @return One mean per group.
#' @noRd
.gcg_gmean <- function(x, g, n_k) {
  ne <- tabulate(g, nbins = n_k)
  if (length(x) == 0L) {
    return(rep(NA_real_, n_k))
  }
  grp <- structure(g, N.groups = n_k, class = "qG")
  gsum <- function(v) {
    collapse::fsum(v, grp, na.rm = FALSE, use.g.names = FALSE)
  }
  s <- gsum(x) / ne
  fin <- is.finite(s)
  s[fin] <- s[fin] + gsum(x - s[g])[fin] / ne[fin]
  s[ne == 0L] <- NA_real_
  s
}


#' Maximal cliques of the per-orthogroup gene graph
#'
#' Builds one graph per hog from the co-expressolog calls in `edges`.
#' Nodes are (species, gene) pairs. Lists every maximal clique of at least
#' `min_size` nodes, one per paralog combination.
#'
#' [find_cliques()] cliques the species graph instead and returns a single
#' best gene assignment per species clique. This function is the
#' computation used by Netotea-style clique workflows.
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
#' The clique `score` adds the edge scores. Log-odds add, so the sum is
#' the log-likelihood ratio of the conserved subnetwork (NetworkBLAST),
#' and it grows with clique size on purpose.
#'
#' A hog of 40 or more genes is enumerated one start gene at a time.
#' The pieces are disjoint, so the cliques are the same. Only the peak
#' memory drops.
#'
#' @section Duplicate rows:
#' Two rows describing the same undirected pair within one HOG collapse
#' to the more significant of the two, broken on `effect_size`
#' descending when the q-values tie -- q-values saturate at the
#' permutation floor, so a tie there must not be settled by input row
#' order. This is a duplicate-row collapse, not a direction combine:
#' [find_coexpressologs()] already merges the two comparison
#' directions with `pval_combine` (default `"max"`, the reciprocal
#' criterion) and emits one row per pair. A caller who row-binds the two
#' directional tables instead gets `"min"` semantics here, so combine
#' upstream if the reciprocal criterion is wanted.
#'
#' @param edges Data frame of co-expressolog calls, as returned by
#'   [find_coexpressologs()]: columns `gene1`, `gene2`,
#'   `species1`, `species2`, `hog` and `q_value`. The `effect_size`
#'   and `score` columns are used when present.
#' @param min_size Minimum number of nodes in a reported clique
#'   (default 3, matching the published workflow).
#' @param alpha_graph Edges with `q_value < alpha_graph` build the
#'   graph. Use the calling threshold (e.g. 0.1) for complete cliques,
#'   a permissive value (0.9) for partially-significant cliques, and a
#'   value **above** 1 (`Inf`) for an unfiltered graph. The comparison
#'   is strict, so `alpha_graph = 1` drops every edge whose `q_value` is
#'   exactly 1 -- not a corner case, since BH q-values cap at 1 and 20 of
#'   511 module q-values sit there on the package's own worked
#'   example. Passing several thresholds and combining
#'   the results is the intended way to feed
#'   [classify_gene_cliques()], since a clique that is maximal
#'   at one threshold need not be maximal at another. Pass the runs to
#'   it as a list.
#'
#' @return A list of class `gene_cliques` with two tables.
#'   `cliques` has one row per clique:
#'   \describe{
#'     \item{clique_id}{Integer clique id, 1 to the number of cliques}
#'     \item{hog}{Ortholog group}
#'     \item{n_members}{Nodes in the clique}
#'     \item{n_species}{Distinct species in the clique (equal to
#'       `n_members` for cross-species edge tables)}
#'     \item{n_edges}{Graph edges inside the clique}
#'     \item{mean_q, max_q}{q-value summaries over those edges}
#'     \item{mean_effect_size}{Present only when `edges` carries
#'       `effect_size`. Prefer it over `mean_q` for ranking: q-values
#'       saturate at the permutation floor and cannot separate cliques
#'       once they get there}
#'     \item{score}{Present only when `edges` carries `score`. The sum
#'       of the clique edges' `score`, in bits}
#'     \item{mean_q_floor, n_cliques_at_q_floor}{Smallest clique
#'       `mean_q` in the run and how many cliques are tied there. A
#'       large tie count means `mean_q` cannot rank those cliques at
#'       all}
#'   }
#'   `members` has one row per clique member: `clique_id`, `species`
#'   (a factor; its levels are the species of `edges` in order of first
#'   appearance) and `gene`.
#'   The two floor diagnostics, plus `alpha_graph`, `min_size` and
#'   `q_floor` (smallest edge q-value in the graph), are also attributes
#'   of `cliques`. Both are always set, empty results included.
#'
#' @examples
#' edges <- data.frame(
#'   gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
#'   species1 = c("SP_A", "SP_A", "SP_B"),
#'   species2 = c("SP_B", "SP_C", "SP_C"),
#'   hog = "HOG1", q_value = c(0.01, 0.02, 0.03)
#' )
#' g <- gene_clique_graph(edges)
#' g$cliques
#' g$members
#'
#' @seealso [classify_gene_cliques()],
#'   [find_cliques()]
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
                                      alpha_graph = 0.1, ...) {
  rlang::check_dots_empty()
  required <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q_value"
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
  .gcg_check_ids(edges)

  has_effect <- "effect_size" %in% names(edges)
  has_score <- "score" %in% names(edges)
  sp_lv <- unique(as.character(c(edges$species1, edges$species2)))
  keep <- !is.na(edges$q_value) & edges$q_value < alpha_graph
  edges <- edges[keep, , drop = FALSE]
  if (nrow(edges) == 0L) {
    return(.gcg_graph_attrs(
      .gcg_empty(has_effect, has_score, sp_lv), alpha_graph, min_size,
      NA_real_, NA_real_, 0L
    ))
  }

  # One integer id per (species, gene) node: the first row end that
  # names it, over both ends of every row.
  n_row <- nrow(edges)
  all_sp <- match(as.character(c(edges$species1, edges$species2)), sp_lv)
  all_gn <- as.character(c(edges$gene1, edges$gene2))
  gid <- match(
    paste(all_sp, all_gn, sep = .gcg_sep),
    paste(all_sp, all_gn, sep = .gcg_sep)
  )
  gid1 <- gid[seq_len(n_row)]
  gid2 <- gid[n_row + seq_len(n_row)]
  qv <- as.numeric(edges$q_value)
  ev <- if (has_effect) {
    as.numeric(edges$effect_size)
  } else {
    rep(NA_real_, n_row)
  }
  sv <- if (has_score) as.numeric(edges$score) else rep(NA_real_, n_row)

  by_hog <- split(seq_len(n_row), as.character(edges$hog))
  blk <- vector("list", length(by_hog))
  # Plain integer cliques skip igraph's vertex-sequence objects.
  old_opt <- igraph::igraph_options(return.vs.es = FALSE)
  on.exit(igraph::igraph_options(old_opt), add = TRUE)
  n_capped <- 0L
  n_dropped <- 0L

  for (k in seq_along(by_hog)) {
    h <- names(by_hog)[k]
    idx <- by_hog[[k]]
    # Both endpoints on the same node is a self-loop, which igraph would
    # happily accept and which would inflate n_edges. A within-species
    # row between two *distinct* paralogs is a real edge and is kept;
    # it surfaces downstream as n_species < n_members.
    idx <- idx[gid1[idx] != gid2[idx]]
    if (length(idx) == 0L) next
    nodes <- unique(c(gid1[idx], gid2[idx]))
    n_v <- length(nodes)
    if (n_v < min_size) next
    node_sp <- all_sp[nodes]

    i1 <- match(gid1[idx], nodes)
    i2 <- match(gid2[idx], nodes)
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
    es <- sv[idx][sel]

    # Paralog copy cap, as in find_cliques(): every paralog combination
    # is its own clique, so c copies per species in a near-complete
    # S-partite graph give up to c^S maximal cliques. A species above
    # the cap keeps its most-connected genes, counted over the collapsed
    # edges that enter the graph; order() is stable, so ties keep the
    # gene that entered the graph first. Dropped genes stay as isolated
    # vertices, which max_cliques(min >= 2) never reports.
    deg <- tabulate(c(i1, i2), nbins = n_v)
    keep_node <- rep(TRUE, n_v)
    for (s in which(tabulate(node_sp) > 10L)) {
      v <- which(node_sp == s)
      top <- v[order(-deg[v])[seq_len(10L)]]
      keep_node[setdiff(v, top)] <- FALSE
    }
    if (!all(keep_node)) {
      n_capped <- n_capped + 1L
      n_dropped <- n_dropped + sum(!keep_node)
      ke <- keep_node[i1] & keep_node[i2]
      i1 <- i1[ke]
      i2 <- i2[ke]
      eq <- eq[ke]
      ee <- ee[ke]
      es <- es[ke]
    }

    g <- igraph::make_graph(as.vector(rbind(i1, i2)),
      n = n_v, directed = FALSE
    )
    em <- matrix(0L, n_v, n_v)
    em[cbind(i1, i2)] <- em[cbind(i2, i1)] <- seq_along(i1)
    sp_h <- match(node_sp, unique(node_sp))
    # A few hogs hold most cliques. One call per start vertex keeps one
    # piece in memory at a time; the pieces are disjoint and together
    # complete only when every vertex id starts one.
    starts <- if (sum(keep_node) >= 40L) seq_len(n_v) else list(NULL)
    out_h <- list()
    for (v in starts) {
      cl <- igraph::max_cliques(g, min = min_size, subset = v)
      for (b in split(seq_along(cl), ceiling(seq_along(cl) / 2e5))) {
        st <- .gcg_clique_stats(cl[b], em, eq, ee, es, sp_h, nodes)
        st$hog <- h
        out_h[[length(out_h) + 1L]] <- st
      }
    }
    blk[[k]] <- out_h
  }

  if (n_capped > 0L) {
    message(
      "gene_clique_graph: ", n_capped, " ortholog groups exceeded ",
      "10 genes in some species; ",
      n_dropped, " genes dropped"
    )
  }

  blk <- unlist(blk, recursive = FALSE, use.names = FALSE)
  if (length(blk) == 0L) {
    return(.gcg_graph_attrs(
      .gcg_empty(has_effect, has_score, sp_lv), alpha_graph, min_size,
      min(qv), NA_real_, 0L
    ))
  }
  hog <- rep.int(
    vapply(blk, `[[`, "", "hog"), vapply(blk, function(b) length(b$sz), 1L)
  )
  # Each column leaves the blocks as it is joined, so the blocks and the
  # result never both hold it: that sets the peak memory.
  take <- function(f) unlist(lapply(blk, `[[`, f), use.names = FALSE)
  gi <- take("gi")
  for (i in seq_along(blk)) blk[[i]]$gi <- NULL
  gene <- all_gn[gi]
  species <- structure(all_sp[gi], levels = sp_lv, class = "factor")
  rm(gi)
  sz <- take("sz")
  n_k <- length(sz)
  members <- data.frame(
    clique_id = rep.int(seq_len(n_k), sz), species = species, gene = gene
  )
  rm(species, gene)
  col <- list()
  for (f in c("ns", "ne", "mq", "xq", "me", "sc")) {
    col[[f]] <- take(f)
    for (i in seq_along(blk)) blk[[i]][[f]] <- NULL
  }
  rm(blk)
  cliques <- data.frame(
    clique_id = seq_len(n_k), hog = hog, n_members = sz,
    n_species = col$ns, n_edges = col$ne, mean_q = col$mq, max_q = col$xq
  )
  if (has_effect) cliques$mean_effect_size <- col$me
  if (has_score) cliques$score <- col$sc
  rm(col)
  # Tie counts go through the tolerant comparison for the same reason
  # pvalue_resolution() does: mean_q is a mean over a different edge
  # subset per clique, so two mathematically equal values need not be
  # bitwise equal, and a tie count that exists to say "these cannot be
  # ranked" must not under-report the tie. The rule is .tol_min_ties()'s;
  # its count of distinct values sorts every clique, and is not needed.
  mq_min <- min(cliques$mean_q)
  n_tied <- sum(cliques$mean_q <= mq_min + .tie_tol() * abs(mq_min))
  .gcg_graph_attrs(
    list(cliques = cliques, members = members), alpha_graph, min_size,
    min(qv), mq_min, n_tied
  )
}


#' Tier order of the classification waterfall
#' @noRd
.gcg_tiers <- c(
  "complete_conserved", "lineage_specific",
  "partial_significant", "partial_present",
  "differentiated", "trait_specific", "underpowered", "unclassified"
)


#' Classify gene-graph cliques into conservation tiers
#'
#' Applies the five-tier taxonomy of Rodriguez et al. (2026), plus a
#' `trait_specific` tier, to the cliques from [gene_clique_graph()].
#' Every threshold follows from the number of species supplied.
#'
#' No threshold is tied to the six species and fifteen species-pairs of
#' the original workflow.
#'
#' The home clade of a clique is the smallest clade that holds all its
#' species. With no such clade, the home is the whole species set. The
#' child clades of the home are the largest clades inside it. A species
#' of the home in no child clade forms its own child clade. The
#' `lineage_specific`, `trait_specific` and `differentiated` tiers, and
#' the within and cross counts, read the home and its child clades. So
#' a clique can be specific to, or differentiated inside, a nested clade.
#'
#' With `S = length(species)` species, `P = choose(S, 2)` species-pairs,
#' top-level clade sizes `n_l`, `W = sum(choose(n_l, 2))` within-clade
#' pairs and `X = P - W` cross-clade pairs, the tiers are:
#' \describe{
#'   \item{complete_conserved}{All `S` species present and all
#'     `choose(S, 2)` pairs significant at `alpha_call`.}
#'   \item{lineage_specific}{A complete clique over one entire clade,
#'     with no species outside it testable against the clique. A species
#'     that *was* compared against every member and came back
#'     non-significant is evidence of a boundary rather than a gap, so
#'     it blocks this tier: that clique is `trait_specific`, or with both
#'     clades conserved within a candidate for `differentiated` on the
#'     cliques of the unfiltered graph.}
#'   \item{partial_significant}{All `S` species present, every clique
#'     edge below `alpha_graph`, and at least `choose(S - 1, 2) + 1`
#'     pairs significant at `alpha_call`.}
#'   \item{partial_present}{`S - g` species present for
#'     `1 <= g <= max_gap`, all `choose(S - g, 2)` pairs significant,
#'     and every absent species an annotation gap or `underpowered`
#'     rather than a rejected test.}
#'   \item{differentiated}{Every species of the home present, every pair
#'     tested, at least one child clade of two or more species, every
#'     such child clade fully significant within itself, and at most
#'     `choose(H - 1, 2) - W_H` cross-clade pairs significant. Here `H`
#'     is the home size and `W_H` its within-child pairs. The child
#'     clades are disjoint, so a nested pair is never `differentiated`.
#'     With the whole species set as home, the bound is `cross_max`.}
#'   \item{trait_specific}{A complete clique over one entire clade, at
#'     least one outside species compared against every member and
#'     rejected at adequate power (`tested_ns`), no outside species
#'     `underpowered`, and no complete clique of a disjoint clade in the
#'     same HOG. Where `lineage_specific` reads the other clade's
#'     absence as a gap, this reads its presence and rejection as a
#'     boundary: one clade conserved, the other present but not
#'     co-conserved. Two clades each conserved within and rejected
#'     across is `differentiated`, scored on the unfiltered graph's
#'     clique, and their one-clade cliques stay `unclassified`. This
#'     tier is rcomplex's addition to the published five.}
#'   \item{underpowered}{A clique that would be `lineage_specific`,
#'     `trait_specific` or `differentiated` but for tests that could not
#'     have succeeded, read from a `power` column in `edges` (see
#'     \code{comparison_to_edges()}). It takes the place of the
#'     two specificity tiers when at least one outside species is
#'     `underpowered` -- reading that species as conserved would extend
#'     the clique, so neither call survives -- and of
#'     `differentiated` when `n_sig_cross + n_underpowered_cross`
#'     exceeds its bound: a specificity or divergence call must survive
#'     treating every underpowered pair as possibly significant. A
#'     low-degree gene cannot reach the call whatever its conservation,
#'     so without this its missing edges read as a clade boundary.}
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
#' `edges`. Without that, a clique whose cross-clade pairs were never
#' tested would be scored as diverged on absent evidence -- the same
#' conflation `partial_present` refuses through `missing_reason`.
#'
#' The classifier works on the sparse gene x clique incidence, in blocks
#' of 100,000 cliques. Sparse products with the edge rows give each
#' missing species its state.
#'
#' @param cliques A [gene_clique_graph()] result, a list of them, or a
#'   table of member rows with `clique_id`, `hog`, `species` and `gene`
#'   columns. Combine runs at several `alpha_graph` values to expose
#'   every tier: a clique complete at `alpha_call` need not be maximal
#'   on a looser graph. Pass the runs as a list. The result then gains a
#'   `run` column, the position of the run in the list.
#' @param edges The full, unfiltered co-expressolog table. It must not
#'   be pre-filtered on `q_value`: the gap tier needs to see rows that
#'   were tested and failed in order to refuse them. An optional `power`
#'   column (from \code{comparison_to_edges()}) enables the
#'   `underpowered` tier; without it, or where it is `NA`, the
#'   classification is unchanged.
#' @param species Character vector of every species in the analysis.
#'   Every species appearing in `cliques` must be listed; a stranger
#'   would be counted into the clique's species total while also being
#'   reported as missing.
#' @param clades Optional named list of species vectors, one per clade.
#'   Clades may nest but must not cross. A species in no clade forms its
#'   own clade. A flat vector `x` of species to clade is
#'   `split(names(x), x)`. Required for `lineage_specific`,
#'   `trait_specific` and `differentiated`. When `NULL`, these tiers are
#'   skipped.
#' @param alpha_call Significance threshold for calling a species pair
#'   conserved (default 0.1).
#' @param alpha_graph Loose threshold defining the `partial_significant`
#'   graph (default 0.9).
#' @param min_power Detection power below which a non-significant pair
#'   is read as uninformative rather than as evidence against
#'   conservation (default 0.8). Only used when `edges` has `power`. For
#'   rank-test edges (`find_coexpressologs(method = "rank")`) use 0.9: under
#'   the default reference rank `p0` their power is at least 0.5 by
#'   construction (0 when a direction of the species pair has no call under
#'   `pval_combine = "max"`, or neither has under `"min"`), and it
#'   overstates detection.
#'
#' @return A data frame with one row per clique:
#'   \describe{
#'     \item{run}{Only for a list of runs: the run of the clique}
#'     \item{clique_id, hog}{Clique identity. The id is the run's own
#'       `clique_id`.}
#'     \item{classification}{Tier, or `"unclassified"`}
#'     \item{hog_class}{Earliest tier reached by any clique of that HOG,
#'       reproducing the HOG-level precedence of the published scripts
#'       without discarding the per-clique detail}
#'     \item{n_members, n_species}{Clique size}
#'     \item{n_pairs, n_present, n_sig}{Member pairs, pairs with a row
#'       in `edges`, and pairs significant at `alpha_call`}
#'     \item{n_sig_within, n_sig_cross}{Significant pairs within one
#'       child clade of the clique's home, and across two (`NA` without
#'       `clades`)}
#'     \item{n_underpowered_cross}{Cross-clade pairs present in
#'       `edges`, not significant, and with `power` below `min_power`
#'       (`NA` without `clades`)}
#'     \item{clade}{Name of the smallest clade that holds every species
#'       of the clique. `NA` when no clade holds them all, or without
#'       `clades`.}
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
#'   hog = "HOG1", q_value = c(0.01, 0.02, 0.03)
#' )
#' cl <- gene_clique_graph(edges)
#' classify_gene_cliques(cl, edges, c("SP_A", "SP_B", "SP_C"))
#'
#' @seealso [gene_clique_graph()]
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
                                          clades = NULL, alpha_call = 0.1,
                                          alpha_graph = 0.9,
                                          min_power = 0.8, ...) {
  rlang::check_dots_empty()
  mem <- .gcg_members(cliques)
  need_ed <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q_value"
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
  # A clique species outside `species` is counted into the clique's own
  # species total but never into choose(S, 2), so the clique would be
  # scored complete while the same row reports a real species missing.
  # A factor is checked on its levels, without a copy per member.
  sp_m <- mem$sp
  stray <- if (is.factor(sp_m)) {
    lv_out <- is.na(match(levels(sp_m), species))
    anyNA(sp_m) ||
      (any(lv_out) && any(tabulate(sp_m, nlevels(sp_m))[lv_out] > 0L))
  } else {
    anyNA(match(as.character(sp_m), species))
  }
  if (stray) {
    s <- as.character(sp_m)
    stop(
      "cliques contain species absent from `species`: ",
      paste(unique(s[is.na(match(s, species))]), collapse = ", ")
    )
  }
  max_gap <- 1L
  .gcg_check_ids(edges)

  n_sp <- length(species)
  n_pair <- choose(n_sp, 2)

  lin <- NULL
  parts <- NULL
  if (!is.null(clades)) {
    clades <- .check_clades(clades, species)
    lin <- .clade_groups(clades, species)
    # Child clades of every clade, and of the whole species set first.
    parts <- lapply(c(list(species), clades), function(h) {
      inner <- vapply(clades, function(v) {
        length(v) < length(h) && all(v %in% h)
      }, logical(1))
      .clade_groups(clades[inner], h)
    })
    names(parts) <- c("", names(clades))
  }
  lin_sizes <- if (is.null(lin)) integer(0) else table(lin)
  w_pairs <- if (is.null(lin)) NA_real_ else sum(choose(lin_sizes, 2))
  x_pairs <- if (is.null(lin)) NA_real_ else n_pair - w_pairs
  # choose(S - 1, 2) - W, the generalisation of the published cut; the
  # original six-species script used a looser hard-coded 6
  cross_max <- if (is.null(lin)) {
    NA_real_
  } else {
    max(0, choose(n_sp - 1L, 2) - w_pairs)
  }

  has_effect <- "effect_size" %in% names(edges)
  const <- list(
    n_sp = n_sp, n_pair = n_pair, w_pairs = w_pairs,
    x_pairs = x_pairs, cross_max = cross_max,
    alpha_call = alpha_call, alpha_graph = alpha_graph,
    max_gap = max_gap
  )
  n_k <- length(mem$ids)
  if (n_k == 0L) {
    out <- .gcg_empty_class(has_effect)
    out$clique_id <- mem$out_ids
    return(.gcg_class_attrs(.gcg_runs(out, mem), const))
  }

  ix <- .gcg_edge_index(edges, species, mem$hog, alpha_call, min_power)
  ix$clades <- clades
  ix$parts <- parts
  ix$alpha_graph <- alpha_graph
  ix$max_gap <- max_gap
  # Blocks of 1e5 cliques; members are sorted by clique.
  sz <- tabulate(mem$id, nbins = mem$id[length(mem$id)])[mem$ids]
  ends <- cumsum(as.numeric(sz))
  starts <- seq.int(1L, n_k, by = 100000L)
  cols <- NULL
  for (b in seq_along(starts)) {
    # R lets garbage grow with the live heap before it collects. With
    # the cliques of a large run held, that garbage alone would add
    # gigabytes to the peak, so it is freed before each block and after
    # the last.
    gc(verbose = FALSE)
    cb <- starts[b]:min(n_k, starts[b] + 99999L)
    rr <- seq(ends[cb[1L]] - sz[cb[1L]] + 1, ends[cb[length(cb)]])
    blk <- .gcg_classify_block(
      ix, mem$id[rr], mem$sp[rr], mem$gene[rr], mem$ids[cb], mem$hog[cb]
    )
    if (is.null(cols)) {
      cols <- lapply(blk, function(v) vector(typeof(v), n_k))
    }
    for (f in names(blk)) cols[[f]][cb] <- blk[[f]]
  }
  rm(blk)
  gc(verbose = FALSE)
  is_core <- cols$is_core
  cols$is_core <- NULL
  if (!has_effect) cols$mean_effect_size <- NULL
  out <- list2DF(c(list(clique_id = mem$out_ids, hog = mem$hog), cols))

  # A trait-specific call needs no disjoint clade to have a complete
  # clique of its own in the HOG: two lineages each conserved within and
  # rejected across is `differentiated`, scored on the unfiltered graph's
  # clique, and their one-lineage cliques stay unclassified.
  if (!is.null(lin)) {
    ts <- which(out$classification == "trait_specific")
    if (length(ts) > 0L) {
      core <- which(is_core)
      cores <- lapply(split(out$clade[core], out$hog[core]), unique)
      other <- vapply(ts, function(i) {
        mine <- clades[[out$clade[i]]]
        any(vapply(cores[[out$hog[i]]], function(k) {
          !any(clades[[k]] %in% mine)
        }, logical(1)))
      }, logical(1))
      out$classification[ts[other]] <- "unclassified"
    }
  }

  # HOG-level precedence: the published scripts removed a whole
  # orthogroup from later tiers once any of its cliques matched.
  hg <- match(out$hog, ix$hog_u)
  best <- collapse::fmin(match(out$classification, .gcg_tiers),
    structure(hg, N.groups = length(ix$hog_u), class = "qG"),
    use.g.names = FALSE
  )
  out$hog_class <- .gcg_tiers[best[hg]]
  .gcg_class_attrs(.gcg_runs(out, mem), const)
}


#' Node index and row relations of the edge table
#'
#' Nodes are (HOG, species, gene) triples over both ends of every row.
#' Row relations are sparse node x node matrices, one entry per row end,
#' duplicates summed: tested, failed, and failed at power below
#' `min_power`, each summed per species, and significant, as a pattern.
#' Duplicate rows of one undirected pair also collapse to a pair lookup.
#'
#' @inheritParams classify_gene_cliques
#' @param hogs HOGs of the cliques.
#' @return A list of the index, the relations and the pair lookup.
#' @noRd
.gcg_edge_index <- function(edges, species, hogs, alpha_call, min_power) {
  n_sp <- length(species)
  n_e <- nrow(edges)
  # Species of the analysis first, so a member's species index is also
  # its index here.
  sp_u <- unique(c(
    species, as.character(edges$species1), as.character(edges$species2)
  ))
  gn_u <- unique(c(as.character(edges$gene1), as.character(edges$gene2)))
  hog_u <- unique(c(as.character(edges$hog), hogs))
  n_s <- length(sp_u)
  n_g <- length(gn_u)
  key <- function(h, s, g) ((h - 1) * n_s + (s - 1)) * n_g + g
  eh <- match(as.character(edges$hog), hog_u)
  es <- match(
    c(as.character(edges$species1), as.character(edges$species2)), sp_u
  )
  ek <- key(
    c(eh, eh), es,
    match(c(as.character(edges$gene1), as.character(edges$gene2)), gn_u)
  )
  nodes <- unique(ek)
  n_nd <- length(nodes)
  ni <- match(ek, nodes)
  i <- ni[seq_len(n_e)]
  j <- ni[n_e + seq_len(n_e)]
  first <- match(nodes, ek)
  node_s <- es[first]
  node_h <- c(eh, eh)[first]

  # Duplicated undirected pairs collapse to their most significant row,
  # matching the graph the cliques were built on. Ties break on
  # effect_size descending, never on row order: q saturates at the
  # permutation floor and the survivor sets mean_effect_size.
  q <- as.numeric(edges$q_value)
  na <- rep(NA_real_, n_e)
  ev <- if ("effect_size" %in% names(edges)) {
    as.numeric(edges$effect_size)
  } else {
    na
  }
  pw <- if ("power" %in% names(edges)) as.numeric(edges$power) else na
  pk <- pmin(i, j) * (n_nd + 1) + pmax(i, j)
  ord <- order(pk, q, -ev, na.last = TRUE)
  u <- ord[!duplicated(pk[ord])]

  # A species is tested against a clique when a row joins one of its
  # genes to a member. A failed row keeps it out of the clique; if every
  # failed row had power below min_power, the failure is no rejection.
  # NA power keeps the old reading.
  sig <- !is.na(q) & q < alpha_call
  sym <- function(w) {
    Matrix::sparseMatrix(c(i[w], j[w]), c(j[w], i[w]),
      x = 1, dims = c(n_nd, n_nd)
    )
  }
  ins <- which(node_s <= n_sp)
  by_sp <- Matrix::sparseMatrix(node_s[ins], ins,
    x = 1, dims = c(n_sp, n_nd)
  )
  sg <- sym(sig)
  sg@x[] <- 1
  present <- matrix(FALSE, n_sp, length(hog_u))
  present[cbind(node_s[ins], node_h[ins])] <- TRUE
  list(
    species = species, hog_u = hog_u, gn_u = gn_u, key = key,
    nodes = nodes, node_s = node_s, n_nd = n_nd,
    lut_k = pk[u], lut_q = q[u], lut_e = ev[u], lut_p = pw[u],
    t_sp = by_sp %*% sym(rep(TRUE, n_e)), f_sp = by_sp %*% sym(!sig),
    u_sp = by_sp %*% sym(!sig & !is.na(pw) & pw < min_power),
    sg = sg, sg_deg = Matrix::colSums(sg), present = present,
    alpha_call = alpha_call,
    min_power = min_power
  )
}


#' Classify one block of cliques on the incidence
#'
#' The block's cliques form a sparse node x clique incidence. Its
#' products with the row relations give each missing species its state.
#' Member pairs, expanded in `combn()` order, give the pair counts. The
#' tier waterfall then runs on the block's vectors.
#'
#' @param ix Context from `.gcg_edge_index()`, plus `clades`, `parts`,
#'   `alpha_graph` and `max_gap`.
#' @param id,sp,gene Clique, species and gene of each member row, sorted
#'   by clique.
#' @param ids,hog The block's clique ids and their HOGs.
#' @return A list of output columns, plus `is_core`: the clique is one
#'   whole clade, fully significant.
#' @noRd
.gcg_classify_block <- function(ix, id, sp, gene, ids, hog) {
  species <- ix$species
  clades <- ix$clades
  n_sp <- length(species)
  nb <- length(ids)
  kb <- match(id, ids)
  sb <- if (is.factor(sp)) {
    match(levels(sp), species)[as.integer(sp)]
  } else {
    match(as.character(sp), species)
  }
  ch <- match(hog, ix$hog_u)
  nd <- match(
    ix$key(ch[kb], sb, match(as.character(gene), ix$gn_u)), ix$nodes
  )
  m <- tabulate(kb, nbins = nb)

  # Missing-species state: 1 absent, 2 untested, 3 tested_ns,
  # 4 underpowered, 5 extendable, 0 in the clique.
  kn <- !is.na(nd)
  mb <- Matrix::sparseMatrix(nd[kn], kb[kn], x = 1, dims = c(ix$n_nd, nb))
  mb@x[] <- 1
  tested <- as.matrix(ix$t_sp %*% mb) > 0
  fails <- as.matrix(ix$f_sp %*% mb)
  ups <- as.matrix(ix$u_sp %*% mb)
  # A gene significant against every member would enlarge the clique:
  # the input clique set was built on a looser graph, so this is a
  # threshold mismatch rather than an annotation gap. Distinct members
  # are counted, not rows. The product has up to the members' summed
  # degree entries per clique, so it runs in chunks of about 1e7, and
  # each chunk's garbage is freed before the next.
  ext <- matrix(FALSE, n_sp, nb)
  cost <- cumsum(as.vector(Matrix::crossprod(mb, ix$sg_deg)))
  for (cc in split(seq_len(nb), cost %/% 1e7)) {
    x <- ix$sg %*% mb[, cc, drop = FALSE]
    hit <- which(x@x >= rep.int(m[cc], diff(x@p)))
    s <- ix$node_s[x@i[hit] + 1L]
    col <- cc[findInterval(hit, x@p, left.open = TRUE)]
    ext[cbind(s, col)[which(s <= n_sp), , drop = FALSE]] <- TRUE
    rm(x)
    gc(verbose = FALSE, full = FALSE)
  }
  code <- ifelse(!ix$present[, ch, drop = FALSE], 1L,
    ifelse(!tested, 2L,
      ifelse(ext, 5L, ifelse(fails > 0 & ups == fails, 4L, 3L))
    )
  )
  code[cbind(sb, kb)] <- 0L
  grp <- collapse::group(lapply(seq_len(n_sp), function(s) code[s, ]))
  pat <- lapply(which(!duplicated(grp)), function(f) code[, f])
  labels <- c("absent", "untested", "tested_ns", "underpowered", "extendable")
  miss_sp <- vapply(pat, function(cc) {
    paste(species[cc > 0L], collapse = ",")
  }, "")[grp]
  miss_re <- vapply(pat, function(cc) {
    paste(labels[cc[cc > 0L]], collapse = ",")
  }, "")[grp]
  home <- rep(NA_character_, nb)
  if (!is.null(clades)) {
    home <- vapply(pat, function(cc) {
      .clade_home(clades, species[cc == 0L])
    }, "")[grp]
  }
  m_sp <- as.integer(colSums(code == 0L))
  n_tns <- colSums(code == 3L)
  n_up <- colSums(code == 4L)
  n_ext <- colSums(code == 5L)

  # Member pairs in combn() order within each clique.
  n_right <- m[kb] - sequence(m)
  left <- rep.int(seq_along(nd), n_right)
  right <- left + sequence(n_right)
  a <- nd[left]
  b <- nd[right]
  hit <- match(pmin(a, b) * (ix$n_nd + 1) + pmax(a, b), ix$lut_k)
  g <- kb[left]
  qp <- ix$lut_q[hit]
  pres <- !is.na(qp)
  sig <- pres & qp < ix$alpha_call
  n_present <- tabulate(g[pres], nbins = nb)
  n_sig <- tabulate(g[sig], nbins = nb)
  max_q <- rep(NA_real_, nb)
  if (any(pres)) {
    gp <- structure(g[pres], N.groups = nb, class = "qG")
    max_q <- collapse::fmax(qp[pres], gp, na.rm = FALSE, use.g.names = FALSE)
    max_q[n_present == 0L] <- NA_real_
  }
  ep <- ix$lut_e[hit]
  ok <- !is.na(ep)
  n_sig_w <- n_sig_x <- n_up_x <- rep(NA_integer_, nb)
  lin_core <- diff_ok <- diff_up <- rep(FALSE, nb)
  n_pairs <- (m * (m - 1L)) %/% 2L
  if (!is.null(clades)) {
    # The clique is scored against the child clades of its home clade,
    # or against the top-level clades when no clade holds it.
    p_lab <- t(vapply(ix$parts, function(l) {
      match(l[species], unique(l))
    }, integer(n_sp)))
    p_w <- vapply(ix$parts, function(l) sum(choose(table(l), 2)), 1)
    p_cmax <- pmax(0, choose(lengths(ix$parts) - 1L, 2) - p_w)
    hi <- match(home, names(clades)) + 1L
    hi[is.na(hi)] <- 1L
    hp <- hi[g]
    within <- p_lab[cbind(hp, sb[left])] == p_lab[cbind(hp, sb[right])]
    pp <- ix$lut_p[hit]
    n_sig_w <- tabulate(g[sig & within], nbins = nb)
    n_sig_x <- tabulate(g[sig & !within], nbins = nb)
    up_x <- pres & !sig & !within & !is.na(pp) & pp < ix$min_power
    n_up_x <- tabulate(g[up_x], nbins = nb)
    n_l <- unname(lengths(clades)[home])
    lin_core <- !is.na(n_l) & m_sp == n_l & n_sig == choose(n_l, 2)
    # Covering the home with one gene per species, every child clade is
    # fully significant within exactly when n_sig_within reaches W_H.
    # n_present == n_pairs: an untested cross pair is absent evidence,
    # not evidence of divergence. W_H > 0 refuses a home of singletons,
    # where the tier would pass a clique with no significant pair.
    diff_ok <- m_sp == lengths(ix$parts)[hi] & n_present == n_pairs &
      p_w[hi] > 0 & n_sig_w == p_w[hi] & n_sig_x <= p_cmax[hi]
    # The call must survive treating every underpowered cross pair as
    # possibly significant.
    diff_up <- diff_ok & n_sig_x + n_up_x > p_cmax[hi]
  }

  n_gone <- n_sp - m_sp
  # Only "absent" and "untested" are annotation gaps. Admitting
  # "tested_ns" would let partial_present claim a clique whose missing
  # species was in fact rejected, and would let lineage_specific claim
  # one whose outside species were compared against every member and
  # diverged -- evidence of a boundary, not a gap, and the case
  # `differentiated` exists to score. The published workflow draws the
  # same line: its dicot- and conifer-specific sets require Cross == 0
  # in a matrix whose 1s mark a *tested* species pair, not a significant
  # one, while its differentiated set requires both lineages present.
  gap_only <- n_tns + n_up + n_ext == 0
  # An underpowered outside species was tested, so it is no gap, but its
  # failure is no boundary either: a lineage- or trait-specific call has
  # to survive reading it as conserved, which it cannot.
  up_only <- n_gone > 0L & n_up > 0 & n_ext == 0
  # At least one outside species was compared against every member and
  # rejected at adequate power, and no rejection is excused by power.
  ts_ok <- n_gone > 0L & n_tns > 0 & n_up == 0 & n_ext == 0
  # partial_present only asks that no missing species was *rejected*; an
  # underpowered one is unknown, not rejected, so it counts as a gap.
  gap_pp <- n_tns + n_ext == 0
  # choose(S - 1, 2) + 1 == choose(S, 2) - (S - 2): the tolerance is
  # S - 2 non-significant edges, the largest that cannot isolate a
  # member, since cutting one loose needs all S - 1 of its edges.
  is_part_sig <- m_sp == n_sp & n_present == choose(n_sp, 2) &
    !is.na(max_q) & max_q < ix$alpha_graph &
    n_sig >= choose(n_sp - 1L, 2) + 1
  is_part_pres <- n_gone >= 1L & n_gone <= ix$max_gap &
    n_sig == choose(m_sp, 2) & gap_pp
  # The waterfall, first match wins. One gene per species is the
  # invariant the pair arithmetic rests on; a within-species edge breaks
  # it, so such a clique is not scored. A one-member "clique" has no
  # pair, which every tier count would vacuously satisfy.
  valid <- m_sp == m & m_sp >= 2L
  rules <- list(
    complete_conserved = m_sp == n_sp & n_sig == choose(n_sp, 2),
    lineage_specific = lin_core & gap_only,
    underpowered = lin_core & up_only,
    partial_significant = is_part_sig,
    partial_present = is_part_pres,
    underpowered = diff_up,
    differentiated = diff_ok,
    trait_specific = lin_core & ts_ok
  )
  cls <- rep("unclassified", nb)
  for (r in rev(seq_along(rules))) {
    cls[which(valid & rules[[r]])] <- names(rules)[r]
  }
  list(
    classification = cls, clade = home, n_members = m, n_species = m_sp,
    n_pairs = n_pairs, n_present = n_present, n_sig = n_sig,
    n_sig_within = n_sig_w, n_sig_cross = n_sig_x,
    n_underpowered_cross = n_up_x, n_missing = n_gone,
    missing_species = miss_sp, missing_reason = miss_re,
    mean_q = .gcg_gmean(qp[pres], g[pres], nb), max_q = max_q,
    mean_effect_size = .gcg_gmean(ep[ok], g[ok], nb),
    is_core = valid & lin_core
  )
}


#' Members of the cliques to classify, sorted by clique
#'
#' Several runs get offset clique ids, so their cliques stay apart. A
#' clique takes the HOG of its first member row. Members sorted by a
#' positive integer clique id are not copied.
#'
#' @param x A `gene_cliques` result, a list of them, or a table of
#'   member rows.
#' @return A list: `id`, `sp` and `gene` of each member, sorted by `id`;
#'   `ids`, the sorted clique ids that `id` holds, and `hog`, their
#'   HOGs; `out_ids`, the clique ids to report, in order of first
#'   appearance. For a list of runs also `run` and `local`, the run and
#'   its own id of each offset id.
#' @noRd
.gcg_members <- function(x) {
  if (is.data.frame(x)) {
    absent <- setdiff(c("clique_id", "hog", "species", "gene"), names(x))
    if (length(absent) > 0L) {
      stop(
        "cliques missing required columns: ",
        paste(absent, collapse = ", ")
      )
    }
    ids <- unique(x$clique_id)
    return(.gcg_sort_members(list(
      id = x$clique_id, sp = x$species, gene = x$gene, ids = ids,
      hog = as.character(x$hog)[match(ids, x$clique_id)]
    )))
  }
  one <- inherits(x, "gene_cliques")
  runs <- if (one) list(x) else x
  ok <- is.list(runs) && length(runs) > 0L &&
    all(vapply(runs, inherits, logical(1), "gene_cliques"))
  if (!ok) {
    stop(
      "cliques must be a gene_clique_graph() result, a list of them, ",
      "or a table of member rows"
    )
  }
  top <- vapply(runs, function(r) max(c(0L, r$cliques$clique_id)), 1L)
  off <- cumsum(c(0L, top))
  mem <- lapply(seq_along(runs), function(i) {
    m <- runs[[i]]$members
    cl <- runs[[i]]$cliques$clique_id
    if (!all(.gcg_unique_ids(m$clique_id) %in% cl)) {
      m <- m[m$clique_id %in% cl, , drop = FALSE]
    }
    if (off[[i]] > 0L) m$clique_id <- m$clique_id + off[[i]]
    m
  })
  col <- function(f) {
    if (one) mem[[1L]][[f]] else unlist(lapply(mem, `[[`, f))
  }
  id <- col("clique_id")
  cl_id <- unlist(lapply(seq_along(runs), function(i) {
    runs[[i]]$cliques$clique_id + off[[i]]
  }))
  ids <- .gcg_unique_ids(id)
  cl_hog <- unlist(lapply(runs, function(r) as.character(r$cliques$hog)))
  .gcg_sort_members(list(
    id = id, sp = col("species"), gene = col("gene"), ids = ids,
    hog = cl_hog[match(ids, cl_id)],
    run = if (!one) rep.int(seq_along(runs), top),
    local = if (!one) sequence(top)
  ))
}


#' Whether clique ids are positive integers in ascending order
#' @noRd
.gcg_sorted_ids <- function(id) {
  is.integer(id) && length(id) > 0L && !anyNA(id) && id[1L] >= 1L &&
    isFALSE(is.unsorted(id))
}


#' Distinct clique ids in order of first appearance
#'
#' Sorted ids are counted rather than hashed: no copy per member.
#' @noRd
.gcg_unique_ids <- function(id) {
  if (.gcg_sorted_ids(id)) {
    which(tabulate(id, nbins = id[length(id)]) > 0L)
  } else {
    unique(id)
  }
}


#' Sort members by positive integer clique unless they already are
#' @noRd
.gcg_sort_members <- function(mem) {
  mem$out_ids <- mem$ids
  if (.gcg_sorted_ids(mem$id)) {
    return(mem)
  }
  k <- match(mem$id, mem$ids)
  o <- order(k)
  mem$id <- k[o]
  mem$sp <- mem$sp[o]
  mem$gene <- mem$gene[o]
  mem$ids <- seq_along(mem$ids)
  mem
}


#' Give each clique its run and its own id back
#' @noRd
.gcg_runs <- function(out, runs) {
  if (is.null(runs$run)) {
    return(out)
  }
  run <- runs$run[out$clique_id]
  out$clique_id <- runs$local[out$clique_id]
  cbind(run = run, out)
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
  fin <- is.finite(out$mean_q)
  floor_q <- if (!any(fin)) {
    NA_real_
  } else if (all(fin)) {
    min(out$mean_q)
  } else {
    min(out$mean_q[fin])
  }
  rm(fin)
  # .tol_min_ties()'s rule, without its count of distinct values, which
  # sorts every clique and is not needed.
  n_tied <- if (is.na(floor_q)) {
    0L
  } else {
    mn <- min(out$mean_q, na.rm = TRUE)
    sum(out$mean_q <= mn + .tie_tol() * abs(mn), na.rm = TRUE)
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
    classification = character(0), clade = character(0),
    n_members = integer(0),
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
