# Tier classification of gene-graph cliques. The builder is in
# clique_gene_graph.R.


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
#' @return A data frame with one row per clique variant:
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
#'     \item{variant, multiplicity}{Variant of the clique, and the
#'       paralog combinations it stands for (see below)}
#'     \item{capped}{`TRUE` when the clique budget of
#'       [gene_clique_graph()] dropped genes from the clique's HOG. For a
#'       list of runs, any run counts}
#'   }
#'   Twin genes share their partners in the graph, but their own tests
#'   can differ at `alpha_call`, `alpha_graph` or `min_power`. A class
#'   splits only where its copies differ in what that clique reads: the
#'   species it misses and its other members. A clique gives one row per
#'   choice of split, numbered by `variant`. Every combination a row
#'   stands for gets its tier, so tiers are exact. `multiplicity` counts
#'   those combinations; `mean_q`, `max_q` and `mean_effect_size`
#'   summarise them as in [gene_clique_graph()]. Most cliques have one variant.
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

  twin <- !is.null(mem$tc)
  ix <- .gcg_edge_index(
    edges, species, mem$hog, alpha_call, min_power, alpha_graph, twin
  )
  if (twin) {
    mem <- .gcg_variants(mem, ix)
    ix$copies <- mem$copies
    n_k <- length(mem$ids)
  }
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
      ix, mem$id[rr], mem$sp[rr], mem$gene[rr], mem$ids[cb], mem$hog[cb],
      mem$sc[rr]
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
  out$variant <- mem$variant %||% rep(1L, n_k)
  out$multiplicity <- mem$mult %||% rep(1L, n_k)
  out$capped <- out$hog %in% mem$capped
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
#' @param rows Whether to keep each row's end nodes and status code, and
#'   each pair's status code, for the twin signatures.
#' @return A list of the index, the relations and the pair lookup.
#' @noRd
.gcg_edge_index <- function(edges, species, hogs, alpha_call, min_power,
                            alpha_graph = Inf, rows = FALSE) {
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
  rel <- NULL
  if (rows) {
    # Status of a row: significant, failed at low power, under
    # alpha_graph. The tiers read nothing else of a row.
    up <- !sig & !is.na(pw) & pw < min_power
    below <- !is.na(q) & q < alpha_graph
    rel <- list(
      row_i = i, row_j = j, row_sig = sig,
      row_up = up, lut_code = (2L * below + up)[u]
    )
  }
  c(rel, list(
    species = species, hog_u = hog_u, gn_u = gn_u, key = key,
    nodes = nodes, node_s = node_s, n_nd = n_nd,
    lut_k = pk[u], lut_q = q[u], lut_e = ev[u], lut_p = pw[u],
    t_sp = by_sp %*% sym(rep(TRUE, n_e)), f_sp = by_sp %*% sym(!sig),
    u_sp = by_sp %*% sym(!sig & !is.na(pw) & pw < min_power),
    sg = sg, sg_deg = Matrix::colSums(sg), present = present,
    alpha_call = alpha_call,
    min_power = min_power
  ))
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
#' @param sc Twin sub-class of each member row (an index into
#'   `ix$copies`), `NA` for a single gene, or `NULL`.
#' @return A list of output columns, plus `is_core`: the clique is one
#'   whole clade, fully significant.
#' @noRd
.gcg_classify_block <- function(ix, id, sp, gene, ids, hog, sc = NULL) {
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
  # A twin pair stands for every gene pair of its two sub-classes: their
  # mean q and effect, their largest q. The status is one for all.
  ep <- ix$lut_e[hit]
  qm <- qx <- qp
  if (!is.null(sc)) {
    tw <- which(!is.na(sc[left]) | !is.na(sc[right]))
    if (length(tw) > 0L) {
      pm <- .gcg_pair_means(ix, a[tw], b[tw], sc[left][tw], sc[right][tw])
      qm[tw] <- pm$q
      qx[tw] <- pm$x
      ep[tw] <- pm$e
    }
  }
  max_q <- rep(NA_real_, nb)
  if (any(pres)) {
    gp <- structure(g[pres], N.groups = nb, class = "qG")
    max_q <- collapse::fmax(qx[pres], gp, na.rm = FALSE, use.g.names = FALSE)
    max_q[n_present == 0L] <- NA_real_
  }
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
    mean_q = .gcg_gmean(qm[pres], g[pres], nb), max_q = max_q,
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
#'   its own id of each offset id. `capped`, the HOGs that any run
#'   dropped genes from. With twin classes also `tc`, the
#'   class of each member (`NA` for a single gene), and `tgenes`,
#'   `thog`, `tsp`: the genes, HOG and species of each class.
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
  cp <- list()
  mem <- lapply(seq_along(runs), function(i) {
    m <- runs[[i]]$members
    cl <- runs[[i]]$cliques$clique_id
    if (!all(.gcg_unique_ids(m$clique_id) %in% cl)) {
      m <- m[m$clique_id %in% cl, , drop = FALSE]
    }
    cp[[i]] <<- .gcg_copies(runs[[i]], m)
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
  out <- list(
    id = id, sp = col("species"), gene = col("gene"), ids = ids,
    hog = cl_hog[match(ids, cl_id)],
    run = if (!one) rep.int(seq_along(runs), top),
    local = if (!one) sequence(top),
    capped = unlist(lapply(runs, function(r) as.character(r$dropped$hog)))
  )
  n_cls <- vapply(cp, function(z) length(z$genes), 1L)
  if (sum(n_cls) > 0L) {
    # Class ids of later runs are offset, so the runs stay apart.
    out$tc <- unlist(Map(`+`, lapply(cp, `[[`, "tc"), cumsum(c(0L, n_cls))[
      seq_along(cp)
    ]))
    for (f in c("genes", "hog", "sp")) {
      out[[paste0("t", f)]] <- unlist(lapply(cp, `[[`, f),
        recursive = FALSE, use.names = FALSE
      )
    }
  }
  .gcg_sort_members(out)
}


#' Twin classes of one run's members
#'
#' @param run A `gene_cliques` result.
#' @param m Its member rows to classify.
#' @return A list: `tc`, the class of each member row (`NA` for a single
#'   gene); `genes`, the genes of each class with the representative
#'   first; `hog` and `sp`, the HOG and species of each class.
#' @noRd
.gcg_copies <- function(run, m) {
  tc <- rep(NA_integer_, nrow(m))
  r <- which(m$n_copies > 1L)
  if (length(r) == 0L) {
    return(list(
      tc = tc, genes = list(), hog = character(0), sp = character(0)
    ))
  }
  cl <- run$cliques
  hog <- as.character(cl$hog)[match(m$clique_id[r], cl$clique_id)]
  sp <- as.character(m$species[r])
  key <- paste(hog, sp, m$gene[r], sep = .gcg_sep)
  u <- unique(key)
  tc[r] <- match(key, u)
  tw <- run$twins
  tk <- paste(tw$hog, as.character(tw$species), tw$representative,
    sep = .gcg_sep
  )
  by <- split(tw$gene, factor(tk, levels = u))
  first <- match(u, key)
  genes <- Map(
    function(g0, t) c(g0, sort(t, method = "radix")),
    m$gene[r][first], by
  )
  list(tc = tc, genes = unname(genes), hog = hog[first], sp = sp[first])
}


#' Split twin classes per clique and expand cliques into variants
#'
#' Copies of one class share their partners at the builder's
#' `alpha_graph`, but their own rows may differ. A clique reads only the
#' rows to the species it misses and the pairs with its other members.
#' In each clique, a class splits into sub-classes of copies that agree
#' on those items, and stays whole where its copies agree. A clique
#' becomes one variant per choice of sub-class in each member; a
#' sub-class's gene is its first copy in radix order. The split is per
#' class, so where classes interact, the product can give more variants
#' than distinct outcomes.
#'
#' @param mem Members from `.gcg_members()`, with twin classes.
#' @param ix Context from `.gcg_edge_index(rows = TRUE)`.
#' @return `mem` with one clique per variant, plus `sc` (sub-class of
#'   each member row, `NA` for one gene), `copies` (nodes of each
#'   sub-class), `variant` and `mult` (per variant); `out_ids` and `hog`
#'   are per variant.
#' @noRd
.gcg_variants <- function(mem, ix) {
  cn <- lengths(mem$tgenes)
  n_c <- length(cn)
  cls <- rep.int(seq_len(n_c), cn)
  cs <- cumsum(c(0L, cn))
  genes <- unlist(mem$tgenes, use.names = FALSE)
  nd <- match(ix$key(
    match(mem$thog, ix$hog_u)[cls], match(mem$tsp, ix$species)[cls],
    match(genes, ix$gn_u)
  ), ix$nodes)
  dis <- .gcg_discord(ix, cls, nd, cn)

  # A member is one unit: its node, or its whole twin class.
  n_sp <- length(ix$species)
  tc <- mem$tc
  k <- match(mem$id, mem$ids)
  sb <- if (is.factor(mem$sp)) {
    match(levels(mem$sp), ix$species)[as.integer(mem$sp)]
  } else {
    match(as.character(mem$sp), ix$species)
  }
  m_nd <- match(ix$key(
    match(mem$hog, ix$hog_u)[k], sb, match(as.character(mem$gene), ix$gn_u)
  ), ix$nodes)
  n_u <- ix$n_nd + n_c + 1
  u_key <- k * n_u + ifelse(is.na(tc), m_nd, ix$n_nd + tc)
  s_key <- k * (n_sp + 1) + sb

  # Discordant items of each twin row's class that its clique reads: a
  # species it misses, or another member.
  tw <- which(!is.na(tc))
  it_n <- tabulate(dis$cls, nbins = n_c)
  it_at <- cumsum(c(0L, it_n))
  r <- rep.int(tw, it_n[tc[tw]])
  i <- it_at[tc[r]] + sequence(it_n[tc[tw]])
  y <- dis$item[i]
  rel <- y < 0
  rel[rel] <- is.na(collapse::fmatch(k[r[rel]] * (n_sp + 1) - y[rel], s_key))
  g <- which(y > 0)
  rel[g] <- !is.na(collapse::fmatch(k[r[g]] * n_u + y[g], u_key))
  # A gene is also a member when it is a copy of another member class.
  o <- order(nd)
  len <- tabulate(nd, ix$n_nd)[y[g]]
  gi <- rep.int(g, len)
  cc <- cls[o][rep.int(match(y[g], nd[o]), len) + sequence(len) - 1L]
  hit <- cc != tc[r[gi]] &
    !is.na(collapse::fmatch(k[r[gi]] * n_u + ix$n_nd + cc, u_key))
  rel[gi[hit]] <- TRUE

  # Copies of a split row with equal values on its items form one part.
  ri <- r[rel]
  sz <- cn[tc[ri]]
  rr <- rep.int(ri, sz)
  pp <- sequence(sz)
  iv <- rep.int(i[rel], sz)
  o <- order(rr, pp, iv, method = "radix")
  rr <- rr[o]
  pp <- pp[o]
  gp <- collapse::group(rr, pp)
  sig <- vapply(split(dis$vals[dis$off[iv[o]] + pp], gp), paste, "",
    collapse = ","
  )
  f <- which(!duplicated(gp))
  rr <- rr[f]
  pp <- pp[f]
  part <- collapse::group(rr, sig)
  pos <- unname(split(pp, part))
  p_row <- rr[!duplicated(part)]
  pk <- paste(tc[p_row], vapply(pos, paste, "", collapse = ","))
  u <- !duplicated(pk)
  # Sub-classes: each whole class, then each distinct part.
  p_sub <- n_c + match(pk, pk[u])
  u_cls <- tc[p_row[u]]
  copies <- c(
    unname(split(nd, cls)), Map(function(c0, p) nd[cs[c0] + p], u_cls, pos[u])
  )
  sub_gene <- c(
    genes[cs[seq_len(n_c)] + 1L],
    genes[cs[u_cls] + vapply(pos[u], `[`, 1L, 1L)]
  )

  n_alt <- rep(1L, length(tc))
  n_alt[unique(p_row)] <- tabulate(match(p_row, unique(p_row)))
  at <- match(seq_along(tc), p_row) - 1L
  nv <- as.integer(collapse::fprod(as.numeric(n_alt),
    structure(k, N.groups = length(mem$ids), class = "qG"),
    use.g.names = FALSE
  ))
  off <- cumsum(c(0L, nv))
  one <- which(nv[k] == 1L)
  src <- list(one)
  sc <- list(tc[one])
  vid <- list(off[k[one]] + 1L)
  # Few cliques hold a split class; each lists its variants.
  for (rs in split(which(nv[k] > 1L), k[nv[k] > 1L])) {
    alt <- lapply(rs, function(x) {
      if (n_alt[x] == 1L) tc[x] else p_sub[at[x] + seq_len(n_alt[x])]
    })
    grid <- as.matrix(expand.grid(alt, KEEP.OUT.ATTRS = FALSE))
    src[[length(src) + 1L]] <- rep(rs, nrow(grid))
    sc[[length(sc) + 1L]] <- as.vector(t(grid))
    vid[[length(vid) + 1L]] <- off[k[rs[1L]]] +
      rep(seq_len(nrow(grid)), each = length(rs))
  }
  vid <- unlist(vid)
  o <- order(vid, method = "radix")
  src <- unlist(src)[o]
  sc <- unlist(sc)[o]
  vid <- vid[o]
  size <- lengths(copies)
  gene <- mem$gene[src]
  gene[!is.na(sc)] <- sub_gene[sc[!is.na(sc)]]
  mult <- collapse::fprod(ifelse(is.na(sc), 1, size[sc]),
    structure(vid, N.groups = off[length(off)], class = "qG"),
    use.g.names = FALSE
  )
  sc[!is.na(sc) & size[sc] == 1L] <- NA_integer_
  mem$id <- vid
  mem$sp <- mem$sp[src]
  mem$gene <- gene
  mem$sc <- sc
  mem$copies <- copies
  mem$ids <- seq_len(off[length(off)])
  mem$hog <- rep.int(mem$hog, nv)
  mem$out_ids <- rep.int(mem$out_ids, nv)
  mem$variant <- sequence(nv)
  mem$mult <- if (max(mult) <= .Machine$integer.max) as.integer(mult) else mult
  mem
}


#' Items on which the copies of each twin class differ
#'
#' A clique reads two things of a copy. Per species it misses: tested,
#' any failed row, any failed row with power, and the significant genes.
#' Per other member: the pair's presence, significance, q under
#' `alpha_graph`, and low power. An item on which every copy of a class
#' agrees never splits it, so only the other items are kept.
#'
#' @param ix Context from `.gcg_edge_index(rows = TRUE)`.
#' @param cls,nd Class and node of each copy, sorted by class.
#' @param cn Copies per class.
#' @return A list sorted by class: `cls`; `item`, minus a species index
#'   or a gene node; `off`, so that copy `p` of the class has the value
#'   `vals[off + p]` on the item.
#' @noRd
.gcg_discord <- function(ix, cls, nd, cn) {
  a <- c(ix$row_i, ix$row_j)
  hit <- which(a %in% nd)
  y <- c(ix$row_j, ix$row_i)[hit]
  t <- ix$node_s[y]
  hit <- hit[t <= length(ix$species)]
  x <- a[hit]
  y <- c(ix$row_j, ix$row_i)[hit]
  t <- ix$node_s[y]
  sig <- rep(ix$row_sig, 2L)[hit]
  nup <- !sig & !rep(ix$row_up, 2L)[hit]
  # Species value: flags negative, significant genes positive.
  sx <- c(x, x[!sig], x[nup], x[sig])
  st <- c(t, t[!sig], t[nup], t[sig])
  sv <- c(
    rep(-1, length(x)), rep(-2, sum(!sig)), rep(-4, sum(nup)), y[sig]
  )
  o <- order(sx, st, sv, method = "radix")
  key <- (sx * (length(ix$species) + 1) + st) * (ix$n_nd + 8) + sv
  o <- o[!duplicated(key[o])]
  gs <- collapse::group(sx[o], st[o])
  s_val <- vapply(split(sv[o], gs), paste, "", collapse = ",")
  f <- o[!duplicated(gs)]
  # Pair value: 1 + significant + 2 * (2 * below alpha_graph + low power).
  h <- match(pmin(x, y) * (ix$n_nd + 1) + pmax(x, y), ix$lut_k)
  ok <- which(!is.na(ix$lut_q[h]) & !duplicated(x * (ix$n_nd + 1) + y))
  nx <- c(sx[f], x[ok])
  item <- c(-st[f], y[ok])
  val <- c(
    match(s_val, unique(s_val)),
    1L + (ix$lut_q[h[ok]] < ix$alpha_call) + 2L * ix$lut_code[h[ok]]
  )

  # Each copy takes the items of its node.
  o <- order(nx)
  len <- tabulate(nx, ix$n_nd)[nd]
  len[is.na(len)] <- 0L
  j <- o[rep.int(match(nd, nx[o]), len) + sequence(len) - 1L]
  cc <- rep.int(cls, len)
  pp <- rep.int(sequence(cn), len)
  item <- item[j]
  val <- val[j]
  g <- collapse::group(cc, item)
  n_g <- attr(g, "N.groups")
  n_v <- collapse::fndistinct(val, g, use.g.names = FALSE)
  f <- which(!duplicated(g))
  g_cls <- cc[f]
  g_item <- item[f]
  bad <- which(tabulate(g, n_g) < cn[g_cls] | n_v > 1L)
  bad <- bad[order(g_cls[bad], g_item[bad])]
  off <- cumsum(c(0, cn[g_cls[bad]]))
  vals <- integer(off[length(off)])
  w <- match(g, bad)
  k <- which(!is.na(w))
  vals[off[w[k]] + pp[k]] <- val[k]
  list(
    cls = g_cls[bad], item = g_item[bad], off = off[seq_along(bad)],
    vals = vals
  )
}


#' Mean and largest q, and mean effect, over the gene pairs of two
#' sub-classes
#'
#' @param ix Context with `copies`.
#' @param a,b Representative nodes of each pair.
#' @param sa,sb Their sub-classes, `NA` for one gene.
#' @return A list of `q`, `x` (largest q) and `e`, one entry per pair.
#' @noRd
.gcg_pair_means <- function(ix, a, b, sa, sb) {
  ka <- ifelse(is.na(sa), -a, sa)
  kb <- ifelse(is.na(sb), -b, sb)
  pk <- paste(ka, kb)
  u <- which(!duplicated(pk))
  side <- function(n, k) {
    lapply(seq_along(u), function(i) {
      if (is.na(k[u[i]])) n[u[i]] else ix$copies[[k[u[i]]]]
    })
  }
  p <- side(a, sa)
  r <- side(b, sb)
  np <- lengths(p)
  nr <- lengths(r)
  n_u <- length(u)
  gid <- rep.int(seq_len(n_u), np * nr)
  pp <- unlist(Map(rep, p, nr), use.names = FALSE)
  rr <- unlist(Map(rep, r, each = np), use.names = FALSE)
  h <- match(pmin(pp, rr) * (ix$n_nd + 1) + pmax(pp, rr), ix$lut_k)
  q <- ix$lut_q[h]
  grp <- structure(gid, N.groups = n_u, class = "qG")
  st <- rbind(
    .gcg_gmean(q, gid, n_u),
    collapse::fmax(q, grp, na.rm = FALSE, use.g.names = FALSE),
    collapse::fmean(ix$lut_e[h], grp, use.g.names = FALSE)
  )
  j <- match(pk, pk[u])
  list(q = st[1L, j], x = st[2L, j], e = st[3L, j])
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
  if (!is.null(mem$tc)) mem$tc <- mem$tc[o]
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
  out$variant <- integer(0)
  out$multiplicity <- integer(0)
  out$capped <- logical(0)
  out
}
