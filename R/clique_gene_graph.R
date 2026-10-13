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
#' @return The zero-row tables with their full column sets.
#' @noRd
.gcg_empty <- function(has_effect, has_score, sp_lv) {
  cliques <- data.frame(
    clique_id = integer(0), hog = character(0),
    n_members = integer(0), n_species = integer(0),
    n_edges = integer(0), mean_q = numeric(0), max_q = numeric(0)
  )
  if (has_effect) cliques$mean_effect_size <- numeric(0)
  if (has_score) cliques$score <- numeric(0)
  cliques$multiplicity <- integer(0)
  members <- data.frame(
    clique_id = integer(0), species = factor(character(0), sp_lv),
    gene = character(0), n_copies = integer(0)
  )
  list(
    cliques = cliques, members = members, twins = .gcg_no_twins(sp_lv),
    dropped = .gcg_no_dropped(sp_lv)
  )
}


#' Empty twin table
#' @noRd
.gcg_no_twins <- function(sp_lv) {
  data.frame(
    hog = character(0), species = factor(character(0), sp_lv),
    gene = character(0), representative = character(0)
  )
}


#' Empty table of genes dropped by the clique budget
#' @noRd
.gcg_no_dropped <- function(sp_lv) {
  data.frame(
    hog = character(0), species = factor(character(0), sp_lv),
    gene = character(0), n_cliques_full = numeric(0),
    n_classes = integer(0), n_classes_kept = integer(0)
  )
}


#' Twin classes to keep in a hog over the clique budget
#'
#' Classes rank by evidence: the sum of their collapsed edge scores, or
#' their edge count without scores. Ties go to the representative gene
#' in radix order, then the species name. The top class of each species
#' always stays. A binary search finds the most classes, in rank order,
#' whose graph has at most `budget` maximal cliques.
#'
#' @param g The hog's twin-collapsed graph.
#' @param i1,i2,es Duplicate-collapsed edges before the twin collapse,
#'   as node indices, and their scores.
#' @param rep_v,nc Representative node and class size of each node.
#' @param sp,gene Species name and gene id of each node.
#' @param budget,min_size Call parameters.
#' @return A list: `keep`, a logical per node; `n_classes` and `n_kept`.
#' @noRd
.gcg_budget <- function(g, i1, i2, es, rep_v, nc, sp, gene, budget,
                        min_size) {
  n_v <- length(rep_v)
  # A collapsed edge carries the mean score of its gene pairs: each pair
  # adds score / (pairs). Sorted sums do not depend on row order.
  w <- if (all(is.na(es))) rep(1, length(i1)) else es
  w <- w / (nc[i1] * nc[i2])
  ev <- vapply(
    split(c(w, w), factor(rep_v[c(i1, i2)], seq_len(n_v))),
    function(x) sum(sort(x)), 0
  )
  cls <- which(rep_v == seq_len(n_v))
  cls <- cls[order(-ev[cls], gene[cls], sp[cls], method = "radix")]
  top <- !duplicated(sp[cls])
  cls <- c(cls[top], cls[!top])
  n_cls <- length(cls)
  fits <- function(m) {
    sub <- igraph::induced_subgraph(g, cls[seq_len(m)])
    igraph::count_max_cliques(sub, min = min_size) <= budget
  }
  lo <- sum(top)
  hi <- n_cls - 1L
  while (lo < hi) {
    mid <- (lo + hi + 1L) %/% 2L
    if (fits(mid)) lo <- mid else hi <- mid - 1L
  }
  list(
    keep = rep_v %in% cls[seq_len(lo)], n_classes = n_cls, n_kept = lo
  )
}


#' False-twin classes of one HOG's gene graph
#'
#' Two genes of one species are false twins when their neighbour sets
#' are equal. They are then never adjacent. The representative of a
#' class is its first gene in radix (C-locale) order.
#'
#' @param i1,i2 Edge end points, as node indices.
#' @param node_sp Species of each node.
#' @param gene Gene id of each node.
#' @return `NULL` when no class has two genes, else a list: `rep`, the
#'   representative node of each node, and `nc`, its class size.
#' @noRd
.gcg_twins <- function(i1, i2, node_sp, gene) {
  n_v <- length(node_sp)
  live <- tabulate(c(i1, i2), nbins = n_v) > 0L
  cnt <- tabulate(node_sp[live], nbins = max(node_sp))
  cand <- which(live & cnt[node_sp] > 1L)
  if (length(cand) == 0L) {
    return(NULL)
  }
  ends <- c(i1, i2)
  hit <- ends %in% cand
  nb <- split(c(i2, i1)[hit], factor(ends[hit], levels = cand))
  key <- paste(node_sp[cand], vapply(nb, function(x) {
    paste(sort(x), collapse = ",")
  }, ""))
  cls <- match(key, key)
  if (anyDuplicated(cls) == 0L) {
    return(NULL)
  }
  o <- order(cls, gene[cand], method = "radix")
  first <- o[!duplicated(cls[o])]
  rep_v <- seq_len(n_v)
  rep_v[cand] <- cand[first][match(cls, cls[first])]
  list(rep = rep_v, nc = tabulate(rep_v, nbins = n_v)[rep_v])
}


#' Attach the floor diagnostics to a [gene_clique_graph()] result
#'
#' Carried as columns as well as attributes: `subset()` drops attributes,
#' so a caller who filters the result would otherwise lose the record of
#' how little resolution `mean_q` has.
#'
#' @param res List of the `cliques`, `members` and `twins` tables.
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
#' @param xq Largest q-value behind each edge.
#' @param nc Twin class size of each node.
#' @return A list of per-clique statistics, the members' node ids and
#'   their class sizes.
#' @noRd
.gcg_clique_stats <- function(cl, em, eq, ee, es, sp_h, nodes, xq, nc) {
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
    xq = collapse::fmax(xq[e_h], grp, na.rm = FALSE, use.g.names = FALSE),
    me = .gcg_gmean(ee[e_h], g, n_k),
    sc = collapse::fsum(es[e_h], grp, na.rm = FALSE, use.g.names = FALSE),
    mu = collapse::fprod(as.numeric(nc[nd]),
      structure(cq, N.groups = n_k, class = "qG"),
      use.g.names = FALSE
    ),
    gi = nodes[nd], nc = nc[nd]
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
#' Two genes of one species are twins when they have the same partners
#' in the graph. Twins are never adjacent, and a clique holds at most one
#' of them. Swapping one twin for another gives another maximal clique,
#' so the twins of a class collapse to one gene before the cliques are
#' listed. Each clique stands for every choice of one gene per class.
#' The collapse is exact: expanding the cliques gives the full set. The
#' representative is the first gene of its class in radix (C-locale)
#' order.
#'
#' Paralogs multiply cliques, so a hog may produce millions of them. A
#' cap on genes per species does not bound that count. One hog with 76
#' twin classes, fewer than 10 per species, gave 2.4 million. The budget
#' counts cliques per hog, after the twin collapse. A hog over
#' `max_cliques` ranks its twin classes by the sum of their edge scores.
#' It keeps the top class of each species. It then keeps the most
#' classes, in rank order, that fit the budget. The genes of the other
#' classes go to `dropped`. Ties go to the gene id, so row order does not
#' matter.
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
#' @param max_cliques Largest number of maximal cliques one hog may
#'   produce. A hog above it keeps its wiring classes with the strongest
#'   edges until it fits. Default 1e5; `Inf` keeps every clique.
#'
#' @return A list of class `gene_cliques` with four tables.
#'   `cliques` has one row per clique:
#'   \describe{
#'     \item{clique_id}{Integer clique id, 1 to the number of cliques}
#'     \item{hog}{Ortholog group}
#'     \item{n_members}{Nodes in the clique}
#'     \item{n_species}{Distinct species in the clique (equal to
#'       `n_members` for cross-species edge tables)}
#'     \item{n_edges}{Graph edges inside the clique}
#'     \item{mean_q, max_q}{q-value summaries over those edges. Twin
#'       classes are joined by all their gene pairs. A clique edge
#'       carries the mean q of those pairs, so `mean_q` is the mean over
#'       the clique's expansions. `max_q` is the largest q of any
#'       expansion}
#'     \item{mean_effect_size}{Present only when `edges` carries
#'       `effect_size`. Prefer it over `mean_q` for ranking: q-values
#'       saturate at the permutation floor and cannot separate cliques
#'       once they get there}
#'     \item{score}{Present only when `edges` carries `score`. The sum
#'       of the clique edges' `score`, in bits. With twins, the mean
#'       over the clique's expansions, as for `mean_effect_size`}
#'     \item{multiplicity}{Cliques this row stands for: the product of
#'       the members' `n_copies`. 1 without twins}
#'     \item{mean_q_floor, n_cliques_at_q_floor}{Smallest clique
#'       `mean_q` in the run and how many cliques are tied there. A
#'       large tie count means `mean_q` cannot rank those cliques at
#'       all}
#'   }
#'   `members` has one row per clique member: `clique_id`, `species`
#'   (a factor; its levels are the species of `edges` in order of first
#'   appearance), `gene` (the representative of its twin class) and
#'   `n_copies` (the genes with the same called partners, itself
#'   included).
#'   `twins` has one row per collapsed gene: `hog`, `species`, `gene`
#'   and `representative`. It has no rows when no two genes are twins.
#'   `dropped` has one row per gene that the clique budget dropped:
#'   `hog`, `species`, `gene`, `n_cliques_full` (the hog's maximal
#'   cliques before the drop), `n_classes` (its twin classes) and
#'   `n_classes_kept`. It has no rows when no hog is over the budget.
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
                                      alpha_graph = 0.1, max_cliques = 1e5,
                                      ...) {
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
  ok_max <- is.numeric(max_cliques) && length(max_cliques) == 1L &&
    !is.na(max_cliques) && max_cliques >= 1
  if (!ok_max) {
    stop(
      "max_cliques must be a single number of at least 1; you gave ",
      deparse1(max_cliques)
    )
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
  twins <- list()
  dropped <- list()

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

    # Twin classes: a clique holds at most one gene of a class, and two
    # classes are joined by all their gene pairs or by none. One
    # representative per class gives every clique once; its gene pairs
    # become one edge with mean q, effect and score, and max q.
    xq <- eq
    nc <- rep(1L, n_v)
    n_q <- n_v
    j1 <- i1
    j2 <- i2
    s0 <- es
    tw <- .gcg_twins(i1, i2, node_sp, all_gn[nodes])
    if (!is.null(tw)) {
      r1 <- tw$rep[i1]
      r2 <- tw$rep[i2]
      pk <- (pmin(r1, r2) - 1) * n_v + pmax(r1, r2)
      u <- which(!duplicated(pk))
      gq <- match(pk, pk[u])
      n_u <- length(u)
      xq <- collapse::fmax(eq, structure(gq, N.groups = n_u, class = "qG"),
        na.rm = FALSE, use.g.names = FALSE
      )
      eq <- .gcg_gmean(eq, gq, n_u)
      ee <- .gcg_gmean(ee, gq, n_u)
      es <- .gcg_gmean(es, gq, n_u)
      i1 <- r1[u]
      i2 <- r2[u]
      nc <- tw$nc
      v <- which(tw$rep != seq_len(n_v))
      n_q <- n_q - length(v)
      twins[[length(twins) + 1L]] <- list(
        hog = h, sp = node_sp[v], gene = all_gn[nodes[v]],
        rep = all_gn[nodes[tw$rep[v]]]
      )
    }

    g <- igraph::make_graph(as.vector(rbind(i1, i2)),
      n = n_v, directed = FALSE
    )
    # Clique budget: a hog over it keeps its strongest twin classes. The
    # genes of a dropped class stay as isolated vertices.
    n_full <- if (is.finite(max_cliques)) {
      igraph::count_max_cliques(g, min = min_size)
    } else {
      0
    }
    if (n_full > max_cliques) {
      rep_v <- if (is.null(tw)) seq_len(n_v) else tw$rep
      bd <- .gcg_budget(g, j1, j2, s0, rep_v, nc, sp_lv[node_sp],
        all_gn[nodes], max_cliques, min_size
      )
      v <- which(!bd$keep)
      v <- v[order(sp_lv[node_sp[v]], all_gn[nodes[v]], method = "radix")]
      dropped[[length(dropped) + 1L]] <- list(
        hog = h, sp = node_sp[v], gene = all_gn[nodes[v]], n_full = n_full,
        n_cls = bd$n_classes, n_kept = bd$n_kept
      )
      ke <- bd$keep[i1] & bd$keep[i2]
      i1 <- i1[ke]
      i2 <- i2[ke]
      eq <- eq[ke]
      ee <- ee[ke]
      es <- es[ke]
      xq <- xq[ke]
      g <- igraph::make_graph(as.vector(rbind(i1, i2)),
        n = n_v, directed = FALSE
      )
    }
    em <- matrix(0L, n_v, n_v)
    em[cbind(i1, i2)] <- em[cbind(i2, i1)] <- seq_along(i1)
    sp_h <- match(node_sp, unique(node_sp))
    # A few hogs hold most cliques. One call per start vertex keeps one
    # piece in memory at a time; the pieces are disjoint and together
    # complete only when every vertex id starts one.
    starts <- if (n_q >= 40L) seq_len(n_v) else list(NULL)
    out_h <- list()
    for (v in starts) {
      cl <- igraph::max_cliques(g, min = min_size, subset = v)
      for (b in split(seq_along(cl), ceiling(seq_along(cl) / 2e5))) {
        st <- .gcg_clique_stats(cl[b], em, eq, ee, es, sp_h, nodes, xq, nc)
        st$hog <- h
        out_h[[length(out_h) + 1L]] <- st
      }
    }
    blk[[k]] <- out_h
  }

  dr_tab <- .gcg_no_dropped(sp_lv)
  if (length(dropped) > 0L) {
    n_g <- lengths(lapply(dropped, `[[`, "gene"))
    per <- function(f) rep.int(vapply(dropped, `[[`, 1, f), n_g)
    dr_tab <- data.frame(
      hog = rep.int(vapply(dropped, `[[`, "", "hog"), n_g),
      species = structure(
        unlist(lapply(dropped, `[[`, "sp"), use.names = FALSE),
        levels = sp_lv, class = "factor"
      ),
      gene = unlist(lapply(dropped, `[[`, "gene"), use.names = FALSE),
      n_cliques_full = per("n_full"),
      n_classes = as.integer(per("n_cls")),
      n_classes_kept = as.integer(per("n_kept"))
    )
    message(
      "gene_clique_graph: ", length(dropped), " hogs exceeded ",
      "`max_cliques` = ", format(max_cliques, big.mark = ",",
        scientific = FALSE
      ), " maximal cliques; ", nrow(dr_tab),
      " genes dropped, see `$dropped`."
    )
    rm(dropped)
  }
  tw_tab <- .gcg_no_twins(sp_lv)
  if (length(twins) > 0L) {
    tw_sp <- unlist(lapply(twins, `[[`, "sp"), use.names = FALSE)
    tw_tab <- data.frame(
      hog = rep.int(
        vapply(twins, `[[`, "", "hog"), lengths(lapply(twins, `[[`, "sp"))
      ),
      species = structure(tw_sp, levels = sp_lv, class = "factor"),
      gene = unlist(lapply(twins, `[[`, "gene"), use.names = FALSE),
      representative = unlist(lapply(twins, `[[`, "rep"), use.names = FALSE)
    )
    rm(twins, tw_sp)
    cls_key <- paste(tw_tab$hog, as.integer(tw_tab$species),
      tw_tab$representative,
      sep = .gcg_sep
    )
    n_cls <- length(unique(cls_key))
    message(
      "gene_clique_graph: ", nrow(tw_tab) + n_cls, " genes collapsed ",
      "into ", n_cls, " twin classes; cliques carry multiplicity"
    )
  }

  blk <- unlist(blk, recursive = FALSE, use.names = FALSE)
  if (length(blk) == 0L) {
    out <- .gcg_empty(has_effect, has_score, sp_lv)
    out$twins <- tw_tab
    out$dropped <- dr_tab
    return(.gcg_graph_attrs(out, alpha_graph, min_size, min(qv), NA_real_, 0L))
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
    clique_id = rep.int(seq_len(n_k), sz), species = species, gene = gene,
    n_copies = take("nc")
  )
  rm(species, gene)
  col <- list()
  for (f in c("ns", "ne", "mq", "xq", "me", "sc", "mu")) {
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
  # Integer while it fits: ten copies in each of ten species overflow.
  fits <- max(col$mu) <= .Machine$integer.max
  cliques$multiplicity <- if (fits) as.integer(col$mu) else col$mu
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
    list(
      cliques = cliques, members = members, twins = tw_tab, dropped = dr_tab
    ),
    alpha_graph, min_size, min(qv), mq_min, n_tied
  )
}


#' Every clique of a twin-reduced result, one row per expansion
#'
#' A clique of [gene_clique_graph()] stands for every choice of one gene
#' per twin class of its members. This lists those choices and
#' recomputes their statistics from `edges`, as an unreduced run would.
#' The tests use it; nothing else does.
#'
#' @param g A [gene_clique_graph()] result, members sorted by clique.
#' @param edges The edge table it was built from.
#' @return A `gene_cliques` list: `cliques`, one row per expansion with
#'   `clique_id` 1 to their number and `quotient_id`, the clique of `g`
#'   it expands; and `members`.
#' @noRd
.gcg_expand <- function(g, edges) {
  cl <- g$cliques
  m <- g$members
  tw <- g$twins
  nc <- m$n_copies
  k <- match(m$clique_id, cl$clique_id)
  stopifnot(!is.unsorted(k))
  hog_m <- as.character(cl$hog)[k]
  sp_m <- as.character(m$species)
  cp <- as.list(m$gene)
  r <- which(nc > 1L)
  if (length(r) > 0L) {
    by <- split(tw$gene, paste(tw$hog, as.character(tw$species),
      tw$representative,
      sep = .gcg_sep
    ))
    key <- paste(hog_m[r], sp_m[r], m$gene[r], sep = .gcg_sep)
    cp[r] <- Map(function(g0, t) c(g0, sort(t, method = "radix")),
      m$gene[r], by[key]
    )
  }
  flat <- unlist(cp, use.names = FALSE)
  at <- cumsum(c(0L, lengths(cp)))[seq_along(cp)]
  fkey <- paste(rep.int(hog_m, lengths(cp)), rep.int(sp_m, lengths(cp)),
    flat,
    sep = .gcg_sep
  )

  # Mixed radix: expansion t of a clique takes copy
  # (t %/% stride) %% n_copies of each member.
  n_k <- nrow(cl)
  sz <- tabulate(k, nbins = n_k)
  f <- match(seq_len(n_k), k)
  stride <- stats::ave(as.numeric(nc), k, FUN = function(x) {
    cumprod(c(1, x))[seq_along(x)]
  })
  mu <- as.numeric(cl$multiplicity)
  ek <- rep.int(seq_len(n_k), mu)
  tt <- sequence(mu) - 1
  er <- rep.int(seq_along(ek), sz[ek])
  row <- f[ek][er] + sequence(sz[ek]) - 1L
  pos <- at[row] + (tt[er] %/% stride[row]) %% nc[row] + 1

  ag <- attr(cl, "alpha_graph")
  e <- edges[!is.na(edges$q_value) & edges$q_value < ag, , drop = FALSE]
  k1 <- paste(e$hog, e$species1, e$gene1, sep = .gcg_sep)
  k2 <- paste(e$hog, e$species2, e$gene2, sep = .gcg_sep)
  nodes <- unique(c(k1, k2))
  i <- match(k1, nodes)
  j <- match(k2, nodes)
  n_nd <- length(nodes)
  pk <- pmin(i, j) * (n_nd + 1) + pmax(i, j)
  na <- rep(NA_real_, nrow(e))
  ev <- if (is.null(e$effect_size)) na else as.numeric(e$effect_size)
  sv <- if (is.null(e$score)) na else as.numeric(e$score)
  u <- order(pk, e$q_value, -ev, na.last = TRUE)
  u <- u[!duplicated(pk[u])]

  nd <- match(fkey, nodes)[pos]
  n_e <- length(ek)
  ms <- sz[ek]
  n_right <- ms[er] - sequence(ms)
  left <- rep.int(seq_along(nd), n_right)
  right <- left + sequence(n_right)
  a <- nd[left]
  b <- nd[right]
  hit <- match(pmin(a, b) * (n_nd + 1) + pmax(a, b), pk[u])
  gp <- er[left]
  grp <- structure(gp, N.groups = n_e, class = "qG")
  q <- as.numeric(e$q_value)[u][hit]
  out <- data.frame(
    clique_id = seq_len(n_e), quotient_id = cl$clique_id[ek],
    hog = cl$hog[ek], n_members = ms, n_species = cl$n_species[ek],
    n_edges = tabulate(gp, nbins = n_e), mean_q = .gcg_gmean(q, gp, n_e),
    max_q = collapse::fmax(q, grp, na.rm = FALSE, use.g.names = FALSE)
  )
  if (!is.null(cl$mean_effect_size)) {
    out$mean_effect_size <- .gcg_gmean(ev[u][hit], gp, n_e)
  }
  if (!is.null(cl$score)) {
    out$score <- collapse::fsum(sv[u][hit], grp,
      na.rm = FALSE,
      use.g.names = FALSE
    )
  }
  attr(out, "alpha_graph") <- ag
  attr(out, "min_size") <- attr(cl, "min_size")
  members <- data.frame(
    clique_id = er, species = m$species[row], gene = flat[pos]
  )
  structure(list(cliques = out, members = members), class = "gene_cliques")
}
