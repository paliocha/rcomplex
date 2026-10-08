# was rcomplex::recurrence_graph print.recurrence_graph recurrence_modules until 0.4.0; see git history (main @ b245ac0)

#' Cross-species recurrence graph of ortholog-group pairs
#'
#' @description
#' Contracts each species' co-expression network to ortholog groups, counts
#' in how many species every pair of groups is co-expressed, and tests that
#' count against independence of the species. Significant pairs form a
#' weighted graph over ortholog groups for [recurrence_modules()].
#'
#' @details
#' **Contraction.** In species \eqn{s} the binary adjacency \eqn{A_s} at
#' `net$threshold` (membership only; edge weights are not read) is
#' contracted to ortholog groups by \eqn{M^T A_s M}, with \eqn{M} the
#' gene-by-group indicator: the pair \eqn{(h_1, h_2)}, \eqn{h_1 \ne h_2},
#' is present in \eqn{s} when any copy of \eqn{h_1} is a neighbour of any
#' copy of \eqn{h_2}. This is MULE's ortholog contraction (Koyuturk et
#' al. 2004). Within-group pairs are not scored.
#'
#' **Null.** With expression shuffled within each species, the species
#' are independent and an edge of species \eqn{s} falls on a gene pair
#' with probability \eqn{d_s}, the realised edge density among the
#' group-mapped genes. A pair whose groups have \eqn{c_1} and \eqn{c_2}
#' copies in \eqn{s} is then present with probability
#' \eqn{p_s = 1 - (1 - d_s)^{c_1 c_2}}, so paralog copy number enters the
#' probability, not a weight. The number of species \eqn{K} in which the
#' pair is present is Poisson-binomial over the \eqn{p_s};
#' `p.val` is the exact upper tail \eqn{P(K \ge k_{obs})} by iterative
#' convolution. Only pairs with \eqn{K \ge} `min_species` are listed,
#' but `q.val` is Benjamini-Hochberg over every pair whose groups share
#' at least `min_species` species (`params$n_tests`), the unlisted ones
#' counted at p = 1. Correcting over the listed pairs alone selects on
#' the outcome: on three shuffled 200-gene networks at d = 0.03 every
#' listed pair has p near \eqn{3 d^2}, and that correction called all 73
#' of them at q < 0.05; the full count calls none. This is CODENSE's
#' summary graph (Hu et al. 2005) with an analytic null. Within-species
#' triangles make edges of one species dependent but leave each edge's
#' marginal alone; on shuffled Pooideae leaf, root and EVOTREE wood
#' networks (six to eight species, BH over the listed pairs) the null
#' gave 0 pairs at q < 0.05 and matched the K distribution within 3 % at
#' every K
#' (`dev/design-notes/module-engine-borrow-scope.md`, Section 10.3).
#'
#' @param nets Named list of at least two network objects from
#'   [compute_network()], one per species.
#' @param hog_map Data frame with columns `species`, `gene` and `hog`, one
#'   row per gene copy (many-to-many). Genes outside a species' network
#'   are ignored.
#' @param alpha FDR level for the edges of `graph` (0.1, the package's
#'   call threshold; the benchmark in the design note used 0.05).
#' @param min_species Smallest \eqn{K} for a pair to be tested.
#'
#' @return An object of class `recurrence_graph`: `edges` (data frame of
#'   the tested pairs: `hog1`, `hog2`, `K`, `profile` = presence per
#'   species in `names(nets)` order as a string such as `"101"`, `p.val`,
#'   `q.val`), `species` (`species`, `density`, `n_hogs`, `n_genes`, over
#'   group-mapped genes in the network), `copies` (group x species integer
#'   matrix of copies in the network), `graph` (undirected igraph of the
#'   pairs with `q.val < alpha`, edge `weight = -log10(p.val)`), `genes`
#'   (named list, the network genes of each species) and `params`
#'   (with `n_tests`, the number of testable pairs BH corrects over).
#'
#' @references
#' Hu, H., Yan, X., Huang, Y., Han, J. & Zhou, X. J. (2005). Mining
#' coherent dense subgraphs across massive biological networks for
#' functional discovery. \emph{Bioinformatics}, 21(Suppl 1), i213--i221.
#' \doi{10.1093/bioinformatics/bti1049}
#'
#' Koyuturk, M., Grama, A. & Szpankowski, W. (2004). An efficient
#' algorithm for detecting frequent subgraphs in biological networks.
#' \emph{Bioinformatics}, 20(Suppl 1), i200--i207.
#' \doi{10.1093/bioinformatics/bth919}
#'
#' @seealso [recurrence_modules()]
#' @export
recurrence_graph <- function(nets, hog_map, alpha = 0.1,
                             min_species = 2L) {
  sp <- names(nets)
  if (!is.list(nets) || length(nets) < 2L || is.null(sp) ||
        any(!nzchar(sp)) || anyDuplicated(sp)) {
    stop("nets must be a list of at least two networks with unique names")
  }
  if (length(sp) > 30L) stop("at most 30 species are supported")
  if (!is.data.frame(hog_map) ||
        !all(c("species", "gene", "hog") %in% names(hog_map))) {
    stop("hog_map must be a data frame with columns species, gene, hog")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha > 1) {
    stop("alpha must be a single number in (0, 1]")
  }
  if (!is.numeric(min_species) || length(min_species) != 1L ||
        is.na(min_species) || min_species < 1 || min_species > length(sp)) {
    stop("min_species must be a single number in 1..length(nets)")
  }
  hog_map <- unique(data.frame(
    species = as.character(hog_map$species),
    gene = as.character(hog_map$gene),
    hog = as.character(hog_map$hog), stringsAsFactors = FALSE
  ))
  genes <- lapply(nets, function(n) rownames(n$network))
  hog_map <- hog_map[hog_map$species %in% sp, , drop = FALSE]
  in_net <- logical(nrow(hog_map))
  for (s in sp) {
    r <- hog_map$species == s
    in_net[r] <- hog_map$gene[r] %in% genes[[s]]
  }
  hog_map <- hog_map[in_net & !is.na(hog_map$hog), , drop = FALSE]
  hogs <- sort(unique(hog_map$hog))
  nh <- length(hogs)
  if (nh < 2L) stop("fewer than two ortholog groups map to the networks")

  con <- lapply(sp, function(s) {
    .recurrence_contract(nets[[s]], hog_map[hog_map$species == s, ], hogs)
  })
  copies <- vapply(con, `[[`, integer(nh), "copies")
  dimnames(copies) <- list(hogs, sp)
  d <- vapply(con, `[[`, 0, "density")

  # species presence as bits: pair present in s adds 2^(s - 1)
  pm <- Reduce(`+`, lapply(seq_along(sp), function(s) {
    Matrix::sparseMatrix(
      i = con[[s]]$i, j = con[[s]]$j, x = 2^(s - 1), dims = c(nh, nh)
    )
  }))
  pm <- methods::as(pm, "TsparseMatrix")
  prof <- as.integer(pm@x)
  bits <- vapply(seq_along(sp) - 1L, function(b) {
    bitwAnd(prof, bitwShiftL(1L, b)) > 0L
  }, logical(length(prof)))
  bits <- matrix(bits, ncol = length(sp))
  k_obs <- rowSums(bits)
  ok <- k_obs >= min_species
  h1 <- pm@i[ok] + 1L
  h2 <- pm@j[ok] + 1L
  bits <- bits[ok, , drop = FALSE]
  k_obs <- k_obs[ok]

  # ponytail: one pairs x (S + 1) matrix in memory (about 2 GB at 3e7
  # pairs and 8 species); chunk over pairs if that binds
  v <- copies[h1, , drop = FALSE] * copies[h2, , drop = FALSE]
  p_s <- v
  for (s in seq_along(sp)) p_s[, s] <- 1 - (1 - d[s])^v[, s]
  p_val <- .poibin_upper(p_s, k_obs)
  n_tests <- .n_testable_pairs(copies, min_species)
  q_val <- stats::p.adjust(p_val, "BH", n = max(n_tests, length(p_val)))

  edges <- data.frame(
    hog1 = hogs[h1], hog2 = hogs[h2], K = as.integer(k_obs),
    profile = do.call(paste0, lapply(seq_along(sp), function(s) {
      ifelse(bits[, s], "1", "0")
    })),
    p.val = p_val, q.val = q_val, stringsAsFactors = FALSE
  )
  sig <- edges[edges$q.val < alpha, , drop = FALSE]
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = sig$hog1, to = sig$hog2,
      weight = -log10(pmax(sig$p.val, .Machine$double.xmin))
    ),
    directed = FALSE
  )
  structure(list(
    edges = edges,
    species = data.frame(
      species = sp, density = d,
      n_hogs = colSums(copies > 0L),
      n_genes = vapply(con, `[[`, 0L, "n_genes"),
      row.names = NULL, stringsAsFactors = FALSE
    ),
    copies = copies, graph = graph, genes = genes,
    params = list(
      alpha = alpha, min_species = min_species, n_tests = n_tests
    )
  ), class = "recurrence_graph")
}


#' @export
print.recurrence_graph <- function(x, ...) {
  e <- x$edges
  cat(
    "Recurrence graph over", nrow(x$species), "species:",
    paste(x$species$species, collapse = ", "), "\n"
  )
  cat("  ortholog groups:", nrow(x$copies), "\n")
  cat("  pairs testable:", x$params$n_tests, "\n")
  cat("  pairs listed (K >=", x$params$min_species, "):", nrow(e), "\n")
  cat(
    "  significant (q <", x$params$alpha, "):",
    igraph::ecount(x$graph), "\n"
  )
  tab <- table(factor(e$K, seq_len(nrow(x$species))))
  tab <- tab[as.integer(names(tab)) >= x$params$min_species]
  cat(
    "  K distribution:",
    paste0("K=", names(tab), ": ", tab, collapse = ", "), "\n"
  )
  invisible(x)
}


#' Contract one species' network to ortholog-group pairs
#'
#' Returns the upper-triangle group pairs present (1-based `i < j` into
#' `hogs`), copies per group, edge density among mapped genes and the
#' number of mapped genes.
#' @noRd
.recurrence_contract <- function(net, map, hogs) {
  net <- .net_as_sparse(net)
  a <- .net_cpp_args(net, net$threshold)
  g <- rownames(net$network)
  n <- length(g)
  col <- rep.int(seq_len(n), diff(a$p))
  keep <- a$x >= a$thr
  adj <- Matrix::sparseMatrix(
    i = a$i[keep] + 1L, j = col[keep], x = 1, dims = c(n, n)
  )
  gi <- match(map$gene, g)
  hi <- match(map$hog, hogs)
  mm <- Matrix::sparseMatrix(
    i = gi, j = hi, x = 1, dims = c(n, length(hogs))
  )
  mapped <- logical(n)
  mapped[gi] <- TRUE
  nm <- sum(mapped)
  n_edges <- sum(mapped[a$i[keep] + 1L] & mapped[col[keep]]) / 2
  pr <- methods::as(Matrix::crossprod(mm, adj %*% mm), "TsparseMatrix")
  up <- pr@i < pr@j & pr@x > 0
  list(
    i = pr@i[up] + 1L, j = pr@j[up] + 1L,
    copies = tabulate(hi, length(hogs)),
    density = if (nm > 1L) n_edges / choose(nm, 2) else 0,
    n_genes = nm
  )
}


#' Number of group pairs that share at least k species
#'
#' Counted over the distinct species-presence patterns of the groups.
#' @noRd
.n_testable_pairs <- function(copies, k) {
  pres <- copies > 0L
  key <- as.vector(pres %*% 2^(seq_len(ncol(pres)) - 1L))
  first <- !duplicated(key)
  n <- tabulate(match(key, key[first]))
  # ponytail: patterns x patterns matrices; fine for the <= 2^S patterns
  # that S <= 10 species allow, quadratic in patterns beyond that
  shared <- tcrossprod(pres[first, , drop = FALSE] * 1)
  both <- outer(n, n)
  diag(both) <- choose(n, 2)
  sum(both[upper.tri(both, diag = TRUE) & shared >= k])
}


#' Poisson-binomial upper tail P(K >= k) per row of a probability matrix
#'
#' Exact by iterative convolution over the columns.
#' @noRd
.poibin_upper <- function(p, k) {
  s <- ncol(p)
  f <- cbind(1, matrix(0, nrow(p), s))
  for (j in seq_len(s)) {
    q <- p[, j]
    f <- cbind(
      f[, 1L] * (1 - q),
      f[, -1L, drop = FALSE] * (1 - q) + f[, -(s + 1L), drop = FALSE] * q
    )
  }
  for (j in s:1L) f[, j] <- f[, j] + f[, j + 1L]
  f[cbind(seq_len(nrow(f)), k + 1L)]
}


#' Modules of ortholog groups from a recurrence graph
#'
#' @description
#' Finds modules of ortholog groups in the significant-pair graph of
#' [recurrence_graph()] and expands each to the gene copies of every
#' species, ready for [module_auroc()].
#'
#' @details
#' `"leiden"` partitions `rg$graph` with `igraph::cluster_leiden()` under
#' the modularity objective, edges weighted by \eqn{-\log_{10} p}, run to
#' convergence. Coarse resolutions (0.5) gave 3-4 giant modules that
#' replicated between sample halves on n = 20 Pooideae networks; fine
#' partitions did not (Section 10.3 of the design note).
#'
#' `"anchored"` returns, for each anchor group, the subgraph of maximum
#' weighted density \eqn{w(E(S)) / |S|} among the sets \eqn{S} that
#' contain the anchor, searched within the anchor's one-hop neighbourhood
#' in `rg$graph`. It is exact by Goldberg's (1984) parametric minimum cut:
#' source to every node with capacity \eqn{m = w(E)}, node \eqn{v} to sink
#' with \eqn{m + 2g - d_v}, edges in both directions at their weight, and
#' the anchor tied to the source; some \eqn{S} has density above \eqn{g}
#' exactly when the cut is below \eqn{m |V|}, and \eqn{g} is found by
#' bisection. There is no resolution parameter. Anchored modules of nearby
#' anchors can share groups.
#'
#' **Species elements.** Each module's gene set in a species is every
#' copy (from `hog_map`) of its groups present in that species' network.
#' When the sets of a species are disjoint the element is an
#' [as_modules()] object. When they overlap (nearby anchors, or a gene in
#' two groups) the element has the same fields, `module_genes` holding
#' every set whole and `modules` assigning each shared gene to its first
#' module only; [module_auroc()] reads `module_genes`, while the
#' preservation tests need a partition and should not be given it.
#'
#' @param rg A [recurrence_graph()] result.
#' @param hog_map The `hog_map` given to [recurrence_graph()].
#' @param method `"leiden"` or `"anchored"`.
#' @param resolution Leiden resolution (modularity).
#' @param anchors Character vector of anchor ortholog groups (required for
#'   `"anchored"`); anchors missing from `rg$graph` are dropped with a
#'   message.
#' @param min_size Smallest module, in ortholog groups.
#' @param seed Seed for Leiden, or `NULL` for the ambient stream (package
#'   RNG contract). `"anchored"` draws nothing.
#'
#' @return A list: `modules` (named list over species, each an
#'   [as_modules()]-shaped object for that species' genes), `hog_modules`
#'   (named list, module -> ortholog groups; Leiden modules are labelled
#'   `M1`, `M2`, ... by decreasing size, anchored modules by their anchor),
#'   `method` and `params`.
#'
#' @references
#' Goldberg, A. V. (1984). Finding a maximum density subgraph. Technical
#' report UCB/CSD-84-171, University of California, Berkeley.
#'
#' @seealso [recurrence_graph()], [module_auroc()]
#' @export
recurrence_modules <- function(rg, hog_map,
                               method = c("leiden", "anchored"),
                               resolution = 0.5, anchors = NULL,
                               min_size = 10L, seed = NULL) {
  if (!inherits(rg, "recurrence_graph")) {
    stop("rg must be a recurrence_graph() result")
  }
  method <- match.arg(method)
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
        min_size < 1) {
    stop("min_size must be a single number >= 1")
  }
  .seed_scope(seed)
  g <- rg$graph
  if (method == "leiden") {
    if (igraph::ecount(g) == 0L) stop("rg$graph has no significant pairs")
    memb <- igraph::membership(igraph::cluster_leiden(
      g,
      objective_function = "modularity", resolution = resolution,
      n_iterations = -1L
    ))
    hog_modules <- split(igraph::V(g)$name, memb)
    hog_modules <- hog_modules[order(-lengths(hog_modules))]
    hog_modules <- hog_modules[lengths(hog_modules) >= min_size]
    names(hog_modules) <- paste0("M", seq_along(hog_modules))
  } else {
    if (!is.character(anchors) || !length(anchors)) {
      stop("method = \"anchored\" needs anchors (ortholog group names)")
    }
    anchors <- unique(anchors)
    found <- anchors %in% igraph::V(g)$name
    if (!all(found)) {
      message(sum(!found), " anchor(s) not in rg$graph dropped")
    }
    hog_modules <- lapply(anchors[found], function(a) {
      ego <- igraph::make_ego_graph(g, 1L, a)[[1L]]
      .anchored_densest(ego, a)
    })
    names(hog_modules) <- anchors[found]
    small <- lengths(hog_modules) < min_size
    if (any(small)) {
      message(sum(small), " anchored module(s) below min_size dropped")
    }
    hog_modules <- hog_modules[!small]
  }

  map <- hog_map[hog_map$hog %in% unlist(hog_modules), , drop = FALSE]
  modules <- lapply(names(rg$genes), function(s) {
    m <- map[map$species == s & map$gene %in% rg$genes[[s]], , drop = FALSE]
    sets <- lapply(hog_modules, function(h) {
      unique(as.character(m$gene[m$hog %in% h]))
    })
    sets <- sets[lengths(sets) > 0L]
    all_g <- unlist(sets, use.names = FALSE)
    if (!anyDuplicated(all_g)) {
      return(as_modules(sets))
    }
    lab <- rep(names(sets), lengths(sets))
    first <- !duplicated(all_g)
    list(
      modules = stats::setNames(lab[first], all_g[first]),
      module_genes = sets, n_modules = length(sets),
      method = "external", params = list(min_size = 1L)
    )
  })
  names(modules) <- names(rg$genes)
  list(
    modules = modules, hog_modules = hog_modules, method = method,
    params = list(
      resolution = resolution, anchors = anchors,
      min_size = min_size
    )
  )
}


#' Densest weighted subgraph containing an anchor (Goldberg 1984)
#'
#' Bisection on the density g; at each g a max flow on the standard
#' construction with the anchor tied to the source. Returns node names.
#' @noRd
.anchored_densest <- function(g, anchor, tol = 1e-3) {
  n <- igraph::vcount(g)
  if (n < 2L) {
    return(anchor)
  }
  w <- igraph::E(g)$weight
  deg <- igraph::strength(g)
  m <- sum(w)
  el <- igraph::as_edgelist(g, names = FALSE)
  ai <- match(anchor, igraph::V(g)$name)
  src <- n + 1L
  snk <- n + 2L
  h <- igraph::make_graph(rbind(
    c(el[, 1L], el[, 2L], rep(src, n), seq_len(n)),
    c(el[, 2L], el[, 1L], seq_len(n), rep(snk, n))
  ), n = n + 2L, directed = TRUE)
  base <- c(w, w, rep(m, n))
  big <- 2 * m * n + sum(deg) + 1
  base[2L * nrow(el) + ai] <- big
  lo <- 0
  hi <- max(deg)
  best <- ai
  while (hi - lo > tol * max(1, hi)) {
    gam <- (lo + hi) / 2
    mf <- igraph::max_flow(h, src, snk, capacity = c(base, m + 2 * gam - deg))
    # cut = m n + 2 min_{S containing anchor} |S| (g - density(S))
    if (mf$value < m * n * (1 - 1e-9)) {
      lo <- gam
      best <- setdiff(as.integer(mf$partition1), src)
    } else {
      hi <- gam
    }
  }
  igraph::V(g)$name[best]
}
