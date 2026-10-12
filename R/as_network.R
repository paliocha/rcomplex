#' Import a network built outside rcomplex
#'
#' Wraps given co-expression weights in the network object that
#' [compute_network()] returns. Every consumer then reads it unchanged.
#'
#' @param x A symmetric numeric matrix, dense or sparse, with gene names
#'   as dimnames. Or an edge list: a data frame with columns `gene1`,
#'   `gene2` and `weight`.
#' @param density Fraction of gene pairs kept as edges, in (0, 1).
#' @param genes Character vector of genes, or `NULL` for all genes in `x`.
#'   It sets the gene order. Genes absent from an edge list get no edges.
#'
#' @return A list as [compute_network()] returns with `sparse = TRUE`:
#'   `network` (a `dgCMatrix` with both triangles), `threshold` (the
#'   weight at the `density` quantile), `n_genes`, `params`,
#'   `store_density` and `store_threshold`.
#'
#' @details
#' Weights are taken as given, for example a TEA-GCN
#' `zScore(Co-exp_Str_MR)` or a WGCNA adjacency. `as_network()` does not
#' recompute mutual ranks; that is the job of [compute_network()] on
#' expression. Larger weights mean stronger co-expression. A dense
#' matrix is stored as [compute_network()] stores it, at the
#' `max(density, 0.05)` quantile. A sparse matrix or an edge list is
#' stored in full, and a pair it omits ranks below every given weight.
#' The threshold is the weight of the k-th strongest of all
#' `n (n - 1) / 2` pairs, with k as in [compute_network()]. It must be a
#' given weight, so `x` must hold at least k pairs. Functions that need
#' expression, such as [null_network()], do not work on the result.
#'
#' @examples
#' x <- matrix(rnorm(400), 40, 10, dimnames = list(paste0("g", 1:40)))
#' net <- as_network(compute_network(x)$network)
#' net$threshold
#'
#' @export
as_network <- function(x, density = 0.03, genes = NULL) {
  ok <- is.numeric(density) && length(density) == 1L && !is.na(density)
  if (!ok || density <= 0 || density >= 1) {
    stop("density must be a single number in (0, 1)")
  }
  if (is.data.frame(x)) {
    x <- .edges_to_dgc(x, genes)
  } else if (methods::is(x, "sparseMatrix")) {
    x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
  } else if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix, a sparse matrix or an edge list")
  }
  dn <- dimnames(x)
  if (is.null(dn[[1L]]) || !identical(dn[[1L]], dn[[2L]])) {
    stop("x must have the same gene names on rows and columns")
  }
  if (!is.null(genes)) {
    miss <- setdiff(genes, dn[[1L]])
    if (length(miss) > 0L) {
      stop("genes not in x: ", paste(utils::head(miss, 5L), collapse = ", "))
    }
    x <- x[genes, genes, drop = FALSE]
  }
  if (anyNA(x)) stop("x must not contain NA")
  if (!Matrix::isSymmetric(x)) stop("x must be symmetric")
  n <- nrow(x)
  if (n < 3L) stop("x must have at least 3 genes")

  if (is.matrix(x)) {
    thr <- density_threshold_cpp(x, density)
    net <- list(network = x, threshold = thr, params = list(density = density))
    net <- as_sparse_network(net, max(density, 0.05))
  } else {
    x <- methods::as(x, "TsparseMatrix")
    off <- x@i != x@j
    x <- Matrix::sparseMatrix(
      i = x@i[off] + 1L, j = x@j[off] + 1L, x = x@x[off],
      dims = c(n, n), dimnames = dimnames(x)
    )
    up <- x@x[x@i < rep.int(seq_len(n) - 1L, diff(x@p))]
    n_pairs <- n * (n - 1) / 2
    k <- min(max(floor(density * n_pairs + 0.5), 1), n_pairs - 1)
    if (length(up) < k) {
      stop("x holds ", length(up), " gene pairs, fewer than the ", k, " kept")
    }
    store_density <- length(up) / n_pairs
    net <- list(
      network = x,
      threshold = sort(up, decreasing = TRUE)[k],
      params = list(density = density, store_density = store_density),
      store_density = store_density,
      store_threshold = min(up)
    )
  }
  append(net, list(n_genes = n), after = 2L)
}


#' Edge list (gene1, gene2, weight) to a symmetric dgCMatrix
#' @noRd
.edges_to_dgc <- function(x, genes) {
  if (!all(c("gene1", "gene2", "weight") %in% names(x))) {
    stop("an edge list must have columns gene1, gene2, weight")
  }
  if (!is.numeric(x$weight)) stop("weight must be numeric")
  g1 <- as.character(x$gene1)
  g2 <- as.character(x$gene2)
  if (is.null(genes)) genes <- unique(c(g1, g2))
  i <- match(g1, genes)
  j <- match(g2, genes)
  keep <- !is.na(i) & !is.na(j) & i != j
  e <- unique(data.frame(
    lo = pmin(i, j)[keep], hi = pmax(i, j)[keep], w = x$weight[keep]
  ))
  if (anyDuplicated(e[c("lo", "hi")])) {
    stop("x gives a gene pair two weights")
  }
  Matrix::sparseMatrix(
    i = c(e$lo, e$hi), j = c(e$hi, e$lo), x = as.numeric(c(e$w, e$w)),
    dims = rep(length(genes), 2L), dimnames = list(genes, genes)
  )
}
