#' Split expression into wiring and deployment parts
#'
#' The function fits each gene on a per-sample block factor, such as time
#' point, tree or zone. The block means are the deployment part, and the
#' residuals are the wiring part. Build the network from the wiring part
#' with [compute_network()].
#'
#' @details
#' **What each part is for.** The *wiring* layer (the residuals) holds
#' co-expression beyond the shared course: two genes that only follow the
#' same time course have no residual correlation. Build the network from
#' it with [compute_network()]. The *deployment* layer (the block means)
#' says where along the course a gene is expressed. `r2` is the share of
#' each gene's sum of squares the block explains (the unadjusted
#' variance fraction; see the methods article).
#'
#' **Two facts not to miss.** The residuals of one gene sum to zero
#' within each block, so within a block of size \eqn{n_b} they are
#' correlated \eqn{-1/(n_b - 1)} (-1/3 for four replicates). A full
#' shuffle would break that structure, so the shuffled-expression null
#' for a wiring network must be `null_network(..., block = block)`. And
#' the wiring part keeps only \eqn{n_b - 1} degrees of freedom per
#' block (`df_residual` in all, 15 on a 5 x 4 design of 20 samples), so
#' each half of a split-half design carries about 5 df: split-half
#' replication is not informative there, and the gate for a wiring
#' network has to be cross-species.
#'
#' **Which part to cluster.** On the Pooideae leaf and root time series,
#' about 75 % of cross-species gene-level conservation (neighbour AUROC
#' above 0.5) survived in the wiring part (leaf 0.78, root 0.71), and the
#' deployment part conserved weakly. Cluster the wiring part and carry
#' `r2` as a per-gene covariate.
#'
#' @param x Expression matrix (genes x samples) or a
#'   `SummarizedExperiment`, whose first assay is used (the default of
#'   [compute_network()]).
#' @param block Per-sample factor (or vector) of length `ncol(x)`, no
#'   `NA`, with at least two levels and at least two samples per level.
#' @return A list of class `split_layers`:
#'   \describe{
#'     \item{wiring}{Residual matrix, dimnames of `x`. Input for
#'       [compute_network()]; pair it with
#'       `null_network(..., block = block)`.}
#'     \item{deployment}{Genes x levels matrix of block means.}
#'     \item{r2}{Named numeric, per-gene share of the sum of squares
#'       explained by `block` (not adjusted).}
#'     \item{df_residual}{`ncol(x) - nlevels(block)`.}
#'     \item{block}{List with `factor` (the block factor) and `sizes`
#'       (samples per level).}
#'     \item{params}{List with `n_genes` and `n_samples`.}
#'   }
#' @examples
#' set.seed(1)
#' block <- factor(rep(c("T1", "T2", "T3"), each = 4))
#' x <- matrix(rnorm(10 * 12), 10, 12,
#'   dimnames = list(paste0("g", 1:10), paste0("s", 1:12))
#' )
#' sl <- split_layers(x, block)
#' sl
#' net <- compute_network(sl$wiring, density = 0.2)
#' @export
split_layers <- function(x, block) {
  if (methods::is(x, "SummarizedExperiment")) {
    x <- SummarizedExperiment::assay(x, 1L)
  }
  x <- as.matrix(x)
  if (length(block) != ncol(x)) {
    stop("block must have one entry per sample (ncol(x) = ", ncol(x), ")")
  }
  if (anyNA(block)) stop("block must not contain NA")
  block <- droplevels(as.factor(block))
  sizes <- table(block)
  if (length(sizes) < 2L) stop("block must have at least two levels")
  if (any(sizes < 2L)) {
    stop(
      "every block level needs at least two samples (singleton: ",
      paste(names(sizes)[sizes < 2L], collapse = ", "), ")"
    )
  }
  # OLS on ~ 0 + block: the coefficients are the block means.
  deployment <- t(rowsum(t(x), block)) / rep(sizes, each = nrow(x))
  colnames(deployment) <- levels(block)
  wiring <- x - deployment[, as.integer(block), drop = FALSE]
  dimnames(wiring) <- dimnames(x)
  r2 <- 1 - rowSums(wiring^2) / rowSums((x - rowMeans(x))^2)
  names(r2) <- rownames(x)
  structure(list(
    wiring = wiring,
    deployment = deployment,
    r2 = r2,
    df_residual = ncol(x) - nlevels(block),
    block = list(
      factor = block,
      sizes = stats::setNames(as.integer(sizes), names(sizes))
    ),
    params = list(n_genes = nrow(x), n_samples = ncol(x))
  ), class = "split_layers")
}


#' @export
print.split_layers <- function(x, ...) {
  s <- x$block$sizes
  cat(
    "Split layers:", x$params$n_genes, "genes x", x$params$n_samples,
    "samples\n"
  )
  cat("  blocks:", paste0(names(s), " (", s, ")", collapse = ", "), "\n")
  cat("  df_residual:", x$df_residual, "\n")
  cat(
    "  median r2:", signif(stats::median(x$r2, na.rm = TRUE), 3), "\n"
  )
  invisible(x)
}
