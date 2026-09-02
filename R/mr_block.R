# Exact local reconstruction of mutual-rank values for a gene subset.
#
# MR_ij = sqrt(R_ij * R_ji) only needs the correlations of genes i and j
# ranked over all n network genes, so a k x k block (including entries a
# sparse network discarded below store_threshold) is computed from a k x n
# correlation slice without rebuilding the dense n x n matrix.


#' Reconstruct mutual-rank values for a gene subset
#'
#' Computes the exact Mutual Rank co-expression values for a subset of
#' genes from the expression matrix, using the network's stored parameters
#' (`cor_method`, `abs_cor`, `mr_log_transform`). A sparse network
#' (`compute_network(sparse = TRUE)`) discards values below its
#' `store_threshold`; this function reconstructs them (for heatmaps or
#' other visualisations of, say, a module) without rebuilding the dense
#' n x n matrix: \eqn{MR_{ij} = \sqrt{R_{ij} R_{ji}}} only needs the
#' correlations of genes i and j ranked over all n network genes, so a
#' k x n correlation slice suffices.
#'
#' The gene universe is `rownames(net$network)`: `x` is subset to those
#' rows first (the `min_var` filter was already applied when the network
#' was built), and ranks span exactly the network genes. Pass the same
#' expression matrix that built the network; the reconstructed block then
#' matches the dense `compute_network(sparse = FALSE)` matrix.
#'
#' Reconstruction is exact up to ties in the correlation values: exactly
#' tied correlations get average ranks, and a mathematical tie can be
#' broken differently by this function (`stats::cor()`) and the network's
#' correlation backend (`Rfast::cora()`), shifting a rank by 0.5. Ties are
#' measure-zero for Pearson on continuous data, but Spearman correlations
#' lie on a finite grid where occasional ties are possible.
#'
#' @param x Expression matrix (genes x samples) the network was built
#'   from. Must contain every network gene as a row; extra rows (e.g.
#'   genes removed by the `min_var` filter) are ignored.
#' @param genes Character vector of gene identifiers (no duplicates), all
#'   present in `rownames(net$network)`.
#' @param net Network object from [compute_network()] (sparse or dense)
#'   with `norm_method = "MR"`.
#' @return A k x k numeric matrix of mutual-rank values with
#'   `dimnames = list(genes, genes)` and zero diagonal, on the same scale
#'   as the network (raw MR, or log-normalized when the network was built
#'   with `mr_log_transform = TRUE`).
#'
#' @examples
#' \dontrun{
#' net <- compute_network(x, density = 0.03)  # sparse
#' # exact values for one module, including sub-threshold entries:
#' blk <- mr_block(x, module_genes, net)
#' heatmap(blk)
#' }
#'
#' @export
mr_block <- function(x, genes, net) {
  universe <- rownames(net$network)
  if (is.null(universe)) {
    stop("net$network must have gene row names")
  }
  params <- net$params
  if (is.null(params$cor_method)) {
    stop("net must carry the parameters of compute_network() ",
         "(net$params$cor_method is missing)")
  }
  if (!identical(params$norm_method, "MR")) {
    stop("mr_block() requires a network built with norm_method = \"MR\"; ",
         "got \"", params$norm_method, "\"")
  }
  if (!is.character(genes) || length(genes) == 0L) {
    stop("genes must be a non-empty character vector")
  }
  if (anyDuplicated(genes) > 0L) {
    stop("genes must not contain duplicates")
  }
  missing_genes <- setdiff(genes, universe)
  if (length(missing_genes) > 0L) {
    stop("genes not in the network: ",
         paste(missing_genes, collapse = ", "))
  }
  if (is.null(rownames(x)) || !all(universe %in% rownames(x))) {
    stop("x must contain every network gene as a row")
  }

  # Same gene universe as the network (min_var filter already applied
  # there); ranks below must span exactly the n network genes.
  x <- x[universe, , drop = FALSE]
  n <- length(universe)

  # k x n correlation slice: genes vs the full universe (self included,
  # as in the full matrix). Clamp and abs as mutual_rank_inplace_cpp().
  cm <- stats::cor(t(x[genes, , drop = FALSE]), t(x),
                   method = params$cor_method)
  cm <- pmin(pmax(cm, -1), 1)  # cm first: pmin/pmax keep its dim
  if (isTRUE(params$abs_cor)) {
    cm <- abs(cm)
  }

  # Row i = average ranks of gene i's correlations over all n network
  # genes: ascending for raw MR, descending (rank the negated values) for
  # log MR - same tie handling as mutual_rank_inplace_cpp().
  rk <- t(apply(if (isTRUE(params$mr_log_transform)) -cm else cm, 1, rank))

  # block[i, j] = sqrt(R[i, genes[j]] * R[j, genes[i]])
  half <- rk[, match(genes, universe), drop = FALSE]
  block <- sqrt(half * t(half))
  if (isTRUE(params$mr_log_transform)) {
    block <- pmin(pmax(1 - log(block) / log(n), 0), 1)
  }
  diag(block) <- 0
  dimnames(block) <- list(genes, genes)
  block
}
