#' Compute co-expression network from an expression matrix
#'
#' Calculates correlation, applies normalization (Mutual Rank or CLR),
#' and determines a density-based co-expression threshold.
#'
#' @param expr_matrix Numeric matrix of expression values (genes x samples).
#'   Row names must be gene identifiers.
#' @param cor_method Correlation method: `"pearson"` (default) or `"spearman"`.
#' @param norm_method Normalization method: `"MR"` (Mutual Rank, default) or
#'   `"CLR"` (Context Likelihood Ratio).
#' @param density Fraction of top edges to keep (default 0.03 = 3%).
#' @param abs_cor If `TRUE`, take absolute value of correlations before
#'   normalization (default `FALSE`).
#' @param mr_log_transform If `FALSE` (default), use the raw MR formula
#'   matching the original RComPlEx R Markdown. If `TRUE`, use Obayashi &
#'   Kinoshita (2009) log-normalized formula (values in \[0,1\]).
#' @param min_var Minimum row variance threshold. Genes with variance below
#'   this value are removed before computing correlations. Default `0` removes
#'   only constant genes (zero variance). Set to a positive value (e.g. `1e-3`)
#'   to also remove near-invariant genes that produce noisy correlations.
#'   Set to `NULL` to disable filtering entirely.
#' @param n_cores Number of threads for parallel computation (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{network}{Named symmetric matrix (genes x genes) with normalized
#'       co-expression values. Diagonal is set to 0.}
#'     \item{threshold}{The co-expression threshold at the given density.}
#'     \item{n_genes}{Number of genes in the network.}
#'     \item{n_removed}{Number of genes removed by variance filter (0 if
#'       `min_var` is `NULL`).}
#'     \item{params}{List of parameters used.}
#'   }
#'
#' @export
compute_network <- function(expr_matrix,
                            cor_method = c("pearson", "spearman"),
                            norm_method = c("MR", "CLR"),
                            density = 0.03,
                            abs_cor = FALSE,
                            mr_log_transform = FALSE,
                            min_var = 0,
                            n_cores = 1L) {
  cor_method <- match.arg(cor_method)
  norm_method <- match.arg(norm_method)
  n_cores <- as.integer(n_cores)

  if (!is.matrix(expr_matrix)) {
    stop("expr_matrix must be a matrix")
  }
  if (is.null(rownames(expr_matrix))) {
    stop("expr_matrix must have row names (gene identifiers)")
  }
  if (density <= 0 || density >= 1) {
    stop("density must be between 0 and 1 (exclusive)")
  }

  # Filter low-variance genes
  n_removed <- 0L
  if (!is.null(min_var)) {
    row_var <- rowSums((expr_matrix - rowMeans(expr_matrix))^2) /
      (ncol(expr_matrix) - 1L)
    keep <- row_var > min_var
    n_removed <- sum(!keep)
    if (n_removed > 0L) {
      expr_matrix <- expr_matrix[keep, , drop = FALSE]
    }
    if (nrow(expr_matrix) < 3L) {
      stop("Fewer than 3 genes remain after variance filtering (min_var = ",
           min_var, ")")
    }
  }

  gene_names <- rownames(expr_matrix)
  n_genes <- nrow(expr_matrix)

  # Correlation
  if (cor_method == "pearson") {
    net <- Rfast::cora(t(expr_matrix))
  } else {
    ranked <- t(apply(expr_matrix, 1, rank))
    net <- Rfast::cora(t(ranked))
  }

  # Clip to [-1, 1]
  net[net > 1] <- 1
  net[net < -1] <- -1

  if (abs_cor) {
    net <- abs(net)
  }

  # Normalization
  if (norm_method == "MR") {
    net <- mutual_rank_transform_cached_cpp(net,
                                            log_transform = mr_log_transform,
                                            n_cores = n_cores)
  } else {
    net <- apply_clr_to_cor_cpp(net, n_cores = n_cores)
  }

  # Set diagonal to 0
  diag(net) <- 0

  # Assign gene names
  dimnames(net) <- list(gene_names, gene_names)

  # Compute density threshold
  thr <- density_threshold_cpp(net, density)

  list(
    network = net,
    threshold = thr,
    n_genes = n_genes,
    n_removed = n_removed,
    params = list(
      cor_method = cor_method,
      norm_method = norm_method,
      density = density,
      abs_cor = abs_cor,
      mr_log_transform = mr_log_transform,
      min_var = min_var
    )
  )
}
