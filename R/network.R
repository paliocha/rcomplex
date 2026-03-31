# ---- Torch helpers --------------------------------------------------------

#' Select torch device and dtype
#'
#' Probes GPU capability with a small matmul smoke test. CUDA prefers float64
#' but falls back to float32 if float64 kernels are missing (e.g. Blackwell
#' sm_100+ with torch R 0.16.3). MPS only supports float32. The FE matrix
#' computation is numerically exact in float32 (binary/integer arithmetic,
#' values < 2^24).
#'
#' @return List with `device` (string) and `dtype` (torch dtype object).
#' @noRd
torch_device_dtype <- function() {
  if (torch::cuda_is_available()) {
    # Try float64 first (full precision), then float32 (still exact for FE).
    # Older torch builds may lack float64 kernels for newer GPU architectures
    # (e.g. Blackwell sm_100+ with torch R 0.16.3).
    for (dt in list(torch::torch_float64(), torch::torch_float32())) {
      ok <- tryCatch({
        m <- matrix(1.0, 64L, 64L)
        t <- torch::torch_tensor(m, dtype = dt, device = "cuda")
        r <- t$mm(t)
        result <- as.matrix(r$cpu())
        (abs(result[1, 1] - 64) < 0.01)
      }, error = function(e) FALSE)
      if (ok) {
        if (as.character(dt) == "Float") {
          message("CUDA float64 unavailable; using CUDA float32")
        }
        return(list(device = "cuda", dtype = dt))
      }
    }
    message("CUDA smoke tests failed; falling back to CPU torch")
  }
  if (torch::backends_mps_is_available()) {
    list(device = "mps", dtype = torch::torch_float32())
  } else {
    list(device = "cpu", dtype = torch::torch_float64())
  }
}

#' Flush stale GPU memory
#'
#' Runs R garbage collection to free dead torch external pointers, then
#' releases the CUDA caching allocator's free blocks. The CUDA empty-cache
#' call is only made when CUDA is available; no-op on MPS/CPU.
#'
#' @noRd
.gpu_gc <- function() {
  gc()
  if (torch::cuda_is_available()) torch::cuda_empty_cache()
}


# ---- Correlation backends ------------------------------------------------

#' Compute gene-gene correlation matrix via torch
#'
#' Uses GPU (CUDA/MPS) when available, otherwise torch on CPU.
#' Pearson: center rows, L2-normalize, GEMM.
#' Spearman: rank rows first, then Pearson on ranks.
#'
#' @param x Numeric matrix (genes x samples).
#' @param method `"pearson"` or `"spearman"`.
#' @return Correlation matrix (genes x genes) as a base R matrix.
#' @noRd
cor_torch <- function(x, method = "pearson") {
  dd <- torch_device_dtype()

  if (method == "spearman") {
    x <- t(apply(x, 1, rank))
  }

  x <- torch::torch_tensor(x, dtype = dd$dtype, device = dd$device)
  x <- x - x$mean(dim = 2L, keepdim = TRUE)
  norms <- x$norm(dim = 2L, keepdim = TRUE)$clamp(min = 1e-30)
  x <- x / norms
  cor_mat <- x$mm(x$t())

  as.matrix(cor_mat$cpu()$to(dtype = torch::torch_float64()))
}

#' Compute gene-gene correlation matrix via Rfast
#'
#' @param x Numeric matrix (genes x samples).
#' @param method `"pearson"` or `"spearman"`.
#' @return Correlation matrix (genes x genes).
#' @noRd
cor_rfast <- function(x, method = "pearson") {
  if (method == "pearson") {
    Rfast::cora(t(x))
  } else {
    Rfast::cora(apply(x, 1, rank))
  }
}

#' Compute co-expression network
#'
#' Calculates correlation, applies normalization (Mutual Rank or CLR),
#' and determines a density-based co-expression threshold. Accepts either
#' a numeric matrix or a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#'
#' @param x Expression data: a numeric matrix (genes x samples) with row
#'   names as gene identifiers, or a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @param ... Arguments passed to methods (see below).
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
#' @param use_torch If `TRUE`, use torch for GPU-accelerated correlation
#'   (CUDA, MPS, or CPU fallback). Requires the
#'   \href{https://torch.mlverse.org/}{torch} package. Default `FALSE`.
#'   On Apple Silicon (MPS), torch uses float32 because Metal does not support
#'   float64. Pearson correlation is accurate to ~1e-7, but Spearman + MR
#'   normalization can exhibit rank-swap artifacts. For exact Spearman + MR on
#'   MPS, keep `use_torch = FALSE` here and use
#'   \code{\link{permutation_hog_test}(use_torch = TRUE)} for the permutation
#'   speedup instead. On CUDA, float64 is used with no precision tradeoff.
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
#' @examples
#' \dontrun{
#' # From a matrix:
#' net <- compute_network(x, cor_method = "spearman",
#'                        norm_method = "mr", density = 0.03)
#'
#' # From a SummarizedExperiment:
#' net <- compute_network(se, assay = "vst", cor_method = "spearman")
#' }
#'
#' @rdname compute_network
#' @export
setGeneric("compute_network", function(x, ...) standardGeneric("compute_network"))

#' @rdname compute_network
#' @export
setMethod("compute_network", "matrix", function(x,
                            cor_method = c("pearson", "spearman"),
                            norm_method = c("MR", "CLR"),
                            density = 0.03,
                            abs_cor = FALSE,
                            mr_log_transform = FALSE,
                            min_var = 0,
                            n_cores = 1L,
                            use_torch = FALSE) {
  cor_method <- match.arg(cor_method)
  norm_method <- match.arg(norm_method)
  if (is.null(rownames(x))) {
    stop("x must have row names (gene identifiers)")
  }
  if (density <= 0 || density >= 1) {
    stop("density must be between 0 and 1 (exclusive)")
  }
  if (use_torch && !requireNamespace("torch", quietly = TRUE)) {
    stop("use_torch = TRUE requires the torch package ",
         "(install.packages('torch'); torch::install_torch())")
  }

  # Filter low-variance genes
  n_removed <- 0L
  if (!is.null(min_var)) {
    row_var <- rowSums((x - rowMeans(x))^2) /
      (ncol(x) - 1L)
    keep <- row_var > min_var
    n_removed <- sum(!keep)
    if (n_removed > 0L) {
      x <- x[keep, , drop = FALSE]
    }
    if (nrow(x) < 3L) {
      stop("Fewer than 3 genes remain after variance filtering (min_var = ",
           min_var, ")")
    }
  }

  gene_names <- rownames(x)
  n_genes <- nrow(x)

  # Correlation
  cor_fn <- if (use_torch) cor_torch else cor_rfast
  net <- cor_fn(x, method = cor_method)
  if (use_torch) .gpu_gc()

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
})

#' @rdname compute_network
#' @param assay Assay name or index to extract from the
#'   SummarizedExperiment (default 1). Requires the
#'   \pkg{SummarizedExperiment} package.
#' @export
setMethod("compute_network", "SummarizedExperiment", function(x,
    assay = 1L, ...) {
  expr <- SummarizedExperiment::assay(x, assay)
  compute_network(expr, ...)
})
