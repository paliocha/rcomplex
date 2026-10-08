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
      ok <- tryCatch(
        {
          m <- matrix(1.0, 64L, 64L)
          t <- torch::torch_tensor(m, dtype = dt, device = "cuda")
          r <- t$mm(t)
          result <- as.matrix(r$cpu())
          (abs(result[1, 1] - 64) < 0.01)
        },
        error = function(e) FALSE
      )
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

# Samples x genes matrix zt with crossprod(zt) equal to the correlation
# matrix cor_rfast() builds (Spearman ranks each gene first).
.standardise_for_cor <- function(x, cor_method) {
  xs <- if (cor_method == "pearson") t(x) else apply(x, 1, rank)
  mat <- t(xs) - Rfast::colmeans(xs)
  t(mat / sqrt(Rfast::rowsums(mat^2)))
}

# Weakest correlation among the gene pairs at or above `thr` in `net`:
# the minimum for sign = "positive", the maximum for "negative". With
# `levels` (the kept partition levels) it is the rectified mean over the
# levels, as compute_network() averages it, with the sign put back.
# Correlations are recomputed for the passing pairs only, in chunks of
# about 2^24 doubles, so dense and blockwise builds report the same value.
.r_threshold <- function(x, net, thr, cor_method, sign, partition = NULL,
                         levels = NULL) {
  e <- .adj_edges(net, thr)
  up <- e$rows < e$cols
  i <- e$rows[up]
  j <- e$cols[up]
  if (length(i) == 0L) {
    return(NA_real_)
  }
  s <- if (sign == "negative") -1 else 1
  groups <- if (is.null(partition)) {
    list(seq_len(ncol(x)))
  } else {
    lapply(levels, function(lv) which(partition == lv))
  }
  r <- 0
  for (g in groups) {
    zt <- .standardise_for_cor(x[, g, drop = FALSE], cor_method)
    chunk <- max(1L, 2^24 %/% nrow(zt))
    rg <- unlist(lapply(
      split(seq_along(i), (seq_along(i) - 1L) %/% chunk),
      function(k) colSums(zt[, i[k], drop = FALSE] * zt[, j[k], drop = FALSE])
    ), use.names = FALSE)
    r <- r + if (is.null(partition)) s * rg else pmax(s * rg, 0, na.rm = TRUE)
  }
  s * min(r / length(groups))
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
#' @param sign `"positive"` (default) ranks the strongest positive
#'   correlations first. `"negative"` negates every correlation between
#'   two genes before normalization, so the strongest anticorrelations
#'   rank first. With `partition`, the function negates the correlation
#'   of each level before it sets negative values to zero.
#' @param mr_log_transform If `FALSE` (default), use the raw MR formula
#'   matching the original RComPlEx R Markdown. If `TRUE`, use Obayashi &
#'   Kinoshita (2009) log-normalized formula (values in \[0,1\]).
#' @param sparse If `TRUE` (default), return the network as a sparse
#'   `Matrix::dgCMatrix` holding both triangles of the thresholded
#'   co-expression matrix (entries below `store_threshold` are discarded;
#'   the dense n x n matrix never leaves this function). If `FALSE`,
#'   return the dense matrix as in rcomplex < 0.2.0.
#' @param store_density Fraction of top edges to keep in the sparse store
#'   (default `NULL` = `max(density, 0.05)`). Must satisfy
#'   `density <= store_density < 1`; the margin above `density` lets
#'   downstream sweeps loosen the threshold (e.g. `density_sweep()`)
#'   without recomputing the network. Only valid with `sparse = TRUE`.
#' @param n_cores Number of threads for parallel computation (default 1).
#' @param block_size `NULL` (default) builds the dense n x n matrix first.
#'   A positive whole number builds the network a block of genes at a
#'   time and never forms the n x n matrix (each block holds
#'   8 * n * `block_size` bytes, so a `block_size` near n saves nothing;
#'   values above n are capped at n with a message; threads rank the
#'   columns of one block, so keep it well above `n_cores`, as the
#'   256-1024 benchmarked below): each gene keeps only its
#'   partners in the top fraction of its correlation ranks, which at the
#'   default `store_density` is 15% of all pairs held at 8 bytes each. On
#'   BDIS leaf data (20 samples, n = 20,000, 8 threads) peak memory was
#'   1.2-1.4 GB against 5.1 GB for the dense build, in 7.6-7.9 s against
#'   8.5 s. The result equals the dense build up to floating-point
#'   near-ties between correlations: the two builds compute correlations
#'   with different BLAS calls, so two correlations of one gene that differ
#'   only in the last bits can be ranked in the other order, shifting an MR
#'   value. With `cor_method = "spearman"` and few samples, correlations
#'   that are equal as rationals are common and may round apart, so the
#'   result can then differ from the dense build in many entries and
#'   between block sizes. Requires `sparse = TRUE` and `norm_method = "MR"`.
#'   With `mr_log_transform = TRUE` a pair is kept
#'   when either gene ranks the other in its top 10% (twice
#'   `store_density`) and a second correlation pass reads the other rank;
#'   on the same data peak memory was 1.7 GB against 5.1 GB dense in
#'   9.2 s against 8.2 s. If the store threshold cannot prove the kept
#'   pairs complete the fraction widens and a message says so; peak memory
#'   grows with it (for log MR it can exceed the dense build), and at all
#'   pairs the build saves no memory.
#' @param partition `NULL` (default), or a per-sample factor of length
#'   `ncol(x)` without `NA`, such as tissue or study. The function
#'   computes the correlation within each level. It sets negative
#'   correlations to zero and takes the mean over the levels. Then it
#'   normalises the result. The function drops levels with fewer than 5
#'   samples and tells you which. At least two levels must remain.
#'   Use the dense build: `block_size` must be `NULL`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{network}{Named symmetric matrix (genes x genes) with normalized
#'       co-expression values. A `dgCMatrix` storing both triangles of the
#'       entries at or above `store_threshold` (diagonal absent) when
#'       `sparse = TRUE`; a dense matrix with zero diagonal when
#'       `sparse = FALSE`.}
#'     \item{threshold}{The co-expression threshold at the given density,
#'       equal to that of the full dense matrix in every mode.}
#'     \item{n_genes}{Number of genes in the network.}
#'     \item{n_removed}{Number of constant genes removed before computing
#'       correlations.}
#'     \item{params}{List of parameters used. It holds `partition` when
#'       you give one. `r_threshold` is the weakest correlation of an edge:
#'       the smallest for `sign = "positive"`, the largest for
#'       `"negative"`. With `partition` it is the rectified mean over the
#'       levels, negated for `"negative"`. It is `NA` when no pair passes.}
#'     \item{store_density, store_threshold}{Sparse networks only: the
#'       stored edge fraction and the value cutoff of the store. Analyses
#'       at thresholds below `store_threshold` are refused (see
#'       [as_sparse_network()]).}
#'   }
#'
#' @details
#' With `sign = "negative"` a gene keeps its own correlation of 1, so it
#' ranks itself first in both signs. Negative correlations are rarer and
#' weaker than positive ones in RNA-seq data. The density threshold still
#' keeps the top fraction of pairs, so a negative network exists at any
#' sample size even when it holds only noise: at n = 20 samples, r >= -0.3
#' is noise. Check the weakest correlation that passed the threshold
#' before you read a negative network.
#'
#' @section Gene universe:
#' Networks are built on all supplied genes and downstream tests use the
#' whole network as the hypergeometric population (Netotea et al. 2014).
#' The canonical ComPlEx implementations restrict the expression tables to
#' genes with an ortholog before building the networks, which changes
#' neighbourhoods, thresholds and calls; rcomplex reproduces them only when
#' the expression matrices are restricted the same way first (see
#' `compare_neighborhoods()`).
#'
#' @section Reconstructing sub-threshold values:
#' A sparse network (`sparse = TRUE`) discards values below
#' `store_threshold`, so visualisations of a gene subset (e.g. module
#' heatmaps) cannot read them back from `net$network`. [mr_block()]
#' reconstructs the exact mutual-rank values for a gene subset from the
#' expression matrix and the network's parameters, without rebuilding the
#' dense n x n matrix.
#'
#' @examples
#' \dontrun{
#' # From a matrix:
#' net <- compute_network(x,
#'   cor_method = "spearman",
#'   norm_method = "mr", density = 0.03
#' )
#'
#' # From a SummarizedExperiment:
#' net <- compute_network(se, assay = "vst", cor_method = "spearman")
#' }
#'
#' @rdname compute_network
#' @export
setGeneric(
  "compute_network",
  function(x, ...) standardGeneric("compute_network")
)

#' @rdname compute_network
#' @export
setMethod("compute_network", "matrix", function(
  x,
  cor_method = c("pearson", "spearman"),
  norm_method = c("MR", "CLR"),
  density = 0.03,
  sign = c("positive", "negative"),
  mr_log_transform = FALSE,
  sparse = TRUE,
  store_density = NULL,
  n_cores = 1L,
  block_size = NULL,
  partition = NULL) {
  cor_method <- match.arg(cor_method)
  norm_method <- match.arg(norm_method)
  sign <- match.arg(sign)
  if (is.null(rownames(x))) {
    stop("x must have row names (gene identifiers)")
  }
  if (density <= 0 || density >= 1) {
    stop("density must be between 0 and 1 (exclusive)")
  }
  if (sparse) {
    if (is.null(store_density)) {
      store_density <- max(density, 0.05)
    }
    if (!is.numeric(store_density) || length(store_density) != 1L ||
          store_density < density || store_density >= 1) {
      stop(
        "store_density must satisfy density <= store_density < 1 (got ",
        store_density, " with density = ", density, ")"
      )
    }
  } else if (!is.null(store_density)) {
    stop("store_density requires sparse = TRUE")
  }
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1L ||
          is.na(block_size) || block_size < 1 ||
          block_size != round(block_size)) {
      stop("block_size must be NULL or a positive whole number")
    }
    if (!sparse) stop("block_size requires sparse = TRUE")
    if (norm_method != "MR") stop("block_size requires norm_method = \"MR\"")
  }
  levels_kept <- NULL
  if (!is.null(partition)) {
    if (!is.null(block_size)) {
      stop("partition needs the dense build; set block_size = NULL")
    }
    if (length(partition) != ncol(x)) {
      stop(
        "partition must have one entry per sample (ncol(x) = ", ncol(x), ")"
      )
    }
    if (anyNA(partition)) stop("partition must not contain NA")
    partition <- droplevels(as.factor(partition))
    min_partition_n <- 5L
    sizes <- table(partition)
    small <- names(sizes)[sizes < min_partition_n]
    if (length(small) > 0L) {
      message(
        "Dropped partition levels with fewer than ", min_partition_n,
        " samples: ", paste(small, collapse = ", ")
      )
    }
    levels_kept <- setdiff(names(sizes), small)
    if (length(levels_kept) < 2L) {
      stop(
        "partition needs at least two levels with ", min_partition_n,
        " or more samples"
      )
    }
  }

  # Filter constant genes
  row_var <- rowSums((x - rowMeans(x))^2) /
    (ncol(x) - 1L)
  # A constant row can come out at ~1e-30 instead of 0 in floating point
  # and pass `> 0`; its correlations are then NaN. Test constancy exactly.
  keep <- row_var > 0 & rowSums(x != x[, 1L]) > 0L
  n_removed <- sum(!keep)
  if (n_removed > 0L) {
    x <- x[keep, , drop = FALSE]
  }
  if (nrow(x) < 3L) {
    stop("Fewer than 3 genes remain after variance filtering")
  }

  gene_names <- rownames(x)
  n_genes <- nrow(x)

  params <- list(
    cor_method = cor_method,
    norm_method = norm_method,
    density = density,
    sign = sign,
    mr_log_transform = mr_log_transform
  )
  # partition applies the sign per level, before the rectification
  negate <- sign == "negative" && is.null(partition)
  params$partition <- partition

  if (!is.null(block_size)) {
    if (block_size >= n_genes) {
      message(
        "block_size >= the ", n_genes, " genes: one block holds the whole ",
        "correlation matrix, so the blockwise build saves no memory"
      )
    }
    slots <- mr_block_network_cpp(
      .standardise_for_cor(x, cor_method), mr_log_transform, negate,
      density, store_density, as.integer(min(block_size, n_genes)), n_cores
    )
    # log MR holds a reverse index and its values, so a wide fraction can
    # need more memory than the dense build
    more <- if (mr_log_transform) {
      "; peak memory may exceed the dense build"
    } else {
      ""
    }
    if (slots$fraction >= 1) {
      message(
        "Blockwise build fell back to all pairs and saved no memory",
        more, " (lower store_density to save memory)"
      )
    } else if (slots$fraction > slots$start_fraction) {
      message(sprintf(
        "Blockwise build widened its candidate fraction from %.3g to %.3g%s",
        slots$start_fraction, slots$fraction, more
      ))
    }
    network <- methods::new(
      "dgCMatrix",
      i = slots$i, p = slots$p, x = slots$x,
      Dim = c(n_genes, n_genes),
      Dimnames = list(gene_names, gene_names)
    )
    params$r_threshold <- .r_threshold(
      x, network, slots$threshold, cor_method, sign
    )
    return(list(
      network = network,
      threshold = slots$threshold,
      n_genes = n_genes,
      n_removed = n_removed,
      params = c(params, list(store_density = store_density)),
      store_density = store_density,
      store_threshold = slots$store_threshold
    ))
  }

  # Correlation
  if (is.null(partition)) {
    net <- cor_rfast(x, method = cor_method)
  } else {
    # Rectified average over levels (TEA-GCN): negative correlations count
    # as zero, and so do the NaN correlations of a gene constant in a level.
    net <- 0
    for (lv in levels_kept) {
      r <- cor_rfast(x[, partition == lv, drop = FALSE], method = cor_method)
      if (sign == "negative") r <- -r
      net <- net + pmax(r, 0, na.rm = TRUE)
    }
    net <- net / length(levels_kept)
    diag(net) <- 1
  }

  # Normalization
  if (norm_method == "MR") {
    # Clip to [-1, 1], negation (if negate), MR ranks and zero diagonal are
    # all done in C++ directly on `net`, which is freshly allocated by cor_fn
    # (refcount 1): in-place mutation is intentional (no n x n temporaries).
    mutual_rank_inplace_cpp(net, mr_log_transform, negate, n_cores)
  } else {
    # Clip to [-1, 1]
    net[net > 1] <- 1
    net[net < -1] <- -1

    if (negate) {
      net <- -net
      diag(net) <- -diag(net)
    }

    net <- apply_clr_to_cor_cpp(net, n_cores = n_cores)

    # Set diagonal to 0
    diag(net) <- 0
  }

  # Assign gene names
  dimnames(net) <- list(gene_names, gene_names)

  # Compute density threshold (always from the full dense matrix)
  thr <- density_threshold_cpp(net, density)

  r_thr <- function(m) {
    .r_threshold(x, m, thr, cor_method, sign, partition, levels_kept)
  }
  if (!sparse) {
    params$r_threshold <- r_thr(net)
    return(list(
      network = net,
      threshold = thr,
      n_genes = n_genes,
      n_removed = n_removed,
      params = params
    ))
  }

  # Sparse store: threshold the dense matrix at store_density, repack the
  # surviving entries as dgCMatrix slots, and free the dense matrix. Field
  # order matches as_sparse_network() so both constructions are identical.
  store_thr <- density_threshold_cpp(net, store_density)
  slots <- extract_sparse_cpp(net, store_thr, n_cores)
  rm(net)
  spnet <- methods::new(
    "dgCMatrix",
    i = slots$i, p = slots$p, x = slots$x,
    Dim = c(n_genes, n_genes),
    Dimnames = list(gene_names, gene_names)
  )
  params$r_threshold <- r_thr(spnet)
  list(
    network = spnet,
    threshold = thr,
    n_genes = n_genes,
    n_removed = n_removed,
    params = c(params, list(store_density = store_density)),
    store_density = store_density,
    store_threshold = store_thr
  )
})

#' @rdname compute_network
#' @param assay Assay name or index to extract from the
#'   SummarizedExperiment (default 1). Requires the
#'   \pkg{SummarizedExperiment} package.
#' @export
setMethod("compute_network", "SummarizedExperiment", function(
  x,
  assay = 1L, ...
) {
  # No requireNamespace guard needed: S4 dispatch to this method
  # requires the SummarizedExperiment class (and package) to be loaded.
  expr <- SummarizedExperiment::assay(x, assay)
  compute_network(expr, ...)
})
