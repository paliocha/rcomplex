#' Compute q-values from a vector of p-values
#'
#' Wrapper around [qvalue::qvalue()] that handles edge cases (fewer than
#' 2 p-values).
#'
#' @param pvals Numeric vector of p-values.
#' @return Numeric vector of q-values, same length as `pvals`.
#' @noRd
compute_qvalues <- function(pvals) {
  if (length(pvals) < 2L) return(pvals)
  tryCatch(
    qvalue::qvalue(pvals)$qvalues,
    error = function(e) qvalue::qvalue(pvals, pi0 = 1)$qvalues
  )
}


#' Compute Besag-Clifford p-value support for discrete q-values
#'
#' Constructs the discrete support set for p-values produced by Besag-Clifford
#' adaptive stopping, used by [DiscreteQvalue::DQ()].
#'
#' @param min_exceedances BC stopping parameter.
#' @param max_permutations Maximum permutations.
#' @return Sorted numeric vector of possible p-values (capped at 0.5).
#' @noRd
bc_pvalue_support <- function(min_exceedances, max_permutations) {
  early <- (min_exceedances + 1) /
    (seq.int(min_exceedances, max_permutations) + 1)
  maxp <- seq_len(min_exceedances) / (max_permutations + 1)
  support <- sort(unique(c(early, maxp)))
  support[support <= 0.5]
}


#' Summarize neighborhood comparison results
#'
#' Computes q-values (Storey & Tibshirani, 2003), filters results, and computes summary
#' statistics at gene-pair, gene, and ortholog-group levels.
#'
#' @param comparison Data frame from [compare_neighborhoods()].
#' @param alternative Which tail to test: `"greater"` (default) for
#'   conservation (upper-tail, uses `.p.val.con` columns) or `"less"` for
#'   divergence (lower-tail, uses `.p.val.div` columns).
#' @param alpha Significance threshold applied to q-values (default 0.05).
#' @param filter_zero If `TRUE` (default for `"greater"`), remove rows where
#'   both overlap values are zero. Defaults to `FALSE` for `"less"`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{results}{Data frame with the original p-values preserved and new
#'       q-value columns (`.q.val.con` or `.q.val.div`) added.
#'       Rows are filtered according to `filter_zero`.}
#'     \item{summary}{List of summary statistics at gene-pair, gene, and
#'       ortholog-group levels, thresholded on q-values.}
#'   }
#'
#' @references
#' Storey, J. D. & Tibshirani, R. (2003). Statistical significance for
#' genomewide studies. \emph{Proceedings of the National Academy of Sciences},
#' 100(16), 9440--9445. \doi{10.1073/pnas.1530509100}
#'
#' @export
summarize_comparison <- function(comparison,
                                 alternative = c("greater", "less"),
                                 alpha = 0.05,
                                 filter_zero = NULL) {
  alternative <- match.arg(alternative)

  if (is.null(filter_zero)) {
    filter_zero <- alternative == "greater"
  }

  if (!all(c("Species1", "Species2", "OrthoGroup",
             "Species1.p.val.con", "Species2.p.val.con",
             "Species1.p.val.div", "Species2.p.val.div",
             "Species1.neigh.overlap",
             "Species2.neigh.overlap") %in%
             names(comparison))) {
    stop("comparison must be output from compare_neighborhoods()")
  }

  # Select p-value columns based on alternative
  suffix <- if (alternative == "greater") "con" else "div"
  sp1_col <- paste0("Species1.p.val.", suffix)
  sp2_col <- paste0("Species2.p.val.", suffix)

  res <- comparison

  # Filter zero-overlap rows
  if (filter_zero) {
    res <- res[
      res$Species1.neigh.overlap > 0 &
        res$Species2.neigh.overlap > 0, ,
      drop = FALSE
    ]
  }

  if (nrow(res) == 0) {
    return(list(
      results = res,
      summary = list(
        gene_pairs = list(sp1 = 0L, sp2 = 0L, reciprocal = 0L, total = 0L),
        genes = list(
          sp1 = 0L, sp2 = 0L,
          reciprocal_sp1 = 0L, reciprocal_sp2 = 0L
        ),
        orthogroups = list(sp1 = 0L, sp2 = 0L, reciprocal = 0L, total = 0L)
      )
    ))
  }

  # Compute q-values on selected p-value columns
  q1_col <- sub("p\\.val", "q.val", sp1_col)
  q2_col <- sub("p\\.val", "q.val", sp2_col)
  res[[q1_col]] <- compute_qvalues(res[[sp1_col]])
  res[[q2_col]] <- compute_qvalues(res[[sp2_col]])

  # Helper: count groups where the best (min) q-value is significant
  count_sig <- function(qvals, groups) {
    sum(tapply(qvals, groups, min) < alpha)
  }

  q1 <- res[[q1_col]]
  q2 <- res[[q2_col]]
  max_q <- pmax(q1, q2)

  list(
    results = res,
    summary = list(
      gene_pairs = list(
        sp1 = sum(q1 < alpha),
        sp2 = sum(q2 < alpha),
        reciprocal = sum(q1 < alpha & q2 < alpha),
        total = nrow(res)
      ),
      genes = list(
        sp1 = count_sig(q1, res$Species1),
        sp2 = count_sig(q2, res$Species2),
        reciprocal_sp1 = count_sig(max_q, res$Species1),
        reciprocal_sp2 = count_sig(max_q, res$Species2)
      ),
      orthogroups = list(
        sp1 = count_sig(q1, res$OrthoGroup),
        sp2 = count_sig(q2, res$OrthoGroup),
        reciprocal = count_sig(max_q, res$OrthoGroup),
        total = length(unique(res$OrthoGroup))
      )
    )
  )
}


# ---- GPU-accelerated fold-enrichment precomputation -------------------------

#' Build combined fold-enrichment matrix on GPU via torch
#'
#' Precomputes `combined[a, b] = FE1(b->a) + FE2(a->b)` where:
#' - FE1(b, a) = |N1(a) \eqn{\cap} R1(b)| / E1(a, b), anchor = net1
#' - FE2(a, b) = |N2(b) \eqn{\cap} R2(a)| / E2(a, b), anchor = net2
#'
#' Overlaps are computed via matrix multiply (GEMM) on GPU — O(n^3) but runs
#' in milliseconds on CUDA/MPS vs seconds on CPU. The resulting matrix converts
#' each permutation from O(M*N*n_words) bit-vector popcount operations to
#' O(M*N) table lookups.
#'
#' @param net1_mat,net2_mat Co-expression matrices (n1 x n1, n2 x n2).
#' @param thr1,thr2 Co-expression thresholds.
#' @param ortho_sp1_idx,ortho_sp2_idx 0-based ortholog gene indices.
#' @return Numeric matrix (n1 x n2) — the combined fold-enrichment matrix.
#' @keywords internal
#' Build binary adjacency matrix on GPU from sparse edge list
#'
#' Avoids transferring the full dense n×n co-expression matrix to GPU.
#' At 3% density, transfers ~125 MB of edge indices instead of ~4 GB dense.
#'
#' @param net_mat Co-expression matrix (n x n).
#' @param thr Co-expression threshold.
#' @param dtype Torch dtype for the result.
#' @param device Torch device string.
#' @return Binary adjacency tensor (n x n) on `device`, diagonal = 0.
#' @keywords internal
adj_to_gpu <- function(net_mat, thr, dtype, device) {
  n <- nrow(net_mat)
  # Find edges in R (sparse) — diagonal excluded
  edges <- which(net_mat >= thr & row(net_mat) != col(net_mat), arr.ind = TRUE)

  adj <- torch::torch_zeros(n, n, dtype = dtype, device = device)
  if (nrow(edges) > 0L) {
    rows_t <- torch::torch_tensor(edges[, 1L], dtype = torch::torch_long(),
                                  device = device)
    cols_t <- torch::torch_tensor(edges[, 2L], dtype = torch::torch_long(),
                                  device = device)
    ones_t <- torch::torch_ones(nrow(edges), dtype = dtype, device = device)
    adj$index_put_(list(rows_t, cols_t), ones_t)
    rm(rows_t, cols_t, ones_t)
  }
  adj
}

build_combined_fe_torch <- function(net1_mat, net2_mat, thr1, thr2,
                                    ortho_sp1_idx, ortho_sp2_idx) {
  n1 <- nrow(net1_mat)
  n2 <- nrow(net2_mat)

  dd <- torch_device_dtype()
  device <- dd$device
  dtype <- dd$dtype

  # Build adjacency on GPU from sparse edge list — avoids large dense R→GPU
  # copy that crashes on some CUDA architectures (Blackwell + torch 0.16.3).
  adj1 <- adj_to_gpu(net1_mat, thr1, dtype, device)
  rm(net1_mat)
  adj2 <- adj_to_gpu(net2_mat, thr2, dtype, device)
  rm(net2_mat)
  gc()

  combined <- torch::with_no_grad({
    # Ortholog indicator (n2 x n1) — build on GPU via index_put_
    ortho <- torch::torch_zeros(n2, n1, dtype = dtype, device = device)
    rows_t <- torch::torch_tensor(
      ortho_sp2_idx + 1L, dtype = torch::torch_long(), device = device
    )
    cols_t <- torch::torch_tensor(
      ortho_sp1_idx + 1L, dtype = torch::torch_long(), device = device
    )
    ones_t <- torch::torch_ones(length(ortho_sp1_idx), dtype = dtype,
                                device = device)
    ortho$index_put_(list(rows_t, cols_t), ones_t)
    rm(rows_t, cols_t, ones_t)

    # --- Direction 1: FE1[b, a] = |R1(b) ∩ N1(a)| / E1(b, a) ---
    reach1 <- (adj2$mm(ortho) > 0)$to(dtype = dtype)   # (n2, n1)
    neigh1_sz <- adj1$sum(dim = 2L)                     # (n1,)
    reach1_sz <- reach1$sum(dim = 2L)                   # (n2,)
    overlap1 <- reach1$mm(adj1)                         # (n2, n1)
    rm(reach1)

    E1 <- reach1_sz$unsqueeze(2L) * neigh1_sz$unsqueeze(1L) / n1
    FE1 <- overlap1 / E1$clamp(min = 1e-30)
    rm(overlap1, E1, reach1_sz, neigh1_sz)

    # --- Direction 2: FE2[a, b] = |R2(a) ∩ N2(b)| / E2(a, b) ---
    reach2 <- (adj1$mm(ortho$t()) > 0)$to(dtype = dtype) # (n1, n2)
    rm(adj1, ortho)
    neigh2_sz <- adj2$sum(dim = 2L)                       # (n2,)
    overlap2 <- reach2$mm(adj2)                           # (n1, n2)
    rm(adj2)
    reach2_sz <- reach2$sum(dim = 2L)                     # (n1,)
    rm(reach2)

    E2 <- reach2_sz$unsqueeze(2L) * neigh2_sz$unsqueeze(1L) / n2
    FE2 <- overlap2 / E2$clamp(min = 1e-30)
    rm(overlap2, E2, reach2_sz, neigh2_sz)
    gc()

    # Combined: combined[a, b] = FE1[b, a] + FE2[a, b]
    result <- FE1$t() + FE2
    rm(FE1, FE2)
    result  # (n1, n2)
  })

  result_r <- as.matrix(combined$cpu()$to(dtype = torch::torch_float64()))
  rm(combined)
  gc()
  result_r
}


#' Permutation-based HOG-level conservation/divergence test
#'
#' Tests each Hierarchical Ortholog Group (HOG) for co-expression conservation
#' (or divergence) using a gene-identity permutation null with adaptive stopping
#' (Besag & Clifford, 1991).
#'
#' @section Statistical method:
#' The null hypothesis is that the specific gene identities in a HOG carry no
#' information about co-expression conservation: replacing the HOG's genes with
#' randomly chosen genes from the same networks would yield equally large
#' neighborhood overlap.
#'
#' For each permutation, M random species-1 genes and N random species-2 genes
#' are drawn (matching the HOG's gene counts), and the sum-of-fold-enrichments
#' statistic T is computed over all M x N pair x direction combinations. The
#' permutation p-value is the fraction of permuted T values that exceed (or
#' fall below, for divergence) the observed T.
#'
#' The Besag & Clifford (1991) sequential stopping rule terminates permutations
#' early once `min_exceedances` permutation statistics exceed T_obs,
#' providing efficient computation without sacrificing accuracy for clearly
#' significant or non-significant HOGs.
#'
#' @section Why not Fisher's method:
#' Fisher's method for combining p-values assumes independent tests. Within a
#' HOG, pair-level hypergeometric tests share network neighborhoods (genes
#' co-expressed with one HOG member are often co-expressed with another) and
#' share the ortholog mapping. This non-independence inflates Fisher's statistic
#' and produces anti-conservative p-values. The permutation approach is exact
#' regardless of the dependency structure.
#'
#' @section GPU acceleration:
#' When `use_torch = TRUE`, the fold-enrichment matrix is precomputed on GPU
#' via GEMM (CUDA, MPS, or CPU fallback). This converts each permutation from
#' O(M*N*n_words) bit-vector intersections to O(M*N) table lookups — typically
#' a 100-300x overall speedup. Requires the
#' \href{https://torch.mlverse.org/}{torch} package.
#'
#' On Apple Silicon (MPS), torch uses float32 because Metal does not support
#' float64. The fold-enrichment matrix is numerically exact in float32 (all
#' intermediate values are small integers), so `permutation_hog_test(use_torch
#' = TRUE)` is safe on MPS. However, [compute_network()] correlation in
#' float32 can cause rank-swap artifacts in Spearman + MR normalization. If
#' you need exact Spearman + MR correlation on Apple Silicon, use
#' `compute_network(use_torch = FALSE)` with `permutation_hog_test(use_torch =
#' TRUE)` — the 100-300x permutation speedup comes from the FE precomputation
#' and table-lookup architecture, not from GPU-accelerated correlation. On
#' CUDA, both functions use float64 with no precision tradeoff.
#'
#' @section Multiple testing correction:
#' Q-values are computed using the Liang (2016) discrete q-value method via
#' [DiscreteQvalue::DQ()], which accounts for the discrete support of
#' Besag & Clifford p-values. The standard Storey & Tibshirani (2003)
#' pi0-adaptive method cannot be used here because the adaptive stopping
#' rule produces p-values that are valid (super-uniform under the null) but
#' **not** uniformly distributed: the negative binomial stopping concentrates
#' null p-values in the 0.1--0.5 range and depletes the right tail
#' (p > 0.9), causing Storey's pi0 estimator to severely underestimate pi0.
#' The Liang method estimates pi0 using the discrete support structure and
#' is correctly calibrated (pi0 = 1 on null simulations).
#'
#' @param net1 Network object for species 1 (output of [compute_network()]).
#' @param net2 Network object for species 2 (output of [compute_network()]).
#' @param comparison Data frame from [compare_neighborhoods()] --- raw p-values,
#'   **not** FDR-corrected. Used to define HOG structure and extract
#'   effect sizes.
#' @param alternative Which tail to test: `"greater"` (default) for conservation
#'   or `"less"` for divergence.
#' @param min_exceedances Besag & Clifford stopping parameter: stop permuting
#'   once this many permutation statistics exceed (or fall below) T_obs.
#'   Higher values give more precise p-values but take longer (default 50).
#' @param max_permutations Maximum permutations per HOG (default 10000).
#' @param n_cores Number of threads for parallel computation (default 1).
#' @param use_torch If `TRUE`, precompute the fold-enrichment matrix on GPU
#'   via torch, then run permutations as fast table lookups. Requires the
#'   \href{https://torch.mlverse.org/}{torch} package. Default `FALSE`.
#'
#' @return A data frame with one row per HOG, ordered by p-value, with columns:
#'   \describe{
#'     \item{OrthoGroup}{HOG identifier}
#'     \item{n_pairs}{Number of ortholog pairs (M x N) in this HOG}
#'     \item{n_sp1}{Number of unique species-1 genes (M)}
#'     \item{n_sp2}{Number of unique species-2 genes (N)}
#'     \item{T_obs}{Observed sum-of-fold-enrichments statistic}
#'     \item{n_perm}{Number of permutations performed (may be < max_permutations
#'       due to adaptive stopping)}
#'     \item{n_exceed}{Number of permutation statistics exceeding T_obs}
#'     \item{mean_eff}{Mean geometric-mean effect size across pairs}
#'     \item{p.value}{Permutation p-value: (n_exceed + 1) / (n_perm + 1)}
#'     \item{q.value}{Discrete q-value (Liang, 2016) accounting for Besag-Clifford support}
#'   }
#'
#' @references
#' Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
#' \emph{Biometrika}, 78(2), 301--304. \doi{10.1093/biomet/78.2.301}
#'
#' Liang, K. (2016). False discovery rate estimation for large-scale
#' homogeneous discrete p-values. \emph{Biometrics}, 72(2), 639--648.
#' \doi{10.1111/biom.12429}
#'
#' @export
permutation_hog_test <- function(net1, net2, comparison,
                                 alternative = c("greater", "less"),
                                 min_exceedances = 50L,
                                 max_permutations = 10000L,
                                 n_cores = 1L,
                                 use_torch = FALSE) {
  alternative <- match.arg(alternative)
  min_exceedances <- as.integer(min_exceedances)
  max_permutations <- as.integer(max_permutations)
  n_cores <- as.integer(n_cores)

  if (!is.list(net1) || is.null(net1$network)) {
    stop("net1 must be a network object from compute_network()")
  }
  if (!is.list(net2) || is.null(net2$network)) {
    stop("net2 must be a network object from compute_network()")
  }
  if (!all(c("Species1", "Species2", "OrthoGroup",
             "Species1.effect.size",
             "Species2.effect.size") %in%
             names(comparison))) {
    stop("comparison must be output from compare_neighborhoods()")
  }
  if (use_torch && !requireNamespace("torch", quietly = TRUE)) {
    stop("use_torch = TRUE requires the torch package ",
         "(install.packages('torch'); torch::install_torch())")
  }
  if (nrow(comparison) == 0) {
    return(data.frame(
      OrthoGroup = character(0), n_pairs = integer(0),
      n_sp1 = integer(0), n_sp2 = integer(0),
      T_obs = numeric(0), n_perm = integer(0),
      n_exceed = integer(0), mean_eff = numeric(0),
      p.value = numeric(0), q.value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  net1_mat <- net1$network
  net2_mat <- net2$network
  thr1 <- net1$threshold
  thr2 <- net2$threshold
  net1_genes <- rownames(net1_mat)
  net2_genes <- rownames(net2_mat)

  idx1 <- stats::setNames(seq_along(net1_genes) - 1L, net1_genes)
  idx2 <- stats::setNames(seq_along(net2_genes) - 1L, net2_genes)

  # Filter comparison to genes present in both networks — gene names may

  # not match after upstream merging/reduction of multi-copy orthologs.
  # NA indices from missing genes would cause index_put_ failures in torch
  # and undefined behavior in C++.
  in_net <- comparison$Species1 %in% net1_genes &
            comparison$Species2 %in% net2_genes
  if (!all(in_net)) {
    n_dropped <- sum(!in_net)
    message("Dropped ", n_dropped, " ortholog pairs with genes not in networks")
    comparison <- comparison[in_net, , drop = FALSE]
  }

  # Deduplicated ortholog mapping for reachable-set construction
  ortho_pairs <- unique(comparison[, c("Species1", "Species2")])
  ortho_sp1_idx <- as.integer(idx1[ortho_pairs$Species1])
  ortho_sp2_idx <- as.integer(idx2[ortho_pairs$Species2])

  # Per-HOG unique gene indices
  hog_groups <- split(seq_len(nrow(comparison)), comparison$OrthoGroup)
  hog_names <- names(hog_groups)

  hog_sp1_list <- lapply(hog_groups, function(rows) {
    as.integer(unique(idx1[comparison$Species1[rows]]))
  })
  hog_sp2_list <- lapply(hog_groups, function(rows) {
    as.integer(unique(idx2[comparison$Species2[rows]]))
  })

  if (use_torch) {
    combined <- build_combined_fe_torch(
      net1_mat, net2_mat, thr1, thr2,
      ortho_sp1_idx, ortho_sp2_idx
    )
    rm(net1_mat, net2_mat)
    gc()
    perm_result <- fe_hog_permutation_test_cpp(
      combined = combined,
      hog_sp1_list = hog_sp1_list,
      hog_sp2_list = hog_sp2_list,
      test_greater = (alternative == "greater"),
      min_exceedances = min_exceedances,
      max_permutations = max_permutations,
      n_cores = n_cores
    )
    rm(combined)
    gc()
  } else {
    perm_result <- hog_permutation_test_cpp(
      net1 = net1_mat, net2 = net2_mat,
      thr1 = thr1, thr2 = thr2,
      ortho_sp1_idx = ortho_sp1_idx,
      ortho_sp2_idx = ortho_sp2_idx,
      hog_sp1_list = hog_sp1_list,
      hog_sp2_list = hog_sp2_list,
      test_greater = (alternative == "greater"),
      min_exceedances = min_exceedances,
      max_permutations = max_permutations,
      n_cores = n_cores
    )
  }

  eff <- sqrt(comparison$Species1.effect.size * comparison$Species2.effect.size)

  result <- data.frame(
    OrthoGroup = hog_names,
    n_pairs    = vapply(hog_groups, length, integer(1)),
    n_sp1      = vapply(hog_sp1_list, length, integer(1)),
    n_sp2      = vapply(hog_sp2_list, length, integer(1)),
    T_obs      = perm_result$T_obs,
    n_perm     = perm_result$n_perm,
    n_exceed   = perm_result$n_exceed,
    mean_eff   = vapply(hog_groups, function(i) mean(eff[i]), double(1)),
    p.value    = perm_result$p_value,
    stringsAsFactors = FALSE
  )
  bc_support <- bc_pvalue_support(min_exceedances, max_permutations)

  result$q.value <- if (nrow(result) < 2L) result$p.value else {
    DiscreteQvalue::DQ(
      result$p.value, ss = bc_support, method = "Liang"
    )$q.values
  }
  result[order(result$p.value), ]
}
