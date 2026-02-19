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
  sp1_col <- if (alternative == "greater") {
    "Species1.p.val.con"
  } else {
    "Species1.p.val.div"
  }
  sp2_col <- if (alternative == "greater") {
    "Species2.p.val.con"
  } else {
    "Species2.p.val.div"
  }

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
#' @section Intersection modes:
#' For networks with max(n1, n2) <= 100,000 genes, bit-vector intersection
#' with popcount is used for maximum throughput. For larger networks, a
#' flag-vector approach (sparse set/count/clear) avoids excessive memory use.
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
                                 n_cores = 1L) {
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
  # Construct Besag-Clifford p-value support for discrete q-values (Liang 2016)
  # Early-stopped HOGs: p = (min_exceedances+1)/(n+1), n = min_exceedances..max_permutations
  early_support <- (min_exceedances + 1) /
    (seq.int(min_exceedances, max_permutations) + 1)
  # Max-perm HOGs: p = (j+1)/(max_permutations+1), j = 0..(min_exceedances-1)
  maxp_support <- seq_len(min_exceedances) / (max_permutations + 1)
  bc_support <- sort(unique(c(early_support, maxp_support)))
  bc_support <- bc_support[bc_support <= 0.5]

  if (nrow(result) < 2L) {
    result$q.value <- result$p.value
  } else {
    result$q.value <- DiscreteQvalue::DQ(
      result$p.value, ss = bc_support, method = "Liang"
    )$q.values
  }
  result[order(result$p.value), ]
}
