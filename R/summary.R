#' Summarize neighborhood comparison results
#'
#' Applies FDR correction, filters results, and computes summary statistics
#' at gene-pair, gene, and ortholog-group levels.
#'
#' @param comparison Data frame from [compare_neighborhoods()].
#' @param alternative Which tail to test: `"greater"` (default) for
#'   conservation (upper-tail, uses `.p.val.con` columns) or `"less"` for
#'   divergence (lower-tail, uses `.p.val.div` columns).
#' @param fdr_method P-value adjustment method passed to [p.adjust()]
#'   (default `"fdr"`).
#' @param alpha Significance threshold (default 0.05).
#' @param filter_zero If `TRUE` (default for `"greater"`), remove rows where
#'   both overlap values are zero. Defaults to `FALSE` for `"less"`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{results}{Data frame with FDR-corrected p-values and filtered rows.
#'       FDR correction is applied to the p-value columns selected by
#'       `alternative`.}
#'     \item{summary}{List of summary statistics at gene-pair, gene, and
#'       ortholog-group levels.}
#'   }
#'
#' @export
summarize_comparison <- function(comparison,
                                 alternative = c("greater", "less"),
                                 fdr_method = "fdr",
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

  # FDR correction on selected p-value columns
  res[[sp1_col]] <- stats::p.adjust(res[[sp1_col]], method = fdr_method)
  res[[sp2_col]] <- stats::p.adjust(res[[sp2_col]], method = fdr_method)

  # Helper: count groups where the best (min) p-value is significant
  count_sig <- function(pvals, groups) {
    sum(tapply(pvals, groups, min) < alpha)
  }

  p1 <- res[[sp1_col]]
  p2 <- res[[sp2_col]]
  max_p <- pmax(p1, p2)

  list(
    results = res,
    summary = list(
      gene_pairs = list(
        sp1 = sum(p1 < alpha),
        sp2 = sum(p2 < alpha),
        reciprocal = sum(p1 < alpha & p2 < alpha),
        total = nrow(res)
      ),
      genes = list(
        sp1 = count_sig(p1, res$Species1),
        sp2 = count_sig(p2, res$Species2),
        reciprocal_sp1 = count_sig(max_p, res$Species1),
        reciprocal_sp2 = count_sig(max_p, res$Species2)
      ),
      orthogroups = list(
        sp1 = count_sig(p1, res$OrthoGroup),
        sp2 = count_sig(p2, res$OrthoGroup),
        reciprocal = count_sig(max_p, res$OrthoGroup),
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
#' @param fdr_method P-value adjustment method passed to [p.adjust()]
#'   (default `"fdr"`).
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
#'     \item{q.value}{FDR-adjusted p-value}
#'   }
#'
#' @references
#' Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
#' \emph{Biometrika}, 78(2), 301--304. \doi{10.1093/biomet/78.2.301}
#'
#' @export
permutation_hog_test <- function(net1, net2, comparison,
                                 alternative = c("greater", "less"),
                                 min_exceedances = 50L,
                                 max_permutations = 10000L,
                                 n_cores = 1L,
                                 fdr_method = "fdr") {
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
  result$q.value <- stats::p.adjust(result$p.value, method = fdr_method)
  result[order(result$p.value), ]
}
