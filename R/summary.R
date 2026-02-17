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
             "Species1.neigh.overlap", "Species2.neigh.overlap") %in%
           names(comparison))) {
    stop("comparison must be output from compare_neighborhoods()")
  }

  # Select p-value columns based on alternative
  sp1_col <- if (alternative == "greater") "Species1.p.val.con" else "Species1.p.val.div"
  sp2_col <- if (alternative == "greater") "Species2.p.val.con" else "Species2.p.val.div"

  res <- comparison

  # Filter zero-overlap rows
  if (filter_zero) {
    res <- res[res$Species1.neigh.overlap > 0 & res$Species2.neigh.overlap > 0, ,
               drop = FALSE]
  }

  if (nrow(res) == 0) {
    return(list(
      results = res,
      summary = list(
        gene_pairs = list(sp1 = 0L, sp2 = 0L, reciprocal = 0L, total = 0L),
        genes = list(sp1 = 0L, sp2 = 0L, reciprocal_sp1 = 0L, reciprocal_sp2 = 0L),
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


#' Summarize neighborhood comparison at the HOG level using Fisher's method
#'
#' Combines per-pair p-values within each Hierarchical Ortholog Group (HOG)
#' using Fisher's method for combining independent tests. This replaces
#' gene-level significance checks with a proper aggregate test that accounts
#' for both the number and strength of significant pairs within a HOG.
#'
#' For each HOG, the combined p-value per pair is `pmax(Species1.p, Species2.p)`
#' (reciprocal criterion), then Fisher's statistic is computed as
#' `T = -2 * sum(log(p_i))` which follows a chi-squared distribution with
#' `df = 2k` under the null (k = number of pairs).
#'
#' For single-pair HOGs (1:1 orthologs), Fisher's test reduces exactly to the
#' original pair-level test.
#'
#' @param comparison Data frame from [compare_neighborhoods()] — raw p-values,
#'   **not** FDR-corrected. FDR correction is applied to the HOG-level p-values.
#' @param alternative Which tail to test: `"greater"` for conservation
#'   (uses `.p.val.con` columns) or `"less"` for divergence
#'   (uses `.p.val.div` columns).
#' @param fdr_method P-value adjustment method passed to [p.adjust()]
#'   (default `"fdr"`).
#'
#' @return A data frame with one row per HOG, ordered by p-value, with columns:
#'   \describe{
#'     \item{OrthoGroup}{HOG identifier}
#'     \item{n_pairs}{Number of ortholog pairs in this HOG}
#'     \item{fisher_stat}{Fisher's combined statistic (-2 * sum(log(p)))}
#'     \item{mean_eff}{Mean geometric-mean effect size across pairs}
#'     \item{min_eff}{Minimum geometric-mean effect size across pairs}
#'     \item{df}{Degrees of freedom (2 * n_pairs)}
#'     \item{p.value}{P-value from chi-squared test}
#'     \item{q.value}{FDR-adjusted p-value}
#'   }
#'
#' @export
summarize_hog_comparison <- function(comparison,
                                     alternative = c("greater", "less"),
                                     fdr_method = "fdr") {
  alternative <- match.arg(alternative)

  if (!all(c("Species1", "Species2", "OrthoGroup",
             "Species1.p.val.con", "Species2.p.val.con",
             "Species1.p.val.div", "Species2.p.val.div",
             "Species1.effect.size", "Species2.effect.size") %in%
           names(comparison))) {
    stop("comparison must be output from compare_neighborhoods()")
  }

  if (nrow(comparison) == 0) {
    return(data.frame(
      OrthoGroup  = character(0),
      n_pairs     = integer(0),
      fisher_stat = numeric(0),
      mean_eff    = numeric(0),
      min_eff     = numeric(0),
      df          = integer(0),
      p.value     = numeric(0),
      q.value     = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  sp1_col <- if (alternative == "greater") "Species1.p.val.con" else "Species1.p.val.div"
  sp2_col <- if (alternative == "greater") "Species2.p.val.con" else "Species2.p.val.div"

  # Combined p-value per pair: max (reciprocal)
  pair_p <- pmax(comparison[[sp1_col]], comparison[[sp2_col]])
  pair_p <- pmax(pair_p, .Machine$double.xmin)  # avoid log(0)

  # Effect sizes per pair (geometric mean of both directions)
  eff <- sqrt(comparison$Species1.effect.size * comparison$Species2.effect.size)

  # Split by HOG
  hogs <- split(seq_len(nrow(comparison)), comparison$OrthoGroup)

  result <- data.frame(
    OrthoGroup  = names(hogs),
    n_pairs     = vapply(hogs, length, integer(1)),
    fisher_stat = vapply(hogs, function(i) -2 * sum(log(pair_p[i])), double(1)),
    mean_eff    = vapply(hogs, function(i) mean(eff[i]), double(1)),
    min_eff     = vapply(hogs, function(i) min(eff[i]), double(1)),
    stringsAsFactors = FALSE
  )
  result$df      <- 2L * result$n_pairs
  result$p.value <- stats::pchisq(result$fisher_stat, df = result$df, lower.tail = FALSE)
  result$q.value <- stats::p.adjust(result$p.value, method = fdr_method)
  result[order(result$p.value), ]
}
