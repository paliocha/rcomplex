# Cross-species module preservation.
#
# Reference-species modules are projected onto a test species through a
# paralog-resolved ortholog map, then scored on how well their topology
# survives. The inferential call comes from permutation p-values; a
# Zsummary-style composite is reported alongside for continuity with the
# WGCNA literature.
#
# Two statistics carry the call, the two NetRep computes when only an
# adjacency matrix is available (Ritchie et al. 2016):
#
#   avg.weight  sum(kIM) / (m^2 - m)  -- module density (density half)
#   cor.degree  cor(kIM_ref, kIM_test) -- intramodular degree concordance
#
# meanClusterCoeff and meanMAR are reported as diagnostics but excluded from
# the call: both are weighted means of the surviving edge weights, and a
# hard-thresholded network leaves those weights nearly constant (a max/min
# ratio of ~1.04 for MR at density 0.03), so they return almost the same value
# for any gene set. Including them in a median-of-three collapses a density
# signal of Z = 124 to Z = 5.3 and misclassifies a perfectly preserved module.
#
# cor.degree is Pearson, matching NetRep's implementation (arma::cor on the
# raw degree vectors; the "relative rank" wording in its documentation
# describes the intent, not the computation).


# Statistic columns returned by the C++ kernel, in order.
.PRES_STATS <- c(
  "avg.weight", "meanClusterCoeff", "meanMAR",
  "cor.degree", "cor.clusterCoeff", "cor.MAR"
)

# The two that carry the call: indices into .PRES_STATS.
.PRES_DENSITY <- 1L
.PRES_CONNECTIVITY <- 4L


#' Test whether co-expression modules are preserved across species
#'
#' Projects the modules of a reference species onto a test species through an
#' ortholog map and tests, per module, whether the module's topology survives
#' in the test species' network. Unlike gene-overlap tests, a module that keeps
#' its gene membership but loses its internal wiring is correctly reported as
#' diverged.
#'
#' @section Statistics:
#' Two statistics carry the call, the pair NetRep computes when only an
#' adjacency matrix is available:
#' \describe{
#'   \item{avg.weight}{`sum(kIM) / (m^2 - m)`, the module's mean within-module
#'     adjacency ("module density"). On a binary network this is the
#'     proportion of realised to possible within-module edges.}
#'   \item{cor.degree}{Pearson correlation between each gene's intramodular
#'     connectivity in the reference network and in the test network. Measures
#'     whether hub identity is conserved.}
#' }
#' `meanClusterCoeff` and `meanMAR` are reported in `observed` as diagnostics
#' but take no part in the call; on a hard-thresholded network the surviving
#' edge weights are nearly constant, which leaves both statistics with almost
#' no dynamic range.
#'
#' @section Null model:
#' Gene identities are shuffled while edges are held constant, and each module
#' is handed a contiguous block of the shuffled genes of its own size. Only
#' ortholog-mappable test-species genes enter the shuffle, matching NetRep's
#' `"overlap"` null model. P-values are one-sided
#' (`(exceedances + 1) / (permutations + 1)`) and combined across the two
#' statistics with `pmax`, so a module is called preserved only when both are
#' significant -- the same reciprocal criterion as
#' `pval_combine = "max"` elsewhere in the package.
#'
#' @section Paralog resolution:
#' Multi-copy HOGs are reduced to one counterpart per gene by
#' [resolve_ortholog_map()] when `map` is not supplied. Resolution only
#' chooses which paralog copy carries a module label; it never changes which
#' genes are mappable, so the tested gene set does not depend on the
#' conservation evidence used to resolve it.
#'
#' @param modules_ref Module detection result for the reference species
#'   (output of [detect_modules()]).
#' @param net_ref,net_test Network objects from [compute_network()] for the
#'   reference and test species.
#' @param orthologs Data frame with columns `Species1`, `Species2`, `hog`.
#'   Ignored when `map` is supplied.
#' @param map Optional pre-built ortholog map from [resolve_ortholog_map()].
#' @param edges,cliques Optional [find_coexpressologs()] and [find_cliques()]
#'   results, passed to [resolve_ortholog_map()] for paralog resolution.
#' @param sp_ref,sp_test Species labels, required only when `edges` or
#'   `cliques` is supplied.
#' @param n_perm Number of permutations (default 10000). The smallest
#'   attainable p-value is `1 / (n_perm + 1)`, so this sets the floor on how
#'   significant any module can be: at 1000 permutations every strongly
#'   preserved module ties at `p = 0.000999` and cannot be ranked by q-value.
#'   The permutation loop costs O(edges) per iteration -- roughly 0.7 s per
#'   1000 permutations on a 6000-gene network at density 0.03 -- so the
#'   default buys a floor of 1e-4 cheaply. Raise it further when many modules
#'   sit at the floor.
#' @param min_module_size Minimum number of mapped test-species genes for a
#'   module to be tested (default 10).
#' @param binary Treat every surviving edge as weight 1 (default `FALSE`).
#'   On a hard-thresholded network this changes little and makes `avg.weight`
#'   exactly the module edge density.
#' @param alpha Significance threshold (default 0.05).
#' @param qvalue_method Multiple-testing correction. `"liang"` (default) uses
#'   [DiscreteQvalue::DQ()] with the exact discrete support of a permutation
#'   p-value, `{1/(n_perm+1), ..., 1}` -- the same treatment
#'   [permutation_hog_test()] gives its permutation p-values. `"bh"` forces
#'   Benjamini-Hochberg. Liang falls back to BH automatically when there are
#'   too few modules to estimate pi0 (fewer than 10) or if the estimator
#'   fails.
#' @param n_cores Number of OpenMP threads (default 1).
#' @param seed Optional RNG seed. Results are independent of `n_cores`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{preservation}{One row per tested module: `module`, `size`,
#'       `size_mapped`, the two headline statistics, their permutation
#'       p-values, the combined `p.value` and `q.value`, `Z.avg.weight`,
#'       `Z.cor.degree`, `Zsummary` and `medianRank`.}
#'     \item{observed}{All six statistics per module, with permutation means
#'       and standard deviations.}
#'     \item{map}{The ortholog map used.}
#'     \item{params}{Call parameters.}
#'   }
#'
#' @references
#' Ritchie, S. C. et al. (2016). A scalable permutation approach reveals
#' replication and preservation patterns of network modules in large datasets.
#' *Cell Systems*, 3(1), 71-82. \doi{10.1016/j.cels.2016.06.012}
#'
#' Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. (2011). Is my network
#' module preserved and reproducible? *PLoS Computational Biology*, 7(1),
#' e1001057. \doi{10.1371/journal.pcbi.1001057}
#'
#' @examples
#' \dontrun{
#' pres <- module_preservation(mods_A, net_A, net_B, ortho,
#'   edges = rcx$edges, cliques = rcx$cliques,
#'   sp_ref = "SP_A", sp_test = "SP_B"
#' )
#' classify_preservation(pres)
#' }
#'
#' @export
module_preservation <- function(modules_ref, net_ref, net_test,
                                orthologs = NULL,
                                map = NULL, edges = NULL, cliques = NULL,
                                sp_ref = NULL, sp_test = NULL,
                                n_perm = 10000L, min_module_size = 10L,
                                binary = FALSE, alpha = 0.05,
                                qvalue_method = c("liang", "bh"),
                                n_cores = 1L, seed = NULL) {
  if (!is.list(modules_ref) || is.null(modules_ref$module_genes) ||
    is.null(modules_ref$modules)) {
    stop("modules_ref must be output from detect_modules()")
  }
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) stop("n_perm must be >= 1")
  qvalue_method <- match.arg(qvalue_method)
  min_module_size <- as.integer(min_module_size)
  if (is.na(min_module_size) || min_module_size < 3L) {
    stop(
      "min_module_size must be >= 3 (the degree correlation is ",
      "undefined below that)"
    )
  }

  mat_ref <- .net_check(net_ref, net_ref$threshold)
  mat_test <- .net_check(net_test, net_test$threshold)
  genes_ref <- rownames(mat_ref)
  genes_test <- rownames(mat_test)

  if (is.null(map)) {
    if (is.null(orthologs)) {
      stop("supply either 'orthologs' or a pre-built 'map'")
    }
    map <- resolve_ortholog_map(orthologs, genes_ref, genes_test,
      sp1 = sp_ref, sp2 = sp_test,
      edges = edges, cliques = cliques,
      alpha = alpha
    )
  }
  if (!is.data.frame(map) ||
    !all(c("gene1", "gene2", "source") %in% names(map))) {
    stop("map must be a data frame from resolve_ortholog_map()")
  }

  # ---- Project reference module labels onto test-species genes ----
  map <- map[map$gene1 %in% genes_ref & map$gene2 %in% genes_test, ,
    drop = FALSE
  ]

  # The permutation pool is every ortholog-mappable test gene, whether or not
  # its reference partner carries a module label -- NetRep's "overlap" null
  # model. Restricting the pool to labelled genes would make the null "genes
  # from the other modules", and since every module is dense that guts the
  # contrast avg.weight is supposed to measure.
  mappable_test <- unique(map$gene2)

  map$module <- as.character(modules_ref$modules[map$gene1])
  map_labelled <- map[!is.na(map$module), , drop = FALSE]
  if (nrow(map_labelled) == 0L) {
    stop("no reference module genes map to the test species")
  }

  proj <- .pres_project(map_labelled)
  if (is.null(proj) || nrow(proj) == 0L) {
    stop("no test-species gene received an unambiguous module label")
  }

  rows_by_mod <- split(seq_len(nrow(proj)), proj$module)
  sizes <- vapply(rows_by_mod, length, integer(1))
  tested <- names(sizes)[sizes >= min_module_size]
  if (length(tested) == 0L) {
    stop(
      "no module has at least min_module_size (", min_module_size,
      ") mapped test-species genes; largest is ", max(sizes)
    )
  }
  rows_by_mod <- rows_by_mod[tested]

  # ---- Local index spaces (ascending, as the C++ entry points require) ----
  keep_test <- sort(match(mappable_test, genes_test))
  loc_test <- match(match(proj$gene2, genes_test), keep_test) - 1L

  ri <- match(proj$gene1, genes_ref)
  keep_ref <- sort(unique(ri))
  loc_ref <- match(ri, keep_ref) - 1L

  # Order each module's genes by test index so the run is reproducible.
  rows_by_mod <- lapply(rows_by_mod, function(rows) rows[order(loc_test[rows])])

  test_members <- lapply(rows_by_mod, function(rows) as.integer(loc_test[rows]))
  ref_members <- lapply(
    rows_by_mod,
    function(rows) as.integer(unique(loc_ref[rows]))
  )

  # ---- Reference-network per-gene statistics (fixed across permutations) ----
  rstats <- .pres_gene_stats(
    net_ref, as.integer(keep_ref - 1L),
    ref_members, binary
  )

  # Index by test gene, so a reference gene mapped to several test paralogs
  # contributes its value once per paralog.
  ref_kIM <- lapply(rows_by_mod, function(rows) rstats$kIM[loc_ref[rows] + 1L])
  ref_cc <- lapply(rows_by_mod, function(rows) rstats$CC[loc_ref[rows] + 1L])
  ref_mar <- lapply(rows_by_mod, function(rows) rstats$MAR[loc_ref[rows] + 1L])

  flat <- vapply(ref_kIM, function(v) stats::var(v) == 0, logical(1))
  if (any(flat)) {
    warning(
      "intramodular connectivity is constant in the reference network ",
      "for module(s) ", paste(tested[flat], collapse = ", "),
      "; cor.degree is not meaningful there (NetRep documents this ",
      "for modules whose nodes are all connected with similar strength)"
    )
  }

  if (!is.null(seed)) set.seed(seed)
  res <- .pres_run(
    net_test, as.integer(keep_test - 1L), test_members,
    ref_kIM, ref_cc, ref_mar, n_perm, as.integer(n_cores),
    binary
  )

  .pres_assemble(
    res, modules_ref, tested, rows_by_mod, proj, map,
    n_perm, min_module_size, binary, alpha, qvalue_method, seed
  )
}


#' Assign one module label per test-species gene (internal)
#'
#' Resolved pairs win outright. A test gene that no resolved pair claims takes
#' the modal label of its unresolved partners, and is dropped on a tie.
#'
#' @noRd
.pres_project <- function(map) {
  res <- map[map$source != "unresolved", , drop = FALSE]
  unres <- map[map$source == "unresolved" &
    !(map$gene2 %in% res$gene2), , drop = FALSE]

  pick <- function(df) {
    if (nrow(df) == 0L) {
      return(NULL)
    }
    do.call(rbind, lapply(split(df, df$gene2), function(d) {
      tab <- table(d$module)
      top <- names(tab)[tab == max(tab)]
      if (length(top) != 1L) {
        return(NULL)
      } # tie: drop the gene
      d <- d[d$module == top, , drop = FALSE]
      d <- d[order(d$gene1), , drop = FALSE]
      data.frame(
        gene1 = d$gene1[1], gene2 = d$gene2[1],
        module = d$module[1], source = d$source[1],
        stringsAsFactors = FALSE
      )
    }))
  }

  out <- rbind(pick(res), pick(unres))
  if (!is.null(out)) rownames(out) <- NULL
  out
}


#' Per-gene intramodular statistics, dense/sparse dispatch (internal)
#' @noRd
.pres_gene_stats <- function(net, keep, members, binary) {
  a <- .net_cpp_args(net, net$threshold)
  if (.net_is_sparse(net)) {
    module_gene_stats_sparse_cpp(a$p, a$i, a$x, a$thr, keep, members, binary)
  } else {
    module_gene_stats_dense_cpp(a$net, a$thr, keep, members, binary)
  }
}


#' Preservation permutation engine, dense/sparse dispatch (internal)
#' @noRd
.pres_run <- function(net, keep, members, ref_kIM, ref_cc, ref_mar,
                      n_perm, n_cores, binary) {
  a <- .net_cpp_args(net, net$threshold)
  if (.net_is_sparse(net)) {
    module_preservation_sparse_cpp(
      a$p, a$i, a$x, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, binary
    )
  } else {
    module_preservation_dense_cpp(
      a$net, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, binary
    )
  }
}


#' Build the result tables from the kernel output (internal)
#' @noRd
.pres_assemble <- function(res, modules_ref, tested, rows_by_mod, proj, map,
                           n_perm, min_module_size, binary, alpha,
                           qvalue_method, seed) {
  colnames(res$observed) <- .PRES_STATS
  colnames(res$perm_mean) <- .PRES_STATS
  colnames(res$perm_sd) <- .PRES_STATS
  colnames(res$p_value) <- .PRES_STATS

  z <- (res$observed - res$perm_mean) / res$perm_sd
  colnames(z) <- .PRES_STATS

  d <- .PRES_DENSITY
  cc <- .PRES_CONNECTIVITY

  # Both statistics must be significant: the reciprocal criterion used by
  # pval_combine = "max" elsewhere in the package.
  p_comb <- pmax(res$p_value[, d], res$p_value[, cc])
  q_comb <- .pres_qvalues(p_comb, n_perm, qvalue_method)

  # medianRank: rank of the observed statistics across modules, 1 = strongest.
  rank_d <- rank(-res$observed[, d], na.last = "keep")
  rank_c <- rank(-res$observed[, cc], na.last = "keep")
  median_rank <- (rank_d + rank_c) / 2

  size_all <- vapply(modules_ref$module_genes, length, integer(1))

  preservation <- data.frame(
    module = tested,
    size = as.integer(size_all[tested]),
    size_mapped = vapply(rows_by_mod, length, integer(1)),
    avg.weight = res$observed[, d],
    cor.degree = res$observed[, cc],
    p.avg.weight = res$p_value[, d],
    p.cor.degree = res$p_value[, cc],
    p.value = p_comb,
    q.value = q_comb,
    Z.avg.weight = z[, d],
    Z.cor.degree = z[, cc],
    Zsummary = (z[, d] + z[, cc]) / 2,
    medianRank = median_rank,
    stringsAsFactors = FALSE
  )
  rownames(preservation) <- NULL

  observed <- data.frame(module = tested, stringsAsFactors = FALSE)
  for (s in .PRES_STATS) {
    observed[[s]] <- res$observed[, s]
    observed[[paste0("perm_mean.", s)]] <- res$perm_mean[, s]
    observed[[paste0("perm_sd.", s)]] <- res$perm_sd[, s]
  }
  rownames(observed) <- NULL

  list(
    preservation = preservation,
    observed = observed,
    projection = proj,
    map = map,
    params = list(
      n_perm = n_perm, min_module_size = min_module_size,
      binary = binary, alpha = alpha, qvalue_method = qvalue_method,
      seed = seed, scale = res$scale,
      n_mapped = nrow(proj)
    )
  )
}


#' Q-values for permutation p-values (internal)
#'
#' A fixed-`n_perm` permutation p-value lives on the discrete support
#' `{1/(n+1), 2/(n+1), ..., 1}`, and `pmax` of two such p-values lives on the
#' same set. Storey's continuous estimator is invalid on that support and
#' Benjamini-Hochberg is valid but assumes pi0 = 1; Liang's discrete method
#' estimates pi0 from the support, which is why the package already uses it
#' for the Besag-Clifford p-values in permutation_hog_test().
#'
#' Falls back to Benjamini-Hochberg when there are too few modules for pi0
#' estimation to mean anything, or if the estimator errors.
#'
#' @noRd
.pres_qvalues <- function(p, n_perm, method = "liang") {
  if (length(p) < 2L) return(p)
  bh <- function() compute_qvalues(p, pi0_method = "none")$qvalues
  if (method != "liang" || length(p) < 10L) return(bh())

  support <- seq_len(n_perm + 1L) / (n_perm + 1L)
  out <- tryCatch(
    DiscreteQvalue::DQ(p, ss = support, method = "Liang")$q.values,
    error = function(e) NULL
  )
  if (is.null(out) || anyNA(out)) bh() else out
}


#' Classify modules as conserved, moderately preserved, or diverged
#'
#' Turns [module_preservation()] output into a per-module call. The call is
#' driven by the combined permutation q-value; the `Zsummary` thresholds of
#' Langfelder et al. (2011) are applied as a secondary split among the
#' significant modules, so the familiar 2 / 10 cut points still appear.
#'
#' @section Criteria:
#' \describe{
#'   \item{conserved}{`q.value < alpha` and `Zsummary >= z_conserved`}
#'   \item{moderate}{`q.value < alpha` and `Zsummary < z_conserved`}
#'   \item{diverged}{`q.value >= alpha`}
#' }
#' Modules whose degree correlation was undefined (constant connectivity)
#' carry `NA` and are reported as `diverged` only if their density statistic
#' also fails.
#'
#' @param pres Output of [module_preservation()].
#' @param alpha Significance threshold for the combined q-value (default 0.05).
#' @param z_conserved `Zsummary` at or above which a significant module is
#'   called conserved rather than moderately preserved (default 10, the
#'   Langfelder et al. 2011 cut point).
#' @param species Optional species label recorded in the `species` column.
#' @param pair_name Optional contrast label recorded in the `pair_name` column.
#'
#' @return A data frame with `module`, `species`, `pair_name`,
#'   `classification`, `Zsummary`, `q.value`, `size` and `size_mapped`.
#'
#' @examples
#' \dontrun{
#' cls <- classify_preservation(pres)
#' table(cls$classification)
#' }
#'
#' @export
classify_preservation <- function(pres, alpha = 0.05, z_conserved = 10,
                                  species = NA_character_,
                                  pair_name = NA_character_) {
  if (!is.list(pres) || is.null(pres$preservation)) {
    stop("pres must be output from module_preservation()")
  }
  p <- pres$preservation

  significant <- !is.na(p$q.value) & p$q.value < alpha
  strong <- !is.na(p$Zsummary) & p$Zsummary >= z_conserved
  classification <- ifelse(!significant, "diverged",
    ifelse(strong, "conserved", "moderate")
  )

  data.frame(
    module = p$module,
    species = species,
    pair_name = pair_name,
    classification = classification,
    Zsummary = p$Zsummary,
    q.value = p$q.value,
    size = p$size,
    size_mapped = p$size_mapped,
    stringsAsFactors = FALSE
  )
}
