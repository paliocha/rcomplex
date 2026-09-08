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
#' @param sensitivity Re-run under a naive ortholog map -- one built from
#'   `orthologs` alone, with no clique or coexpressolog resolution -- and
#'   report both results side by side (default `FALSE`). Doubles the runtime.
#'   Because resolution may only choose which paralog copy carries a label,
#'   the two runs must map the identical set of test-species genes; the
#'   returned table records whether they did. A large `Zsummary` gap with
#'   matching gene sets means the copy choice mattered; a mismatched gene set
#'   means the resolution layer is filtering rather than choosing, which is a
#'   bug.
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
#'     \item{projection}{One row per test-species gene that received a module
#'       label: the reference gene it came from, the module, and which
#'       resolution layer chose the pair.}
#'     \item{map}{The ortholog map used.}
#'     \item{sensitivity}{Only when `sensitivity = TRUE`: per-module
#'       `Zsummary` and `q.value` under the resolved and naive maps, plus the
#'       `same_gene_set` attribute recording whether both mapped the identical
#'       test-species genes.}
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
                                sensitivity = FALSE,
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

  supplied_map <- !is.null(map)
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

  out <- .pres_assemble(
    res, modules_ref, tested, rows_by_mod, proj, map,
    n_perm, min_module_size, binary, alpha, qvalue_method, seed
  )

  if (isTRUE(sensitivity)) {
    if (is.null(orthologs)) {
      warning("sensitivity = TRUE needs 'orthologs' to build the naive map; ",
              "skipping the comparison")
    } else if (!supplied_map && is.null(edges) && is.null(cliques)) {
      # Nothing resolved any copy, so the map already IS the naive map and a
      # second run would spend a full n_perm to report a delta of exactly zero.
      warning("sensitivity = TRUE has nothing to compare: without 'edges', ",
              "'cliques' or a supplied 'map' the resolved map is already the ",
              "naive map; skipping the second run")
    } else {
      naive_map <- resolve_ortholog_map(orthologs, genes_ref, genes_test)
      # The naive run must never discard the primary result: it can legitimately
      # stop(), e.g. when all its modules fall below min_module_size.
      naive <- tryCatch(
        module_preservation(
          modules_ref, net_ref, net_test, orthologs = orthologs,
          map = naive_map, n_perm = n_perm,
          min_module_size = min_module_size, binary = binary, alpha = alpha,
          qvalue_method = qvalue_method, sensitivity = FALSE,
          n_cores = n_cores, seed = seed
        ),
        error = function(e) {
          warning("the naive-map run failed (", conditionMessage(e),
                  "); skipping the sensitivity comparison")
          NULL
        }
      )
      if (!is.null(naive)) {
        out$sensitivity <- .pres_sensitivity(out, naive, map, naive_map)
      }
    }
  }

  out
}


#' Compare a resolved run against its naive-map counterpart (internal)
#'
#' Resolution may only change which paralog copy carries a module label, never
#' which genes are mappable, so both runs must cover the identical set of
#' test-species genes. That invariant is what keeps the circularity in check:
#' coexpressologs are defined by conserved neighbourhoods and preservation
#' measures conserved topology, so a resolution layer that also filtered the
#' mapped set would be selecting the tested genes on the statistic being
#' tested.
#'
#' @noRd
.pres_sensitivity <- function(resolved, naive, map, naive_map) {
  a <- resolved$preservation
  b <- naive$preservation
  idx <- match(a$module, b$module)
  if (anyNA(idx)) {
    warning(sum(is.na(idx)), " module(s) tested under the resolved map were ",
            "not tested under the naive map; their naive columns are NA")
  }

  out <- data.frame(
    module = a$module,
    Zsummary = a$Zsummary,
    Zsummary_naive = b$Zsummary[idx],
    q.value = a$q.value,
    q.value_naive = b$q.value[idx],
    stringsAsFactors = FALSE
  )
  out$Zsummary_delta <- out$Zsummary - out$Zsummary_naive

  same <- setequal(unique(map$gene2), unique(naive_map$gene2))
  attr(out, "same_gene_set") <- same
  if (!same) {
    warning("the resolved and naive ortholog maps cover different ",
            "test-species genes; paralog resolution should only choose which ",
            "copy carries a label, so this indicates a filtering bug and the ",
            "Zsummary comparison is not interpretable")
  }
  out
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
  # A statistic that could not be computed carries an NA p-value; correct the
  # rest and leave those NA, rather than letting them error here or be scored
  # as significant downstream.
  ok <- !is.na(p)
  out <- rep(NA_real_, length(p))
  if (sum(ok) < 2L) {
    out[ok] <- p[ok]
    return(out)
  }
  pv <- p[ok]

  bh <- function() compute_qvalues(pv, pi0_method = "none")$qvalues
  q <- if (method != "liang" || length(pv) < 10L) {
    bh()
  } else {
    support <- seq_len(n_perm + 1L) / (n_perm + 1L)
    liang <- tryCatch(
      DiscreteQvalue::DQ(pv, ss = support, method = "Liang")$q.values,
      error = function(e) NULL
    )
    if (is.null(liang) || anyNA(liang)) bh() else liang
  }
  out[ok] <- q
  out
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
#' A module whose degree correlation could not be computed -- `cor.degree` is
#' undefined when intramodular connectivity is constant, or when fewer than
#' three genes map -- gets `NA` for that statistic, so the `pmax` combination
#' and hence `q.value` are `NA` too and the module is reported `diverged`. It
#' is never called preserved on the density statistic alone.
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


#' Match modules across species by ortholog overlap
#'
#' Cross-tabulates the modules of two species over a paralog-resolved ortholog
#' map and tests each module pair for excess overlap. This answers "which
#' module corresponds to which", a different question from whether a module's
#' topology is preserved ([module_preservation()]); [classify_hub_conservation()]
#' needs the correspondence, not the preservation call.
#'
#' The map assigns one reference gene to each test-species gene, so the
#' hypergeometric's independence assumption holds. The multi-copy expansion
#' that makes the same test anti-conservative in `compare_modules()` is gone
#' here: a HOG with three paralogs no longer contributes three correlated
#' draws to the same urn.
#'
#' @param modules_ref,modules_test Module detection results
#'   (output of [detect_modules()]) for the two species.
#' @param map Ortholog map from [resolve_ortholog_map()], with `gene1` in the
#'   reference species and `gene2` in the test species.
#' @param qvalue_method Passed to `compute_qvalues()`; `"randomized"`
#'   (default) estimates pi0 on randomized p-values, which is what the package
#'   uses elsewhere for discrete hypergeometric p-values.
#'
#' @return A list with one element, `pairs`: a data frame of `module_sp1`,
#'   `module_sp2`, `size_sp1`, `size_sp2`, `overlap`, `jaccard`, `p.value` and
#'   `q.value`, one row per module pair. The column names match what
#'   [classify_hub_conservation()] expects.
#'
#' @examples
#' \dontrun{
#' map <- resolve_ortholog_map(ortho, rownames(net_a$network),
#'   rownames(net_b$network))
#' corr <- module_correspondence(mods_a, mods_b, map)
#' subset(corr$pairs, q.value < 0.05)
#' }
#'
#' @export
module_correspondence <- function(modules_ref, modules_test, map,
                                  qvalue_method = "randomized") {
  for (nm in c("modules_ref", "modules_test")) {
    m <- get(nm)
    if (!is.list(m) || is.null(m$modules) || is.null(m$module_genes)) {
      stop(nm, " must be output from detect_modules()")
    }
  }
  if (!is.data.frame(map) ||
    !all(c("gene1", "gene2", "source") %in% names(map))) {
    stop("map must be a data frame from resolve_ortholog_map()")
  }

  map$module <- as.character(modules_ref$modules[map$gene1])
  map <- map[!is.na(map$module), , drop = FALSE]
  proj <- .pres_project(map)
  if (is.null(proj) || nrow(proj) == 0L) {
    stop("no test-species gene received an unambiguous module label")
  }

  proj$module_test <- as.character(modules_test$modules[proj$gene2])
  proj <- proj[!is.na(proj$module_test), , drop = FALSE]
  if (nrow(proj) == 0L) {
    stop("no mapped gene falls in a module of the test species")
  }

  tab <- table(proj$module, proj$module_test)
  n_total <- nrow(proj)
  ref_n <- rowSums(tab)
  test_n <- colSums(tab)

  i <- rep(seq_len(nrow(tab)), ncol(tab))
  j <- rep(seq_len(ncol(tab)), each = nrow(tab))
  overlap <- as.integer(tab)
  m <- as.integer(ref_n[i])
  k <- as.integer(test_n[j])

  p_gt <- stats::phyper(overlap, m, n_total - m, k, lower.tail = FALSE)
  p_eq <- stats::dhyper(overlap, m, n_total - m, k)

  pairs <- data.frame(
    module_sp1 = rownames(tab)[i],
    module_sp2 = colnames(tab)[j],
    size_sp1 = m,
    size_sp2 = k,
    overlap = overlap,
    jaccard = overlap / (m + k - overlap),
    p.value = p_gt + p_eq,
    stringsAsFactors = FALSE
  )

  pairs$q.value <- if (nrow(pairs) < 2L) {
    pairs$p.value
  } else {
    compute_qvalues(
      pairs$p.value,
      p_rand_fn = function() p_gt + stats::runif(length(p_gt)) * p_eq,
      pi0_method = qvalue_method
    )$qvalues
  }
  rownames(pairs) <- NULL

  list(pairs = pairs)
}


#' Run module preservation across many species pairs
#'
#' Applies [module_preservation()] to each contrast in `pairs`. Preservation is
#' directional -- whether species A's modules survive in B is a different
#' question from the reverse -- so both directions are run and reported
#' separately.
#'
#' @param modules Named list of [detect_modules()] results, keyed by species.
#' @param networks Named list of [compute_network()] results, keyed by species.
#' @param orthologs Data frame with columns `Species1`, `Species2`, `hog`.
#' @param pairs Data frame with columns `sp1`, `sp2` and optionally
#'   `pair_name`.
#' @param group Optional named vector mapping species to a trait group. When
#'   supplied, a `group` column records `"conserved"` for preserved modules and
#'   the owning species' group for diverged ones.
#' @param edges,cliques Optional [find_coexpressologs()] and [find_cliques()]
#'   results, used for paralog resolution.
#' @param alpha,z_conserved Passed to [classify_preservation()].
#' @param ... Further arguments passed to [module_preservation()].
#'
#' @return A list with `classification` (one row per module per direction,
#'   carrying `pair_name`, `module`, `species`, `reference`, `test` and
#'   `classification`), `summary` (counts per contrast and direction) and `raw`
#'   (the [module_preservation()] results, keyed by `"<reference>.<test>"`).
#'
#'   Note that [tag_permutation()] cannot consume this table yet: it selects
#'   rows on `classification == "species_specific"` and `species %in%
#'   c("sp1", "sp2")`, the vocabulary of the overlap engine, whereas this
#'   function emits `"conserved"` / `"moderate"` / `"diverged"` and real
#'   species names. Feeding it straight in yields an empty result with no
#'   error. [tag_permutation()] is updated when the overlap engine is removed.
#'
#' @examples
#' \dontrun{
#' res <- preservation_paired(mods, nets, ortho,
#'   pairs = data.frame(sp1 = "BDIS", sp2 = "BSYL"),
#'   group = c(BDIS = "annual", BSYL = "perennial")
#' )
#' res$summary
#' }
#'
#' @export
preservation_paired <- function(modules, networks, orthologs, pairs,
                                group = NULL, edges = NULL, cliques = NULL,
                                alpha = 0.05, z_conserved = 10, ...) {
  if (!is.list(modules) || is.null(names(modules))) {
    stop("modules must be a named list keyed by species")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  if (!is.data.frame(pairs) || !all(c("sp1", "sp2") %in% names(pairs))) {
    stop("pairs must have columns 'sp1' and 'sp2'")
  }
  if (any(pairs$sp1 == pairs$sp2)) {
    stop("pairs must not compare a species with itself")
  }
  # Both directions of every contrast are run, so results are keyed by
  # "<reference>.<test>". A repeated contrast -- in either orientation -- would
  # collide on that key and silently overwrite the earlier one.
  unordered <- paste(pmin(pairs$sp1, pairs$sp2), pmax(pairs$sp1, pairs$sp2))
  if (anyDuplicated(unordered) > 0L) {
    stop("pairs lists the same species pair more than once (both directions ",
         "of each contrast are run, so (A, B) and (B, A) are the same row): ",
         paste(unique(unordered[duplicated(unordered)]), collapse = ", "))
  }
  species <- unique(c(pairs$sp1, pairs$sp2))
  missing_sp <- setdiff(species, intersect(names(modules), names(networks)))
  if (length(missing_sp) > 0L) {
    stop("modules and networks must both cover: ",
         paste(missing_sp, collapse = ", "))
  }
  if (!is.null(group)) {
    missing_grp <- setdiff(species, names(group))
    if (length(missing_grp) > 0L) {
      stop("group missing entries for: ",
           paste(missing_grp, collapse = ", "))
    }
  }
  if (!"pair_name" %in% names(pairs)) {
    pairs$pair_name <- paste(pairs$sp1, pairs$sp2, sep = ".")
  }

  raw <- list()
  class_list <- list()

  for (p in seq_len(nrow(pairs))) {
    for (direction in list(c(pairs$sp1[p], pairs$sp2[p]),
                           c(pairs$sp2[p], pairs$sp1[p]))) {
      ref <- direction[1]
      test <- direction[2]
      key <- paste(ref, test, sep = ".")

      pres <- module_preservation(
        modules[[ref]], networks[[ref]], networks[[test]],
        .orient_orthologs(
          orthologs, rownames(networks[[ref]]$network),
          rownames(networks[[test]]$network)
        ),
        edges = edges, cliques = cliques, sp_ref = ref, sp_test = test,
        alpha = alpha, ...
      )
      raw[[key]] <- pres

      cls <- classify_preservation(pres, alpha = alpha,
                                   z_conserved = z_conserved,
                                   species = ref,
                                   pair_name = pairs$pair_name[p])
      cls$reference <- ref
      cls$test <- test
      if (!is.null(group)) {
        cls$group <- ifelse(cls$classification != "diverged",
                            "conserved", as.character(group[ref]))
      }
      class_list[[key]] <- cls
    }
  }

  classification <- do.call(rbind, class_list)
  rownames(classification) <- NULL

  count_col <- if (is.null(group)) "classification" else "group"
  summary_df <- stats::aggregate(
    stats::as.formula(paste("module ~ pair_name + reference +", count_col)),
    data = classification, FUN = length
  )
  names(summary_df)[ncol(summary_df)] <- "n"

  list(classification = classification, summary = summary_df, raw = raw)
}


#' Put an ortholog table in the orientation a direction needs (internal)
#'
#' `Species1` / `Species2` hold gene identifiers, and the rest of the package
#' selects rows by testing them against each network's gene names. A table
#' written for the A -> B direction therefore yields nothing when B is the
#' reference, so the columns are swapped whenever that recovers more pairs.
#'
#' @noRd
.orient_orthologs <- function(orthologs, genes_ref, genes_test) {
  as_is <- sum(orthologs$Species1 %in% genes_ref &
    orthologs$Species2 %in% genes_test)
  swapped <- sum(orthologs$Species2 %in% genes_ref &
    orthologs$Species1 %in% genes_test)
  if (swapped > as_is) {
    orthologs[c("Species1", "Species2")] <-
      orthologs[c("Species2", "Species1")]
  }
  orthologs
}
