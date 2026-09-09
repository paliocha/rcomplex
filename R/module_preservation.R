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
#' Multi-copy HOGs are reduced toward one counterpart per gene by
#' [resolve_ortholog_map()] when `map` is not supplied. Resolution never
#' changes which genes are *mappable* -- that is the invariant
#' [resolve_ortholog_map()] enforces. It can still change which genes end up
#' *tested*: an unresolved gene whose candidate module labels tie is dropped
#' by the majority vote, and resolving its copy rescues it. The tested set
#' therefore does depend on the resolution, so the circularity defence rests
#' on the `p_copy` columns of `sensitivity`, not on the mappable-set
#' invariant.
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
#' @param calibrate How the combined p-value is calibrated before FDR
#'   correction. `"mixture"` (default) recalibrates `pmax` toward the
#'   empirical joint null of the two statistics; `"none"` uses `pmax`
#'   unchanged.
#'
#'   `pmax` is a valid p-value for the intersection-union null, but it is
#'   calibrated against a bound rather than the actual joint null and runs
#'   roughly `1/t` conservative when the two statistics are near-independent,
#'   as they are here -- measured on this engine, a realised false discovery
#'   rate of 1.2e-4 against a nominal 0.05, with the smallest attainable
#'   q-value 0.10. The mixture blends `pmax` with the permutation joint null
#'   in proportion to the estimated fraction of modules null on *both*
#'   statistics, which is a super-uniform bound for any dependence structure.
#'   The rejection region stays `max(p1, p2) <= c`, so a module still cannot
#'   be called preserved on density alone.
#'
#'   One caveat, measured rather than argued: when neither statistic shows
#'   signal anywhere the estimated both-null fraction approaches 1 and the
#'   procedure calibrates against the intersection null for that contrast.
#'   Realised FDR there was 0.031 against a nominal 0.05.
#' @param qvalue_method Deprecated and ignored; the Liang discrete path was
#'   measured bit-identical to Benjamini-Hochberg on this engine. Use
#'   `calibrate = "none"` for the uncalibrated result.
#' @param sensitivity Re-run under a naive ortholog map -- one built from
#'   `orthologs` alone, with no clique or coexpressolog resolution -- and
#'   report both results side by side (default `FALSE`). Costs the naive run
#'   plus `copy_draws` nested analyses; see `copy_draws`.
#'   Because resolution may only choose which paralog copy carries a label,
#'   the two runs must map the identical set of test-species genes; the
#'   returned table records whether they did. A large `Zsummary` gap with
#'   matching gene sets means the copy choice mattered; a mismatched gene set
#'   means the resolution layer is filtering rather than choosing, which is a
#'   bug.
#' @param copy_draws Number of random copy choices drawn for the
#'   `sensitivity` comparison (default 200). Each draw resolves the same
#'   candidate map by picking one species-1 partner per species-2 gene
#'   uniformly at random, so only the copy choice varies.
#'
#'   Each draw is a nested analysis of the reference network, so
#'   `sensitivity = TRUE` costs roughly `copy_draws` extra projections on top
#'   of the naive-map run. Set `copy_draws = 0L` to skip the copy null and
#'   keep only the naive-map comparison.
#' @param n_cores Number of OpenMP threads (default 1).
#' @param seed Optional RNG seed. Results are independent of `n_cores`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{preservation}{One row per tested module: `module`, `size`,
#'       `size_mapped`, the two headline statistics, their permutation
#'       p-values, the combined `p.value` (raw `pmax`), `p.calibrated` (the
#'       mixture recalibration described under `calibrate`), `q.value`
#'       (Benjamini-Hochberg on `p.calibrated`, NOT on `p.value`),
#'       `Z.avg.weight`, `Z.cor.degree`, `Zsummary`, `Zsummary_null_sd`,
#'       `Zsummary_std` and `medianRank`.}
#'     \item{observed}{All six statistics per module, with permutation means
#'       and standard deviations, plus `n_perm.<stat>` -- the number of
#'       permutations each statistic was actually scored over, which can be
#'       fewer than `n_perm` when a statistic was undefined -- and `n_joint`,
#'       the number scorable for both headline statistics at once.}
#'     \item{coverage}{One row per reference module -- `size`,
#'       `size_mapped`, `tested`, and the `reason` it was not -- so the tested
#'       set reconciles against the partition. Modules below
#'       `min_module_size`, and modules whose genes never reach the test
#'       species, leave the analysis entirely; without this the preservation
#'       table looks like a complete accounting when it is not.}
#'     \item{projection}{One row per test-species gene that received a module
#'       label: the reference gene it came from, the module, and which
#'       resolution layer chose the pair.}
#'     \item{map}{The ortholog map used.}
#'     \item{sensitivity}{Only when `sensitivity = TRUE`: per-module
#'       `size_mapped`, `Zsummary` and `q.value` under the resolved and naive
#'       maps, plus `p_copy.avg.weight` and `p_copy.cor.degree` -- where each
#'       statistic sits in a null over `copy_draws` random copy choices of the
#'       same candidate map. Attributes: `same_candidate_set` (always `TRUE`
#'       -- the mappable set is invariant by construction, so this is a
#'       structural check only), `same_projected_set` (whether the two runs
#'       tested the same genes, which resolution CAN change by rescuing genes
#'       from a tied majority vote), `n_rescued` and `n_lost` counting that
#'       difference, and `n_multi_copy` / `n_copy_draws` for the copy null
#'       (skipped, and its columns absent, when nothing is multi-copy).
#'       `size_mapped` is the deterministic consequence of the copy choice;
#'       `Zsummary_delta` also absorbs permutation-stream drift when the two
#'       runs have different block sizes.}
#'     \item{params}{Call parameters, including `calibrate`, the estimated
#'       both-null fraction `w00`, and `n_joint` per module.}
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
                                calibrate = c("mixture", "none"),
                                qvalue_method = NULL,
                                sensitivity = FALSE, copy_draws = 200L,
                                n_cores = 1L, seed = NULL) {
  if (!is.list(modules_ref) || is.null(modules_ref$module_genes) ||
    is.null(modules_ref$modules)) {
    stop("modules_ref must be output from detect_modules()")
  }
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) stop("n_perm must be >= 1")
  calibrate <- match.arg(calibrate)
  if (!is.null(qvalue_method)) {
    warning("qvalue_method is deprecated and ignored: the Liang path was ",
            "measured bit-identical to Benjamini-Hochberg on this engine. ",
            "Use calibrate = \"none\" for the uncalibrated pmax + BH result.")
  }
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

  # Every module the reference has, and what became of it. Modules below
  # min_module_size leave the analysis entirely, and modules whose genes never
  # reach the test species never appear in `sizes` at all -- reporting only
  # the tested ones makes a partial analysis look complete.
  all_mods <- names(modules_ref$module_genes)
  mapped_n <- integer(length(all_mods))
  names(mapped_n) <- all_mods
  mapped_n[names(sizes)] <- sizes
  coverage <- data.frame(
    module = all_mods,
    size = as.integer(vapply(modules_ref$module_genes, length, integer(1))),
    size_mapped = as.integer(mapped_n),
    tested = all_mods %in% tested,
    stringsAsFactors = FALSE
  )
  coverage$reason <- ifelse(
    coverage$tested, NA_character_,
    ifelse(coverage$size_mapped == 0L, "no mapped gene",
           paste0("fewer than min_module_size (", min_module_size,
                  ") mapped genes"))
  )
  rownames(coverage) <- NULL

  n_dropped <- sum(!coverage$tested)
  if (n_dropped > 0L) {
    message(n_dropped, " of ", nrow(coverage), " reference modules were not ",
            "tested (", sum(coverage$size_mapped == 0L), " with no mapped ",
            "gene, ", sum(!coverage$tested & coverage$size_mapped > 0L),
            " below min_module_size = ", min_module_size,
            "); see $coverage")
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
    n_perm, min_module_size, binary, alpha, calibrate, seed
  )
  out$coverage <- coverage

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
          calibrate = calibrate, sensitivity = FALSE,
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
        cn <- .pres_copy_null(
          out$preservation, modules_ref, net_ref, net_test, orthologs,
          genes_ref, genes_test, unique(out$projection$gene2), copy_draws,
          min_module_size, binary, n_cores
        )
        if (!is.null(cn$p_copy.avg.weight)) {
          out$sensitivity$p_copy.avg.weight <- cn$p_copy.avg.weight
          out$sensitivity$p_copy.cor.degree <- cn$p_copy.cor.degree
          attr(out$sensitivity, "n_multi_copy") <- cn$n_multi
          attr(out$sensitivity, "n_copy_draws") <- cn$n_draws
        } else {
          attr(out$sensitivity, "n_multi_copy") <- cn$n_multi
        }
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
  # The reverse direction produces no NA and would otherwise pass silently,
  # yet it is the more alarming one: resolution lost a module the naive map
  # could test.
  lost <- setdiff(b$module, a$module)
  if (length(lost) > 0L) {
    warning(length(lost), " module(s) tested under the naive map were not ",
            "tested under the resolved map (", paste(lost, collapse = ", "),
            "); paralog resolution concentrated their genes")
  }

  out <- data.frame(
    module = a$module,
    size_mapped = a$size_mapped,
    size_mapped_naive = b$size_mapped[idx],
    Zsummary = a$Zsummary,
    Zsummary_naive = b$Zsummary[idx],
    q.value = a$q.value,
    q.value_naive = b$q.value[idx],
    stringsAsFactors = FALSE
  )
  # size_mapped is the deterministic consequence of the copy choice; the
  # Zsummary difference also absorbs permutation-stream drift, because
  # differing block sizes make the two runs consume the RNG differently.
  out$Zsummary_delta <- out$Zsummary - out$Zsummary_naive

  # Compare the PROJECTED sets, not the candidate sets. resolve_ortholog_map()
  # guarantees both maps carry every candidate gene2 -- that is its
  # preserved-gene-set invariant -- so comparing map$gene2 is a tautology that
  # can never fire. What actually varies is which genes survive projection:
  # an unresolved gene whose candidate labels tie is dropped by the majority
  # vote, and resolving its copy rescues it.
  same_cand <- setequal(unique(map$gene2), unique(naive_map$gene2))
  res_proj <- unique(resolved$projection$gene2)
  nai_proj <- unique(naive$projection$gene2)
  same <- setequal(res_proj, nai_proj)

  attr(out, "same_candidate_set") <- same_cand
  attr(out, "same_projected_set") <- same
  attr(out, "n_rescued") <- length(setdiff(res_proj, nai_proj))
  attr(out, "n_lost") <- length(setdiff(nai_proj, res_proj))

  if (!same) {
    warning("paralog resolution changed which test-species genes carry a ",
            "module label: ", length(setdiff(res_proj, nai_proj)),
            " rescued from a tied majority vote, ",
            length(setdiff(nai_proj, res_proj)), " lost. The resolved and ",
            "naive Zsummary are therefore not measured on the same gene set, ",
            "so read Zsummary_delta with that in mind and prefer the p_copy ",
            "columns, which hold the set fixed.")
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

  # Vectorised: one sort plus linear passes. Splitting per gene2 and building a
  # one-row data frame each costs O(mappable genes) allocations, and on the
  # default path every row is unresolved, so the whole map goes through here --
  # once per module_preservation() and again per module_correspondence(), twice
  # more under sensitivity, and twice per contrast in preservation_paired().
  pick <- function(df) {
    if (nrow(df) == 0L) {
      return(NULL)
    }
    # gene2, then module, then gene1: the first row of each (gene2, module)
    # run is that module's smallest gene1, which is the deterministic pick.
    d <- df[order(df$gene2, df$module, df$gene1), , drop = FALSE]
    runs <- rle(paste(d$gene2, d$module, sep = "\x01"))
    # Sorted by (gene2, module), so each run start is that module's smallest
    # gene1 and the run lengths are already the per-cell counts.
    starts <- cumsum(c(1L, utils::head(runs$lengths, -1L)))
    fd <- d[starts, , drop = FALSE]
    fc <- runs$lengths

    # Modal module per gene2, dropped when two modules tie for the maximum.
    mx <- stats::ave(fc, fd$gene2, FUN = max)
    n_top <- stats::ave(as.integer(fc == mx), fd$gene2, FUN = sum)
    keep <- fc == mx & n_top == 1L

    out <- fd[keep, c("gene1", "gene2", "module", "source"), drop = FALSE]
    if (nrow(out) == 0L) {
      return(NULL)
    }
    rownames(out) <- NULL
    out
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
                      n_perm, n_cores, binary, store_perm = FALSE) {
  a <- .net_cpp_args(net, net$threshold)
  if (.net_is_sparse(net)) {
    module_preservation_sparse_cpp(
      a$p, a$i, a$x, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, binary, store_perm
    )
  } else {
    module_preservation_dense_cpp(
      a$net, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, binary, store_perm
    )
  }
}


#' Build the result tables from the kernel output (internal)
#' @noRd
.pres_assemble <- function(res, modules_ref, tested, rows_by_mod, proj, map,
                           n_perm, min_module_size, binary, alpha,
                           calibrate, seed) {
  colnames(res$observed) <- .PRES_STATS
  colnames(res$perm_mean) <- .PRES_STATS
  colnames(res$perm_sd) <- .PRES_STATS
  colnames(res$p_value) <- .PRES_STATS
  colnames(res$n_perm_used) <- .PRES_STATS

  z <- (res$observed - res$perm_mean) / res$perm_sd
  colnames(z) <- .PRES_STATS

  d <- .PRES_DENSITY
  cc <- .PRES_CONNECTIVITY

  # Both statistics must be significant: the reciprocal criterion used by
  # pval_combine = "max" elsewhere in the package. Computed on the JOINTLY
  # scorable permutations so that pmax and the NPC combination share a
  # denominator -- mismatched denominators would break the identity the
  # calibration below rests on.
  p_comb <- pmax(res$p_joint[, 1L], res$p_joint[, 2L])

  # pmax is a valid intersection-union p-value (Berger) but is calibrated
  # against a bound, not the joint null, and runs about 1/t conservative when
  # the two statistics are near-independent -- measured here as a realised FDR
  # of 1.2e-4 against a nominal 0.05. Recalibrate toward the empirical joint
  # null by the estimated fraction of modules null on BOTH statistics:
  #   F(t) <= w00 * C(t, t) + (1 - w00) * t
  # is a super-uniform bound for any dependence and any alternative, because
  # each partial-null term has one exactly-uniform marginal. p_npc is C(t, t)
  # at the observed value and p_comb is t, so p_cal is that bound evaluated
  # where it matters. Under-estimating w00 moves the result toward plain pmax,
  # i.e. toward the conservative side.
  w00 <- if (identical(calibrate, "mixture")) {
    .pres_w00(res$p_joint[, 1L], res$p_joint[, 2L])
  } else {
    0
  }
  p_cal <- ifelse(is.na(res$p_npc), p_comb,
                  w00 * res$p_npc + (1 - w00) * p_comb)
  q_comb <- .pres_qvalues(p_cal)

  # medianRank: rank of the observed statistics across modules, 1 = strongest.
  rank_d <- rank(-res$observed[, d], na.last = "keep")
  rank_c <- rank(-res$observed[, cc], na.last = "keep")
  median_rank <- (rank_d + rank_c) / 2

  # Null scale of Zsummary. Each Z is standardized by its own permutation
  # mean and sd, so it has unit null variance by construction; their mean does
  # not, and its variance depends on how correlated the two statistics are
  # under the null: sd = sqrt(2 + 2*rho) / 2. Langfelder et al.'s 10 / 2 cut
  # points were set for a Zsummary built from medians over several statistics,
  # a quantity with a different null spread, so reading them against a raw
  # mean of two is reading them on an unknown scale.
  # Computed in the kernel over the jointly scorable draws, so the covariance
  # and both marginal spreads come from the same set and the same convention.
  rho <- res$rho_null
  # rho is a sample correlation over n_perm draws, so a negative value is an
  # ordinary sampling outcome at small n_perm -- and sqrt(2 + 2*rho)/2 shrinks
  # toward 0 as it goes negative, which would inflate Zsummary_std without
  # bound and turn a diverged module into a conserved one. Clamp at 0: the two
  # statistics are near-independent by construction (cor.degree is
  # scale-invariant, avg.weight is pure scale), so a negative estimate is
  # noise, and clamping keeps the divisor in the honest range [1/sqrt(2), 1].
  rho <- pmin(pmax(rho, 0), 1)
  z_null_sd <- sqrt(2 + 2 * rho) / 2

  size_all <- vapply(modules_ref$module_genes, length, integer(1))

  preservation <- data.frame(
    module = tested,
    size = as.integer(size_all[tested]),
    size_mapped = vapply(rows_by_mod, length, integer(1)),
    avg.weight = res$observed[, d],
    cor.degree = res$observed[, cc],
    p.avg.weight = res$p_joint[, 1L],
    p.cor.degree = res$p_joint[, 2L],
    p.value = p_comb,
    p.calibrated = p_cal,
    q.value = q_comb,
    Z.avg.weight = z[, d],
    Z.cor.degree = z[, cc],
    Zsummary = (z[, d] + z[, cc]) / 2,
    Zsummary_null_sd = z_null_sd,
    Zsummary_std = ((z[, d] + z[, cc]) / 2) / z_null_sd,
    medianRank = median_rank,
    stringsAsFactors = FALSE
  )
  rownames(preservation) <- NULL

  observed <- data.frame(module = tested, stringsAsFactors = FALSE)
  for (s in .PRES_STATS) {
    observed[[s]] <- res$observed[, s]
    observed[[paste0("perm_mean.", s)]] <- res$perm_mean[, s]
    observed[[paste0("perm_sd.", s)]] <- res$perm_sd[, s]
    # The denominator each p-value was actually computed over. The kernel has
    # counted these since day one and the R layer discarded them; they are the
    # audit trail for the support.
    observed[[paste0("n_perm.", s)]] <- res$n_perm_used[, s]
  }
  observed$n_joint <- res$n_joint
  rownames(observed) <- NULL

  list(
    preservation = preservation,
    observed = observed,
    projection = proj,
    map = map,
    params = list(
      n_perm = n_perm, min_module_size = min_module_size,
      binary = binary, alpha = alpha, calibrate = calibrate,
      w00 = w00, n_joint = res$n_joint,
      seed = seed, scale = res$scale,
      n_mapped = nrow(proj)
    )
  )
}


#' Fraction of modules null on BOTH statistics (internal)
#'
#' The Frechet lower bound on the both-null fraction: pi0 for each margin,
#' summed and shifted. It under-estimates by construction, and it estimates the
#' unconditional both-null fraction where the bound wants the conditional one,
#' so it errs conservative twice over. Storey's estimator at lambda = 0.5 is
#' high-variance on a handful of p-values, so below `min_m` modules it returns
#' 0, which degrades the calibration to plain pmax rather than guessing.
#'
#' @noRd
.pres_w00 <- function(p1, p2, lambda = 0.5, min_m = 10L) {
  ok <- !is.na(p1) & !is.na(p2)
  if (sum(ok) < min_m) {
    return(0)
  }
  st0 <- function(p) min(1, mean(p > lambda) / (1 - lambda))
  max(0, st0(p1[ok]) + st0(p2[ok]) - 1)
}


#' Benjamini-Hochberg on the calibrated p-value (internal)
#'
#' After calibration the p-value lives on a grid of a thousand points or more,
#' where a discrete correction buys about 1% of rejections. DiscreteQvalue's
#' Liang path was measured bit-identical to BH on this engine, and
#' qvalue::qvalue() errors outright on the p-value shapes produced at these
#' module counts, so BH is the whole of it.
#'
#' @noRd
.pres_qvalues <- function(p) {
  # BH on the calibrated p-value. After calibration the p-value lives on a
  # grid of a thousand points or more, where the discrete correction buys
  # about 1% of rejections; DiscreteQvalue's Liang path was measured
  # bit-identical to BH on this engine, and qvalue::qvalue errors outright on
  # the p-value shapes produced at these module counts.
  ok <- !is.na(p)
  out <- rep(NA_real_, length(p))
  if (sum(ok) < 2L) {
    out[ok] <- p[ok]
    return(out)
  }
  out[ok] <- compute_qvalues(p[ok], pi0_method = "none")$qvalues
  out
}



#' Null distribution over paralog copy choices (internal)
#'
#' The resolved map picks copies using coexpressolog evidence, and a
#' coexpressolog is a gene pair whose co-expression neighbourhoods overlap --
#' the same signal `avg.weight` and `cor.degree` measure. The permutation null
#' cannot see that selection, because it draws random gene blocks with no copy
#' choice at all, so the observed statistics are inflated relative to it.
#'
#' This draws `n_draws` alternative resolutions restricted to the species-2
#' genes the run under test actually projected, each picking one species-1
#' partner per gene uniformly at random, and reports where the real map's
#' statistics sit in that distribution. Restricting to the projected set is
#' what makes only the copy choice vary: an unrestricted draw projects a
#' different number of genes, because resolving a copy rescues genes whose
#' candidate labels would otherwise tie in the majority vote, and that
#' set-size difference would be read as a copy-choice effect.
#'
#' A large `p_copy` means the result is typical of an arbitrary copy choice and
#' the resolution did not manufacture it. A small `p_copy` means the finding
#' depends on having picked those particular copies, which is exactly the
#' circularity to distrust.
#'
#' @noRd
.pres_copy_null <- function(observed, modules_ref, net_ref, net_test,
                            orthologs, genes_ref, genes_test, projected,
                            n_draws, min_module_size, binary,
                            n_cores = 1L) {
  cand <- resolve_ortholog_map(orthologs, genes_ref, genes_test)
  # Only candidates whose reference partner carries a module label can project,
  # and only the genes the run under test actually projected may enter -- a
  # draw that projected a different gene set would confound a set-size
  # difference with the copy-choice effect this null exists to isolate.
  cand <- cand[!is.na(modules_ref$modules[cand$gene1]), , drop = FALSE]
  cand <- cand[cand$gene2 %in% projected, , drop = FALSE]
  by_g2 <- split(seq_len(nrow(cand)), cand$gene2)
  # With a caller-supplied map the candidates need not cover the projected
  # genes, and then every draw would score fewer genes -- reintroducing the
  # set-size artefact this null exists to remove.
  if (!setequal(names(by_g2), projected)) {
    warning("the ortholog table does not cover every projected gene (",
            length(setdiff(projected, names(by_g2))), " missing), so the ",
            "copy-choice null cannot hold the gene set fixed; skipping it")
    return(list(draws = NULL, n_multi = 0L))
  }
  multi <- sum(lengths(by_g2) > 1L)
  if (multi == 0L) {
    return(list(draws = NULL, n_multi = 0L))
  }

  draws <- vector("list", n_draws)
  for (d in seq_len(n_draws)) {
    idx <- vapply(by_g2, function(ix) {
      if (length(ix) == 1L) ix else ix[sample.int(length(ix), 1L)]
    }, integer(1))
    m <- cand[idx, , drop = FALSE]
    # One row per gene2, so the majority vote is unanimous and the projected
    # set is exactly `projected` -- the set the observed run used.
    m$source <- "random"
    # Each draw is a nested analysis, so mute its reporting: otherwise a run
    # with many untested modules prints the coverage message and the
    # flat-connectivity warning once per draw.
    draws[[d]] <- tryCatch(
      suppressMessages(suppressWarnings(
        module_preservation(
          modules_ref, net_ref, net_test, map = m, n_perm = 1L,
          min_module_size = min_module_size, binary = binary,
          sensitivity = FALSE, n_cores = n_cores
        )$preservation
      )),
      error = function(e) NULL
    )
  }
  draws <- Filter(Negate(is.null), draws)
  if (length(draws) == 0L) {
    return(list(draws = NULL, n_multi = multi))
  }

  stat_p <- function(col) {
    vapply(seq_len(nrow(observed)), function(i) {
      vals <- vapply(draws, function(d) {
        j <- match(observed$module[i], d$module)
        if (is.na(j)) NA_real_ else d[[col]][j]
      }, numeric(1))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0L || is.na(observed[[col]][i])) {
        return(NA_real_)
      }
      (sum(vals >= observed[[col]][i]) + 1) / (length(vals) + 1)
    }, numeric(1))
  }
  list(
    n_multi = multi,
    n_draws = length(draws),
    p_copy.avg.weight = stat_p("avg.weight"),
    p_copy.cor.degree = stat_p("cor.degree")
  )
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
#'   \item{conserved}{`q.value < alpha` and the `Zsummary` scale chosen by
#'     `z_scale` is at or above `z_conserved`}
#'   \item{moderate}{`q.value < alpha` and it is below `z_conserved`}
#'   \item{diverged}{`q.value >= alpha`}
#'   \item{untested}{`q.value` is `NA` -- a statistic could not be computed,
#'     so neither preservation nor divergence was measured}
#' }
#' A module whose degree correlation could not be computed -- `cor.degree` is
#' undefined when intramodular connectivity is constant in either network --
#' gets `NA` for that statistic, so the `pmax` combination and hence `q.value`
#' are `NA` too. Such a module is reported `"untested"` rather than
#' `"diverged"`: nothing was measured, so divergence would be a positive claim
#' the data does not support. It is never called preserved on the density
#' statistic alone.
#'
#' @param pres Output of [module_preservation()].
#' @param alpha Significance threshold for the combined q-value (default 0.05).
#' @param z_conserved Cut point at or above which a significant module is
#'   called conserved rather than moderately preserved (default 10).
#' @param z_scale Which `Zsummary` the cut point is read against.
#'   `"standardized"` (default) uses `Zsummary_std`, which has unit variance
#'   under the permutation null; `"raw"` uses `Zsummary` itself.
#'
#'   This matters because the familiar 10 / 2 cut points are not
#'   scale-free. Langfelder et al. (2011) set them for a `Zsummary` built as
#'   the mean of two medians over several density and connectivity statistics.
#'   Here only two statistics are available from an adjacency matrix, so
#'   `Zsummary` is the mean of two standardized values and its null spread is
#'   `sqrt(2 + 2 * rho) / 2`, where `rho` is the null correlation between
#'   `avg.weight` and `cor.degree` -- between 0.71 and 1 rather than the
#'   smaller spread a median-of-several carries. Reading 10 against the raw
#'   mean therefore imports a number calibrated on a different quantity.
#'   Dividing by the null sd restores its intended meaning, "this many null
#'   standard deviations". Use `"raw"` only to reproduce older output.
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
                                  z_scale = c("standardized", "raw"),
                                  species = NA_character_,
                                  pair_name = NA_character_) {
  z_scale <- match.arg(z_scale)
  if (!is.list(pres) || is.null(pres$preservation)) {
    stop("pres must be output from module_preservation()")
  }
  p <- pres$preservation

  testable <- !is.na(p$q.value)
  significant <- testable & p$q.value < alpha
  # Read the cut point on a scale where it means what it is meant to mean:
  # one null standard deviation. Zsummary is a mean of two standardized
  # statistics, so its own null sd is sqrt(2 + 2*rho)/2, not 1.
  z_used <- if (z_scale == "standardized" && "Zsummary_std" %in% names(p)) {
    p$Zsummary_std
  } else {
    p$Zsummary
  }
  strong <- !is.na(z_used) & z_used >= z_conserved
  classification <- ifelse(!testable, "untested",
    ifelse(!significant, "diverged",
      ifelse(strong, "conserved", "moderate")
    )
  )
  if (any(!testable)) {
    where <- paste(stats::na.omit(c(species, pair_name)), collapse = " / ")
    warning(sum(!testable), " module(s) could not be tested (a statistic was ",
            "undefined); reported as \"untested\"",
            if (nzchar(where)) paste0(" [", where, "]") else "",
            ": ", paste(p$module[!testable], collapse = ", "))
  }

  # rep() rather than recycling: with zero modules the scalar species and
  # pair_name would otherwise clash with the 0-length columns.
  n <- nrow(p)
  data.frame(
    module = p$module,
    species = rep(species, length.out = n),
    pair_name = rep(pair_name, length.out = n),
    classification = classification,
    Zsummary = p$Zsummary,
    Zsummary_std = if ("Zsummary_std" %in% names(p)) {
      p$Zsummary_std
    } else {
      rep(NA_real_, n)
    },
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
#' topology is preserved ([module_preservation()]);
#' [classify_hub_conservation()] needs the correspondence, not the
#' preservation call.
#'
#' The map assigns one reference gene to each test-species gene, so the
#' hypergeometric's independence assumption holds. The multi-copy expansion
#' that made the same test anti-conservative under the retired
#' gene-overlap engine is gone here: a HOG with three paralogs no longer
#' contributes three correlated draws to the same urn.
#'
#' @param modules_ref,modules_test Module detection results
#'   (output of [detect_modules()]) for the two species.
#' @param map Ortholog map from [resolve_ortholog_map()], with `gene1` in the
#'   reference species and `gene2` in the test species.
#' @param qvalue_method Passed to `compute_qvalues()`; `"randomized"`
#'   (default) estimates pi0 on randomized p-values, which is what the package
#'   uses elsewhere for discrete hypergeometric p-values.
#' @param sp_ref,sp_test Optional species labels recorded on the result.
#'   `module_sp1` belongs to `modules_ref`, and that orientation cannot be
#'   recovered from the table, so supplying these lets
#'   [classify_hub_conservation()] catch a transposed call.
#'
#' @return A list with `pairs` -- a data frame of `module_sp1`,
#'   `module_sp2`, `size_sp1`, `size_sp2`, `overlap`, `jaccard`, `p.value` and
#'   `q.value`, one row per module pair -- plus the `sp_ref` / `sp_test`
#'   labels. The column names match what [classify_hub_conservation()]
#'   expects.
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
                                  qvalue_method = "randomized",
                                  sp_ref = NULL, sp_test = NULL) {
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

  # Orientation is not recoverable from the table, so record it: module_sp1
  # belongs to modules_ref. classify_hub_conservation() checks it against the
  # key when present.
  list(pairs = pairs, sp_ref = sp_ref, sp_test = sp_test)
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
#'   the owning species' group for diverged ones. Modules reported
#'   `"untested"` are counted under `"untested"` rather than a trait group:
#'   nothing was measured, so attributing them to one would overstate the
#'   evidence.
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
#'   This table feeds [tag_permutation()] directly: it selects the
#'   `"diverged"` rows and uses `reference` / `test` to pick each contrast's
#'   two sides.
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
#' @param ... Additional arguments passed to the default method.
#' @export
preservation_paired <- function(modules, ...) {
  UseMethod("preservation_paired")
}

#' @rdname preservation_paired
#' @export
preservation_paired.default <- function(modules, networks, orthologs, pairs,
                                        group = NULL, edges = NULL,
                                        cliques = NULL, alpha = 0.05,
                                        z_conserved = 10, ...) {
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
        # "untested" is not evidence of species-specific divergence, so it
        # earns no trait-group attribution -- but it is carried as its own
        # level rather than NA, which stats::aggregate() would silently drop
        # from the summary under its default na.omit.
        cls$group <- ifelse(
          cls$classification == "untested", "untested",
          ifelse(cls$classification != "diverged",
                 "conserved", as.character(group[ref]))
        )
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
