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
.PRES_STATS <- c(  # nolint
  "avg.weight", "meanClusterCoeff", "meanMAR",
  "cor.degree", "cor.clusterCoeff", "cor.MAR"
)

# The two that carry the call: indices into .PRES_STATS.
.PRES_DENSITY <- 1L  # nolint
.PRES_CONNECTIVITY <- 4L  # nolint


#' Test whether co-expression modules are preserved across species
#'
#' Projects the modules of a reference species onto a test species through
#' an ortholog map. Tests, per module, whether its topology survives in the
#' test species network.
#'
#' Unlike gene-overlap tests, this test reports a module as diverged when
#' it keeps its gene membership but loses its internal wiring.
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
#' [resolve_ortholog_map()]. Resolution never
#' changes which genes are *mappable* -- that is the invariant
#' [resolve_ortholog_map()] enforces. It can still change which genes end up
#' *tested*: an unresolved gene whose candidate module labels tie is dropped
#' by the majority vote, and resolving its copy rescues it. The tested set
#' therefore does depend on the resolution.
#'
#' @param modules_ref Module assignment for the reference species, from
#'   [detect_modules()] or [as_modules()].
#' @param net_ref,net_test Network objects from [compute_network()] for the
#'   reference and test species.
#' @param orthologs Data frame with columns `gene1`, `gene2`, `hog`.
#' @param edges Optional [find_coexpressologs()] result, passed to
#'   [resolve_ortholog_map()] for paralog resolution.
#' @param species_ref,species_test Species labels, required only when `edges`
#'   is supplied.
#' @param n_perm Number of permutations (default 10000). The smallest
#'   attainable p-value is `1 / (n_perm + 1)`, so this sets the floor on how
#'   significant any module can be: at 1000 permutations every strongly
#'   preserved module ties at `p = 0.000999` and cannot be ranked by q-value.
#'   The permutation loop costs O(edges) per iteration -- roughly 0.7 s per
#'   1000 permutations on a 6000-gene network at density 0.03 -- so the
#'   default buys a floor of 1e-4 cheaply. Raise it further when many modules
#'   sit at the floor.
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
#' @param n_cores Number of OpenMP threads (default 1).
#' @param seed Optional RNG seed. Results are independent of `n_cores`.
#'   `NULL` (default) draws from the ambient stream and leaves it advanced;
#'   a seed draws from a private stream and restores the caller's on exit,
#'   the package-wide contract described under [detect_modules()].
#'
#' @return A list with components:
#'   \describe{
#'     \item{preservation}{One row per tested module: `module`, `size`,
#'       `size_mapped`, the two headline statistics, their permutation
#'       p-values, the combined `p_value` (raw `pmax`), `p_calibrated` (the
#'       mixture recalibration described under `calibrate`), `q_value`
#'       (Benjamini-Hochberg on `p_calibrated`, NOT on `p_value`),
#'       `evalue` (`p_calibrated` times the number of modules tested: the
#'       expected count of modules this preserved by chance),
#'       `Z.avg.weight`, `Z.cor.degree`, `Zsummary`, `Zsummary_null_sd`,
#'       `Zsummary_std`, and `medianRank` -- the mean of the `avg.weight`
#'       and `cor.degree` ranks across the tested modules, 1 = strongest.
#'       No permutation moments enter it, so it avoids the null-sd
#'       module-size dependence of `Zsummary` -- it is a rank *within one
#'       run* and so is not comparable between runs that tested different
#'       numbers of modules.
#'       [classify_preservation()] carries it through.}
#'     \item{observed}{All six statistics per module, with permutation means
#'       and standard deviations, plus `n_perm.<stat>` -- the number of
#'       permutations each statistic was actually scored over, which can be
#'       fewer than `n_perm` when a statistic was undefined -- and `n_joint`,
#'       the number scorable for both headline statistics at once.}
#'     \item{coverage}{One row per reference module -- `size`,
#'       `size_mapped`, `tested`, and the `reason` it was not -- so the tested
#'       set reconciles against the partition. Modules below
#'       10 mapped genes, and modules whose genes never reach the test
#'       species, leave the analysis entirely; without this the preservation
#'       table looks like a complete accounting when it is not.}
#'     \item{projection}{One row per test-species gene that received a module
#'       label: the reference gene it came from, the module, and which
#'       resolution layer chose the pair.}
#'     \item{map}{The ortholog map used.}
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
#'   edges = rcx$edges,
#'   species_ref = "SP_A", species_test = "SP_B"
#' )
#' classify_preservation(pres)
#' }
#'
#' @export
module_preservation <- function(modules_ref, net_ref, net_test,
                                orthologs, edges = NULL,
                                species_ref = NULL, species_test = NULL,
                                n_perm = 10000L,
                                calibrate = c("mixture", "none"),
                                n_cores = 1L, seed = NULL) {
  if (!is.list(modules_ref) || is.null(modules_ref$module_genes) ||
        is.null(modules_ref$modules)) {
    stop(
      "modules_ref must be a module assignment from detect_modules() ",
      "or as_modules()"
    )
  }
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) stop("n_perm must be >= 1")
  calibrate <- match.arg(calibrate)
  # Covers everything that draws: the C++ kernel's per-thread seeds. Nothing
  # between here and the kernel call consumes the stream, so a seeded run is
  # bit-identical to the pre-0.3.0 set.seed() that sat just above it. See
  # .seed_scope() in R/rng.R.
  .seed_scope(seed)

  mat_ref <- .net_check(net_ref, net_ref$threshold)
  mat_test <- .net_check(net_test, net_test$threshold)
  genes_ref <- rownames(mat_ref)
  genes_test <- rownames(mat_test)

  map <- resolve_ortholog_map(orthologs, genes_ref, genes_test,
    species1 = species_ref, species2 = species_test, edges = edges
  )
  # character, not factor: a factor would index modules by its integer
  # codes, and radix-order by codes that carry another session's collation
  map$gene1 <- as.character(map$gene1)
  map$gene2 <- as.character(map$gene2)

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

  # radix (C-locale) order, so a seeded null is drawn in the same module
  # order on every machine whatever the labels' case or alphabet
  rows_by_mod <- split(seq_len(nrow(proj)), .radix_factor(proj$module))
  sizes <- vapply(rows_by_mod, length, integer(1))
  tested <- names(sizes)[sizes >= 10L]
  if (length(tested) == 0L) {
    stop(
      "no module has at least 10 mapped test-species genes; largest is ",
      max(sizes)
    )
  }

  # Every module the reference has, and what became of it. Modules below
  # 10 mapped genes leave the analysis entirely, and modules whose genes never
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
      "fewer than 10 mapped genes"
    )
  )
  rownames(coverage) <- NULL

  n_dropped <- sum(!coverage$tested)
  if (n_dropped > 0L) {
    message(
      n_dropped, " of ", nrow(coverage), " reference modules were not ",
      "tested (", sum(coverage$size_mapped == 0L), " with no mapped ",
      "gene, ", sum(!coverage$tested & coverage$size_mapped > 0L),
      " below 10 mapped genes); see $coverage"
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
    ref_members
  )

  # Index by test gene, so a reference gene mapped to several test paralogs
  # contributes its value once per paralog.
  ref_kIM <- lapply(rows_by_mod, function(rows) rstats$kIM[loc_ref[rows] + 1L])  # nolint
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

  res <- .pres_run(
    net_test, as.integer(keep_test - 1L), test_members,
    ref_kIM, ref_cc, ref_mar, n_perm, as.integer(n_cores)
  )

  out <- .pres_assemble(
    res, modules_ref, tested, rows_by_mod, proj, map,
    n_perm, calibrate, seed
  )
  out$coverage <- coverage


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
  # once per module_preservation() and again per module_correspondence(), and
  # twice per contrast in preservation_paired().
  pick <- function(df) {
    if (nrow(df) == 0L) {
      return(NULL)
    }
    # gene2, then module, then gene1: the first row of each (gene2, module)
    # run is that module's smallest gene1, which is the deterministic pick.
    # Radix (C-locale) order, so the pick is the same on every machine.
    d <- df[order(df$gene2, df$module, df$gene1, method = "radix"), ,
      drop = FALSE
    ]
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
.pres_gene_stats <- function(net, keep, members) {
  a <- .net_cpp_args(net, net$threshold)
  if (.net_is_sparse(net)) {
    module_gene_stats_sparse_cpp(a$p, a$i, a$x, a$thr, keep, members, FALSE)
  } else {
    module_gene_stats_dense_cpp(a$net, a$thr, keep, members, FALSE)
  }
}


#' Preservation permutation engine, dense/sparse dispatch (internal)
#' @noRd
.pres_run <- function(net, keep, members, ref_kIM, ref_cc, ref_mar,  # nolint
                      n_perm, n_cores, store_perm = FALSE) {
  a <- .net_cpp_args(net, net$threshold)
  if (.net_is_sparse(net)) {
    module_preservation_sparse_cpp(
      a$p, a$i, a$x, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, FALSE, store_perm
    )
  } else {
    module_preservation_dense_cpp(
      a$net, a$thr, keep, members,
      ref_kIM, ref_cc, ref_mar, n_perm,
      n_cores, FALSE, store_perm
    )
  }
}


#' Build the result tables from the kernel output (internal)
#' @noRd
.pres_assemble <- function(res, modules_ref, tested, rows_by_mod, proj, map,
                           n_perm, calibrate, seed) {
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
  #   F(t) <= w00 * C(t, t) + (1 - w00) * t  # nolint
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
    w00 * res$p_npc + (1 - w00) * p_comb
  )
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
    p_value = p_comb,
    p_calibrated = p_cal,
    q_value = q_comb,
    evalue = sum(!is.na(p_cal)) * p_cal,
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
      n_perm = n_perm, calibrate = calibrate,
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



#' Classify modules as conserved, moderately preserved, or diverged
#'
#' Turns [module_preservation()] output into a call per module. The
#' combined permutation q-value drives the call.
#'
#' The `Zsummary` thresholds of Langfelder et al. (2011) split the
#' significant modules a second time, so the familiar 2 / 10 cut points
#' still appear.
#'
#' @section Criteria:
#' \describe{
#'   \item{conserved}{`q_value < 0.1` and `Zsummary_std` is at or above 10}
#'   \item{moderate}{`q_value < 0.1` and it is below 10}
#'   \item{diverged}{`q_value >= 0.1`}
#'   \item{untested}{`q_value` is `NA` -- a statistic could not be computed,
#'     so neither preservation nor divergence was measured}
#' }
#' A module whose degree correlation could not be computed -- `cor.degree` is
#' undefined when intramodular connectivity is constant in either network --
#' gets `NA` for that statistic, so the `pmax` combination and hence `q_value`
#' are `NA` too. Such a module is reported `"untested"` rather than
#' `"diverged"`: nothing was measured, so divergence would be a positive claim
#' the data does not support. It is never called preserved on the density
#' statistic alone.
#'
#' @section Ordering modules:
#' `Zsummary_std` decides the call, but it is not the only ordering worth
#' reading. It is a Z-score, so its denominator is a permutation standard
#' deviation, and that shrinks as a module grows: over 110 modules from 12
#' Pooideae contrasts the null sd of `cor.degree` tracks mapped module size
#' at Spearman -0.996 and that of `avg.weight` at -0.873. A large, weakly
#' preserved module can therefore outrank a small, strongly preserved one.
#'
#' `medianRank` is the size-independent complement Langfelder et al. (2011)
#' report next to `Zsummary` for exactly that reason. It ranks the observed
#' `avg.weight` and `cor.degree` across the tested modules, 1 = strongest,
#' and averages the two ranks; no permutation moment enters it, so module
#' size cannot set its scale. Prefer `medianRank` when asking which modules
#' are the best preserved *relative to each other*, and `Zsummary_std` when
#' asking how far from its own null any one module sits. The two are not
#' redundant: on those same 12 contrasts they ordered 11% of module pairs
#' differently.
#'
#' `medianRank` ranks across the modules of *one run*. Its scale is set by
#' how many modules that run tested, so it is not comparable between runs --
#' between the two directions of a contrast, or between contrasts in a
#' [preservation_paired()] table -- unless those runs tested equally many
#' modules. It takes no part in the classification.
#'
#' @param pres Output of [module_preservation()].
#' @param species Optional species label recorded in the `species` column.
#' @param pair_name Optional contrast label recorded in the `pair_name` column.
#'
#' @return A data frame with `module`, `species`, `pair_name`,
#'   `classification`, `Zsummary`, `Zsummary_std`, `medianRank`, `q_value`,
#'   `size` and `size_mapped`. `Zsummary_std` is the column the default
#'   criterion is read against, so reproducing
#'   the call from `Zsummary` alone will not match. `medianRank` is carried
#'   through from [module_preservation()] and takes no part in the call; see
#'   *Ordering modules* for when to read it instead.
#'
#' @examples
#' \dontrun{
#' cls <- classify_preservation(pres)
#' table(cls$classification)
#' }
#'
#' @export
classify_preservation <- function(pres, species = NA_character_,
                                  pair_name = NA_character_) {
  if (!is.list(pres) || is.null(pres$preservation)) {
    stop("pres must be output from module_preservation()")
  }
  p <- pres$preservation

  testable <- !is.na(p$q_value)
  significant <- testable & p$q_value < 0.1
  # Read the cut point on a scale where it means what it is meant to mean:
  # one null standard deviation. Zsummary is a mean of two standardized
  # statistics, so its own null sd is sqrt(2 + 2*rho)/2, not 1.
  z_used <- if ("Zsummary_std" %in% names(p)) {
    p$Zsummary_std
  } else {
    p$Zsummary
  }
  # Zsummary_std is NA when the null correlation of the two statistics
  # could not be estimated -- the kernel guards on `nb > 2`, so fewer
  # than three jointly scorable permutations, or one statistic constant
  # across them. Left alone that makes `strong`
  # FALSE and quietly demotes an otherwise significant module to
  # "moderate", which reads as a measurement rather than a missing
  # normaliser. Fall back to the raw scale for those rows and say so.
  # Both warnings below label themselves with the same context.
  where <- paste(stats::na.omit(c(species, pair_name)), collapse = " / ")
  # Only significant rows can be affected: a row with q >= 0.1 is
  # "diverged" whatever z_used says, so warning about it reports a
  # threshold-scale hazard where the fallback is inert.
  fell_back <- is.na(z_used) & !is.na(p$Zsummary) & significant
  if (any(fell_back)) {
    z_used[fell_back] <- p$Zsummary[fell_back]

    warning(
      sum(fell_back), " module(s) have no null correlation for ",
      "Zsummary_std, so the raw Zsummary was used against ",
      "the cut point of 10 for them",
      if (nzchar(where)) paste0(" [", where, "]") else "",
      "; the cut point means a different number of null standard ",
      "deviations there. Raise n_perm."
    )
  }
  strong <- !is.na(z_used) & z_used >= 10
  classification <- ifelse(!testable, "untested",
    ifelse(!significant, "diverged",
      ifelse(strong, "conserved", "moderate")
    )
  )
  if (any(!testable)) {
    warning(
      sum(!testable), " module(s) could not be tested (a statistic was ",
      "undefined); reported as \"untested\"",
      if (nzchar(where)) paste0(" [", where, "]") else "",
      ": ", paste(p$module[!testable], collapse = ", ")
    )
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
    # Carried, not recomputed: the ranks are over the tested modules, which
    # is exactly this table's rows. Guarded because on a preservation object
    # saved before medianRank existed the column is NULL, which data.frame()
    # drops without a word -- the caller would get the silent absence this
    # column was added to end.
    medianRank = if ("medianRank" %in% names(p)) {
      p$medianRank
    } else {
      rep(NA_real_, n)
    },
    q_value = p$q_value,
    size = p$size,
    size_mapped = p$size_mapped,
    stringsAsFactors = FALSE
  )
}


#' Match modules across species by ortholog overlap
#'
#' Cross-tabulates the modules of two species over a paralog-resolved
#' ortholog map. Tests each module pair for excess overlap.
#'
#' This answers "which module corresponds to which". That differs from
#' whether a module's topology is preserved ([module_preservation()]).
#' [classify_hub_conservation()] needs the correspondence, not the
#' preservation call.
#'
#' The map assigns one reference gene to each test-species gene, so the
#' hypergeometric's independence assumption holds. The multi-copy expansion
#' that made the same test anti-conservative under the retired
#' gene-overlap engine is gone here: a HOG with three paralogs no longer
#' contributes three correlated draws to the same urn.
#'
#' @param modules_ref,modules_test Module assignments for the two species,
#'   from [detect_modules()] or [as_modules()].
#' @param map Ortholog map from [resolve_ortholog_map()], with `gene1` in the
#'   reference species and `gene2` in the test species.
#' @param species_ref,species_test Optional species labels recorded on the
#'   result.
#'   `module1` belongs to `modules_ref`, and that orientation cannot be
#'   recovered from the table, so supplying these lets
#'   [classify_hub_conservation()] catch a transposed call.
#' @param seed Integer seed for the randomized-p draws behind pi0, or `NULL`
#'   (default) to draw from the ambient stream and leave it advanced. A seed
#'   draws from a private stream and restores the caller's on exit, the
#'   package-wide contract described under [detect_modules()].
#'
#' @return A list with `pairs` -- a data frame of `module1`,
#'   `module2`, `size1`, `size2`, `overlap`, `jaccard`, `p_value` and
#'   `q_value`, one row per module pair -- plus the `species_ref` /
#'   `species_test` labels. The column names match what
#'   [classify_hub_conservation()] expects.
#'
#' @examples
#' \dontrun{
#' map <- resolve_ortholog_map(
#'   ortho, rownames(net_a$network),
#'   rownames(net_b$network)
#' )
#' corr <- module_correspondence(mods_a, mods_b, map)
#' subset(corr$pairs, q_value < 0.1)
#' }
#'
#' @export
module_correspondence <- function(modules_ref, modules_test, map,
                                  species_ref = NULL, species_test = NULL,
                                  seed = NULL) {
  # The randomized pi0 draws uniforms. See R/rng.R.
  .seed_scope(seed)

  for (nm in c("modules_ref", "modules_test")) {
    m <- get(nm)
    if (!is.list(m) || is.null(m$modules) || is.null(m$module_genes)) {
      stop(
        nm, " must be a module assignment from detect_modules() ",
        "or as_modules()"
      )
    }
  }
  if (!is.data.frame(map) ||
        !all(c("gene1", "gene2", "source") %in% names(map))) {
    stop("map must be a data frame from resolve_ortholog_map()")
  }
  # character, not factor: a factor would index modules by its integer
  # codes, and radix-order by codes that carry another session's collation
  map$gene1 <- as.character(map$gene1)
  map$gene2 <- as.character(map$gene2)

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

  tab <- table(.radix_factor(proj$module), .radix_factor(proj$module_test))
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
    module1 = rownames(tab)[i],
    module2 = colnames(tab)[j],
    size1 = m,
    size2 = k,
    overlap = overlap,
    jaccard = overlap / (m + k - overlap),
    p_value = p_gt + p_eq,
    stringsAsFactors = FALSE
  )

  pairs$q_value <- if (nrow(pairs) < 2L) {
    pairs$p_value
  } else {
    compute_qvalues(
      pairs$p_value,
      p_rand_fn = function() p_gt + stats::runif(length(p_gt)) * p_eq,
      pi0_method = "randomized"
    )$qvalues
  }
  rownames(pairs) <- NULL

  # Orientation is not recoverable from the table, so record it: module1
  # belongs to modules_ref. classify_hub_conservation() checks it against the
  # key when present.
  list(pairs = pairs, species_ref = species_ref, species_test = species_test)
}


#' Run module preservation across many species pairs
#'
#' Applies [module_preservation()] to each contrast in `pairs`. Runs both
#' directions and reports them separately.
#'
#' Preservation is directional. Whether the modules of species A survive
#' in B differs from the reverse.
#'
#' @param modules Named list of module assignments ([detect_modules()] or
#'   [as_modules()]), keyed by species.
#' @param networks Named list of [compute_network()] results, keyed by species.
#' @param orthologs Data frame with columns `gene1`, `gene2`, `hog`.
#' @param pairs Data frame with columns `species1`, `species2` and optionally
#'   `pair_name`.
#' @param group Optional named vector mapping species to a clade. When
#'   supplied, a `group` column records `"conserved"` for preserved modules and
#'   the owning species' group for diverged ones. Modules reported
#'   `"untested"` are counted under `"untested"` rather than a clade:
#'   nothing was measured, so attributing them to one would overstate the
#'   evidence.
#' @param edges Optional [find_coexpressologs()] results, used for paralog
#'   resolution.
#' @param seed Integer seed for the whole run, or `NULL` (default) to draw
#'   from the ambient stream and leave it advanced. The seed is applied once
#'   here and the per-direction [module_preservation()] calls are left
#'   unseeded, so the directions draw in sequence from one stream instead of
#'   every one of them reusing the same permutations. A seeded call restores
#'   the caller's stream on exit, the package-wide contract described under
#'   [detect_modules()].
#'
#'   Before 0.3.0 `seed` reached [module_preservation()] through `...` and
#'   handed every direction the identical seed.
#' @param ... Further arguments passed to [module_preservation()].
#'
#' @return A list with `classification` (one row per module per direction,
#'   carrying `pair_name`, `module`, `species`, `reference`, `test`,
#'   `classification` and both ordering columns, `Zsummary_std` and
#'   `medianRank` -- the latter ranks within a direction, so it is not
#'   comparable across the rows of this table), `summary` (counts per
#'   contrast and direction) and `raw`
#'   (the [module_preservation()] results, keyed by `"<reference>.<test>"`).
#'
#' @examples
#' \dontrun{
#' res <- preservation_paired(mods, nets, ortho,
#'   pairs = data.frame(species1 = "BDIS", species2 = "BSYL"),
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
                                        seed = NULL, ...) {
  # Seeded once for the whole run; the module_preservation() calls below
  # leave seed at NULL and continue this stream, so the directions do not
  # all reuse the same permutations. See .seed_scope() in R/rng.R.
  .seed_scope(seed)

  if (!is.list(modules) || is.null(names(modules))) {
    stop("modules must be a named list keyed by species")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  if (!is.data.frame(pairs) ||
        !all(c("species1", "species2") %in% names(pairs))) {
    stop("pairs must have columns 'species1' and 'species2'")
  }
  if (any(pairs$species1 == pairs$species2)) {
    stop("pairs must not compare a species with itself")
  }
  # Both directions of every contrast are run, so results are keyed by
  # "<reference>.<test>". A repeated contrast -- in either orientation -- would
  # collide on that key and silently overwrite the earlier one.
  unordered <- paste(
    pmin(pairs$species1, pairs$species2), pmax(pairs$species1, pairs$species2)
  )
  if (anyDuplicated(unordered) > 0L) {
    stop(
      "pairs lists the same species pair more than once (both directions ",
      "of each contrast are run, so (A, B) and (B, A) are the same row): ",
      paste(unique(unordered[duplicated(unordered)]), collapse = ", ")
    )
  }
  species <- unique(c(pairs$species1, pairs$species2))
  missing_sp <- setdiff(species, intersect(names(modules), names(networks)))
  if (length(missing_sp) > 0L) {
    stop(
      "modules and networks must both cover: ",
      paste(missing_sp, collapse = ", ")
    )
  }
  if (!is.null(group)) {
    missing_grp <- setdiff(species, names(group))
    if (length(missing_grp) > 0L) {
      stop(
        "group missing entries for: ",
        paste(missing_grp, collapse = ", ")
      )
    }
  }
  if (!"pair_name" %in% names(pairs)) {
    pairs$pair_name <- paste(pairs$species1, pairs$species2, sep = ".")
  }

  # `raw` is keyed "<reference>.<test>" -- the public format this function
  # documents, and the one callers outside the package use (the vignette
  # reads raw[[paste(sp, collapse = ".")]]). A "." inside a species name
  # could make two different contrasts key the same: ref "A" / test "B.C"
  # and ref "A.B" / test "C" both give "A.B.C". Refuse that up front rather
  # than silently overwriting an entry -- the previous defence keyed on
  # "\x01" instead, which avoided the collision but made every key
  # unreadable to callers.
  keys <- c(
    paste(pairs$species1, pairs$species2, sep = "."),
    paste(pairs$species2, pairs$species1, sep = ".")
  )
  if (anyDuplicated(keys) > 0L) {
    stop(
      "species names give colliding contrast keys: ",
      paste(unique(keys[duplicated(keys)]), collapse = ", "),
      ". Rename the species so that '<reference>.<test>' is unique."
    )
  }

  raw <- list()
  class_list <- list()

  for (p in seq_len(nrow(pairs))) {
    for (direction in list(
      c(pairs$species1[p], pairs$species2[p]),
      c(pairs$species2[p], pairs$species1[p])
    )) {
      ref <- direction[1]
      test <- direction[2]
      # The public key, matching this function's documented return and
      # what callers outside the package index `raw` with. Collisions are
      # impossible here: they were rejected up front.
      key <- paste(ref, test, sep = ".")

      pres <- module_preservation(
        modules[[ref]], networks[[ref]], networks[[test]],
        .orient_orthologs(
          orthologs, rownames(networks[[ref]]$network),
          rownames(networks[[test]]$network)
        ),
        edges = edges,
        species_ref = ref, species_test = test, ...
      )
      raw[[key]] <- pres

      cls <- classify_preservation(pres,
        species = ref,
        pair_name = pairs$pair_name[p]
      )
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
            "conserved", as.character(group[ref])
          )
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
#' `gene1` / `gene2` hold gene identifiers, and the rest of the package
#' selects rows by testing them against each network's gene names. A table
#' written for the A -> B direction therefore yields nothing when B is the
#' reference, so the columns are swapped whenever that recovers more pairs.
#'
#' @noRd
.orient_orthologs <- function(orthologs, genes_ref, genes_test) {
  as_is <- sum(orthologs$gene1 %in% genes_ref &
                 orthologs$gene2 %in% genes_test)
  swapped <- sum(orthologs$gene2 %in% genes_ref &
                   orthologs$gene1 %in% genes_test)
  if (swapped > as_is) {
    orthologs[c("gene1", "gene2")] <-
      orthologs[c("gene2", "gene1")]
  }
  orthologs
}


# Factor with levels in radix (C-locale) order: split() and table() would
# otherwise sort character labels and gene IDs by the session's collation,
# which differs between machines for mixed-case or non-ASCII values
# (as_modules() allows any label). Digit labels sort the same either way.
.radix_factor <- function(x) {
  factor(x, levels = sort(unique(x), method = "radix"))
}
