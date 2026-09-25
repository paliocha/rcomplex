# Density-integrated coexpressolog strength — design notes
#
# Problem: compare_neighborhoods() / summarize_comparison() answer a
# binary question at one analysis density ("is this ortholog pair's
# neighbourhood overlap significant at this threshold?"). Comparable
# methods do the same: Netotea et al. (2014, the base method this
# package extends), Curci et al. (2022, Plant Physiology, the DS1-DS5
# density-subnetwork design that motivated the density grid below),
# EVOTREE, Aro et al. and Van Zalen et al. (all cross-species
# co-expressolog/regulatory-core papers) report a call (or a single
# q-value) per edge at one chosen network density. None assigns a
# continuous, threshold-integrated confidence to an individual edge.
#
# This file implements that missing piece: re-examine one already
# estimated network at several matched densities (never recomputing
# correlation/normalization — see .net_check() / density_threshold()),
# rank each edge's evidence within each density, and integrate the
# ranks into a bounded, descriptive "strength" score. This is
# deliberately NOT combined with q-values as extra statistical
# evidence: densities share the same underlying estimate, so their
# p-values are not independent (see the package-level RNG/permutation
# notes on dependent p-values). Formal inference stays anchored at one
# explicit `reference_density`.
#
# The design is grounded in, rather than invented from scratch:
#
# - The "integrate a property over a sweep" aggregation is the
#   standard integrated-robustness construction from network science:
#   Albert, Jeong & Barabasi, "Error and attack tolerance of complex
#   networks", Nature 406, 2000 (robustness curves); Schneider,
#   Moreira, Andrade, Havlin & Herrmann, "Mitigation of malicious
#   attacks on networks", PNAS 108, 2011 (the R = mean(s(Q)) integrated
#   robustness metric this package's AUC-of-rank-vs-log-density
#   aggregation mirrors, applied per-edge instead of to the giant
#   component).
# - Anchoring a *single* reference density in data, rather than an
#   arbitrary round number, is inspired by (but is not a re-run of) the
#   WGCNA soft-thresholding stopping rule of Zhang & Horvath, "A
#   general framework for weighted gene coexpression network
#   analysis", Stat. Appl. Genet. Mol. Biol. 4, 2005, and Langfelder &
#   Horvath, "WGCNA: an R package for weighted correlation network
#   analysis", BMC Bioinformatics 9, 2008. That criterion picks a soft
#   power on a *fully connected, continuously weighted* adjacency
#   (a_ij = |cor_ij|^beta over every pair; nothing is ever thresholded
#   to zero), so its degree is a continuous row sum, not an edge
#   count. rcomplex's networks are hard-thresholded by construction
#   (a density cutoff, a sparse dgCMatrix store, a store_threshold
#   guard below which analysis is refused) — grafting true soft
#   thresholding onto that architecture would mean abandoning the
#   sparse-store contract that `.net_check()` / `store_threshold()`
#   exist to enforce, so it is not adopted here. What
#   suggest_reference_density() actually computes is a *hard-threshold*
#   degree-distribution scale-free check: at each candidate density,
#   count edges retained per gene (a binary cutoff, matching how
#   rcomplex networks are built and stored) and test that count
#   sequence for power-law-like decay. That is the older, harder-cut
#   lineage of scale-free-network diagnostics (Barabasi & Albert,
#   "Emergence of scaling in random networks", Science 286, 1999),
#   applied to the same MR/CLR-normalized, density-thresholded matrix
#   this package already builds — not a reproduction of Zhang &
#   Horvath's power-selection procedure. Only the "smallest cutoff
#   whose degree distribution clears an R-squared bar" *stopping rule*
#   is carried over from WGCNA; the underlying quantity it is applied
#   to is different, and the roxygen docs below say so explicitly.
# - The choice class ("small grid of matched densities" vs. a single
#   fixed cutoff) is benchmarked in Borate, Chesler, Langston, Saxton &
#   Voy, "Comparison of threshold selection methods for microarray
#   gene co-expression matrices", BMC Research Notes 2, 2009, and
#   updated in Bleker et al., "A Comparative Study of Gene
#   Co-Expression Thresholding Algorithms", J. Comput. Biol. 31, 2024.
# - Reporting a sweep across sparsification levels as the evaluation
#   itself, rather than tuning toward one level, follows Hassan &
#   Arifuzzaman, "Learning on Incomplete Graphs: Benchmarking GNN
#   Robustness to Sparsification", IEEE ICMLA, 2021, and the
#   structure-preserving sparsification-level guidance in Su, Meyer,
#   Kurths, Marwan & Meyer, "Generic network sparsification via
#   degree- and subgraph-based edge sampling", Information Sciences,
#   2024.
# - Species leave-k-out (clique_stability()) is an exhaustive jackknife
#   over a finite set of units, not a sampled replicate count; sampled
#   procedures elsewhere in the package (permutation, bootstrap) follow
#   the standard replicate-count guidance in Efron & Tibshirani, "An
#   Introduction to the Bootstrap", Chapman & Hall, 1993. The two
#   should not be conflated when deciding how many times to repeat a
#   randomized step.
#
# On mr_log_transform (checked empirically, not assumed): raw MR
# (default, `mr_log_transform = FALSE`) ranks ascending and the
# Obayashi & Kinoshita (2009) log-MR transform (`mr_log_transform =
# TRUE`) ranks descending before taking sqrt(rank_ij * rank_ji), so the
# two are NOT an exact monotonic reparameterization of each other
# despite both producing a "bigger value = stronger coexpression"
# network; on a 80-gene/50-sample simulated matrix at density = 0.05,
# they agreed on 98% of retained edges (Jaccard) and had value-ranking
# Spearman rho = 0.9999, with disagreement concentrated at the
# threshold boundary. Practical upshot: hold mr_log_transform fixed
# across every network compared in one density sweep or strength score
# (it already is, since it is one compute_network() argument shared by
# every species); do not assume the two modes are interchangeable. The
# residual boundary sensitivity is itself a small instance of exactly
# the threshold-sensitivity phenomenon this feature exists to quantify.
#
# On sample-level robustness (a still-missing axis): ATTED-II's own
# history is a worked example of the same "which axis handles which
# uncertainty" separation these notes insist on. It moved from a
# single-run PCC->MR (v4-v8, Obayashi et al. 2007-2016) to a 1000-fold
# bootstrap ensemble average of PCC->MR->logit over all samples (v9,
# Obayashi et al. 2018), to subagging over principal components -
# repeatedly subsampling a fraction of the informative PCs (~1000-2000
# kept out of ~15000, contribution decaying as a power law; 30% per
# draw in v10, 5% in v11) of a gene-centered sample-space PCA,
# reconstructing a pseudo-expression matrix from each draw, and
# averaging PCC->MR->logit results across draws (v10-v11, Obayashi,
# Hibara, Kagaya, Aoki & Kinoshita, "ATTED-II v11: a plant gene
# coexpression database using a sample balancing technique by
# subagging of principal components", Plant Cell Physiol. 63, 2022) -
# then, separately, standardizing the resulting coexpression as a
# z-score specifically to make it comparable across species and
# platforms (also v10-v11; carried into v12/13, atted.jp, 2024-2026,
# which further moved from Pearson to cosine similarity in PC-score
# space). Two points carry over directly: (1) subagging targets *which
# samples/conditions* dominate the correlation estimate - the
# sample-bootstrap axis the package design already keeps distinct from
# species leave-k-out and density integration, and rcomplex does not
# yet implement any version of it; (2) ATTED-II independently arrived
# at standardizing/matching before cross-species comparison (their
# z-score) for the same reason this file's density grid uses matched
# densities rather than a shared threshold multiplier - external,
# independent convergence on that specific design choice.
#
# Adding subagging to rcomplex is feasible but is a separate, upstream
# piece of work, not part of this file: it only touches correlation
# estimation (a new randomized cor_* backend in R/network.R, producing
# an ordinary correlation matrix that then flows through the existing
# MR/CLR/threshold/sparse pipeline unchanged), so it composes with
# density integration rather than competing with it. It would need: a
# seeded entry point following .seed_scope() (an ensemble over random
# PC subsamples is a new randomized step); a justified subsample
# fraction (ATTED-II itself changed 30% -> 5% between v10 and v11 with
# no closed-form derivation, so this is an empirically-tuned
# hyperparameter, not settled theory, and should be flagged as such
# rather than copied uncritically); and a real compute-cost budget (an
# ensemble of O(100-1000) correlation-matrix recomputations per
# network is a large multiplier on compute_network(), and would likely
# need the existing torch/OpenMP paths to be practical at genome
# scale). Track as a future roadmap item, e.g. `compute_network(...,
# robust = "subagging")`, rather than folding into
# coexpressolog_strength(), which assumes a single fixed correlation
# estimate and asks a different (threshold-robustness) question.
#
#' Reciprocal mid-p value from a compare_neighborhoods() row (internal)
#'
#' Uses the *exact* self-excluded hypergeometric quantities
#' `compare_neighborhoods()` already computes -- `Species{1,2}.p.val.gt`
#' (\eqn{P(X > x)}) and `Species{1,2}.p.val.eq` (\eqn{P(X = x)}) -- rather
#' than a normal approximation to the hypergeometric tail. This is the
#' same exact-hypergeometric formalism Tumminello et al. (2011, PLoS ONE
#' 6(3):e17994, the "statistically validated networks" method of Mantegna
#' and coauthors) use to validate one-mode-projected edges, and the same
#' exact upper-tail statistic already underlying every q-value this
#' package reports (see `compare_neighborhoods()`'s own documentation and
#' `compute_qvalues()`). A plain exact p-value `P(X >= x) = p.val.gt +
#' p.val.eq` is heavily tied at small counts (the same "piles up" problem
#' documented for pair-level q-values), which is unusable for a
#' *continuous* rank; the fix used elsewhere in this package for exactly
#' that problem is the randomized p-value `p.val.gt + U * p.val.eq` --
#' but that draws from the RNG, which the density profile must not do
#' (see file header: it must stay fully deterministic). The **mid-p
#' value** `p.val.gt + 0.5 * p.val.eq` is the standard deterministic
#' continuity correction for exact discrete tests (Lancaster 1961) --
#' the expectation of the randomized p-value, without drawing anything --
#' and is used here purely as a fine-grained, deterministic tie-breaker
#' underneath the formal `type` call (see `.coexpr_density_profile()`),
#' never as the significance call itself. Reciprocal mid-p is the *worse*
#' (larger, i.e. less significant) of the two directions, `pmax(p1, p2)`,
#' matching the package's reciprocal `pval_combine = "max"` convention:
#' both directions must independently look good for the edge to rank
#' well, and a small mid-p in only one direction cannot substitute for
#' the other.
#'
#' @noRd
.coexpr_reciprocal_midp <- function(comparison) {
  midp1 <- comparison$Species1.p.val.gt + 0.5 * comparison$Species1.p.val.eq
  midp2 <- comparison$Species2.p.val.gt + 0.5 * comparison$Species2.p.val.eq
  pmax(midp1, midp2)
}


#' Deterministic within-(pair, density) percentile rank (internal)
#'
#' Primary key is the reciprocal mid-p value (ascending: smaller/more
#' significant ranks stronger); ties are broken, in order, by combined
#' (geometric-mean) effect size, combined Jaccard, minimum directional
#' overlap, then gene identifiers -- so the result does not depend on row
#' order in `df`, only on its content. `normalized_rank` is `1` for the
#' strongest edge, `0` for the weakest, defined as `1` for every row when
#' there is exactly one (no spread to normalize against).
#'
#' @noRd
.coexpr_normalized_rank <- function(df) {
  n <- nrow(df)
  if (n == 0L) {
    df$normalized_rank <- numeric(0)
    return(df)
  }
  ord <- order(
    df$reciprocal_midp, -df$effect_size, -df$jaccard, -df$min_overlap,
    df$gene1, df$gene2,
    na.last = TRUE
  )
  rank <- integer(n)
  rank[ord] <- seq_len(n)
  df$normalized_rank <- if (n == 1L) 1 else 1 - (rank - 1) / (n - 1)
  df
}


#' One density's fully-validated comparison for a species pair (internal)
#'
#' Runs the *complete* statistically-validated-network procedure at this
#' single density -- [compare_neighborhoods()] for the exact hypergeometric
#' quantities, then [summarize_comparison()] for Benjamini-Hochberg
#' multiple-testing correction (the same correction family Tumminello et
#' al. 2011 apply as their FDR option), then [comparison_to_edges()] for
#' the alpha-thresholded `type` call -- exactly the machinery already used
#' for `reference`, not a cheaper substitute. `pi0_method = "storey"` (not
#' the package default `"randomized"`) is used here deliberately: Storey's
#' estimator is deterministic, so the density profile still draws nothing
#' from the RNG (see file header), while `reference` below uses the
#' package's usual randomized-p pi0 for its one formal, citable estimate.
#' Every tested ortholog pair is kept, including zero-overlap and `type
#' == "ns"` rows, so ranks and support fractions share one comparison
#' universe across the whole grid rather than a shrinking one.
#'
#' @noRd
.coexpr_validate_at_density <- function(net1, net2, orthologs, sp1, sp2,
                                        alpha, alternative, pval_combine,
                                        n_cores, pi0_method = "storey",
                                        seed = NULL) {
  comparison <- compare_neighborhoods(net1, net2, orthologs, n_cores)
  # filter_zero = FALSE: keep every tested ortholog pair, including
  # zero-overlap rows, so ranks and support fractions share one
  # comparison universe across the whole density grid (see file header
  # and .coexpr_density_profile() docs) rather than a shrinking one.
  summary_res <- summarize_comparison(comparison, alternative, alpha,
    filter_zero = FALSE, sp1 = sp1, sp2 = sp2,
    pi0_method = pi0_method, seed = seed
  )
  edges <- comparison_to_edges(summary_res$results, sp1, sp2, alternative,
    alpha,
    pval_combine = pval_combine
  )
  edges$reciprocal_midp <- .coexpr_reciprocal_midp(summary_res$results)
  edges$min_overlap <- pmin(
    summary_res$results$Species1.neigh.overlap,
    summary_res$results$Species2.neigh.overlap
  )
  edges
}


#' Density profile of one species pair's coexpressolog edges (internal)
#'
#' Re-derives the analysis threshold for each requested density from the
#' networks already built by [compute_network()] (never recomputing
#' correlation or normalization -- see `.net_density_threshold()`), fully
#' statistically validates the comparison at each density via
#' `.coexpr_validate_at_density()`, and ranks the result. Every tested
#' ortholog pair is kept at every density, including zero-overlap and
#' non-significant rows, so ranks share one comparison universe across the
#' grid rather than a shrinking one.
#'
#' @noRd
.coexpr_density_profile <- function(net1, net2, orthologs, sp1, sp2,
                                    densities, alpha, alternative,
                                    pval_combine, n_cores) {
  rows <- lapply(densities, function(d) {
    thr1 <- .net_density_threshold(net1, d)
    thr2 <- .net_density_threshold(net2, d)
    net1_d <- net1
    net1_d$threshold <- thr1
    net2_d <- net2
    net2_d$threshold <- thr2
    edges <- .coexpr_validate_at_density(
      net1_d, net2_d, orthologs, sp1, sp2, alpha, alternative, pval_combine,
      n_cores
    )
    edges <- .coexpr_normalized_rank(edges)
    edges$density <- d
    edges$threshold1 <- thr1
    edges$threshold2 <- thr2
    # A density is "supported" when its own independently BH-corrected
    # hypergeometric test calls the edge significant at this density
    # alone (Tumminello et al.'s validated-network criterion, run fresh
    # per density) -- never a reused reference-density decision. Recurrence
    # of this per-density call across the grid (not a pooled p-value) is
    # exactly the robustness measure of Curci et al. (2022); see file
    # header for both citations.
    edges$supported <- edges$type != "ns"
    edges
  })
  out <- do.call(rbind, rows)
  out[, c(
    "species1", "species2", "gene1", "gene2", "hog", "density",
    "threshold1", "threshold2", "q.value", "effect_size", "jaccard",
    "type", "reciprocal_midp", "min_overlap", "normalized_rank", "supported"
  ), drop = FALSE]
}


#' Area under normalized rank vs. log-density, scaled to [0, 1] (internal)
#'
#' Trapezoidal integration of `normalized_rank` against `log(density)`
#' over the ordered density grid, divided by the total log-density span so
#' the result is comparable across calls using different grids -- the same
#' "integrate a property over a sweep" construction as the robustness
#' curves of Albert, Jeong & Barabasi (2000) and the `R = mean(s(Q))`
#' metric of Schneider et al. (2011); see the file header. A single
#' density has no span to integrate over, so its own `normalized_rank` is
#' returned directly.
#'
#' @noRd
.coexpr_auc_log_density <- function(density, normalized_rank) {
  ord <- order(density)
  d <- density[ord]
  r <- normalized_rank[ord]
  n <- length(d)
  if (n == 1L) {
    return(r[1L])
  }
  ld <- log(d)
  span <- ld[n] - ld[1L]
  if (span <= 0) {
    return(mean(r))
  }
  trap <- sum((r[-1L] + r[-n]) / 2 * diff(ld))
  trap / span
}


#' Summarize one edge's density profile into the strength table (internal)
#' @noRd
.coexpr_summarize_edge <- function(density, normalized_rank, supported) {
  ord <- order(density)
  density <- density[ord]
  normalized_rank <- normalized_rank[ord]
  supported <- supported[ord]
  strength <- .coexpr_auc_log_density(density, normalized_rank)
  q <- stats::quantile(normalized_rank, c(0.25, 0.5, 0.75),
    na.rm = TRUE,
    names = FALSE, type = 7
  )
  strictest <- if (any(supported)) min(density[supported]) else NA_real_
  transitions <- if (length(supported) > 1L) {
    sum(diff(supported) != 0L)
  } else {
    0L
  }
  data.frame(
    strength = strength,
    median_rank = q[2L],
    rank_iqr = q[3L] - q[1L],
    density_support = mean(supported),
    strictest_density = strictest,
    n_transitions = transitions
  )
}


#' Density-integrated coexpressolog edge strength
#'
#' Re-examines each requested species pair's already-estimated networks at
#' several matched densities (`densities`; never recomputing correlation
#' or normalization), fully statistically validates the comparison at
#' each density -- exact hypergeometric p-values plus Benjamini-Hochberg
#' multiple-testing correction, the same procedure Tumminello et al.
#' (2011, PLoS ONE 6(3):e17994) use for "statistically validated
#' networks" -- ranks every tested ortholog pair within each density by a
#' deterministic exact-hypergeometric mid-p statistic, and integrates the
#' ranks into a bounded, descriptive strength score per edge. The
#' per-density recurrence of independently-validated significance calls
#' is the robustness measure of Curci et al. (2022); see the file header
#' (`R/coexpressolog-strength.R`) for both citations and the full design
#' rationale. This is a complement to, not a replacement for, the
#' pair-level q-values, effect sizes and Jaccard overlap already reported
#' by [find_coexpressologs()], the [density_sweep()] threshold sweep, and
#' the clique-level robustness analyses ([clique_persistence()],
#' [clique_threshold_sweep()], [clique_perturbation_test()]).
#'
#' @section Not circular:
#' Every density, including `reference_density`, is validated by the
#' *same* complete procedure ([compare_neighborhoods()] +
#' [summarize_comparison()] + [comparison_to_edges()]) -- there is no
#' separate, cheaper approximation for the profile grid. `strength` never
#' reads any density's `q.value`; it is built entirely from the
#' deterministic mid-p rank (`.coexpr_reciprocal_midp()`), a quantity
#' `reference`'s q-values do not use and do not feed. `density_support` at
#' each density is that density's own independently BH-corrected
#' significance call (`type != "ns"`), never a reused decision from
#' `reference_density` or from any other density in the grid; the fraction
#' of the grid supporting a call is reported purely descriptively (as
#' recurrence, following Curci et al. 2022), never pooled into a combined
#' p-value across the (mutually dependent) densities. `reference` is
#' additionally computed with the package's usual randomized-p pi0
#' estimator (`pi0_method = "randomized"`, needs `seed`); every other
#' density in the grid uses the deterministic Storey estimator
#' (`pi0_method = "storey"`) so the profile itself never draws from the
#' RNG.
#'
#' @param networks Named list of network objects (output of
#'   [compute_network()]), keyed by species.
#' @param orthologs Data frame with columns `Species1`, `Species2`, `hog`
#'   (output of [parse_orthologs()]).
#' @param densities Numeric vector of two or more matched densities in
#'   `(0, 1)`, e.g. `c(0.01, 0.02, 0.03, 0.05, 0.075)` (see
#'   `suggest_reference_density()`). Equivalent densities are comparable
#'   across species of different sizes; a shared threshold multiplier is
#'   not. Every value must be representable by every supplied network
#'   (see `.net_density_threshold()`/`store_density`).
#' @param reference_density A single density in `densities` at which
#'   ordinary q-value inference is run (via [compare_neighborhoods()] +
#'   [summarize_comparison()] + [comparison_to_edges()]). Not combined
#'   with the other densities -- dependent p-values across a density grid
#'   cannot be pooled into one q-value; see the file header.
#' @param species_pairs List of length-2 character vectors naming species
#'   pairs to compare, or `NULL` (default) for every pairwise combination
#'   of `names(networks)`, matching [find_coexpressologs()].
#' @param alpha Significance threshold for `reference`'s `type` column
#'   (default `0.05`).
#' @param alternative `"greater"` (conservation, default) or `"less"`
#'   (divergence); passed to `reference`'s q-value computation.
#' @param pval_combine `"max"` (default, reciprocal criterion) or `"min"`;
#'   passed to `reference`'s [comparison_to_edges()] call.
#' @param n_cores Number of threads for [compare_neighborhoods()] (default
#'   `1`).
#' @param seed Optional RNG seed for `reference`'s randomized-p pi0
#'   estimation; see [summarize_comparison()] and `.seed_scope()` in
#'   `R/rng.R`. The density profile itself draws nothing.
#'
#' @return A list with:
#'   \describe{
#'     \item{strength}{One row per species-pair coexpressolog edge:
#'       `species1`, `species2`, `gene1`, `gene2`, `hog`, `strength`
#'       (bounded descriptive ranking score, not a probability),
#'       `median_rank`, `rank_iqr`, `density_support`,
#'       `strictest_density`, `n_transitions`.}
#'     \item{profiles}{Long-form per-edge, per-density measurements:
#'       `density`, `threshold1`, `threshold2`, `q.value`, `effect_size`,
#'       `jaccard`, `type` (this density's own BH-corrected significance
#'       call), `reciprocal_midp` (deterministic exact-hypergeometric
#'       rank statistic), `min_overlap`, `normalized_rank`, `supported`
#'       (`type != "ns"`), plus the edge identity columns above.}
#'     \item{reference}{[comparison_to_edges()]-format results at
#'       `reference_density` alone, using the package's usual
#'       randomized-p pi0 estimator: `species1`, `species2`, `gene1`,
#'       `gene2`, `hog`, `q.value`, `effect_size`, `jaccard`, `type`.}
#'     \item{params}{List echoing `densities`, `reference_density`,
#'       `alpha`, `alternative`, `pval_combine` and `method =
#'       "hypergeometric"` (the only density-profile method implemented so
#'       far; HOG permutation inference stays reserved for `reference`).}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- coexpressolog_strength(networks, orthologs,
#'   densities = c(0.01, 0.02, 0.03, 0.05, 0.075),
#'   reference_density = 0.03
#' )
#' head(res$strength[order(-res$strength$strength), ])
#' }
#'
#' @param ... Passed to methods.
#'
#' @seealso [suggest_reference_density()] for a data-driven
#'   `reference_density`; [density_sweep()] for the pair/HOG-level q-value
#'   threshold sweep this complements.
#' @export
coexpressolog_strength <- function(networks, ...) {
  UseMethod("coexpressolog_strength")
}

#' @rdname coexpressolog_strength
#' @export
coexpressolog_strength.default <- function(networks, orthologs, densities,
                                           reference_density,
                                           species_pairs = NULL,
                                           alpha = 0.05,
                                           alternative = c("greater", "less"),
                                           pval_combine = c("max", "min"),
                                           n_cores = 1L, seed = NULL, ...) {
  alternative <- match.arg(alternative)
  pval_combine <- match.arg(pval_combine)
  .seed_scope(seed)

  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  if (length(networks) < 2) {
    stop("networks must contain at least 2 species")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  if (!is.numeric(densities) || length(densities) < 2L ||
        any(densities <= 0) || any(densities >= 1)) {
    stop("densities must be two or more numbers in (0, 1) (exclusive)")
  }
  densities <- sort(unique(densities))
  if (!is.numeric(reference_density) || length(reference_density) != 1L ||
        is.na(reference_density)) {
    stop("reference_density must equal one of the values in densities")
  }
  ref_tol <- .Machine$double.eps^0.5
  if (!isTRUE(any(abs(densities - reference_density) < ref_tol))) {
    stop("reference_density must equal one of the values in densities")
  }

  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(names(networks), 2, simplify = FALSE)
  }

  profile_list <- vector("list", length(species_pairs))
  reference_list <- vector("list", length(species_pairs))

  for (i in seq_along(species_pairs)) {
    sp_a <- species_pairs[[i]][1]
    sp_b <- species_pairs[[i]][2]
    if (!sp_a %in% names(networks)) {
      stop("species '", sp_a, "' not found in networks")
    }
    if (!sp_b %in% names(networks)) {
      stop("species '", sp_b, "' not found in networks")
    }

    profile_list[[i]] <- .coexpr_density_profile(
      networks[[sp_a]], networks[[sp_b]], orthologs, sp_a, sp_b,
      densities, alpha, alternative, pval_combine, n_cores
    )

    ref_thr_a <- .net_density_threshold(networks[[sp_a]], reference_density)
    ref_thr_b <- .net_density_threshold(networks[[sp_b]], reference_density)
    net_a_ref <- networks[[sp_a]]
    net_a_ref$threshold <- ref_thr_a
    net_b_ref <- networks[[sp_b]]
    net_b_ref$threshold <- ref_thr_b
    reference_list[[i]] <- .coexpr_validate_at_density(
      net_a_ref, net_b_ref, orthologs, sp_a, sp_b, alpha, alternative,
      pval_combine, n_cores,
      pi0_method = "randomized", seed = NULL
    )[, c(
      "species1", "species2", "gene1", "gene2", "hog",
      "q.value", "effect_size", "jaccard", "type"
    )]
  }

  profiles <- do.call(rbind, profile_list)
  reference <- do.call(rbind, reference_list)

  key <- paste(profiles$species1, profiles$species2, profiles$gene1,
    profiles$gene2, profiles$hog,
    sep = "\x01"
  )
  strength_rows <- lapply(split(seq_len(nrow(profiles)), key), function(idx) {
    p <- profiles[idx, , drop = FALSE]
    cbind(
      p[1L, c("species1", "species2", "gene1", "gene2", "hog"),
        drop = FALSE
      ],
      .coexpr_summarize_edge(p$density, p$normalized_rank, p$supported)
    )
  })
  strength <- do.call(rbind, strength_rows)
  rownames(strength) <- NULL

  list(
    strength = strength,
    profiles = profiles,
    reference = reference,
    params = list(
      densities = densities, reference_density = reference_density,
      method = "hypergeometric", alpha = alpha, alternative = alternative,
      pval_combine = pval_combine
    )
  )
}


#' Suggest a reference density from a hard-threshold scale-free diagnostic
#'
#' Borrows the WGCNA "smallest cutoff whose fit clears an R-squared bar"
#' *stopping rule* (Zhang & Horvath 2005), but applies it to a different
#' quantity: WGCNA's own scale-free fit index is defined on a fully
#' connected, continuously weighted adjacency (a soft-thresholding power
#' applied to every gene pair, nothing set to zero), which does not match
#' rcomplex's networks — those are hard-thresholded by construction (a
#' density cutoff, a sparse store, a store-density guard). This function
#' instead checks, at each candidate density, whether the *hard-cutoff*
#' degree distribution (edges retained per gene, a count, not a
#' continuous weighted sum) follows a power law — the older, harder-cut
#' lineage of scale-free-network diagnostics (Barabasi & Albert 1999),
#' evaluated on the same MR/CLR-normalized matrix this package already
#' builds. The output (a scale-free fit index, `r_squared`, the signed
#' R-squared of `log10(frequency) ~ log10(degree)` over binned degrees)
#' gives a data-driven way to anchor the single `reference_density` used
#' for formal inference (q-values, permutation tests), instead of picking
#' an arbitrary round number, while a wider density grid can still be
#' used for descriptive threshold-sensitivity summaries elsewhere.
#'
#' This is a diagnostic, not a guarantee: biological co-expression
#' networks only approximately follow scale-free topology, and a small
#' gene universe can leave every candidate density below the cutoff (see
#' `$recommended` below). Hold `mr_log_transform` fixed across every
#' network being compared: raw MR and log-MR rank genes in opposite
#' directions before combining them, so they are not an exact monotonic
#' reparameterization of one another (empirically close on simulated
#' data — Spearman rho > 0.999, ~98% edge-set agreement at one tested
#' density — but not identical, and the residual disagreement sits at
#' the threshold boundary this function is measuring sensitivity to).
#'
#' @param net A network object from [compute_network()] (dense or
#'   sparse). For a sparse network, every requested density must be no
#'   greater than `store_density` (the store cannot represent a lower
#'   threshold); see [as_sparse_network()].
#' @param densities Numeric vector of candidate densities in (0, 1),
#'   e.g. `c(0.01, 0.02, 0.03, 0.05, 0.075)`. Matched, explicit density
#'   values are preferred over threshold multipliers so candidates are
#'   comparable across networks of different sizes.
#' @param r_squared_cutoff Minimum signed scale-free fit index to accept
#'   a density as adequately scale-free (default `0.80`, the criterion
#'   used in Zhang & Horvath 2005).
#' @param n_bins Number of equal-width degree bins used to estimate the
#'   degree-frequency distribution at each density (default `10`, as in
#'   `WGCNA::pickSoftThreshold()`).
#'
#' @return A list with:
#'   \describe{
#'     \item{fit}{Data frame with one row per requested density:
#'       `density`, `threshold`, `mean_connectivity`, `r_squared`
#'       (signed scale-free fit index) and `slope`.}
#'     \item{recommended}{The smallest density whose `r_squared` is at
#'       least `r_squared_cutoff`, or `NA_real_` with a warning if none
#'       qualifies (all requested densities are still returned in
#'       `fit`, so the caller can inspect and decide manually).}
#'   }
#'
#' @examples
#' \dontrun{
#' net <- compute_network(x, density = 0.03, sparse = FALSE)
#' sft <- suggest_reference_density(net, c(0.01, 0.02, 0.03, 0.05, 0.075))
#' sft$recommended
#' }
#'
#' @export
suggest_reference_density <- function(net, densities,
                                      r_squared_cutoff = 0.80,
                                      n_bins = 10L) {
  if (!is.list(net) || is.null(net$network)) {
    stop("net must be a network object from compute_network()")
  }
  if (!is.numeric(densities) || length(densities) == 0L ||
        any(densities <= 0) || any(densities >= 1)) {
    stop("densities must be one or more numbers in (0, 1) (exclusive)")
  }
  densities <- sort(unique(densities))

  rows <- lapply(densities, function(d) {
    thr <- .net_density_threshold(net, d)
    k <- .net_degree_at_threshold(net, thr)
    sft <- .scale_free_fit_index(k, n_bins = n_bins)
    data.frame(
      density = d, threshold = thr, mean_connectivity = mean(k),
      r_squared = sft$r_squared, slope = sft$slope
    )
  })
  fit <- do.call(rbind, rows)

  ok <- which(fit$r_squared >= r_squared_cutoff)
  recommended <- if (length(ok) > 0L) {
    fit$density[ok[which.min(fit$density[ok])]]
  } else {
    warning(
      "no candidate density reached r_squared_cutoff = ",
      r_squared_cutoff, "; inspect $fit and choose manually"
    )
    NA_real_
  }

  list(fit = fit, recommended = recommended)
}

#' Analysis threshold for one density, from a network's stored matrix
#'
#' Mirrors `density_threshold_cpp()` exactly (same upper-triangle
#' population size and rounding), so dense and sparse networks built from
#' the same expression data must return identical thresholds. A sparse
#' network's `dgCMatrix` stores both triangles as exact duplicate pairs
#' (`m[i, j] == m[j, i]` bit-for-bit, since both entries are written from
#' the same symmetric source value in [compute_network()]), so the k-th
#' largest upper-triangle value is the `2k`-th largest entry of `m@x`.
#' Errors if the store does not hold enough entries for the requested
#' density (i.e. `density > store_density`); see `.net_check()`.
#' @noRd
.net_density_threshold <- function(net, density) {
  if (.net_is_sparse(net)) {
    m <- net$network
    n <- nrow(m)
    tri_size <- n * (n - 1) / 2
    # R's round() is round-half-to-even; density_threshold_cpp uses C++
    # std::round() (round-half-away-from-zero). They disagree exactly at
    # x.5 boundaries (e.g. density * tri_size == 148.5), so replicate
    # std::round()'s convention here rather than R's.
    k <- floor(density * tri_size + 0.5)
    if (k == 0) k <- 1
    if (k >= tri_size) k <- tri_size - 1
    n_top <- 2L * k
    if (n_top > length(m@x)) {
      stop(
        "store does not hold enough entries for density ", density,
        "; rebuild with a larger store_density"
      )
    }
    # kit::topn is a partial (quickselect-based) top-n, avoiding a full
    # O(n log n) sort of m@x when only the n_top-th largest value is
    # needed -- the same partial-selection idea as nth_element() in
    # density_threshold_cpp().
    thr <- kit::topn(m@x,
      n = n_top, decreasing = TRUE,
      hasna = FALSE, index = FALSE
    )[n_top]
    .net_check(net, thr)
    thr
  } else {
    density_threshold_cpp(net$network, density)
  }
}

#' Per-gene degree (neighbour count) at a given edge-weight threshold
#' @noRd
.net_degree_at_threshold <- function(net, thr) {
  m <- net$network
  if (.net_is_sparse(net)) {
    if (thr <= 0) {
      stop("threshold must be > 0 for degree counting on a sparse network")
    }
    Matrix::colSums(m >= thr)
  } else {
    mm <- m
    diag(mm) <- -Inf
    unname(colSums(mm >= thr))
  }
}

#' Signed scale-free fit index for a hard-threshold degree sequence
#'
#' Bins a degree vector into `n_bins` equal-width bins, regresses
#' `log10(frequency)` on `log10(bin midpoint)`, and returns the R-squared
#' signed by the slope (a proper power-law degree distribution has a
#' negative slope; the sign flags networks that trend the wrong way).
#' The binning/regression follows the same recipe as
#' `WGCNA::pickSoftThreshold()`'s `scaleFreeFitIndex()`, but `k` here is
#' a hard-cutoff edge count, not a continuous soft-thresholded weighted
#' degree; see `suggest_reference_density()`.
#'
#' @noRd
.scale_free_fit_index <- function(k, n_bins = 10L) {
  k <- k[k > 0]
  if (length(unique(k)) < 2L) {
    return(list(r_squared = NA_real_, slope = NA_real_))
  }
  breaks <- seq(min(k), max(k), length.out = n_bins + 1L)
  h <- graphics::hist(k, breaks = breaks, plot = FALSE)
  freq <- h$counts / sum(h$counts)
  mid <- h$mids
  keep <- freq > 0
  if (sum(keep) < 3L) {
    return(list(r_squared = NA_real_, slope = NA_real_))
  }
  df <- data.frame(x = log10(mid[keep]), y = log10(freq[keep]))
  fit <- stats::lm(y ~ x, df)
  s <- summary(fit)
  slope <- unname(stats::coef(fit)[2L])
  list(r_squared = sign(slope) * s$r.squared, slope = slope)
}
