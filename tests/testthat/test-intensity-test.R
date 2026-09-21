# Tests for clique_intensity_test()
# Fixtures: make_clique_fixture() and make_clique_fixture_3sp()
# from helper-clique-fixtures.R


test_that("clique_intensity_test output has correct structure", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L, seed = 42L
  )

  expect_true(is.data.frame(result))
  expected_cols <- c(
    "clique_idx", "hog", "observed_intensity",
    "null_mean", "null_sd", "z_score", "p_value",
    "n_perm", "n_matched"
  )
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true(all(result$clique_idx >= 1L))
  expect_true(all(result$n_perm == 3L))
  expect_true(all(result$n_matched >= 0L & result$n_matched <= 3L))
  expect_true(is.numeric(result$observed_intensity))
  expect_true(is.numeric(result$null_mean))
  expect_true(is.numeric(result$null_sd))
})


test_that("clique_intensity_test seed produces reproducible results", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  r1 <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 5L, seed = 123L
  )
  r2 <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 5L, seed = 123L
  )

  expect_equal(r1$null_mean, r2$null_mean)
  expect_equal(r1$null_sd, r2$null_sd)
  expect_equal(r1$z_score, r2$z_score)
  expect_equal(r1$p_value, r2$p_value)
  expect_equal(r1$n_matched, r2$n_matched)
})


test_that("clique_intensity_test empty cliques returns 0-row dataframe", {
  setup <- make_clique_fixture()
  empty_cliques <- setup$cliques[0, , drop = FALSE]
  result <- clique_intensity_test(
    empty_cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L
  )

  expect_equal(nrow(result), 0)
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test n_perm = 1 gives valid output", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 1L, seed = 42L
  )

  expect_equal(nrow(result), nrow(setup$cliques))
  expect_equal(result$n_perm[1], 1L)
  # With n_perm=1, sd is NA (need >= 2 values), so z_score should be NA
  expect_true(all(is.na(result$z_score) | is.finite(result$z_score)))
})


test_that("clique_intensity_test n_perm = 0 returns empty", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 0L
  )
  expect_equal(nrow(result), 0)
})


test_that(
  "clique_intensity_test observed_intensity matches compute_clique_edge_stats",
  {
    setup <- make_clique_fixture()
    if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

    result <- clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 2L, seed = 42L
    )
    stats <- rcomplex:::compute_clique_edge_stats(
      setup$cliques, setup$edges, setup$target_species
    )
    expect_equal(result$observed_intensity, stats$intensity)
  }
)


test_that(
  "alternative = 'less' runs without error and returns valid p-values",
  {
    setup <- make_clique_fixture()
    if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

    r_less <- clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 5L, seed = 42L, alternative = "less"
    )

    expect_true(is.data.frame(r_less))
    expect_true("p_value" %in% names(r_less))
    # p-values should be in [0, 1] or NA
    valid <- !is.na(r_less$p_value)
    if (any(valid)) {
      expect_true(all(r_less$p_value[valid] >= 0 &
                        r_less$p_value[valid] <= 1))
    }
  }
)


test_that("clique_intensity_test works with 3+ species", {
  setup <- make_clique_fixture_3sp()
  if (nrow(setup$cliques) == 0) skip("No 3-species cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L, seed = 42L, min_species = 2L
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test validates inputs", {
  setup <- make_clique_fixture()

  expect_error(
    clique_intensity_test(
      "not_df", setup$target_species,
      setup$networks, setup$orthologs
    ),
    "data frame"
  )
  expect_error(
    clique_intensity_test(
      data.frame(x = 1), setup$target_species,
      setup$networks, setup$orthologs
    ),
    "data frame"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, "SP_A",
      setup$networks, setup$orthologs
    ),
    "at least 2"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, setup$target_species,
      list(), setup$orthologs
    ),
    "named list"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, setup$target_species,
      setup$networks, data.frame(x = 1)
    ),
    "Species1, Species2, hog"
  )
})


test_that("max_missing_edges is forwarded to find_cliques", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 2L, seed = 1L,
    max_missing_edges = 1L
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
})


test_that("pval_combine/pi0_method reach the baseline and null reruns", {
  # Baseline cliques built with pval_combine = "min" (one conserved edge).
  # Observed intensity is a Jaccard-index percentile (#11), and Jaccard
  # does not depend on how directional q-values are combined, so it cannot
  # show whether pval_combine was forwarded: it matches under either value.
  # The null reruns can. Under "min" a permutation that keeps the A1-B1
  # hub mapping rebuilds the clique and counts as a match; under "max" the
  # diluted B-side direction never reaches alpha, so no permutation can
  # produce a conserved edge and nothing matches.
  setup <- make_asym_clique_fixture()
  expect_equal(nrow(setup$cliques), 1L)

  run <- function(pc) {
    clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 20L, seed = 7L,
      pval_combine = pc, pi0_method = "none"
    )
  }
  r_min <- run("min")
  r_max <- run("max")

  stats_min <- rcomplex:::compute_clique_edge_stats(
    setup$cliques, setup$edges_min, setup$target_species
  )
  expect_equal(r_min$observed_intensity, stats_min$intensity)
  expect_gt(r_min$n_matched, 0L)
  expect_identical(r_max$n_matched, 0L)
})


test_that("a pre-filtered edges argument warns", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")
  conserved <- setup$edges[setup$edges$type == "conserved", , drop = FALSE]
  expect_gt(nrow(conserved), 0L)

  expect_warning(
    clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 1L, seed = 1L, edges = conserved
    ),
    class = "rcomplex_prefiltered_edges"
  )
})



test_that("null intensities are the permuted runs' own Jaccard weights", {
  # A calibration test (strong fixture rejects, null fixture does not) is
  # not feasible at this scale: the null reshuffles every ortholog, so
  # permuted runs rebuild a baseline clique only by chance and small
  # fixtures match none. What can be pinned exactly is the null itself.
  # With pi0_method = "none" the only draw in a permutation is the
  # ortholog shuffle, so replaying the loop under the same seed rebuilds
  # every permuted table; null_mean must equal the mean intensity of the
  # matched permuted cliques, weighted over each permuted run's own table.
  setup <- make_asym_clique_fixture()
  sp <- setup$target_species
  n_perm <- 20L
  res <- clique_intensity_test(
    setup$cliques, sp, setup$networks, setup$orthologs,
    n_perm = n_perm, seed = 7L,
    pval_combine = "min", pi0_method = "none"
  )
  expect_gt(res$n_matched, 0L)

  set.seed(7L)
  null_int <- numeric(0)
  for (b in seq_len(n_perm)) {
    sh <- setup$orthologs
    sh$Species2 <- sample(sh$Species2)
    e_p <- find_coexpressologs(setup$networks, sh,
      method = "analytical", pval_combine = "min", pi0_method = "none"
    )
    if (nrow(e_p) == 0L) next
    cl_p <- find_cliques(e_p, sp, min_species = length(sp))
    hit <- cl_p[cl_p$hog == setup$cliques$hog[1L], , drop = FALSE]
    if (nrow(hit) == 0L) next
    st <- rcomplex:::compute_clique_edge_stats(hit[1L, ], e_p, sp)
    if (!is.na(st$intensity)) null_int <- c(null_int, st$intensity)
  }
  expect_length(null_int, res$n_matched)
  expect_equal(res$null_mean, mean(null_int), tolerance = 1e-12)
})


test_that("within_hog null matches cliques where global cannot", {
  # The global shuffle moves genes out of their own HOG, so a permuted
  # run almost never rebuilds the observed clique's HOG and nothing can
  # be matched -- measured as 0 of 204 cliques on the eight-species
  # Pooideae run, which left every z undefined. The within-HOG shuffle
  # keeps the grouping, so the same cliques do match.
  #
  # make_clique_fixture() cannot show this: it puts exactly one gene in
  # every HOG, which makes the within-HOG shuffle a no-op and the null a
  # point mass. This fixture pairs genes two to a HOG instead, so the
  # shuffle has something to permute.
  n <- 20L
  ga <- paste0("A", seq_len(n))
  gb <- paste0("B", seq_len(n))
  make_net <- function(genes) {
    m <- matrix(0, n, n, dimnames = list(genes, genes))
    for (i in 2:10) m[1, i] <- m[i, 1] <- 10
    m[1, 11] <- m[11, 1] <- 3
    for (i in 13:15) m[12, i] <- m[i, 12] <- 4
    m
  }
  networks <- list(
    SP_A = list(network = make_net(ga), threshold = 2),
    SP_B = list(network = make_net(gb), threshold = 2)
  )
  sp <- c("SP_A", "SP_B")
  orthologs <- data.frame(
    Species1 = ga, Species2 = gb,
    hog = paste0("HOG", ceiling(seq_len(n) / 2)),
    stringsAsFactors = FALSE
  )
  edges <- find_coexpressologs(networks, orthologs,
    method = "analytical", pi0_method = "storey"
  )
  cliques <- find_cliques(edges, sp, min_species = 2L)
  expect_gt(nrow(cliques), 0L)

  run <- function(nm) {
    clique_intensity_test(cliques, sp, networks, orthologs,
      n_perm = 30L, seed = 11L, null_model = nm, pi0_method = "storey"
    )
  }
  global <- run("global")
  within <- run("within_hog")

  # global destroys the HOGs: nothing matches, so no z is computable
  expect_true(all(global$n_matched == 0L))
  expect_true(all(is.na(global$z_score)))
  # within_hog keeps them: every clique matches and gets a usable z
  expect_true(all(within$n_matched > 0L))
  expect_true(all(!is.na(within$z_score)))
})


test_that("null_model is validated and defaults to global", {
  fx <- formals(rcomplex:::clique_intensity_test.default)
  expect_identical(eval(fx$null_model)[1], "global")
  expect_setequal(eval(fx$null_model), c("global", "within_hog"))
})
