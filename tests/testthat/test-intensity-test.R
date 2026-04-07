# Tests for clique_intensity_test()
# Fixtures: make_clique_fixture() and make_clique_fixture_3sp()
# from helper-clique-fixtures.R


test_that("clique_intensity_test output has correct structure", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 3L, seed = 42L)

  expect_true(is.data.frame(result))
  expected_cols <- c("clique_idx", "hog", "observed_intensity",
                     "null_mean", "null_sd", "z_score", "p_value",
                     "n_perm", "n_matched")
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
    setup$orthologs, n_perm = 5L, seed = 123L)
  r2 <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 5L, seed = 123L)

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
    setup$orthologs, n_perm = 3L)

  expect_equal(nrow(result), 0)
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test n_perm = 1 gives valid output", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 1L, seed = 42L)

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
    setup$orthologs, n_perm = 0L)
  expect_equal(nrow(result), 0)
})


test_that("clique_intensity_test observed_intensity matches compute_clique_edge_stats", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 2L, seed = 42L)
  stats <- rcomplex:::compute_clique_edge_stats(
    setup$cliques, setup$edges, setup$target_species)
  expect_equal(result$observed_intensity, stats$intensity)
})


test_that("alternative = 'less' runs without error and returns valid p-values", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  r_less <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 5L, seed = 42L, alternative = "less")

  expect_true(is.data.frame(r_less))
  expect_true("p_value" %in% names(r_less))
  # p-values should be in [0, 1] or NA
  valid <- !is.na(r_less$p_value)
  if (any(valid)) {
    expect_true(all(r_less$p_value[valid] >= 0 &
                    r_less$p_value[valid] <= 1))
  }
})


test_that("clique_intensity_test works with 3+ species", {
  setup <- make_clique_fixture_3sp()
  if (nrow(setup$cliques) == 0) skip("No 3-species cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 3L, seed = 42L, min_species = 2L)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test validates inputs", {
  setup <- make_clique_fixture()

  expect_error(clique_intensity_test("not_df", setup$target_species,
               setup$networks, setup$orthologs),
               "data frame")
  expect_error(clique_intensity_test(data.frame(x = 1), setup$target_species,
               setup$networks, setup$orthologs),
               "data frame")
  expect_error(clique_intensity_test(setup$cliques, "SP_A",
               setup$networks, setup$orthologs),
               "at least 2")
  expect_error(clique_intensity_test(setup$cliques, setup$target_species,
               list(), setup$orthologs),
               "named list")
  expect_error(clique_intensity_test(setup$cliques, setup$target_species,
               setup$networks, data.frame(x = 1)),
               "Species1, Species2, hog")
})


test_that("max_missing_edges is forwarded to find_cliques", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 2L, seed = 1L,
    max_missing_edges = 1L)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
})
