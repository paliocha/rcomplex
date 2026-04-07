# Tests for clique_perturbation_test()
# Fixture: make_clique_fixture() from helper-clique-fixtures.R


test_that("zero noise gives 100% survival and jaccard = 1.0", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 5L, noise_sd = 0, seed = 42L)

  expect_true(all(result$survival_rate == 1.0))
  expect_true(all(result$mean_jaccard == 1.0))
  expect_true(all(result$n_matched == result$n_boot))
})


test_that("large noise reduces survival rate", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  # noise_sd=100 on networks with values ~3-10 destroys signal;
  # may trigger "no matching cliques" warning (expected)
  result <- suppressWarnings(clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 10L, noise_sd = 100, seed = 42L))

  expect_true(any(result$survival_rate < 1.0))
})


test_that("output has correct structure", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 3L, noise_sd = 0.1, seed = 1L)

  expect_true(is.data.frame(result))
  expect_true(all(c("clique_idx", "hog", "survival_rate",
                     "mean_jaccard", "n_boot", "n_matched") %in%
                    names(result)))

  # One row per baseline clique

  expect_equal(nrow(result), nrow(setup$cliques))

  # Types
  expect_type(result$clique_idx, "integer")
  expect_type(result$hog, "character")
  expect_type(result$survival_rate, "double")
  expect_type(result$mean_jaccard, "double")
  expect_type(result$n_boot, "integer")
  expect_type(result$n_matched, "integer")

  # clique_idx is 1-based sequential
  expect_equal(result$clique_idx, seq_len(nrow(setup$cliques)))

  # survival_rate in [0, 1]
  expect_true(all(result$survival_rate >= 0 & result$survival_rate <= 1))

  # mean_jaccard in [0, 1] or NA (when n_matched == 0)
  expect_true(all(is.na(result$mean_jaccard) |
                  (result$mean_jaccard >= 0 & result$mean_jaccard <= 1)))

  # n_boot matches input; n_matched <= n_boot
  expect_true(all(result$n_boot == 3L))
  expect_true(all(result$n_matched <= result$n_boot))
})


test_that("seed gives reproducible results", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  r1 <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 5L, noise_sd = 0.5, seed = 123L)

  r2 <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 5L, noise_sd = 0.5, seed = 123L)

  expect_equal(r1$survival_rate, r2$survival_rate)
  expect_equal(r1$mean_jaccard, r2$mean_jaccard)
})


test_that("empty cliques returns 0-row dataframe", {
  setup <- make_clique_fixture()
  empty <- setup$cliques[0, , drop = FALSE]

  result <- clique_perturbation_test(
    empty, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 5L, noise_sd = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
  expect_true(all(c("clique_idx", "hog", "survival_rate",
                     "mean_jaccard", "n_boot", "n_matched") %in%
                    names(result)))
})


test_that("n_boot = 0 returns 0-row dataframe", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 0L, noise_sd = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})


test_that("input validation: bad cliques", {
  setup <- make_clique_fixture()

  expect_error(
    clique_perturbation_test("bad", setup$target_species, setup$networks,
                              setup$orthologs),
    "cliques must be")

  expect_error(
    clique_perturbation_test(data.frame(x = 1), setup$target_species,
                              setup$networks, setup$orthologs),
    "cliques must be")
})


test_that("input validation: bad target_species", {
  setup <- make_clique_fixture()

  expect_error(
    clique_perturbation_test(setup$cliques, c("SP_A"), setup$networks,
                              setup$orthologs),
    "at least 2 species")
})


test_that("input validation: bad networks", {
  setup <- make_clique_fixture()

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species, list(),
                              setup$orthologs),
    "networks must be a named list")

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              list(SP_A = list(network = matrix(0), threshold = 1)),
                              setup$orthologs),
    "networks missing")

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              list(SP_A = list(), SP_B = list()),
                              setup$orthologs),
    "network.*threshold")
})


test_that("input validation: bad orthologs", {
  setup <- make_clique_fixture()

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              setup$networks, data.frame(x = 1)),
    "orthologs must have columns")
})


test_that("input validation: bad noise_sd", {
  setup <- make_clique_fixture()

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              setup$networks, setup$orthologs,
                              noise_sd = -1),
    "noise_sd must be a non-negative scalar")

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              setup$networks, setup$orthologs,
                              noise_sd = c(0.1, 0.2)),
    "noise_sd must be a non-negative scalar")

  expect_error(
    clique_perturbation_test(setup$cliques, setup$target_species,
                              setup$networks, setup$orthologs,
                              noise_sd = "bad"),
    "noise_sd must be a non-negative scalar")
})


test_that("hog column matches baseline cliques", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 3L, noise_sd = 0.1, seed = 1L)

  expect_equal(result$hog, setup$cliques$hog)
})


test_that("max_missing_edges is forwarded to find_cliques", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_perturbation_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_boot = 2L, noise_sd = 0, seed = 1L,
    max_missing_edges = 1L)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
})
