# Tests for clique_intensity_test()

# Helper: build 2-species test setup with enough genes for significant results.
make_intensity_setup <- function() {
  n <- 20
  ga <- paste0("A", seq_len(n))
  gb <- paste0("B", seq_len(n))

  # SP_A network: A1 strongly connected to A2..A10
  mat_a <- matrix(0, n, n, dimnames = list(ga, ga))
  for (i in 2:10) {
    mat_a[1, i] <- mat_a[i, 1] <- 10  # strong
  }
  mat_a[1, 11] <- mat_a[11, 1] <- 3   # moderate
  for (i in 13:15) {
    mat_a[12, i] <- mat_a[i, 12] <- 4
  }

  # SP_B network: mirror structure
  mat_b <- matrix(0, n, n, dimnames = list(gb, gb))
  for (i in 2:10) {
    mat_b[1, i] <- mat_b[i, 1] <- 10
  }
  mat_b[1, 11] <- mat_b[11, 1] <- 3
  for (i in 13:15) {
    mat_b[12, i] <- mat_b[i, 12] <- 4
  }

  networks <- list(
    SP_A = list(network = mat_a, threshold = 2.0),
    SP_B = list(network = mat_b, threshold = 2.0)
  )

  # 1:1 orthologs
  orthologs <- data.frame(
    Species1 = ga, Species2 = gb,
    hog = paste0("HOG", seq_len(n)),
    stringsAsFactors = FALSE
  )

  edges <- find_coexpressologs(networks, orthologs, method = "analytical")
  baseline <- find_cliques(edges, c("SP_A", "SP_B"), min_species = 2L)

  list(networks = networks, orthologs = orthologs,
       cliques = baseline, target_species = c("SP_A", "SP_B"),
       edges = edges)
}

# Helper: 3-species setup
make_intensity_setup_3sp <- function() {
  n <- 15
  ga <- paste0("A", seq_len(n))
  gb <- paste0("B", seq_len(n))
  gc <- paste0("C", seq_len(n))

  make_net <- function(genes) {
    mat <- matrix(0, n, n, dimnames = list(genes, genes))
    for (i in 2:8) mat[1, i] <- mat[i, 1] <- 10
    mat[1, 9] <- mat[9, 1] <- 3
    mat
  }

  networks <- list(
    SP_A = list(network = make_net(ga), threshold = 2.0),
    SP_B = list(network = make_net(gb), threshold = 2.0),
    SP_C = list(network = make_net(gc), threshold = 2.0)
  )

  orthologs <- rbind(
    data.frame(Species1 = ga, Species2 = gb,
               hog = paste0("HOG", seq_len(n)), stringsAsFactors = FALSE),
    data.frame(Species1 = ga, Species2 = gc,
               hog = paste0("HOG", seq_len(n)), stringsAsFactors = FALSE),
    data.frame(Species1 = gb, Species2 = gc,
               hog = paste0("HOG", seq_len(n)), stringsAsFactors = FALSE)
  )

  target <- c("SP_A", "SP_B", "SP_C")
  edges <- find_coexpressologs(networks, orthologs, method = "analytical")
  baseline <- find_cliques(edges, target, min_species = 2L)

  list(networks = networks, orthologs = orthologs,
       cliques = baseline, target_species = target, edges = edges)
}


test_that("clique_intensity_test output has correct structure", {
  setup <- make_intensity_setup()
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
  setup <- make_intensity_setup()
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
  setup <- make_intensity_setup()
  empty_cliques <- setup$cliques[0, , drop = FALSE]
  result <- clique_intensity_test(
    empty_cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 3L)

  expect_equal(nrow(result), 0)
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test n_perm = 1 gives valid output", {
  setup <- make_intensity_setup()
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
  setup <- make_intensity_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 0L)
  expect_equal(nrow(result), 0)
})


test_that("clique_intensity_test observed_intensity matches compute_clique_edge_stats", {
  setup <- make_intensity_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 2L, seed = 42L)
  stats <- rcomplex:::compute_clique_edge_stats(
    setup$cliques, setup$edges, setup$target_species)
  expect_equal(result$observed_intensity, stats$intensity)
})


test_that("alternative = 'less' runs without error and returns valid p-values", {
  setup <- make_intensity_setup()
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
  setup <- make_intensity_setup_3sp()
  if (nrow(setup$cliques) == 0) skip("No 3-species cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, n_perm = 3L, seed = 42L, min_species = 2L)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test validates inputs", {
  setup <- make_intensity_setup()

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
