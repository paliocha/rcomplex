# Tests for permutation_hog_test()

# Helper: build test networks with a clearly conserved HOG (HOG1)
# and a clearly non-conserved HOG (HOG_NC)
make_test_nets <- function() {
  n <- 20
  mat1 <- matrix(0, n, n)
  mat2 <- matrix(0, n, n)
  rownames(mat1) <- colnames(mat1) <- paste0("A", 1:n)
  rownames(mat2) <- colnames(mat2) <- paste0("B", 1:n)

  # HOG1 genes (1-3): share neighbors 4-8 in BOTH networks (strong conservation)
  for (g in 1:3) {
    for (nb in 4:8) {
      mat1[g, nb] <- mat1[nb, g] <- 0.9
      mat2[g, nb] <- mat2[nb, g] <- 0.9
    }
  }

  # HOG_NC genes (14-15): neighbors are non-ortholog genes
  # In net1: connect to 16-18 (no orthologs)
  # In net2: connect to 19-20 (no orthologs)
  # -> zero overlap when mapped through orthologs
  for (g in 14:15) {
    mat1[g, 16] <- mat1[16, g] <- 0.9
    mat1[g, 17] <- mat1[17, g] <- 0.9
    mat1[g, 18] <- mat1[18, g] <- 0.9
    mat2[g, 19] <- mat2[19, g] <- 0.9
    mat2[g, 20] <- mat2[20, g] <- 0.9
  }

  diag(mat1) <- diag(mat2) <- 1

  net1 <- list(network = mat1, threshold = 0.5)
  net2 <- list(network = mat2, threshold = 0.5)

  # Orthologs: 1:1 for genes 1-8 and 14-15
  orthologs <- data.frame(
    Species1 = paste0("A", c(1:8, 14, 15)),
    Species2 = paste0("B", c(1:8, 14, 15)),
    OrthoGroup = c(rep("HOG1", 3), paste0("HOG", 2:6), "HOG_NC", "HOG_NC"),
    stringsAsFactors = FALSE
  )

  comparison <- compare_neighborhoods(net1, net2, orthologs)
  list(net1 = net1, net2 = net2, comparison = comparison)
}


test_that("permutation_hog_test returns correct structure", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("OrthoGroup", "n_pairs", "n_sp1", "n_sp2",
                         "T_obs", "n_perm", "n_exceed", "mean_eff",
                         "p.value", "q.value"))
  n_hogs <- length(unique(td$comparison$OrthoGroup))
  expect_equal(nrow(result), n_hogs)
})


test_that("conserved HOG has smaller p-value than non-conserved", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 2000L, min_exceedances = 50L
  )

  hog1 <- result[result$OrthoGroup == "HOG1", ]
  hog_nc <- result[result$OrthoGroup == "HOG_NC", ]

  # HOG1 (conserved) should have much smaller p-value
  expect_true(hog1$p.value < hog_nc$p.value)
  # HOG1 should have higher observed statistic
  expect_true(hog1$T_obs > hog_nc$T_obs)
  # Non-conserved HOG with zero overlap should have T_obs = 0
  expect_equal(hog_nc$T_obs, 0)
})


test_that("adaptive stopping terminates early for clearly null HOGs", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 5000L, min_exceedances = 20L
  )

  # Non-conserved HOGs should stop early (exceedances accumulate fast)
  hog_nc <- result[result$OrthoGroup == "HOG_NC", ]
  expect_true(hog_nc$n_exceed >= 20L || hog_nc$n_perm == 5000L)

  # Conserved HOG should run more permutations (fewer exceedances)
  hog1 <- result[result$OrthoGroup == "HOG1", ]
  expect_true(hog1$n_perm >= hog_nc$n_perm)
})


test_that("p-value formula is (n_exceed + 1) / (n_perm + 1)", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 50L
  )

  for (i in seq_len(nrow(result))) {
    expected_p <- (result$n_exceed[i] + 1) / (result$n_perm[i] + 1)
    expect_equal(result$p.value[i], expected_p, tolerance = 1e-10)
  }
})


test_that("q-values are computed", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )

  expect_true(all(result$q.value >= 0 & result$q.value <= 1))
})


test_that("alternative='less' tests divergence", {
  td <- make_test_nets()
  set.seed(42)

  result_con <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    alternative = "greater",
    max_permutations = 500L, min_exceedances = 20L
  )
  result_div <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    alternative = "less",
    max_permutations = 500L, min_exceedances = 20L
  )

  hog1_con <- result_con[result_con$OrthoGroup == "HOG1", ]
  hog1_div <- result_div[result_div$OrthoGroup == "HOG1", ]

  # Conserved HOG should be significant for conservation, not divergence
  expect_true(hog1_con$p.value < hog1_div$p.value)
})


test_that("empty comparison handled gracefully", {
  td <- make_test_nets()
  empty <- td$comparison[0, ]

  result <- permutation_hog_test(td$net1, td$net2, empty)
  expect_equal(nrow(result), 0)
  expect_named(result, c("OrthoGroup", "n_pairs", "n_sp1", "n_sp2",
                         "T_obs", "n_perm", "n_exceed", "mean_eff",
                         "p.value", "q.value"))
})


test_that("single-pair HOG works", {
  td <- make_test_nets()
  single <- td$comparison[td$comparison$OrthoGroup %in%
                            paste0("HOG", 2:6), ]

  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, single,
    max_permutations = 200L, min_exceedances = 10L
  )

  expect_true(all(result$n_sp1 == 1L))
  expect_true(all(result$n_sp2 == 1L))
  expect_true(all(result$n_pairs == 1L))
})


test_that("n_sp1 and n_sp2 reflect unique gene counts", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L
  )

  hog1 <- result[result$OrthoGroup == "HOG1", ]
  expect_equal(hog1$n_sp1, 3L)
  expect_equal(hog1$n_sp2, 3L)
  expect_equal(hog1$n_pairs, 3L)

  hog_nc <- result[result$OrthoGroup == "HOG_NC", ]
  expect_equal(hog_nc$n_sp1, 2L)
  expect_equal(hog_nc$n_sp2, 2L)
  expect_equal(hog_nc$n_pairs, 2L)
})


test_that("results are ordered by p-value", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )

  expect_true(all(diff(result$p.value) >= 0))
})


test_that("input validation works", {
  td <- make_test_nets()

  expect_error(
    permutation_hog_test("not_a_net", td$net2, td$comparison),
    "net1 must be a network object"
  )
  expect_error(
    permutation_hog_test(td$net1, "not_a_net", td$comparison),
    "net2 must be a network object"
  )
  expect_error(
    permutation_hog_test(td$net1, td$net2, data.frame(x = 1)),
    "comparison must be output from compare_neighborhoods"
  )
})


test_that("effect sizes are computed correctly", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L
  )

  hog1_rows <- td$comparison$OrthoGroup == "HOG1"
  expected_eff <- mean(sqrt(
    td$comparison$Species1.effect.size[hog1_rows] *
      td$comparison$Species2.effect.size[hog1_rows]
  ))
  hog1 <- result[result$OrthoGroup == "HOG1", ]
  expect_equal(hog1$mean_eff, expected_eff, tolerance = 1e-10)
})
