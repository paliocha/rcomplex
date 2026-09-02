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
    hog = c(rep("HOG1", 3), paste0("HOG", 2:6), "HOG_NC", "HOG_NC"),
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
  expect_named(result, c("hog", "n_pairs", "n_sp1", "n_sp2",
                         "T_obs", "n_perm", "n_exceed", "mean_eff",
                         "p.value", "q.value"))
  n_hogs <- length(unique(td$comparison$hog))
  expect_equal(nrow(result), n_hogs)
})


test_that("conserved HOG has smaller p-value than non-conserved", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 2000L, min_exceedances = 50L
  )

  hog1 <- result[result$hog == "HOG1", ]
  hog_nc <- result[result$hog == "HOG_NC", ]

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

  # Non-conserved HOGs: T_obs=0 skips permutations entirely, or
  # exceedances accumulate fast and trigger early stopping
  hog_nc <- result[result$hog == "HOG_NC", ]
  expect_true(hog_nc$n_perm == 0L || hog_nc$n_exceed >= 20L)

  # Conserved HOG should run more permutations (fewer exceedances)
  hog1 <- result[result$hog == "HOG1", ]
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

  hog1_con <- result_con[result_con$hog == "HOG1", ]
  hog1_div <- result_div[result_div$hog == "HOG1", ]

  # Conserved HOG should be significant for conservation, not divergence
  expect_true(hog1_con$p.value < hog1_div$p.value)

  # Non-conserved HOG with T_obs=0: under "less", T_obs=0 is extreme
  # divergence — permutations must actually run (not be short-circuited)
  hog_nc_div <- result_div[result_div$hog == "HOG_NC", ]
  expect_true(hog_nc_div$n_perm > 0L)
})


test_that("empty comparison handled gracefully", {
  td <- make_test_nets()
  empty <- td$comparison[0, ]

  result <- permutation_hog_test(td$net1, td$net2, empty)
  expect_equal(nrow(result), 0)
  expect_named(result, c("hog", "n_pairs", "n_sp1", "n_sp2",
                         "T_obs", "n_perm", "n_exceed", "mean_eff",
                         "p.value", "q.value"))
})


test_that("single-pair HOG works", {
  td <- make_test_nets()
  single <- td$comparison[td$comparison$hog %in%
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

  hog1 <- result[result$hog == "HOG1", ]
  expect_equal(hog1$n_sp1, 3L)
  expect_equal(hog1$n_sp2, 3L)
  expect_equal(hog1$n_pairs, 3L)

  hog_nc <- result[result$hog == "HOG_NC", ]
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
    "comparison missing required columns"
  )
})


test_that("effect sizes are computed correctly", {
  td <- make_test_nets()
  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L
  )

  hog1_rows <- td$comparison$hog == "HOG1"
  expected_eff <- mean(sqrt(
    td$comparison$Species1.effect.size[hog1_rows] *
      td$comparison$Species2.effect.size[hog1_rows]
  ))
  hog1 <- result[result$hog == "HOG1", ]
  expect_equal(hog1$mean_eff, expected_eff, tolerance = 1e-10)
})


# ---- torch backend tests ----

test_that("use_torch errors when torch not installed", {
  skip_if(requireNamespace("torch", quietly = TRUE),
          "torch is installed — cannot test missing-package error")
  td <- make_test_nets()
  expect_error(
    permutation_hog_test(td$net1, td$net2, td$comparison,
                         use_torch = TRUE),
    "requires the torch package"
  )
})


test_that("torch backend T_obs matches bit-vector backend", {
  skip_if_not_installed("torch")
  skip_if_not(tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE),
              "torch backend (Lantern) not available")
  td <- make_test_nets()

  set.seed(42)
  result_bv <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L
  )
  set.seed(42)
  result_fe <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L,
    use_torch = TRUE
  )

  # T_obs should match closely — MPS uses float32, so ~1e-6 differences
  bv_order <- order(result_bv$hog)
  fe_order <- order(result_fe$hog)
  expect_equal(result_fe$T_obs[fe_order], result_bv$T_obs[bv_order],
               tolerance = 1e-5)
})


test_that("torch backend conserved vs non-conserved HOGs", {
  skip_if_not_installed("torch")
  skip_if_not(tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE),
              "torch backend (Lantern) not available")
  td <- make_test_nets()

  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 2000L, min_exceedances = 50L,
    use_torch = TRUE
  )

  hog1 <- result[result$hog == "HOG1", ]
  hog_nc <- result[result$hog == "HOG_NC", ]

  expect_true(hog1$p.value < hog_nc$p.value)
  expect_true(hog1$T_obs > hog_nc$T_obs)
  expect_equal(hog_nc$T_obs, 0)
})


test_that("torch backend p-value formula is correct", {
  skip_if_not_installed("torch")
  skip_if_not(tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE),
              "torch backend (Lantern) not available")
  td <- make_test_nets()

  set.seed(42)
  result <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 50L,
    use_torch = TRUE
  )

  for (i in seq_len(nrow(result))) {
    expected_p <- (result$n_exceed[i] + 1) / (result$n_perm[i] + 1)
    expect_equal(result$p.value[i], expected_p, tolerance = 1e-10)
  }
})


# ---- sparse (dgCMatrix) network dispatch ----

test_that("sparse permutation_hog_test equals dense (seeded)", {
  td <- make_test_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)

  set.seed(42)
  res_d <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )
  set.seed(42)
  res_s <- permutation_hog_test(
    net1_s, net2_s, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )
  expect_equal(res_s, res_d)

  # divergence direction too
  set.seed(7)
  res_d <- permutation_hog_test(
    td$net1, td$net2, td$comparison, alternative = "less",
    max_permutations = 200L, min_exceedances = 10L
  )
  set.seed(7)
  res_s <- permutation_hog_test(
    net1_s, net2_s, td$comparison, alternative = "less",
    max_permutations = 200L, min_exceedances = 10L
  )
  expect_equal(res_s, res_d)
})


test_that("sparse permutation_hog_test with tighter threshold equals dense", {
  # graded weights 10 / 7 / 4 (see make_graded_nets()): store threshold 5,
  # analysis threshold 8 keeps tier 10 and drops the stored tier 7, so
  # T_obs > 0 for the conserved hub HOG is structural, not seed-dependent
  td <- make_graded_nets()
  comparison <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  thr <- 8
  net1_s <- modifyList(sparse_net(td$net1), list(threshold = thr))
  net2_s <- modifyList(sparse_net(td$net2), list(threshold = thr))
  net1_d <- modifyList(td$net1, list(threshold = thr))
  net2_d <- modifyList(td$net2, list(threshold = thr))

  expect_gt(sum(net1_s$network@x >= thr), 0L)
  expect_lt(sum(net1_s$network@x >= thr), length(net1_s$network@x))

  set.seed(11)
  res_d <- permutation_hog_test(
    net1_d, net2_d, comparison,
    max_permutations = 200L, min_exceedances = 10L
  )
  set.seed(11)
  res_s <- permutation_hog_test(
    net1_s, net2_s, comparison,
    max_permutations = 200L, min_exceedances = 10L
  )
  expect_equal(res_s, res_d)

  # partially filtered store must actually be exercised
  expect_gt(res_d$T_obs[res_d$hog == "HOG1"], 0)
  expect_equal(res_d$T_obs[res_d$hog == "HOG4"], 0)
  expect_gt(sum(res_d$n_perm), 0L)
})


test_that("permutation_hog_test rejects non-dgCMatrix Matrix classes", {
  td <- make_test_nets()
  net1_bad <- modifyList(td$net1, list(
    network = Matrix::Matrix(td$net1$network, sparse = TRUE)
  ))
  expect_s4_class(net1_bad$network, "dsCMatrix")
  expect_error(
    permutation_hog_test(net1_bad, sparse_net(td$net2), td$comparison),
    "network must be a dgCMatrix"
  )
  expect_error(
    permutation_hog_test(sparse_net(td$net1), net1_bad, td$comparison),
    "network must be a dgCMatrix"
  )
})


test_that("torch backend accepts sparse networks and matches dense", {
  skip_if_not_installed("torch")
  skip_if_not(tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE),
              "torch backend (Lantern) not available")
  td <- make_test_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)

  set.seed(42)
  res_d <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 100L, min_exceedances = 10L, use_torch = TRUE
  )
  set.seed(42)
  res_s <- permutation_hog_test(
    net1_s, net2_s, td$comparison,
    max_permutations = 100L, min_exceedances = 10L, use_torch = TRUE
  )
  expect_equal(res_s, res_d)
})


test_that(".adj_edges gives identical edges for dense and sparse input", {
  td <- make_graded_nets()
  m <- td$net1$network
  n <- nrow(m)
  thr <- 8

  # every non-zero entry stored, including the diagonal (10) and the tiers
  # below thr (7, 4): both the rows != cols and the @x >= thr filters bite
  ij <- which(m != 0, arr.ind = TRUE)
  s <- Matrix::sparseMatrix(i = ij[, 1L], j = ij[, 2L], x = m[ij],
                            dims = dim(m), dimnames = dimnames(m))
  expect_s4_class(s, "dgCMatrix")
  expect_equal(unname(Matrix::diag(s)), rep(10, n))
  expect_gt(sum(s@x < thr), 0L)

  e_d <- .adj_edges(m, thr)
  e_s <- .adj_edges(s, thr)
  ref <- which((m >= thr) & (row(m) != col(m)), arr.ind = TRUE)
  expect_equal(e_d, list(rows = unname(ref[, 1L]), cols = unname(ref[, 2L])))
  expect_equal(e_s, e_d)
  expect_type(e_s$rows, "integer")
  expect_type(e_s$cols, "integer")
  expect_gt(length(e_d$rows), 0L)

  # empty store
  expect_equal(.adj_edges(dense_to_dgc(m, 100), thr),
               list(rows = integer(0), cols = integer(0)))
})


test_that("permutation_hog_test rejects mixed dense/sparse inputs", {
  td <- make_test_nets()
  expect_error(
    permutation_hog_test(sparse_net(td$net1), td$net2, td$comparison),
    "both dense or both sparse"
  )
  expect_error(
    permutation_hog_test(td$net1, sparse_net(td$net2), td$comparison),
    "both dense or both sparse"
  )
})


test_that("torch backend with tighter threshold: sparse equals dense", {
  skip_if_not_installed("torch")
  skip_if_not(tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE),
              "torch backend (Lantern) not available")
  td <- make_graded_nets()
  comparison <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  thr <- 8
  net1_s <- modifyList(sparse_net(td$net1), list(threshold = thr))
  net2_s <- modifyList(sparse_net(td$net2), list(threshold = thr))
  net1_d <- modifyList(td$net1, list(threshold = thr))
  net2_d <- modifyList(td$net2, list(threshold = thr))

  set.seed(11)
  res_d <- permutation_hog_test(
    net1_d, net2_d, comparison,
    max_permutations = 200L, min_exceedances = 10L
  )
  set.seed(11)
  res_td <- permutation_hog_test(
    net1_d, net2_d, comparison,
    max_permutations = 200L, min_exceedances = 10L, use_torch = TRUE
  )
  set.seed(11)
  res_ts <- permutation_hog_test(
    net1_s, net2_s, comparison,
    max_permutations = 200L, min_exceedances = 10L, use_torch = TRUE
  )

  # sparse store filtered on the GPU path (adj_to_gpu) equals the dense path
  expect_equal(res_ts, res_td)
  # and the torch statistic agrees with the C++ backend (float32 on MPS)
  expect_equal(res_ts$T_obs[order(res_ts$hog)], res_d$T_obs[order(res_d$hog)],
               tolerance = 1e-5)
  expect_gt(res_ts$T_obs[res_ts$hog == "HOG1"], 0)
})


# ---- self-excluded urn (D5) ----

test_that("T_obs uses the self-excluded urn and matches the R oracle", {
  # Fixture: make_self_excluded_nets() (helper-reference.R) -- the anchor
  # is itself ortholog-reachable from a HOG-mate's neighbours and must
  # leave the reachable set; the population is N - 1 in both directions.
  td <- make_self_excluded_nets()
  sp1 <- paste0("A", 1:3)
  sp2 <- paste0("B", 1:3)

  set.seed(3)
  res <- permutation_hog_test(td$net1, td$net2, td$comparison,
                              max_permutations = 50L, min_exceedances = 5L)
  t_hog1 <- res$T_obs[res$hog == "HOG1"]
  expected <- reference_T_obs(td$net1$network, td$net2$network, 0.5, 0.5,
                              td$ortho, sp1, sp2)
  expect_equal(t_hog1, expected, tolerance = 1e-12)

  # the exclusion is exercised: the pre-0.2.0 urn (k, N) gives a
  # different statistic
  old <- reference_T_obs(td$net1$network, td$net2$network, 0.5, 0.5,
                         td$ortho, sp1, sp2, self_exclude = FALSE)
  expect_false(isTRUE(all.equal(t_hog1, old)))

  # sparse path identical
  set.seed(3)
  res_s <- permutation_hog_test(sparse_net(td$net1), sparse_net(td$net2),
                                td$comparison,
                                max_permutations = 50L, min_exceedances = 5L)
  expect_equal(res_s, res)
})


# ---- forced flag-vector mode (internal hook) ----

# The flag-vector engine (the max(n1, n2) > 100000 path) shares loop order,
# skip conditions, and RNG consumption with the bit-vector engine, so a
# seeded forced run must be IDENTICAL to the bit-vector run — not merely
# close. The forced engine announces itself on stderr; asserting on that
# marker guards these tests against passing vacuously on the bit-vector
# path.

test_that("forced flag-vector engine equals bit-vector engine (seeded)", {
  skip_if_not_installed("withr")
  td <- make_test_nets()

  set.seed(42)
  res_bv <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 500L, min_exceedances = 20L
  )
  set.seed(7)
  res_bv_less <- permutation_hog_test(
    td$net1, td$net2, td$comparison, alternative = "less",
    max_permutations = 500L, min_exceedances = 20L
  )

  withr::local_options(rcomplex.force_flag_vector = TRUE)
  set.seed(42)
  msgs <- capture.output(
    res_fl <- permutation_hog_test(
      td$net1, td$net2, td$comparison,
      max_permutations = 500L, min_exceedances = 20L
    ),
    type = "message"
  )
  expect_true(any(grepl("flag-vector", msgs)))
  expect_equal(res_fl, res_bv)

  # divergence direction too
  set.seed(7)
  msgs <- capture.output(
    res_fl_less <- permutation_hog_test(
      td$net1, td$net2, td$comparison, alternative = "less",
      max_permutations = 500L, min_exceedances = 20L
    ),
    type = "message"
  )
  expect_true(any(grepl("flag-vector", msgs)))
  expect_equal(res_fl_less, res_bv_less)
})


test_that("forced flag-vector T_obs matches the self-excluded R oracle", {
  skip_if_not_installed("withr")
  # Same fixture as the bit-vector D5 test above
  # (make_self_excluded_nets(), helper-reference.R): the anchor is itself
  # ortholog-reachable and compute_T_flags must drop it from the urn
  # (k - 1, population N - 1).
  td <- make_self_excluded_nets()

  withr::local_options(rcomplex.force_flag_vector = TRUE)
  set.seed(3)
  msgs <- capture.output(
    res <- permutation_hog_test(td$net1, td$net2, td$comparison,
                                max_permutations = 50L,
                                min_exceedances = 5L),
    type = "message"
  )
  expect_true(any(grepl("flag-vector", msgs)))

  t_hog1 <- res$T_obs[res$hog == "HOG1"]
  expected <- reference_T_obs(td$net1$network, td$net2$network, 0.5, 0.5,
                              td$ortho, paste0("A", 1:3), paste0("B", 1:3))
  expect_equal(t_hog1, expected, tolerance = 1e-12)

  # the exclusion is exercised: the pre-0.2.0 urn (k, N) differs
  old <- reference_T_obs(td$net1$network, td$net2$network, 0.5, 0.5,
                         td$ortho, paste0("A", 1:3), paste0("B", 1:3),
                         self_exclude = FALSE)
  expect_false(isTRUE(all.equal(t_hog1, old)))
})


test_that("forced flag-vector mode: sparse equals dense (seeded)", {
  skip_if_not_installed("withr")
  td <- make_test_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)

  withr::local_options(rcomplex.force_flag_vector = TRUE)
  set.seed(42)
  msgs_d <- capture.output(
    res_d <- permutation_hog_test(
      td$net1, td$net2, td$comparison,
      max_permutations = 200L, min_exceedances = 10L
    ),
    type = "message"
  )
  set.seed(42)
  msgs_s <- capture.output(
    res_s <- permutation_hog_test(
      net1_s, net2_s, td$comparison,
      max_permutations = 200L, min_exceedances = 10L
    ),
    type = "message"
  )
  expect_true(any(grepl("flag-vector", msgs_d)))
  expect_true(any(grepl("flag-vector", msgs_s)))
  expect_equal(res_s, res_d)
})


test_that("forced flag-vector engine: n_cores = 2 matches serial T_obs", {
  skip_if_not_installed("withr")
  td <- make_test_nets()

  withr::local_options(rcomplex.force_flag_vector = TRUE)
  set.seed(42)
  msgs_1 <- capture.output(
    res_1 <- permutation_hog_test(
      td$net1, td$net2, td$comparison,
      max_permutations = 200L, min_exceedances = 10L, n_cores = 1L
    ),
    type = "message"
  )
  set.seed(42)
  msgs_2 <- capture.output(
    res_2 <- permutation_hog_test(
      td$net1, td$net2, td$comparison,
      max_permutations = 200L, min_exceedances = 10L, n_cores = 2L
    ),
    type = "message"
  )
  expect_true(any(grepl("flag-vector", msgs_1)))
  expect_true(any(grepl("flag-vector", msgs_2)))
  # T_obs is permutation-free, so the parallel run must reproduce the
  # serial statistic exactly (without OpenMP n_cores = 2 degrades to
  # serial and the assertion stays valid).
  expect_equal(res_2$T_obs[order(res_2$hog)], res_1$T_obs[order(res_1$hog)])
})


test_that("option set to FALSE keeps the bit-vector path unchanged", {
  skip_if_not_installed("withr")
  td <- make_test_nets()

  set.seed(42)
  res_plain <- permutation_hog_test(
    td$net1, td$net2, td$comparison,
    max_permutations = 200L, min_exceedances = 10L
  )

  withr::local_options(rcomplex.force_flag_vector = FALSE)
  set.seed(42)
  msgs <- capture.output(
    res_false <- permutation_hog_test(
      td$net1, td$net2, td$comparison,
      max_permutations = 200L, min_exceedances = 10L
    ),
    type = "message"
  )
  expect_false(any(grepl("flag-vector", msgs)))
  expect_equal(res_false, res_plain)
})
