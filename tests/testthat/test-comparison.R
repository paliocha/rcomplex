test_that("compare_neighborhoods returns correct structure", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:30)),
    Species2 = paste0("B_", sprintf("%03d", 1:30)),
    OrthoGroup = 1:30,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 30)
  expected_cols <- c("Species1", "Species2", "OrthoGroup",
                     "Species1.neigh", "Species1.ortho.neigh",
                     "Species1.neigh.overlap", "Species1.p.val.con",
                     "Species1.p.val.div", "Species1.effect.size",
                     "Species2.neigh", "Species2.ortho.neigh",
                     "Species2.neigh.overlap", "Species2.p.val.con",
                     "Species2.p.val.div", "Species2.effect.size")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("p-values are in [0,1] and effect sizes are positive", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:30)),
    Species2 = paste0("B_", sprintf("%03d", 1:30)),
    OrthoGroup = 1:30,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_true(all(result$Species1.p.val.con >= 0 & result$Species1.p.val.con <= 1))
  expect_true(all(result$Species2.p.val.con >= 0 & result$Species2.p.val.con <= 1))
  expect_true(all(result$Species1.p.val.div >= 0 & result$Species1.p.val.div <= 1))
  expect_true(all(result$Species2.p.val.div >= 0 & result$Species2.p.val.div <= 1))
  expect_true(all(result$Species1.effect.size >= 0))
  expect_true(all(result$Species2.effect.size >= 0))
})

test_that("C++ comparison matches R reference", {
  set.seed(123)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:30)),
    Species2 = paste0("B_", sprintf("%03d", 1:30)),
    OrthoGroup = 1:30,
    stringsAsFactors = FALSE
  )

  cpp_result <- compare_neighborhoods(net1, net2, ortho)

  # Reference R implementation for each pair
  for (i in 1:5) {
    ref <- reference_compare_pair(
      net1$network, net2$network,
      net1$threshold, net2$threshold,
      ortho,
      cpp_result$Species1[i], cpp_result$Species2[i]
    )
    expect_equal(cpp_result$Species1.neigh[i], ref$Species1.neigh)
    expect_equal(cpp_result$Species1.ortho.neigh[i], ref$Species1.ortho.neigh)
    expect_equal(cpp_result$Species1.neigh.overlap[i], ref$Species1.neigh.overlap)
    expect_equal(cpp_result$Species1.p.val.con[i], ref$Species1.p.val.con, tolerance = 1e-12)
    expect_equal(cpp_result$Species1.p.val.div[i], ref$Species1.p.val.div, tolerance = 1e-12)
    expect_equal(cpp_result$Species1.effect.size[i], ref$Species1.effect.size, tolerance = 1e-12)
    expect_equal(cpp_result$Species2.neigh[i], ref$Species2.neigh)
    expect_equal(cpp_result$Species2.ortho.neigh[i], ref$Species2.ortho.neigh)
    expect_equal(cpp_result$Species2.neigh.overlap[i], ref$Species2.neigh.overlap)
    expect_equal(cpp_result$Species2.p.val.con[i], ref$Species2.p.val.con, tolerance = 1e-12)
    expect_equal(cpp_result$Species2.p.val.div[i], ref$Species2.p.val.div, tolerance = 1e-12)
    expect_equal(cpp_result$Species2.effect.size[i], ref$Species2.effect.size, tolerance = 1e-12)
  }
})

test_that("multicopy orthologs handled correctly", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  # 1:N mapping: A_001 maps to B_001 and B_031
  ortho <- data.frame(
    Species1 = c(paste0("A_", sprintf("%03d", 1:30)), "A_001"),
    Species2 = c(paste0("B_", sprintf("%03d", 1:30)), "B_031"),
    OrthoGroup = c(1:30, 1),
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_equal(nrow(result), 31)
  # No errors
  expect_true(all(!is.na(result$Species1.p.val.con)))
  expect_true(all(!is.na(result$Species1.p.val.div)))
})

test_that("orthologs not in network are filtered", {
  set.seed(42)
  expr1 <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:10))
  expr2 <- matrix(rnorm(80), nrow = 8, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:8))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  # Include orthologs that aren't in the networks
  ortho <- data.frame(
    Species1 = c("A_001", "A_002", "A_999"),
    Species2 = c("B_001", "B_002", "B_999"),
    OrthoGroup = 1:3,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  # Only the first two should be in the result
  expect_equal(nrow(result), 2)
})

test_that("input validation works for compare_neighborhoods", {
  expect_error(compare_neighborhoods(list(), list(), data.frame()),
               "network object")
})

test_that("divergence p-values detect disjoint neighborhoods", {
  # Build controlled networks large enough for statistical power
  # Need n large relative to neighborhood size for phyper to detect depletion
  n <- 100

  net1_mat <- matrix(0, n, n)
  rownames(net1_mat) <- colnames(net1_mat) <- paste0("A_", sprintf("%03d", 1:n))
  # Gene 1 neighbors: 2-21 (20 neighbors)
  for (j in 2:21) { net1_mat[1, j] <- 10; net1_mat[j, 1] <- 10 }
  # Gene 51 neighbors: 52-71 (completely disjoint from gene 1's neighbors)
  for (j in 52:71) { net1_mat[51, j] <- 10; net1_mat[j, 51] <- 10 }

  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  # Gene 1 neighbors in net2: 2-21 (same as net1 = conserved with A_001)
  for (j in 2:21) { net2_mat[1, j] <- 10; net2_mat[j, 1] <- 10 }
  # Gene 51 neighbors in net2: 2-21 (same as gene 1, disjoint from A_051's)
  for (j in 2:21) { net2_mat[51, j] <- 10; net2_mat[j, 51] <- 10 }

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    OrthoGroup = 1:n,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  # Conserved pair (A_001, B_001): identical neighborhoods -> low con p-val
  row1 <- result[result$Species1 == "A_001" & result$Species2 == "B_001", ]
  expect_true(row1$Species1.p.val.con < 0.05)
  expect_true(row1$Species1.effect.size > 1)

  # Diverged pair (A_051, B_051): disjoint neighborhoods -> low div p-val
  row51 <- result[result$Species1 == "A_051" & result$Species2 == "B_051", ]
  expect_true(row51$Species1.p.val.div < 0.05)
  expect_true(row51$Species1.effect.size < 1)
})

test_that("identical neighborhoods give high divergence p-value", {
  n <- 20
  # Both networks have gene 1 connected to genes 2-6
  net_mat <- matrix(0, n, n)
  rownames(net_mat) <- colnames(net_mat) <- paste0("G_", sprintf("%03d", 1:n))
  for (j in 2:6) { net_mat[1, j] <- 10; net_mat[j, 1] <- 10 }

  net1 <- list(network = net_mat, threshold = 5)
  # Use same structure for net2 with different names
  net2_mat <- net_mat
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("H_", sprintf("%03d", 1:n))
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("G_", sprintf("%03d", 1:n)),
    Species2 = paste0("H_", sprintf("%03d", 1:n)),
    OrthoGroup = 1:n,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  row1 <- result[result$Species1 == "G_001" & result$Species2 == "H_001", ]

  # Conserved: low p-value
  expect_true(row1$Species1.p.val.con < 0.05)
  # Divergence: high p-value (neighborhoods are identical, not diverged)
  expect_true(row1$Species1.p.val.div > 0.5)
})

test_that("effect size < 1 when overlap is less than expected", {
  n <- 30
  net1_mat <- matrix(0, n, n)
  rownames(net1_mat) <- colnames(net1_mat) <- paste0("A_", sprintf("%03d", 1:n))
  # Gene 1 connected to 2-11 (10 neighbors)
  for (j in 2:11) { net1_mat[1, j] <- 10; net1_mat[j, 1] <- 10 }

  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  # Gene 1 connected to 21-30 (10 neighbors, completely disjoint from net1's)
  for (j in 21:30) { net2_mat[1, j] <- 10; net2_mat[j, 1] <- 10 }

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    OrthoGroup = 1:n,
    stringsAsFactors = FALSE
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  row1 <- result[result$Species1 == "A_001" & result$Species2 == "B_001", ]

  # 0 overlap, so effect size = 0
  expect_equal(row1$Species1.neigh.overlap, 0)
  expect_equal(row1$Species1.effect.size, 0)
  expect_true(row1$Species1.p.val.div < 0.05)
})
