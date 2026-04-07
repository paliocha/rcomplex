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
    hog = 1:30
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 30)
  expected_cols <- c("Species1", "Species2", "hog",
                     "Species1.neigh", "Species1.ortho.neigh",
                     "Species1.neigh.overlap", "Species1.p.val.con",
                     "Species1.p.val.div", "Species1.effect.size",
                     "Species1.jaccard",
                     "Species2.neigh", "Species2.ortho.neigh",
                     "Species2.neigh.overlap", "Species2.p.val.con",
                     "Species2.p.val.div", "Species2.effect.size",
                     "Species2.jaccard")
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
    hog = 1:30
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_true(all(
    result$Species1.p.val.con >= 0 & result$Species1.p.val.con <= 1
  ))
  expect_true(all(
    result$Species2.p.val.con >= 0 & result$Species2.p.val.con <= 1
  ))
  expect_true(all(
    result$Species1.p.val.div >= 0 & result$Species1.p.val.div <= 1
  ))
  expect_true(all(
    result$Species2.p.val.div >= 0 & result$Species2.p.val.div <= 1
  ))
  expect_true(all(result$Species1.effect.size >= 0))
  expect_true(all(result$Species2.effect.size >= 0))

  # Jaccard in [0, 1]
  expect_true(all(result$Species1.jaccard >= 0 & result$Species1.jaccard <= 1))
  expect_true(all(result$Species2.jaccard >= 0 & result$Species2.jaccard <= 1))
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
    hog = 1:30
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
    expect_equal(
      cpp_result$Species1.neigh.overlap[i],
      ref$Species1.neigh.overlap
    )
    expect_equal(
      cpp_result$Species1.p.val.con[i],
      ref$Species1.p.val.con, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.p.val.div[i],
      ref$Species1.p.val.div, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.effect.size[i],
      ref$Species1.effect.size, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.neigh[i],
      ref$Species2.neigh
    )
    expect_equal(
      cpp_result$Species2.ortho.neigh[i],
      ref$Species2.ortho.neigh
    )
    expect_equal(
      cpp_result$Species2.neigh.overlap[i],
      ref$Species2.neigh.overlap
    )
    expect_equal(
      cpp_result$Species2.p.val.con[i],
      ref$Species2.p.val.con, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.p.val.div[i],
      ref$Species2.p.val.div, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.effect.size[i],
      ref$Species2.effect.size, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.jaccard[i],
      ref$Species1.jaccard, tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.jaccard[i],
      ref$Species2.jaccard, tolerance = 1e-12
    )
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
    hog = c(1:30, 1)
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
    hog = 1:3
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
  for (j in 2:21) {
    net1_mat[1, j] <- 10
    net1_mat[j, 1] <- 10
  }
  # Gene 51 neighbors: 52-71 (completely disjoint from gene 1's neighbors)
  for (j in 52:71) {
    net1_mat[51, j] <- 10
    net1_mat[j, 51] <- 10
  }

  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  # Gene 1 neighbors in net2: 2-21 (same as net1 = conserved with A_001)
  for (j in 2:21) {
    net2_mat[1, j] <- 10
    net2_mat[j, 1] <- 10
  }
  # Gene 51 neighbors in net2: 2-21 (same as gene 1, disjoint from A_051's)
  for (j in 2:21) {
    net2_mat[51, j] <- 10
    net2_mat[j, 51] <- 10
  }

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    hog = 1:n
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
  for (j in 2:6) {
    net_mat[1, j] <- 10
    net_mat[j, 1] <- 10
  }

  net1 <- list(network = net_mat, threshold = 5)
  # Use same structure for net2 with different names
  net2_mat <- net_mat
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("H_", sprintf("%03d", 1:n))
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("G_", sprintf("%03d", 1:n)),
    Species2 = paste0("H_", sprintf("%03d", 1:n)),
    hog = 1:n
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
  for (j in 2:11) {
    net1_mat[1, j] <- 10
    net1_mat[j, 1] <- 10
  }

  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  # Gene 1 connected to 21-30 (10 neighbors, completely disjoint from net1's)
  for (j in 21:30) {
    net2_mat[1, j] <- 10
    net2_mat[j, 1] <- 10
  }

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    hog = 1:n
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  row1 <- result[result$Species1 == "A_001" & result$Species2 == "B_001", ]

  # 0 overlap, so effect size = 0
  expect_equal(row1$Species1.neigh.overlap, 0)
  expect_equal(row1$Species1.effect.size, 0)
  expect_true(row1$Species1.p.val.div < 0.05)
})


# --- Tests for Jaccard overlap ---

test_that("Jaccard = 1 for identical neighborhoods", {
  n <- 20
  net_mat <- matrix(0, n, n)
  rownames(net_mat) <- colnames(net_mat) <- paste0("G_", sprintf("%03d", 1:n))
  for (j in 2:6) {
    net_mat[1, j] <- 10
    net_mat[j, 1] <- 10
  }

  net1 <- list(network = net_mat, threshold = 5)
  net2_mat <- net_mat
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("H_", sprintf("%03d", 1:n))
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("G_", sprintf("%03d", 1:n)),
    Species2 = paste0("H_", sprintf("%03d", 1:n)),
    hog = 1:n
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  row1 <- result[result$Species1 == "G_001" & result$Species2 == "H_001", ]

  # Identical neighborhoods -> Jaccard = 1
  expect_equal(row1$Species1.jaccard, 1.0)
  expect_equal(row1$Species2.jaccard, 1.0)
})

test_that("Jaccard = 0 for disjoint neighborhoods", {
  n <- 30
  net1_mat <- matrix(0, n, n)
  rownames(net1_mat) <- colnames(net1_mat) <- paste0("A_", sprintf("%03d", 1:n))
  for (j in 2:11) {
    net1_mat[1, j] <- 10
    net1_mat[j, 1] <- 10
  }

  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  for (j in 21:30) {
    net2_mat[1, j] <- 10
    net2_mat[j, 1] <- 10
  }

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    hog = 1:n
  )

  result <- compare_neighborhoods(net1, net2, ortho)
  row1 <- result[result$Species1 == "A_001" & result$Species2 == "B_001", ]

  # Disjoint neighborhoods -> Jaccard = 0
  expect_equal(row1$Species1.jaccard, 0.0)
})

test_that("Jaccard = 0 (not NaN) when both neighborhoods are empty", {
  n <- 10
  # All-zero networks: no edges above threshold
  net1_mat <- matrix(0, n, n)
  rownames(net1_mat) <- colnames(net1_mat) <- paste0("A_", sprintf("%03d", 1:n))
  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    hog = 1:n
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  # Both directions: empty neighborhoods -> neigh=0, ortho_neigh=0, union=0 -> jaccard=0

  expect_true(all(result$Species1.jaccard == 0.0))
  expect_true(all(result$Species2.jaccard == 0.0))
  expect_true(all(!is.nan(result$Species1.jaccard)))
  expect_true(all(!is.nan(result$Species2.jaccard)))
})


# --- Tests for comparison_to_edges() ---

test_that("comparison_to_edges produces correct edge format", {
  comp <- data.frame(
    Species1 = c("A1", "A2"),
    Species2 = c("B1", "B2"),
    hog = c(1L, 2L),
    Species1.effect.size = c(4.0, 1.0),
    Species2.effect.size = c(9.0, 1.0),
    Species1.jaccard = c(0.8, 0.0),
    Species2.jaccard = c(0.5, 0.0),
    Species1.q.val.con = c(0.01, 0.80),
    Species2.q.val.con = c(0.03, 0.90)
  )

  edges <- comparison_to_edges(comp, "SP_A", "SP_B")

  expect_equal(names(edges), c("gene1", "gene2", "species1", "species2",
                                "hog", "q.value", "effect_size", "jaccard",
                                "type"))
  expect_equal(edges$gene1, c("A1", "A2"))
  expect_equal(edges$species1, c("SP_A", "SP_A"))
  expect_equal(edges$species2, c("SP_B", "SP_B"))
  # q.value = min of two directions
  expect_equal(edges$q.value, c(0.01, 0.80))
  # effect_size = geometric mean
  expect_equal(edges$effect_size, c(sqrt(4 * 9), sqrt(1 * 1)))
  # jaccard = geometric mean of directional Jaccards
  expect_equal(edges$jaccard, c(sqrt(0.8 * 0.5), sqrt(0.0 * 0.0)))
  # type classification
  expect_equal(edges$type, c("conserved", "ns"))
})


test_that("comparison_to_edges handles alternative='less'", {
  comp <- data.frame(
    Species1 = "A1", Species2 = "B1", hog = 1L,
    Species1.effect.size = 0.2, Species2.effect.size = 0.3,
    Species1.jaccard = 0.1, Species2.jaccard = 0.2,
    Species1.q.val.div = 0.01, Species2.q.val.div = 0.02
  )

  edges <- comparison_to_edges(comp, "SP_A", "SP_B", alternative = "less")

  expect_equal(edges$q.value, 0.01)
  expect_equal(edges$type, "diverged")
})


test_that("comparison_to_edges validates missing columns", {
  expect_error(
    comparison_to_edges(data.frame(x = 1), "A", "B"),
    "comparison missing required columns"
  )
})


# --- Tests for run_pairwise_comparisons() ---

test_that("run_pairwise_comparisons returns combined edges for 2 species", {
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
    hog = paste0("HOG", 1:30)
  )

  result <- run_pairwise_comparisons(
    networks = list(SP_A = net1, SP_B = net2),
    orthologs = ortho
  )

  expect_true(is.data.frame(result))
  expect_true(all(c("gene1", "gene2", "species1", "species2", "hog",
                     "q.value", "effect_size", "type") %in% names(result)))
  expect_true(all(result$species1 == "SP_A"))
  expect_true(all(result$species2 == "SP_B"))
})


test_that("run_pairwise_comparisons handles 3 species (all pairs)", {
  set.seed(42)
  make_net <- function(prefix, n = 30) {
    expr <- matrix(rnorm(n * 10), nrow = n)
    rownames(expr) <- paste0(prefix, "_", sprintf("%03d", seq_len(n)))
    compute_network(expr, density = 0.1, mr_log_transform = FALSE)
  }

  nets <- list(SP_A = make_net("A"), SP_B = make_net("B"), SP_C = make_net("C"))

  ortho <- data.frame(
    Species1 = c(paste0("A_", sprintf("%03d", 1:20)),
                 paste0("A_", sprintf("%03d", 1:20)),
                 paste0("B_", sprintf("%03d", 1:20))),
    Species2 = c(paste0("B_", sprintf("%03d", 1:20)),
                 paste0("C_", sprintf("%03d", 1:20)),
                 paste0("C_", sprintf("%03d", 1:20))),
    hog = rep(paste0("HOG", 1:20), 3)
  )

  result <- run_pairwise_comparisons(nets, ortho)

  # Should have edges from multiple pairs
  sp_pairs <- unique(paste(result$species1, result$species2))
  expect_true(length(sp_pairs) >= 2)
  expect_true("SP_A SP_B" %in% sp_pairs)
})


test_that("run_pairwise_comparisons validates inputs", {
  expect_error(
    run_pairwise_comparisons(list(A = list(network = matrix(0))),
                              data.frame(Species1 = "a",
                                         Species2 = "b", hog = 1)),
    "at least 2 species"
  )
  expect_error(
    run_pairwise_comparisons(c(a = 1, b = 2), data.frame(x = 1)),
    "networks must be a named list"
  )
  expect_error(
    run_pairwise_comparisons(list(A = list(), B = list()), data.frame(x = 1)),
    "orthologs must have columns"
  )
})


test_that("run_pairwise_comparisons with custom species_pairs", {
  set.seed(42)
  make_net <- function(prefix, n = 30) {
    expr <- matrix(rnorm(n * 10), nrow = n)
    rownames(expr) <- paste0(prefix, "_", sprintf("%03d", seq_len(n)))
    compute_network(expr, density = 0.1, mr_log_transform = FALSE)
  }

  nets <- list(SP_A = make_net("A"), SP_B = make_net("B"), SP_C = make_net("C"))

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:20)),
    Species2 = paste0("B_", sprintf("%03d", 1:20)),
    hog = paste0("HOG", 1:20)
  )

  # Only compare A vs B, skip A-C and B-C
  result <- run_pairwise_comparisons(nets, ortho,
                                      species_pairs = list(c("SP_A", "SP_B")))

  if (nrow(result) > 0) {
    expect_true(all(result$species1 == "SP_A" & result$species2 == "SP_B"))
  }
})


# --- Shared fixtures for find_coexpressologs / density_sweep ---

make_coexpr_fixtures <- function(n1 = 50, n2 = 40, n_ortho = 30,
                                  density = 0.1) {
  set.seed(42)
  expr1 <- matrix(rnorm(n1 * 10), nrow = n1, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", seq_len(n1)))
  expr2 <- matrix(rnorm(n2 * 10), nrow = n2, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", seq_len(n2)))

  net1 <- compute_network(expr1, density = density, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = density, mr_log_transform = FALSE)

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", seq_len(n_ortho))),
    Species2 = paste0("B_", sprintf("%03d", seq_len(n_ortho))),
    hog = paste0("HOG", seq_len(n_ortho))
  )

  list(nets = list(SP_A = net1, SP_B = net2), ortho = ortho)
}


test_that("find_coexpressologs default method is analytical", {
  skip_on_cran()
  fix <- make_coexpr_fixtures()
  result_default <- find_coexpressologs(fix$nets, fix$ortho)
  result_explicit <- find_coexpressologs(fix$nets, fix$ortho,
                                          method = "analytical")
  expect_identical(result_default, result_explicit)
})


test_that("find_coexpressologs with method='permutation' returns correct structure", {
  skip_on_cran()
  fix <- make_coexpr_fixtures()
  result <- find_coexpressologs(fix$nets, fix$ortho, method = "permutation")

  expect_true(is.data.frame(result))
  expect_true(all(c("gene1", "gene2", "species1", "species2",
                     "hog", "q.value", "effect_size", "type")
                  %in% names(result)))
  expect_true(all(result$type %in% c("conserved", "ns")))
})


test_that("find_coexpressologs alternative='less' produces 'diverged' labels", {
  skip_on_cran()
  fix <- make_coexpr_fixtures()

  result_perm <- find_coexpressologs(fix$nets, fix$ortho,
                                      method = "permutation",
                                      alternative = "less")
  expect_true(all(result_perm$type %in% c("diverged", "ns")))

  result_anal <- find_coexpressologs(fix$nets, fix$ortho,
                                      method = "analytical",
                                      alternative = "less")
  expect_true(all(result_anal$type %in% c("diverged", "ns")))
})


# --- Tests for density_sweep() ---

test_that("density_sweep returns correct structure", {
  fix <- make_coexpr_fixtures()
  mults <- c(0.98, 1.0, 1.02)

  result <- suppressMessages(density_sweep(
    networks = fix$nets, orthologs = fix$ortho, multipliers = mults
  ))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(c("multiplier", "eff_density", "n_significant", "edges",
                     "species_densities") %in% names(result)))
  expect_equal(result$multiplier, mults)
  expect_true(is.numeric(result$eff_density))
  expect_true(is.integer(result$n_significant))
  expect_true(is.list(result$edges))
  for (j in seq_len(nrow(result))) {
    expect_s3_class(result$edges[[j]], "data.frame")
  }
  # Per-species densities
  expect_true(is.list(result$species_densities))
  expect_true(is.numeric(result$species_densities[[1]]))
  expect_equal(length(result$species_densities[[1]]), length(fix$nets))
})


test_that("density_sweep at multiplier=1 matches find_coexpressologs", {
  fix <- make_coexpr_fixtures()

  result <- suppressMessages(density_sweep(
    networks = fix$nets, orthologs = fix$ortho,
    multipliers = 1.0, method = "analytical"
  ))

  expect_equal(nrow(result), 1)
  expect_equal(result$multiplier, 1.0)

  direct <- find_coexpressologs(
    networks = fix$nets, orthologs = fix$ortho, method = "analytical"
  )

  sweep_edges <- result$edges[[1]]
  sweep_edges <- sweep_edges[order(sweep_edges$gene1, sweep_edges$gene2), ]
  direct <- direct[order(direct$gene1, direct$gene2), ]
  rownames(sweep_edges) <- rownames(direct) <- NULL
  expect_identical(sweep_edges, direct)
})


test_that("density_sweep validates inputs", {
  dummy_ortho <- data.frame(Species1 = "A", Species2 = "B", hog = "H1")

  expect_error(
    density_sweep(c(a = 1, b = 2), dummy_ortho),
    "networks must be a named list"
  )
  expect_error(
    density_sweep(list(A = list(network = matrix(0), threshold = 1)),
                   dummy_ortho),
    "networks must contain at least 2 species"
  )
  expect_error(
    density_sweep(list(A = list(network = matrix(0)),
                        B = list(network = matrix(0), threshold = 1)),
                   dummy_ortho),
    "each network must have 'network' and 'threshold' elements"
  )

  net_a <- list(network = matrix(0, 2, 2), threshold = 1)
  net_b <- list(network = matrix(0, 2, 2), threshold = 1)
  expect_error(
    density_sweep(list(A = net_a, B = net_b), dummy_ortho,
                   multipliers = character(0)),
    "multipliers must be a non-empty numeric vector"
  )
  expect_error(
    density_sweep(list(A = net_a, B = net_b), dummy_ortho,
                   multipliers = -1),
    "all multipliers must be positive"
  )
  expect_error(
    density_sweep(list(A = net_a, B = net_b), data.frame(x = 1)),
    "orthologs must have columns"
  )
})
