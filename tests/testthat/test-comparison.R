test_that("compare_neighborhoods returns correct structure", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:30)),
    Species2 = paste0("B_", sprintf("%03d", 1:30)),
    hog = 1:30
  )

  result <- compare_neighborhoods(net1, net2, ortho)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 30)
  expected_cols <- c(
    "Species1", "Species2", "hog",
    "Species1.neigh", "Species1.ortho.neigh",
    "Species1.neigh.overlap", "Species1.p.val.con",
    "Species1.p.val.div", "Species1.effect.size",
    "Species1.jaccard",
    "Species2.neigh", "Species2.ortho.neigh",
    "Species2.neigh.overlap", "Species2.p.val.con",
    "Species2.p.val.div", "Species2.effect.size",
    "Species2.jaccard"
  )
  expect_true(all(expected_cols %in% names(result)))
})

test_that("p-values are in [0,1] and effect sizes are positive", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

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

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

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
      ref$Species1.p.val.con,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.p.val.div[i],
      ref$Species1.p.val.div,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.effect.size[i],
      ref$Species1.effect.size,
      tolerance = 1e-12
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
      ref$Species2.p.val.con,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.p.val.div[i],
      ref$Species2.p.val.div,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.effect.size[i],
      ref$Species2.effect.size,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.jaccard[i],
      ref$Species1.jaccard,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.jaccard[i],
      ref$Species2.jaccard,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.p.val.gt[i],
      ref$Species1.p.val.gt,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species1.p.val.eq[i],
      ref$Species1.p.val.eq,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.p.val.gt[i],
      ref$Species2.p.val.gt,
      tolerance = 1e-12
    )
    expect_equal(
      cpp_result$Species2.p.val.eq[i],
      ref$Species2.p.val.eq,
      tolerance = 1e-12
    )
  }
})


test_that(
  "self-excluded urn: anchor leaves the mapped set and the population",
  {
    # A1 has two orthologs (B1, B2). For the pair (A1, B1), B1's neighbour B2
    # maps back to the anchor A1, which must be dropped from the mapped set;
    # the population is the other N - 1 genes in both directions.
    n <- 10
    mk <- function(prefix) {
      m <- matrix(0, n, n)
      rownames(m) <- colnames(m) <- paste0(prefix, 1:n)
      m
    }
    link <- function(m, i, j) {
      m[i, j] <- m[j, i] <- 1
      m
    }
    a <- mk("A")
    for (j in 2:4) a <- link(a, 1, j) # N1(A1) = {A2, A3, A4}
    b <- mk("B")
    for (j in c(2:4, 7)) b <- link(b, 1, j) # N2(B1) = {B2, B3, B4, B7}
    net1 <- list(network = a, threshold = 0.5)
    net2 <- list(network = b, threshold = 0.5)
    ortho <- data.frame(
      Species1 = c("A1", "A1", "A2", "A3", "A4", "A5"),
      Species2 = c("B1", "B2", "B3", "B4", "B5", "B6"),
      hog = c("H1", "H1", "H2", "H3", "H4", "H5"),
      stringsAsFactors = FALSE
    )

    res <- compare_neighborhoods(net1, net2, ortho)
    r <- res[res$Species1 == "A1" & res$Species2 == "B1", ]

    # direction 1: m = 3; mapped {A1 (anchor, dropped), A2, A3} -> k = 2;
    # x = 2; population 9
    expect_equal(r$Species1.neigh, 3L)
    expect_equal(r$Species1.ortho.neigh, 2L)
    expect_equal(r$Species1.neigh.overlap, 2L)
    expect_equal(r$Species1.p.val.con, phyper(1, 3, 6, 2, lower.tail = FALSE))
    expect_equal(r$Species1.p.val.gt, phyper(2, 3, 6, 2, lower.tail = FALSE))
    expect_equal(r$Species1.p.val.eq, dhyper(2, 3, 6, 2))
    expect_equal(r$Species1.p.val.div, phyper(2, 3, 6, 2))
    expect_equal(r$Species1.effect.size, (2 / 2) / (3 / 9))
    expect_equal(r$Species1.jaccard, 2 / 3)

    # direction 2: anchor B1 is not in the mapped set {B3, B4, B5}; m = 4,
    # k = 3, x = 2; population 9 (the old urn would use 10)
    expect_equal(r$Species2.neigh, 4L)
    expect_equal(r$Species2.ortho.neigh, 3L)
    expect_equal(r$Species2.neigh.overlap, 2L)
    expect_equal(r$Species2.p.val.con, phyper(1, 4, 5, 3, lower.tail = FALSE))
    expect_equal(r$Species2.p.val.gt, phyper(2, 4, 5, 3, lower.tail = FALSE))
    expect_equal(r$Species2.p.val.eq, dhyper(2, 4, 5, 3))
    expect_equal(r$Species2.p.val.div, phyper(2, 4, 5, 3))
    expect_equal(r$Species2.effect.size, (2 / 3) / (4 / 9))

    # (A1, B2): B2's only neighbour maps back to the anchor -> empty mapped
    # set; p.val.gt / p.val.eq still form a proper randomized p-value
    r2 <- res[res$Species1 == "A1" & res$Species2 == "B2", ]
    expect_equal(r2$Species1.ortho.neigh, 0L)
    expect_equal(r2$Species1.neigh.overlap, 0L)
    expect_equal(r2$Species1.p.val.con, 1)
    expect_equal(r2$Species1.p.val.div, 1)
    expect_equal(r2$Species1.effect.size, 1)
    expect_equal(r2$Species1.p.val.gt, 0)
    expect_equal(r2$Species1.p.val.eq, 1)

    # pure-R oracle and sparse path agree
    ref <- reference_compare_pair(a, b, 0.5, 0.5, ortho, "A1", "B1")
    got <- r[names(ref)]
    rownames(got) <- NULL
    expect_equal(got, ref)
    res_s <- compare_neighborhoods(sparse_net(net1), sparse_net(net2), ortho)
    expect_equal(res_s, res)
  }
)


test_that("p.val.gt / p.val.eq decompose the exact hypergeometric tail", {
  td <- make_graded_nets()
  res <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  for (sp in c("Species1", "Species2")) {
    gt <- res[[paste0(sp, ".p.val.gt")]]
    eq <- res[[paste0(sp, ".p.val.eq")]]
    con <- res[[paste0(sp, ".p.val.con")]]
    div <- res[[paste0(sp, ".p.val.div")]]
    x <- res[[paste0(sp, ".neigh.overlap")]]
    expect_gt(sum(x > 1), 0L)
    expect_true(all(gt >= 0 & eq >= 0 & gt + eq <= 1 + 1e-12))
    # reported conservation p-value keeps the x > 1 gate
    expect_equal(con[x > 1], (gt + eq)[x > 1])
    expect_equal(con[x <= 1], rep(1, sum(x <= 1)))
    # lower tail contains the point mass
    expect_true(all(div - eq >= -1e-12))
  }
})

test_that("multicopy orthologs handled correctly", {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

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

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

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
  expect_error(
    compare_neighborhoods(list(), list(), data.frame()),
    "network object"
  )
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

  # Both directions: empty neighborhoods -> neigh=0, ortho_neigh=0,
  # union=0 -> jaccard=0  # nolint

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

  expect_equal(names(edges), c(
    "gene1", "gene2", "species1", "species2",
    "hog", "q.value", "effect_size", "jaccard",
    "power", "type"
  ))
  # The fixture carries no neighbourhood sizes, so power is undefined.
  expect_true(all(is.na(edges$power)))
  expect_equal(edges$gene1, c("A1", "A2"))
  expect_equal(edges$species1, c("SP_A", "SP_A"))
  expect_equal(edges$species2, c("SP_B", "SP_B"))
  # q.value = max of two directions (default reciprocal combine, D2)
  expect_equal(edges$q.value, c(0.03, 0.90))
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

  expect_equal(edges$q.value, 0.02)
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

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE, sparse = FALSE
  )

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
  expect_true(all(c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value", "effect_size", "type"
  ) %in% names(result)))
  expect_true(all(result$species1 == "SP_A"))
  expect_true(all(result$species2 == "SP_B"))
})


test_that("run_pairwise_comparisons handles 3 species (all pairs)", {
  set.seed(42)
  make_net <- function(prefix, n = 30) {
    expr <- matrix(rnorm(n * 10), nrow = n)
    rownames(expr) <- paste0(prefix, "_", sprintf("%03d", seq_len(n)))
    compute_network(expr,
      density = 0.1, mr_log_transform = FALSE, sparse = FALSE
    )
  }

  nets <- list(SP_A = make_net("A"), SP_B = make_net("B"), SP_C = make_net("C"))

  ortho <- data.frame(
    Species1 = c(
      paste0("A_", sprintf("%03d", 1:20)),
      paste0("A_", sprintf("%03d", 1:20)),
      paste0("B_", sprintf("%03d", 1:20))
    ),
    Species2 = c(
      paste0("B_", sprintf("%03d", 1:20)),
      paste0("C_", sprintf("%03d", 1:20)),
      paste0("C_", sprintf("%03d", 1:20))
    ),
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
    run_pairwise_comparisons(
      list(A = list(network = matrix(0))),
      data.frame(
        Species1 = "a",
        Species2 = "b", hog = 1
      )
    ),
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
    compute_network(expr,
      density = 0.1, mr_log_transform = FALSE, sparse = FALSE
    )
  }

  nets <- list(SP_A = make_net("A"), SP_B = make_net("B"), SP_C = make_net("C"))

  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:20)),
    Species2 = paste0("B_", sprintf("%03d", 1:20)),
    hog = paste0("HOG", 1:20)
  )

  # Only compare A vs B, skip A-C and B-C
  result <- run_pairwise_comparisons(nets, ortho,
    species_pairs = list(c("SP_A", "SP_B"))
  )

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

  net1 <- compute_network(expr1,
    density = density, mr_log_transform = FALSE, sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = density, mr_log_transform = FALSE, sparse = FALSE
  )

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
  # randomized pi0 (the default) draws from the global RNG
  set.seed(1)
  result_default <- find_coexpressologs(fix$nets, fix$ortho)
  set.seed(1)
  result_explicit <- find_coexpressologs(fix$nets, fix$ortho,
    method = "analytical"
  )
  expect_identical(result_default, result_explicit)
})


test_that(
  "find_coexpressologs with method='permutation' returns correct structure",
  {
    skip_on_cran()
    fix <- make_coexpr_fixtures()
    result <- find_coexpressologs(fix$nets, fix$ortho, method = "permutation")

    expect_true(is.data.frame(result))
    expect_true(all(c(
      "gene1", "gene2", "species1", "species2",
      "hog", "q.value", "effect_size", "type"
    )
    %in% names(result)))
    expect_true(all(result$type %in% c("conserved", "ns")))
  }
)


test_that("find_coexpressologs alternative='less' produces 'diverged' labels", {
  skip_on_cran()
  fix <- make_coexpr_fixtures()

  result_perm <- find_coexpressologs(fix$nets, fix$ortho,
    method = "permutation",
    alternative = "less"
  )
  expect_true(all(result_perm$type %in% c("diverged", "ns")))

  result_anal <- find_coexpressologs(fix$nets, fix$ortho,
    method = "analytical",
    alternative = "less"
  )
  expect_true(all(result_anal$type %in% c("diverged", "ns")))
})


test_that(
  "find_coexpressologs(out_file = ) streams edges identical to in-memory",
  {
    skip_on_cran()
    fix <- make_coexpr_fixtures()
    out <- withr::local_tempfile(fileext = ".csv")

    set.seed(1)
    in_memory <- find_coexpressologs(fix$nets, fix$ortho, method = "analytical")

    set.seed(1)
    ret <- find_coexpressologs(fix$nets, fix$ortho,
      method = "analytical", out_file = out
    )

    expect_identical(ret, out)
    expect_true(file.exists(out))
    from_file <- as.data.frame(data.table::fread(out))
    # No pair is called on these random networks, so `power` is all NA
    # and fread() reads the empty column back as logical.
    from_file$power <- as.numeric(from_file$power)
    # fread() infers types from text; compare on values, not attributes
    expect_equal(from_file, in_memory, ignore_attr = TRUE)
  }
)


test_that("find_coexpressologs(out_file = ) overwrites a stale file", {
  skip_on_cran()
  fix <- make_coexpr_fixtures()
  out <- withr::local_tempfile(fileext = ".csv")
  writeLines("stale,content", out)

  find_coexpressologs(fix$nets, fix$ortho,
    method = "analytical", out_file = out
  )
  from_file <- data.table::fread(out)
  expect_false("stale" %in% names(from_file))
  expect_true(all(c("gene1", "gene2", "q.value") %in% names(from_file)))
})


test_that("find_coexpressologs(out_file = ) validates its argument", {
  fix <- make_coexpr_fixtures()
  expect_error(
    find_coexpressologs(fix$nets, fix$ortho, out_file = character(0)),
    "out_file must be"
  )
  expect_error(
    find_coexpressologs(fix$nets, fix$ortho, out_file = c("a", "b")),
    "out_file must be"
  )
  expect_error(
    find_coexpressologs(fix$nets, fix$ortho, out_file = NA_character_),
    "out_file must be"
  )
})


test_that("find_coexpressologs out_file: empty file if no pair succeeds", {
  fix <- make_coexpr_fixtures()
  out <- withr::local_tempfile(fileext = ".csv")
  # An ortholog table with no rows means no pair produces edges
  empty_ortho <- fix$ortho[0, ]

  result <- suppressWarnings(
    find_coexpressologs(fix$nets, empty_ortho, out_file = out)
  )
  # Documented contract: out_file is always created and always what is
  # returned, even when no pair produces edges.
  expect_identical(result, out)
  expect_true(file.exists(out))
  from_file <- data.table::fread(out)
  expect_equal(nrow(from_file), 0)
  expect_true(all(c("gene1", "gene2", "q.value") %in% names(from_file)))
})


test_that("find_coexpressologs out_file: errors if stale file removal fails", {
  fix <- make_coexpr_fixtures()
  out <- withr::local_tempfile()
  # A non-empty directory at out_file cannot be removed by file.remove(),
  # so this must abort rather than silently proceed to append onto it.
  dir.create(out)
  file.create(file.path(out, "child"))

  expect_error(
    find_coexpressologs(fix$nets, fix$ortho,
      method = "analytical", out_file = out
    ),
    "could not remove existing out_file"
  )
})


test_that("find_coexpressologs out_file: one header across pairs", {
  skip_on_cran()
  set.seed(42)
  n <- 30
  expr1 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", seq_len(n)))
  expr2 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", seq_len(n)))
  expr3 <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  rownames(expr3) <- paste0("C_", sprintf("%03d", seq_len(n)))

  net1 <- compute_network(expr1,
    density = 0.1, mr_log_transform = FALSE,
    sparse = FALSE
  )
  net2 <- compute_network(expr2,
    density = 0.1, mr_log_transform = FALSE,
    sparse = FALSE
  )
  net3 <- compute_network(expr3,
    density = 0.1, mr_log_transform = FALSE,
    sparse = FALSE
  )
  nets <- list(SP_A = net1, SP_B = net2, SP_C = net3)

  n_ortho <- 20
  ortho <- rbind(
    data.frame(
      Species1 = paste0("A_", sprintf("%03d", seq_len(n_ortho))),
      Species2 = paste0("B_", sprintf("%03d", seq_len(n_ortho))),
      hog = paste0("HOG_AB", seq_len(n_ortho))
    ),
    data.frame(
      Species1 = paste0("A_", sprintf("%03d", seq_len(n_ortho))),
      Species2 = paste0("C_", sprintf("%03d", seq_len(n_ortho))),
      hog = paste0("HOG_AC", seq_len(n_ortho))
    ),
    data.frame(
      Species1 = paste0("B_", sprintf("%03d", seq_len(n_ortho))),
      Species2 = paste0("C_", sprintf("%03d", seq_len(n_ortho))),
      hog = paste0("HOG_BC", seq_len(n_ortho))
    )
  )
  out <- withr::local_tempfile(fileext = ".csv")

  set.seed(1)
  in_memory <- find_coexpressologs(nets, ortho, method = "analytical")
  set.seed(1)
  ret <- find_coexpressologs(nets, ortho,
    method = "analytical", out_file = out
  )

  expect_identical(ret, out)
  lines <- readLines(out)
  header_lines <- grep(
    "^gene1,gene2,species1,species2,hog,q\\.value",
    lines
  )
  # Exactly one header row, no matter how many pairs were written
  expect_equal(length(header_lines), 1)
  expect_equal(header_lines, 1)

  from_file <- as.data.frame(data.table::fread(out))
  # No pair is called on these random networks, so `power` is all NA and
  # fread() reads the empty column back as logical.
  from_file$power <- as.numeric(from_file$power)
  # Row order can differ from the in-memory rbind() order across pairs;
  # compare as sets keyed on the edge identity columns.
  key_cols <- c("gene1", "gene2", "species1", "species2", "hog")
  from_file <- from_file[do.call(order, from_file[key_cols]), ]
  in_memory <- in_memory[do.call(order, in_memory[key_cols]), ]
  rownames(from_file) <- NULL
  rownames(in_memory) <- NULL
  expect_equal(from_file, in_memory, ignore_attr = TRUE)
})


test_that("find_coexpressologs.rcomplex rejects out_file", {
  fix <- make_coexpr_fixtures()
  species <- names(fix$nets)
  traits <- stats::setNames(rep("trait_a", length(species)), species)
  rcx <- rcomplex(species, traits, fix$nets, fix$ortho)
  out <- withr::local_tempfile(fileext = ".csv")

  expect_error(
    find_coexpressologs(rcx, out_file = out),
    "out_file is not supported on the rcomplex method"
  )
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
  expect_true(all(c(
    "multiplier", "eff_density", "n_significant", "edges",
    "species_densities"
  ) %in% names(result)))
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

  set.seed(1)
  result <- suppressMessages(density_sweep(
    networks = fix$nets, orthologs = fix$ortho,
    multipliers = 1.0, method = "analytical"
  ))

  expect_equal(nrow(result), 1)
  expect_equal(result$multiplier, 1.0)

  set.seed(1)
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
    density_sweep(
      list(A = list(network = matrix(0), threshold = 1)),
      dummy_ortho
    ),
    "networks must contain at least 2 species"
  )
  expect_error(
    density_sweep(
      list(
        A = list(network = matrix(0)),
        B = list(network = matrix(0), threshold = 1)
      ),
      dummy_ortho
    ),
    "each network must have 'network' and 'threshold' elements"
  )

  net_a <- list(network = matrix(0, 2, 2), threshold = 1)
  net_b <- list(network = matrix(0, 2, 2), threshold = 1)
  expect_error(
    density_sweep(list(A = net_a, B = net_b), dummy_ortho,
      multipliers = character(0)
    ),
    "multipliers must be a non-empty numeric vector"
  )
  expect_error(
    density_sweep(list(A = net_a, B = net_b), dummy_ortho,
      multipliers = -1
    ),
    "all multipliers must be positive"
  )
  expect_error(
    density_sweep(list(A = net_a, B = net_b), data.frame(x = 1)),
    "orthologs must have columns"
  )
})


# --- Tests for sparse (dgCMatrix) network dispatch ---

test_that("dense_to_dgc builds a dgCMatrix with off-diagonal entries >= thr", {
  m <- matrix(0, 4, 4)
  dimnames(m) <- list(letters[1:4], letters[1:4])
  m[1, 2] <- m[2, 1] <- 10
  m[3, 4] <- m[4, 3] <- 2
  diag(m) <- 10

  s <- dense_to_dgc(m, 5)
  expect_s4_class(s, "dgCMatrix")
  expect_equal(dim(s), c(4L, 4L))
  expect_equal(dimnames(s), dimnames(m))
  # diagonal dropped, sub-threshold pair dropped
  expect_equal(length(s@x), 2L)
  expect_equal(as.matrix(s)[1, 2], 10)
  expect_equal(as.matrix(s)[2, 1], 10)
  expect_equal(as.matrix(s)[3, 4], 0)
  expect_equal(as.matrix(s)[1, 1], 0)

  # all-zero input -> empty dgCMatrix, no error
  e <- dense_to_dgc(matrix(0, 3, 3), 1)
  expect_s4_class(e, "dgCMatrix")
  expect_equal(length(e@x), 0L)
})


test_that("sparse compare_neighborhoods equals dense (compute_network nets)", {
  td <- make_cmp_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)

  expect_true(.net_is_sparse(net1_s))
  expect_false(.net_is_sparse(td$net1))

  res_d <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  res_s <- compare_neighborhoods(net1_s, net2_s, td$ortho)
  expect_equal(res_s, res_d)
})


test_that("sparse compare_neighborhoods equals dense (hand-built nets)", {
  n <- 100
  net1_mat <- matrix(0, n, n)
  rownames(net1_mat) <- colnames(net1_mat) <- paste0("A_", sprintf("%03d", 1:n))
  for (j in 2:21) {
    net1_mat[1, j] <- 10
    net1_mat[j, 1] <- 10
  }
  for (j in 52:71) {
    net1_mat[51, j] <- 10
    net1_mat[j, 51] <- 10
  }
  net2_mat <- matrix(0, n, n)
  rownames(net2_mat) <- colnames(net2_mat) <- paste0("B_", sprintf("%03d", 1:n))
  for (j in 2:21) {
    net2_mat[1, j] <- 10
    net2_mat[j, 1] <- 10
    net2_mat[51, j] <- 10
    net2_mat[j, 51] <- 10
  }
  diag(net1_mat) <- diag(net2_mat) <- 10

  net1 <- list(network = net1_mat, threshold = 5)
  net2 <- list(network = net2_mat, threshold = 5)
  ortho <- data.frame(
    Species1 = paste0("A_", sprintf("%03d", 1:n)),
    Species2 = paste0("B_", sprintf("%03d", 1:n)),
    hog = 1:n
  )

  res_d <- compare_neighborhoods(net1, net2, ortho)
  res_s <- compare_neighborhoods(
    sparse_net(net1), sparse_net(net2),
    ortho
  )
  expect_equal(res_s, res_d)
})


test_that("sparse compare_neighborhoods with tighter threshold equals dense", {
  # graded weights 10 / 7 / 4 (see make_graded_nets()): store threshold 5,
  # analysis threshold 8 keeps tier 10 and drops the stored tier 7, so the
  # cut provably removes some overlapping edges and keeps others
  td <- make_graded_nets()
  thr <- 8
  net1_s <- modifyList(sparse_net(td$net1), list(threshold = thr))
  net2_s <- modifyList(sparse_net(td$net2), list(threshold = thr))
  net1_d <- modifyList(td$net1, list(threshold = thr))
  net2_d <- modifyList(td$net2, list(threshold = thr))

  expect_setequal(unique(net1_s$network@x), c(10, 7))
  expect_gt(sum(net1_s$network@x >= thr), 0L)
  expect_lt(sum(net1_s$network@x >= thr), length(net1_s$network@x))

  res_d <- compare_neighborhoods(net1_d, net2_d, td$ortho)
  res_s <- compare_neighborhoods(net1_s, net2_s, td$ortho)
  expect_equal(res_s, res_d)

  # overlapping edges survive the cut, and the cut removed some
  res_base <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  expect_gt(sum(res_d$Species1.neigh.overlap), 0L)
  expect_gt(sum(res_d$Species2.neigh.overlap), 0L)
  expect_lt(
    sum(res_d$Species1.neigh.overlap),
    sum(res_base$Species1.neigh.overlap)
  )
  expect_lt(
    sum(res_d$Species2.neigh.overlap),
    sum(res_base$Species2.neigh.overlap)
  )
})


test_that("compare_neighborhoods rejects non-dgCMatrix Matrix classes", {
  td <- make_cmp_nets()
  n <- 20
  m <- matrix(0, n, n)
  rownames(m) <- colnames(m) <- paste0("A_", sprintf("%03d", 1:n))
  m[1, 2] <- m[2, 1] <- 10

  # symmetric input -> dsCMatrix; general triplet form -> dgTMatrix
  net_dsc <- list(network = Matrix::Matrix(m, sparse = TRUE), threshold = 5)
  expect_s4_class(net_dsc$network, "dsCMatrix")
  net_dgt <- list(
    network = methods::as(methods::as(m, "generalMatrix"), "TsparseMatrix"),
    threshold = 5
  )
  expect_s4_class(net_dgt$network, "dgTMatrix")

  net2_s <- sparse_net(td$net2)
  ortho <- data.frame(Species1 = "A_001", Species2 = "B_001", hog = 1)
  ortho_rev <- data.frame(Species1 = "B_001", Species2 = "A_001", hog = 1)

  expect_error(
    compare_neighborhoods(net_dsc, net2_s, ortho),
    "network must be a dgCMatrix"
  )
  expect_error(
    compare_neighborhoods(net_dgt, net2_s, ortho),
    "network must be a dgCMatrix"
  )
  expect_error(
    compare_neighborhoods(net2_s, net_dsc, ortho_rev),
    "network must be a dgCMatrix"
  )
})


test_that("compare_neighborhoods rejects mixed dense/sparse inputs", {
  td <- make_cmp_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)
  expect_error(
    compare_neighborhoods(net1_s, td$net2, td$ortho),
    "both dense or both sparse"
  )
  expect_error(
    compare_neighborhoods(td$net1, net2_s, td$ortho),
    "both dense or both sparse"
  )
})


test_that("store guard errors when threshold is below store_threshold", {
  td <- make_cmp_nets()
  net1_s <- sparse_net(td$net1)
  net2_s <- sparse_net(td$net2)

  # no store_threshold -> fallback guard on the smallest stored value
  net1_s$store_threshold <- NULL
  net1_loose <- modifyList(net1_s, list(threshold = td$net1$threshold * 0.5))
  expect_error(
    compare_neighborhoods(net1_loose, net2_s, td$ortho),
    "below the smallest stored value"
  )

  # exact equality with the smallest stored value passes
  thr_min <- min(net1_s$network@x)
  net1_min_s <- modifyList(net1_s, list(threshold = thr_min))
  net1_min_d <- modifyList(td$net1, list(threshold = thr_min))
  expect_equal(
    compare_neighborhoods(net1_min_s, net2_s, td$ortho),
    compare_neighborhoods(net1_min_d, td$net2, td$ortho)
  )

  net1_guarded <- modifyList(net1_s, list(
    store_threshold = td$net1$threshold, store_density = 0.1,
    threshold = td$net1$threshold * 0.5
  ))
  expect_error(
    compare_neighborhoods(net1_guarded, net2_s, td$ortho),
    "\\(store_density = 0.1\\); recompute"
  )

  # store_density absent -> no parenthetical in the message
  net1_guarded$store_density <- NULL
  expect_error(
    compare_neighborhoods(net1_guarded, net2_s, td$ortho),
    "stored threshold [^(]*; recompute"
  )

  # threshold at or above store_threshold passes
  net1_ok <- modifyList(net1_guarded, list(threshold = td$net1$threshold))
  expect_equal(
    compare_neighborhoods(net1_ok, net2_s, td$ortho),
    compare_neighborhoods(td$net1, td$net2, td$ortho)
  )
})


test_that("sparse compare_neighborhoods rejects malformed dgCMatrix slots", {
  td <- make_cmp_nets()
  net2_s <- sparse_net(td$net2)
  good <- sparse_net(td$net1)$network
  n <- ncol(good)
  with_net <- function(m) modifyList(td$net1, list(network = m))

  # empty p
  m <- good
  m@p <- integer(0)
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "slot p must have length >= 1"
  )

  # p[n] != length(x)  # nolint
  m <- good
  m@p[n + 1L] <- m@p[n + 1L] - 1L
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "p\\[n\\] = .* length\\(x\\) = "
  )

  # length(i) != length(x)  # nolint
  m <- good
  m@i <- m@i[-1L]
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "length\\(i\\) = "
  )

  # duplicate row indices within a column (sparseMatrix(check = FALSE) sums
  # duplicates in Matrix >= 1.5, so build them by slot assignment)
  m <- good
  first_col <- which(diff(m@p) >= 2L)[1L]
  k <- m@p[first_col] + 1L
  m@i[k + 1L] <- m@i[k]
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "strictly increasing"
  )

  # row index out of range
  m <- good
  m@i[length(m@i)] <- n
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "strictly increasing"
  )

  # a p entry above length(x) with intact endpoints: p must be verified
  # non-decreasing for the WHOLE vector before any row-index scan, or the
  # validator itself reads i[] out of bounds for the inflated column
  m <- good
  m@p[2L] <- length(m@x) + 50L
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "non-decreasing"
  )

  # a well-formed network still passes
  expect_equal(
    compare_neighborhoods(with_net(good), net2_s, td$ortho),
    compare_neighborhoods(td$net1, td$net2, td$ortho)
  )
})


test_that("sparse compare_neighborhoods edge cases equal dense", {
  # diagonal = 1, values exactly at thr (>= boundary), a stored entry below
  # thr (store threshold < analysis threshold), an isolated gene, and an
  # entry below the store threshold
  thr <- 0.5
  thr_store <- 0.3
  n <- 8
  build <- function(prefix, extra = NULL) {
    m <- matrix(0, n, n)
    rownames(m) <- colnames(m) <- paste0(prefix, 1:n)
    set <- function(a, b, v) m[a, b] <<- m[b, a] <<- v
    set(1, 2, 0.9)
    set(1, 3, 0.5) # exactly at thr
    set(1, 4, 0.4) # stored (>= thr_store) but below thr
    set(2, 3, 0.9)
    set(5, 6, 0.5) # exactly at thr
    set(5, 8, 0.2) # below thr_store: never stored
    if (!is.null(extra)) set(extra[1], extra[2], extra[3])
    diag(m) <- 1 # gene 7 isolated
    m
  }
  net1 <- list(network = build("A"), threshold = thr)
  net2 <- list(network = build("B", extra = c(2, 4, 0.5)), threshold = thr)
  ortho <- data.frame(
    Species1 = paste0("A", 1:n), Species2 = paste0("B", 1:n), hog = 1:n
  )

  net1_s <- sparse_net(net1, thr_store)
  net2_s <- sparse_net(net2, thr_store)

  # store exercised: sub-threshold entry kept, diagonal dropped, empty column
  expect_gt(length(net1_s$network@x), sum(net1_s$network@x >= thr))
  expect_equal(min(net1_s$network@x), 0.4)
  expect_equal(diff(net1_s$network@p)[7L], 0L)
  expect_equal(unname(Matrix::diag(net1_s$network)), rep(0, n))

  res_d <- compare_neighborhoods(net1, net2, ortho)
  res_s <- compare_neighborhoods(net1_s, net2_s, ortho)
  expect_equal(res_s, res_d)

  expect_equal(res_d$Species1.neigh, c(2, 2, 2, 0, 1, 1, 0, 0))
  expect_equal(res_d$Species2.neigh, c(2, 3, 2, 1, 1, 1, 0, 0))
})


test_that("sparse networks must be square with identical row/col names", {
  td <- make_cmp_nets()
  net2_s <- sparse_net(td$net2)
  good <- sparse_net(td$net1)$network
  with_net <- function(m) modifyList(td$net1, list(network = m))

  # non-square (extra column), row names intact
  wide <- Matrix::sparseMatrix(
    i = good@i + 1L, j = rep.int(seq_len(ncol(good)), diff(good@p)),
    x = good@x, dims = c(nrow(good), ncol(good) + 1L),
    dimnames = list(rownames(good), c(colnames(good), "extra"))
  )
  expect_s4_class(wide, "dgCMatrix")
  expect_error(
    compare_neighborhoods(with_net(wide), net2_s, td$ortho),
    "square"
  )

  # missing column names
  m <- good
  dimnames(m) <- list(rownames(good), NULL)
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "row and column names"
  )

  # row/col names differ
  m <- good
  dimnames(m) <- list(rownames(good), rev(colnames(good)))
  expect_error(
    compare_neighborhoods(with_net(m), net2_s, td$ortho),
    "row and column names"
  )
})


test_that(
  "find_coexpressologs (analytical) equals dense with sparse networks",
  {
    td <- make_cmp_nets()
    nets_d <- list(A = td$net1, B = td$net2)
    nets_s <- list(A = sparse_net(td$net1), B = sparse_net(td$net2))

    set.seed(1)
    res_d <- find_coexpressologs(nets_d, td$ortho, method = "analytical")
    set.seed(1)
    res_s <- find_coexpressologs(nets_s, td$ortho, method = "analytical")
    expect_equal(res_s, res_d)
    expect_gt(nrow(res_d), 0L)
  }
)


# ---- pi0_method / pval_combine pass-through (D4, D2) ----

test_that(
  "find_coexpressologs passes pi0_method through to summarize_comparison",
  {
    td <- make_graded_nets()
    nets <- list(A = td$net1, B = td$net2)
    cmp <- compare_neighborhoods(td$net1, td$net2, td$ortho)
    pool <- cmp$Species1.neigh.overlap > 0 & cmp$Species2.neigh.overlap > 0
    cmp <- cmp[pool, ]
    # default combine is "max" (D2, reciprocal criterion)
    bh <- pmax(
      p.adjust(cmp$Species1.p.val.con, "BH"),
      p.adjust(cmp$Species2.p.val.con, "BH")
    )

    edges <- find_coexpressologs(nets, td$ortho, pi0_method = "none")
    idx <- match(
      paste(cmp$Species1, cmp$Species2),
      paste(edges$gene1, edges$gene2)
    )
    expect_false(anyNA(idx))
    expect_equal(edges$q.value[idx], bh)

    # pi0_method = "none" leaves the RNG untouched
    set.seed(9)
    u <- runif(1)
    set.seed(9)
    invisible(find_coexpressologs(nets, td$ortho, pi0_method = "none"))
    expect_identical(runif(1), u)
  }
)


test_that("comparison_to_edges combines directional q-values by min or max", {
  comp <- data.frame(
    Species1 = c("A1", "A2", "A3"), Species2 = c("B1", "B2", "B3"),
    hog = 1:3,
    Species1.effect.size = c(4, 1, 2), Species2.effect.size = c(9, 1, 2),
    Species1.jaccard = c(0.8, 0, 0.5), Species2.jaccard = c(0.5, 0, 0.5),
    Species1.q.val.con = c(0.01, 0.80, 0.03),
    Species2.q.val.con = c(0.03, 0.90, 0.20)
  )
  e_min <- comparison_to_edges(comp, "SP_A", "SP_B", pval_combine = "min")
  e_def <- comparison_to_edges(comp, "SP_A", "SP_B")
  e_max <- comparison_to_edges(comp, "SP_A", "SP_B", pval_combine = "max")
  expect_identical(e_max, e_def)
  expect_equal(e_min$q.value, c(0.01, 0.80, 0.03))
  expect_equal(e_max$q.value, c(0.03, 0.90, 0.20))
  expect_equal(e_min$type, c("conserved", "ns", "conserved"))
  expect_equal(e_max$type, c("conserved", "ns", "ns"))
  expect_error(comparison_to_edges(comp, "SP_A", "SP_B",
    pval_combine = "mean"
  ))

  # NA in one direction: the other direction's value is used either way
  comp$Species2.q.val.con[1] <- NA
  expect_equal(comparison_to_edges(comp, "SP_A", "SP_B",
                 pval_combine = "max"
               )$q.value[1], 0.01)
  expect_equal(comparison_to_edges(comp, "SP_A", "SP_B",
                 pval_combine = "min"
               )$q.value[1], 0.01)
})


test_that(
  "summarize_comparison and find_coexpressologs pass pval_combine through",
  {
    td <- make_graded_nets()
    cmp <- compare_neighborhoods(td$net1, td$net2, td$ortho)
    s <- summarize_comparison(cmp,
      sp1 = "A", sp2 = "B", pi0_method = "none",
      pval_combine = "max"
    )
    expect_equal(
      s$edges,
      comparison_to_edges(s$results, "A", "B", pval_combine = "max")
    )
    expect_equal(
      s$edges$q.value,
      pmax(s$results$Species1.q.val.con, s$results$Species2.q.val.con)
    )

    nets <- list(A = td$net1, B = td$net2)
    e_min <- find_coexpressologs(nets, td$ortho,
      pi0_method = "none",
      pval_combine = "min"
    )
    e_max <- find_coexpressologs(nets, td$ortho,
      pi0_method = "none",
      pval_combine = "max"
    )
    expect_equal(e_max$q.value, s$edges$q.value)
    expect_true(all(e_max$q.value >= e_min$q.value))
    called_min <- paste(e_min$gene1, e_min$gene2)[e_min$type == "conserved"]
    called_max <- paste(e_max$gene1, e_max$gene2)[e_max$type == "conserved"]
    expect_gt(length(called_max), 0L)
    expect_true(all(called_max %in% called_min))
  }
)


test_that("density_sweep forwards pi0_method and pval_combine (D2)", {
  # Asymmetric network sizes: N differs per direction, so the two
  # directional q-values (and hence min vs max combine) differ.
  n_a <- 40
  n_b <- 25
  ga <- paste0("A", seq_len(n_a))
  gb <- paste0("B", seq_len(n_b))
  mat_a <- matrix(0, n_a, n_a, dimnames = list(ga, ga))
  mat_a[1:12, 1:12] <- 10
  diag(mat_a) <- 0
  mat_b <- matrix(0, n_b, n_b, dimnames = list(gb, gb))
  mat_b[1:8, 1:8] <- 10
  diag(mat_b) <- 0
  nets <- list(
    SP_A = list(network = mat_a, threshold = 5),
    SP_B = list(network = mat_b, threshold = 5)
  )
  ortho <- data.frame(
    Species1 = paste0("A", 1:25),
    Species2 = paste0("B", 1:25),
    hog = paste0("HOG", 1:25),
    stringsAsFactors = FALSE
  )

  # unknown values must error, not vanish into `...`
  expect_error(suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, pi0_method = "bogus"
  )))
  expect_error(suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, pval_combine = "mean"
  )))

  sw_min <- suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, method = "analytical",
    pi0_method = "none", pval_combine = "min"
  ))
  sw_max <- suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, method = "analytical",
    pi0_method = "none", pval_combine = "max"
  ))
  e_min <- find_coexpressologs(nets, ortho,
    pi0_method = "none",
    pval_combine = "min"
  )
  e_max <- find_coexpressologs(nets, ortho,
    pi0_method = "none",
    pval_combine = "max"
  )
  expect_equal(sw_min$edges[[1]], e_min)
  expect_equal(sw_max$edges[[1]], e_max)
  expect_true(any(sw_max$edges[[1]]$q.value != sw_min$edges[[1]]$q.value))
  expect_equal(sw_min$n_significant, sum(e_min$type == "conserved"))
  expect_equal(sw_max$n_significant, sum(e_max$type == "conserved"))

  # pi0_method = "none" must leave the global RNG untouched
  set.seed(9)
  u <- runif(1)
  set.seed(9)
  invisible(suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, method = "analytical",
    pi0_method = "none"
  )))
  expect_identical(runif(1), u)
})


test_that("default pval_combine is 'max' (D2, Netotea reciprocal criterion)", {
  comp <- data.frame(
    Species1 = c("A1", "A2"), Species2 = c("B1", "B2"), hog = 1:2,
    Species1.effect.size = c(4, 2), Species2.effect.size = c(9, 2),
    Species1.jaccard = c(0.8, 0.5), Species2.jaccard = c(0.5, 0.5),
    Species1.q.val.con = c(0.01, 0.03),
    Species2.q.val.con = c(0.03, 0.20)
  )
  expect_identical(
    comparison_to_edges(comp, "SP_A", "SP_B"),
    comparison_to_edges(comp, "SP_A", "SP_B", pval_combine = "max")
  )

  # asymmetric network sizes so the two directional q-values differ
  n_a <- 40
  n_b <- 25
  ga <- paste0("A", seq_len(n_a))
  gb <- paste0("B", seq_len(n_b))
  mat_a <- matrix(0, n_a, n_a, dimnames = list(ga, ga))
  mat_a[1:12, 1:12] <- 10
  diag(mat_a) <- 0
  mat_b <- matrix(0, n_b, n_b, dimnames = list(gb, gb))
  mat_b[1:8, 1:8] <- 10
  diag(mat_b) <- 0
  nets <- list(
    SP_A = list(network = mat_a, threshold = 5),
    SP_B = list(network = mat_b, threshold = 5)
  )
  ortho <- data.frame(
    Species1 = paste0("A", 1:25),
    Species2 = paste0("B", 1:25),
    hog = paste0("HOG", 1:25),
    stringsAsFactors = FALSE
  )

  cmp <- compare_neighborhoods(nets$SP_A, nets$SP_B, ortho)
  s_def <- summarize_comparison(cmp,
    sp1 = "SP_A", sp2 = "SP_B",
    pi0_method = "none"
  )
  s_max <- summarize_comparison(cmp,
    sp1 = "SP_A", sp2 = "SP_B",
    pi0_method = "none", pval_combine = "max"
  )
  expect_identical(s_def$edges, s_max$edges)

  e_def <- find_coexpressologs(nets, ortho, pi0_method = "none")
  e_max <- find_coexpressologs(nets, ortho,
    pi0_method = "none",
    pval_combine = "max"
  )
  e_min <- find_coexpressologs(nets, ortho,
    pi0_method = "none",
    pval_combine = "min"
  )
  expect_identical(e_def, e_max)
  expect_true(any(e_def$q.value != e_min$q.value))

  sw_def <- suppressMessages(density_sweep(
    nets, ortho,
    multipliers = 1.0, method = "analytical",
    pi0_method = "none"
  ))
  expect_equal(sw_def$edges[[1]], e_max)
})


# ---- seed threading (analytical path) -----------------------------------
# find_coexpressologs() and density_sweep() seed once for the whole loop
# and leave `seed` at NULL in the calls below them, so pairs and
# multipliers draw in sequence from one stream instead of every level
# reusing the same uniforms.

test_that("a seed on find_coexpressologs reaches summarize_comparison", {
  td <- make_graded_nets()
  nets <- list(SP_A = td$net1, SP_B = td$net2)
  cmp <- compare_neighborhoods(td$net1, td$net2, td$ortho)

  s <- summarize_comparison(cmp,
    alternative = "greater", alpha = 0.05,
    seed = 11
  )
  direct <- comparison_to_edges(s$results, "SP_A", "SP_B", "greater",
    0.05,
    pval_combine = "max"
  )
  wrapped <- find_coexpressologs(nets, td$ortho,
    method = "analytical", seed = 11
  )
  expect_equal(wrapped$q.value, direct$q.value)
  expect_identical(wrapped$type, direct$type)
})


test_that("find_coexpressologs(seed = ) is reproducible; seeds differ", {
  td <- make_graded_nets()
  nets <- list(SP_A = td$net1, SP_B = td$net2)

  a <- find_coexpressologs(nets, td$ortho, seed = 5)
  expect_identical(find_coexpressologs(nets, td$ortho, seed = 5), a)
  b <- find_coexpressologs(nets, td$ortho, seed = 6)
  expect_false(identical(a$q.value, b$q.value))

  # seed = NULL reproduces the old behaviour exactly
  set.seed(5)
  expect_identical(find_coexpressologs(nets, td$ortho), a)
})


test_that("find_coexpressologs restores the caller's stream", {
  td <- make_graded_nets()
  nets <- list(SP_A = td$net1, SP_B = td$net2)
  set.seed(7)
  before <- get(".Random.seed", envir = globalenv())
  invisible(find_coexpressologs(nets, td$ortho, seed = 42))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  # and it is genuinely restoration, not a coincidence of pinning: the
  # exit state must not be the one set.seed(42) would give.
  set.seed(42)
  expect_false(identical(before, get(".Random.seed", envir = globalenv())))
})


test_that("density_sweep(seed = ) is reproducible and pins level one", {
  td <- make_graded_nets()
  nets <- list(SP_A = td$net1, SP_B = td$net2)
  run <- function(s) {
    suppressMessages(density_sweep(nets, td$ortho,
      multipliers = c(1, 1.02), method = "analytical", seed = s
    ))
  }
  a <- run(3)
  expect_identical(run(3), a)
  expect_false(identical(a$edges[[1]]$q.value, run(4)$edges[[1]]$q.value))

  # nothing is drawn before the first multiplier, so level one reproduces
  # a bare find_coexpressologs() at the same seed
  direct <- find_coexpressologs(nets, td$ortho,
    method = "analytical", seed = 3
  )
  expect_equal(a$edges[[1]]$q.value, direct$q.value)

  set.seed(7)
  before <- get(".Random.seed", envir = globalenv())
  invisible(run(42))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})


# --- Detection power of each edge (#12) ---

# Six co-expression modules shared by two species, gene loadings graded
# so neighbourhood size (and with it detection power) varies by gene.
make_power_comparison <- function() {
  withr::local_seed(7)
  n <- 120
  s <- 30
  mod <- rep(1:6, each = 20)
  load <- rep(seq(0.3, 1.5, length.out = 20), 6)
  sim <- function(prefix) {
    f <- matrix(stats::rnorm(6 * s), 6, s)
    e <- f[mod, ] * load + matrix(stats::rnorm(n * s), n, s)
    rownames(e) <- paste0(prefix, seq_len(n))
    e
  }
  n1 <- compute_network(sim("A"), density = 0.08, sparse = FALSE)
  n2 <- compute_network(sim("B"), density = 0.08, sparse = FALSE)
  ortho <- data.frame(
    Species1 = paste0("A", seq_len(n)), Species2 = paste0("B", seq_len(n)),
    hog = paste0("H", seq_len(n))
  )
  cmp <- compare_neighborhoods(n1, n2, ortho)
  list(
    res = summarize_comparison(cmp, pi0_method = "none")$results,
    np = n - 1
  )
}

# Gated P(X >= x): the reported conservation p-value, 1 for x <= 1.
power_gated_p <- function(x, m, np, k) {
  ifelse(x <= 1, 1, stats::phyper(x - 1, m, np - m, k, lower.tail = FALSE))
}

# One direction's power by exhaustive search over every overlap.
power_brute <- function(res, d, np, alpha, called, f0 = NULL) {
  g <- function(s) res[[paste0(d, s)]]
  m <- g(".neigh")
  k <- g(".ortho.neigh")
  x <- g(".neigh.overlap")
  pcut <- max(g(".p.val.con")[g(".q.val.con") < alpha])
  if (is.null(f0)) f0 <- stats::median((x / m)[called & m > 0])
  vapply(seq_along(m), function(i) {
    nn <- min(k[i], m[i])
    if (nn == 0) {
      return(0)
    }
    xs <- 0:nn
    ok <- which(power_gated_p(xs, m[i], np, k[i]) <= pcut)
    if (length(ok) == 0L) {
      return(0)
    }
    sum(stats::dbinom(xs[min(ok)]:nn, nn, f0))
  }, numeric(1))
}


test_that(".power_xstar matches brute force, tiny p-values included", {
  np <- 999
  for (m in c(5, 30, 50)) {
    for (k in c(4, 40)) {
      n <- min(k, m)
      for (x0 in 2:n) {
        p0 <- power_gated_p(x0, m, np, k)
        for (pcut in c(p0, p0 * (1 - 1e-9), p0 * (1 + 1e-9))) {
          hits <- which(power_gated_p(0:n, m, np, k) <= pcut)
          want <- if (length(hits) == 0L) n + 1 else (0:n)[min(hits)]
          expect_equal(rcomplex:::.power_xstar(pcut, m, k, np), want)
        }
      }
    }
  }
  # qhyper(lower.tail = FALSE) returns 17 here; the true cut is 30.
  p30 <- power_gated_p(30, 50, np, 40)
  expect_lt(p30, 1e-30)
  expect_equal(rcomplex:::.power_xstar(p30, 50, 40, np), 30)
})


test_that("edge power matches a brute-force computation", {
  fx <- make_power_comparison()
  res <- fx$res
  q1 <- res$Species1.q.val.con
  q2 <- res$Species2.q.val.con
  called_max <- pmax(q1, q2) < 0.05
  called_min <- pmin(q1, q2) < 0.05
  expect_gt(sum(called_max), 10L)

  dir_pw <- function(called, f0 = NULL) {
    lapply(c("Species1", "Species2"), power_brute,
      res = res, np = fx$np, alpha = 0.05, called = called, f0 = f0
    )
  }
  b <- dir_pw(called_max)
  pw <- rcomplex:::.edge_power(res, 0.05, "greater", "max")
  expect_equal(pw, pmin(b[[1]], b[[2]]), tolerance = 1e-12)
  # Power has to vary for it to separate anything.
  expect_gt(diff(range(pw)), 0.3)

  # f0 overrides the median in both directions.
  b3 <- dir_pw(called_max, f0 = 0.3)
  expect_equal(
    rcomplex:::.edge_power(res, 0.05, "greater", "max", f0 = 0.3),
    pmin(b3[[1]], b3[[2]]),
    tolerance = 1e-12
  )
  # "min" takes the stronger direction, with f0 from the pairs "min" calls.
  expect_equal(
    rcomplex:::.edge_power(res, 0.05, "greater", "min", f0 = 0.3),
    pmax(b3[[1]], b3[[2]]),
    tolerance = 1e-12
  )
  bm <- dir_pw(called_min)
  expect_equal(
    rcomplex:::.edge_power(res, 0.05, "greater", "min"),
    pmax(bm[[1]], bm[[2]]),
    tolerance = 1e-12
  )

  edges <- comparison_to_edges(res, "SP_A", "SP_B", f0 = 0.3)
  expect_equal(edges$power, pmin(b3[[1]], b3[[2]]), tolerance = 1e-12)
})


test_that("edge power is NA where it is undefined", {
  fx <- make_power_comparison()
  res <- fx$res
  nas <- rep(NA_real_, nrow(res))

  expect_equal(rcomplex:::.edge_power(res, 0.05, "less", "max"), nas)
  # No called pair: no cut to measure power against.
  none <- res
  none$Species1.q.val.con <- 1
  none$Species2.q.val.con <- 1
  expect_equal(rcomplex:::.edge_power(none, 0.05, "greater", "max"), nas)
  # Effect sizes that disagree on the urn size.
  bad <- res
  hit <- which(bad$Species1.neigh.overlap > 0)[1]
  bad$Species1.effect.size[hit] <- bad$Species1.effect.size[hit] * 2
  expect_warning(
    pw <- rcomplex:::.edge_power(bad, 0.05, "greater", "max"),
    "urn size for Species1"
  )
  expect_equal(pw, nas)

  expect_error(
    comparison_to_edges(res, "SP_A", "SP_B", f0 = 0),
    "f0 must be NULL or a single number"
  )
  expect_error(
    comparison_to_edges(res, "SP_A", "SP_B", f0 = c(0.2, 0.3)),
    "f0 must be NULL"
  )
})


test_that("find_coexpressologs carries power on both paths", {
  fx <- make_clique_fixture()
  an <- find_coexpressologs(fx$networks, fx$orthologs,
    pi0_method = "none", f0 = 0.5
  )
  expect_equal(
    names(an)[7:10], c("effect_size", "jaccard", "power", "type")
  )
  expect_true(any(!is.na(an$power)))

  sweep <- suppressMessages(density_sweep(fx$networks, fx$orthologs,
    multipliers = 1, method = "analytical", pi0_method = "none", f0 = 0.5
  ))
  expect_equal(sweep$edges[[1]]$power, an$power)

  perm <- find_coexpressologs(fx$networks, fx$orthologs,
    method = "permutation", min_exceedances = 5L,
    max_permutations = 200L, seed = 1L
  )
  expect_equal(names(perm), names(an))
  expect_true(all(is.na(perm$power)))

  expect_error(
    find_coexpressologs(fx$networks, fx$orthologs, f0 = 2),
    "f0 must be NULL"
  )
})
