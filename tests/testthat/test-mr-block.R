# Tests for mr_block() (P4): exact local reconstruction of mutual-rank
# values for a gene subset, including entries below a sparse network's
# store_threshold. Oracle: the same block of the dense
# compute_network(sparse = FALSE) matrix.


# 40 genes x 25 samples; seed chosen so the query genes' correlation
# ranks are tie-free for both Pearson and Spearman (Spearman rho lies on
# a finite grid: mathematically tied values are broken differently by
# Rfast::cora and stats::cor, which would shift ranks by 0.5)
make_mr_block_expr <- function() {
  set.seed(107)
  x <- matrix(rnorm(40 * 25), nrow = 40, ncol = 25)
  rownames(x) <- paste0("g", sprintf("%02d", 1:40))
  colnames(x) <- paste0("s", seq_len(25))
  x
}

mr_block_genes <- paste0("g", sprintf("%02d", c(3, 7, 11, 25, 40)))


test_that("mr_block matches the dense block for raw MR", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, sparse = FALSE)

  blk <- mr_block(x, mr_block_genes, net)
  expect_equal(blk, as.matrix(net$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)
  expect_identical(dimnames(blk), list(mr_block_genes, mr_block_genes))
})


test_that("mr_block matches the dense block for log MR", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, mr_log_transform = TRUE,
                         sparse = FALSE)

  expect_equal(mr_block(x, mr_block_genes, net),
               as.matrix(net$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)
})


test_that("mr_block matches the dense block with abs_cor", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, abs_cor = TRUE, sparse = FALSE)

  expect_equal(mr_block(x, mr_block_genes, net),
               as.matrix(net$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)
})


test_that("mr_block matches the dense block for Spearman", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, cor_method = "spearman",
                         sparse = FALSE)

  expect_equal(mr_block(x, mr_block_genes, net),
               as.matrix(net$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)
})


test_that("mr_block reconstructs sub-threshold values from a sparse network", {
  x <- make_mr_block_expr()
  dense <- compute_network(x, density = 0.05, sparse = FALSE)
  sp <- compute_network(x, density = 0.05)
  expect_s4_class(sp$network, "dgCMatrix")

  blk <- mr_block(x, mr_block_genes, sp)
  expect_equal(blk, as.matrix(dense$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)

  # the reconstruction must cover values the sparse store discarded
  off <- blk[row(blk) != col(blk)]
  expect_true(any(off > 0 & off < sp$store_threshold))
})


test_that("mr_block uses the network's gene universe, not x's rows", {
  x <- make_mr_block_expr()
  # a constant gene is removed by the min_var filter -> the network's
  # universe is smaller than x's rows; ranks must span network genes only
  x2 <- rbind(x, gz = rep(1, ncol(x)))
  net <- compute_network(x2, density = 0.05, sparse = FALSE)
  expect_equal(net$n_removed, 1L)
  expect_false("gz" %in% rownames(net$network))

  expect_equal(mr_block(x2, mr_block_genes, net),
               as.matrix(net$network[mr_block_genes, mr_block_genes]),
               tolerance = 1e-8)
})


test_that("mr_block handles a single gene", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, sparse = FALSE)

  blk <- mr_block(x, "g05", net)
  expect_equal(blk, matrix(0, 1, 1, dimnames = list("g05", "g05")))
})


test_that("mr_block validates its inputs", {
  x <- make_mr_block_expr()
  net <- compute_network(x, density = 0.05, sparse = FALSE)

  # gene not in the network
  expect_error(mr_block(x, c("g03", "nope"), net), "nope")

  # duplicated genes
  expect_error(mr_block(x, c("g03", "g03"), net), "duplicate")

  # x lacks a network gene (g01 is in the universe, not in the subset)
  expect_error(mr_block(x[-1, , drop = FALSE], mr_block_genes, net),
               "network gene")

  # hand-built network without params
  bare <- list(network = net$network, threshold = net$threshold)
  expect_error(mr_block(x, mr_block_genes, bare), "params")

  # CLR networks cannot be reconstructed by mutual ranks
  clr <- compute_network(x, density = 0.05, norm_method = "CLR",
                         sparse = FALSE)
  expect_error(mr_block(x, mr_block_genes, clr), "MR")
})
