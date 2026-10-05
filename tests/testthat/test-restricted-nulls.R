# null_network(block =) and module_replication().

block_fixture <- function() {
  x <- withr::with_seed(1L, matrix(stats::rnorm(40L * 20L), 40L))
  rownames(x) <- paste0("g", seq_len(40L))
  list(
    x = x, block = rep(paste0("t", 1:5), each = 4L),
    net = compute_network(x, density = 0.1)
  )
}


test_that("the block shuffle keeps block means and moves within blocks", {
  f <- block_fixture()
  xp <- withr::with_seed(1L, .shuffle_genes(f$x, f$block))
  means <- function(m) t(rowsum(t(m), f$block)) / 4
  expect_equal(means(xp), means(f$x))
  for (b in unique(f$block)) {
    j <- f$block == b
    expect_identical(
      t(apply(xp[, j], 1L, sort)), t(apply(f$x[, j], 1L, sort))
    )
  }
  expect_false(identical(xp, f$x))
  # each gene gets its own permutation, so the correlations move
  expect_false(isTRUE(all.equal(stats::cor(t(xp)), stats::cor(t(f$x)))))
})


test_that("null_network(block =) records the permutation space", {
  f <- block_fixture()
  nn <- null_network(f$x, f$net, seed = 1L, block = f$block)
  expect_identical(
    nn$params$block, c(t1 = 4L, t2 = 4L, t3 = 4L, t4 = 4L, t5 = 4L)
  )
  expect_equal(nn$params$log10_perms_per_gene, 5 * log10(24))
  full <- null_network(f$x, f$net, seed = 1L)
  expect_null(full$params$block)
  expect_equal(full$params$log10_perms_per_gene, log10(factorial(20)))
  expect_false(identical(nn$network, full$network))
})


test_that("null_network(block =) validates the grouping", {
  f <- block_fixture()
  expect_error(null_network(f$x, f$net, block = 1:3), "one entry per sample")
  expect_error(
    null_network(f$x, f$net, block = c(NA, f$block[-1])), "NA"
  )
})


# eight planted 50-gene modules over 400 genes; the halves are samples
# 1-10 and 11-20
replication_fixture <- function() {
  x <- withr::with_seed(1L, {
    f <- matrix(stats::rnorm(8L * 20L), 8L)
    matrix(stats::rnorm(400L * 20L), 400L) + 1.5 * f[rep(1:8, each = 50L), ]
  })
  rownames(x) <- paste0("g", seq_len(400L))
  list(
    net_a = compute_network(x[, 1:10], density = 0.05),
    net_b = compute_network(x[, 11:20], density = 0.05)
  )
}


test_that("modules from half A replicate on half B, random sets do not", {
  fx <- replication_fixture()
  m <- detect_modules(fx$net_a, objective_function = "modularity", seed = 1L)
  r <- module_replication(m, fx$net_b, max_draws = 300L, seed = 1L)
  expect_identical(attr(r, "kind"), "replication")
  expect_true(all(r$z > 3))

  g <- rownames(fx$net_b$network)
  rnd <- as_modules(stats::setNames(
    rep(paste0("m", 1:40), each = 10L),
    withr::with_seed(2L, sample(g))
  ))
  r0 <- module_replication(rnd, fx$net_b, max_draws = 300L, seed = 1L)
  expect_lt(abs(mean(r0$z)), 0.5)
})


test_that("the identity table is the default ortholog table", {
  fx <- replication_fixture()
  g <- rownames(fx$net_b$network)
  m <- as_modules(stats::setNames(rep(paste0("m", 1:8), each = 50L), g))
  id <- data.frame(Species1 = g, Species2 = g, hog = g)
  expect_identical(
    module_replication(m, fx$net_b,
      max_draws = 50L, n_null = 50L, seed = 1L
    ),
    module_replication(m, fx$net_b, id,
      max_draws = 50L, n_null = 50L, seed = 1L
    )
  )
})
