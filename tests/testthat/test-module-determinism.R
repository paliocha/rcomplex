# Regression test: consensus module detection must be bit-reproducible,
# and identical across core counts (see R/modules.R .task_seed()).

# Deliberately AMBIGUOUS fixture: 24 small blocks nested in 6 super-blocks
# plus a dense background, so Leiden has no unique optimum. A clean planted
# network is useless here -- every seed recovers the same answer and the test
# would pass vacuously against broken code.
make_ambiguous_net <- function(n = 300L, seed = 7L) {
  withr::with_seed(seed, {
    blk <- rep(seq_len(24L), length.out = n)
    super <- ((blk - 1L) %/% 4L) + 1L
    p <- matrix(0.02, n, n)
    p[outer(super, super, "==")] <- 0.10
    p[outer(blk, blk, "==")] <- 0.45
    a <- matrix(stats::runif(n * n), n, n)
    a <- (a + t(a)) / 2
    w <- matrix(stats::runif(n * n, 0.3, 1.0), n, n)
    w <- (w + t(w)) / 2
    m <- (a < p) * w
    diag(m) <- 0
    g <- paste0("G", seq_len(n))
    dimnames(m) <- list(g, g)
    list(network = m, threshold = 0.01)
  })
}

test_that("consensus modules are bit-reproducible across core counts", {
  skip_on_cran()
  skip_on_os("windows")  # mclapply falls back to serial

  net <- make_ambiguous_net()
  res <- seq(0.25, 2.5, by = 0.25)
  run <- function(nc) {
    detect_modules(net, resolution = res, objective_function = "modularity",
                   seed = 42L, n_cores = nc, max_consensus_iter = 10L,
                   test_k1 = TRUE, n_perm_k1 = 100L, alpha_k1 = 0.05)
  }

  r1  <- run(1L)
  r2a <- run(2L)
  r2b <- run(2L)
  r3  <- run(3L)

  # (a) run-to-run stability at a fixed core count
  expect_identical(r2a$modules, r2b$modules)
  # (b) core-count invariance of the partition
  expect_identical(r1$modules, r2a$modules)
  expect_identical(r1$modules, r3$modules)
  # (c) core-count invariance of the K = 1 test, including the stopping grid
  expect_identical(r1$k1_test$n_perm_completed, r3$k1_test$n_perm_completed)
  expect_identical(r1$k1_test$p_value, r3$k1_test$p_value)
  expect_identical(r1$k1_test$lambda_null, r3$k1_test$lambda_null)
})

test_that("detect_modules leaves the ambient RNG stream core-count invariant", {
  skip_on_cran()
  skip_on_os("windows")

  net <- make_ambiguous_net(n = 150L)
  res <- c(0.5, 1.0, 1.5)
  after <- function(nc) {
    set.seed(99L)
    detect_modules(net, resolution = res, objective_function = "modularity",
                   seed = 42L, n_cores = nc, test_k1 = FALSE)
    runif(1L)
  }
  expect_identical(after(1L), after(2L))
})
