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
  skip_on_os("windows") # mclapply falls back to serial
  # r-lib's check-r-package action sets NOT_CRAN=true (so skip_on_cran()
  # does not fire here) but still runs R CMD check --as-cran, which sets
  # _R_CHECK_LIMIT_CORES_ and makes parallel:::.check_ncores() error on
  # any mclapply() call requesting more than 2 cores -- exactly what the
  # n_cores = 3 run below does.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  skip_if(
    nzchar(chk) && chk != "false",
    "R CMD check --as-cran limits mclapply() to 2 cores"
  )

  net <- make_ambiguous_net()
  res <- seq(0.25, 2.5, by = 0.25)
  run <- function(nc) {
    detect_modules(net,
      resolution = res, objective_function = "modularity",
      seed = 42L, n_cores = nc, max_consensus_iter = 10L,
      test_k1 = TRUE, n_perm_k1 = 100L, alpha_k1 = 0.05
    )
  }

  r1 <- run(1L)
  r2a <- run(2L)
  r2b <- run(2L)
  r3 <- run(3L)

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

test_that("forked K = 1 workers survive the parent's OpenMP threads", {
  skip_on_cran()
  skip_on_os("windows") # mclapply falls back to serial
  # Regression for a Linux-only deadlock. The parent runs the
  # co-classification scan with n_cores = 2, which starts a libgomp thread
  # pool; the K = 1 test then forks workers, and a worker that enters any
  # OpenMP region inherits that pool without its threads and blocks for
  # good. The package's own kernels skip OpenMP at n_cores = 1, but
  # Armadillo parallelised the dense x sparse product inside eigs_sym()
  # on its own, until src/Makevars set ARMA_DONT_USE_OPENMP. macOS (LLVM
  # libomp) never hung, so only Linux CI can fail this. Unlike the
  # core-count test above, n_cores = 2 stays inside the limit
  # R CMD check --as-cran imposes, so R-CMD-check runs it too. These are
  # that test's run(2L) arguments, the call that hung on ubuntu CI.
  net <- make_ambiguous_net()
  res <- detect_modules(net,
    resolution = seq(0.25, 2.5, by = 0.25),
    objective_function = "modularity",
    seed = 42L, n_cores = 2L, max_consensus_iter = 10L,
    test_k1 = TRUE, n_perm_k1 = 100L, alpha_k1 = 0.05
  )
  expect_gt(res$k1_test$n_perm_completed, 0L)
})

test_that("detect_modules leaves the ambient RNG stream core-count invariant", {
  skip_on_cran()
  skip_on_os("windows")

  net <- make_ambiguous_net(n = 150L)
  res <- c(0.5, 1.0, 1.5)
  after <- function(nc) {
    set.seed(99L)
    detect_modules(net,
      resolution = res, objective_function = "modularity",
      seed = 42L, n_cores = nc, test_k1 = FALSE
    )
    runif(1L)
  }
  expect_identical(after(1L), after(2L))
})


test_that("a seeded call restores the caller's stream", {
  # The single-resolution path used to leave the global stream wherever
  # the clustering backend stopped, which is backend- and build-
  # dependent. Anything drawn afterwards without its own set.seed() --
  # summarize_comparison()'s randomized-p pi0, for one -- then started
  # from an unpredictable position. Pinning the exit state at
  # set.seed(seed) fixed that but replaced it with a second problem: the
  # downstream draw then depended on THIS call's seed rather than on the
  # caller's. Restoring fixes both, and both paths agree.
  set.seed(1)
  e <- matrix(stats::rnorm(120 * 8), 120, 8)
  rownames(e) <- paste0("g", seq_len(120))
  net <- compute_network(e, density = 0.1, sparse = FALSE)

  restored <- function(f) {
    set.seed(7)
    before <- get(".Random.seed", envir = globalenv())
    f()
    identical(before, get(".Random.seed", envir = globalenv()))
  }

  expect_true(restored(function() {
    detect_modules(net,
      resolution = 1.0, seed = 42,
      objective_function = "modularity"
    )
  }))
  expect_true(restored(function() {
    detect_modules(net,
      resolution = c(0.8, 1.0), seed = 42,
      objective_function = "modularity",
      n_iterations = 1L, max_consensus_iter = 1L
    )
  }))

  # seed = NULL must still advance the stream, or consecutive unseeded
  # calls would return the same partition.
  set.seed(7)
  before <- get(".Random.seed", envir = globalenv())
  invisible(detect_modules(net,
    resolution = 1.0, seed = NULL,
    objective_function = "modularity"
  ))
  expect_false(identical(before, get(".Random.seed", envir = globalenv())))
})


# ---- K = 1 stopping rule ------------------------------------------------
# The rule these replace stopped after the first batch in every run: with
# ceiling(1 / alpha) permutations done, zero exceedances ended the test and
# so did one, for any alpha < 0.5. n_perm_k1 above that grid point changed
# nothing, and one early exceedance sank a network that the full budget
# would have called structured.

test_that("a stop needs the outstanding permutations to be irrelevant", {
  settled <- rcomplex:::.k1_settled

  # one exceedance out of the first 20 still leaves p = 2/101 in reach
  expect_false(settled(1L, 20L, 100L, 0.05))
  expect_false(settled(4L, 20L, 100L, 0.05))
  # five cannot: even a perfect remaining run finishes at 6/101 > 0.05
  expect_true(settled(5L, 20L, 100L, 0.05))
  # a clean first batch is not yet a decision either
  expect_false(settled(0L, 20L, 100L, 0.05))
  expect_true(settled(0L, 100L, 100L, 0.05))
  # the decision the rule protects is the one the full budget would make
  expect_true(settled(2L, 60L, 100L, 0.02))
  expect_false(settled(0L, 60L, 100L, 0.02))
})


test_that("n_perm_k1 sets the resolution of the K = 1 p-value", {
  skip_on_cran()

  net <- make_ambiguous_net(n = 200L)
  run <- function(np) {
    detect_modules(net,
      resolution = c(0.5, 1.0, 1.5, 2.0),
      objective_function = "modularity", seed = 42L,
      test_k1 = TRUE, n_perm_k1 = np, alpha_k1 = 0.05,
      max_consensus_iter = 5L
    )$k1_test
  }
  k20 <- run(20L)
  k100 <- run(100L)

  expect_identical(k20$n_perm_completed, 20L)
  expect_identical(k100$n_perm_completed, 100L)
  expect_equal(k20$p_value, 1 / 21)
  expect_equal(k100$p_value, 1 / 101)
  expect_true(k100$has_structure)
})


test_that("the K = 1 null is optimised as hard as the observed sweep", {
  skip_on_cran()

  # A null partition found with fewer Leiden iterations than the observed
  # one carries less structure for the same graph, so lambda_null lands
  # low and the test leans toward calling structure that is not there.
  seen <- new.env(parent = emptyenv())
  seen$n_iterations <- integer(0)
  real_leiden <- igraph::cluster_leiden
  recorder <- function(graph, ..., n_iterations = 2L) {
    seen$n_iterations <- c(seen$n_iterations, as.integer(n_iterations))
    real_leiden(graph, ..., n_iterations = n_iterations)
  }

  g <- withr::with_seed(3L, igraph::sample_gnp(60L, 0.15))
  igraph::V(g)$name <- paste0("G", seq_len(60L))
  igraph::E(g)$weight <- withr::with_seed(
    4L, stats::runif(igraph::ecount(g), 0.3, 1.0)
  )
  edge_list_0 <- igraph::as_edgelist(g, names = FALSE) - 1L
  storage.mode(edge_list_0) <- "integer"
  resolutions <- c(0.5, 1.0)
  memberships <- lapply(resolutions, function(r) {
    mem <- igraph::membership(igraph::cluster_leiden(
      g,
      resolution = r, objective_function = "modularity",
      n_iterations = 3L
    ))
    names(mem) <- igraph::V(g)$name
    mem
  })

  testthat::with_mocked_bindings(
    rcomplex:::test_community_structure(
      g, igraph::V(g)$name, resolutions, "modularity", 3L, memberships,
      edge_list_0,
      n_perm = 4L, n_cores = 1L, alpha = 0.05, seed_root = 1L
    ),
    cluster_leiden = recorder, .package = "igraph"
  )

  expect_gt(length(seen$n_iterations), 0L)
  expect_identical(unique(seen$n_iterations), 3L)
})


# ---- consensus convergence ---------------------------------------------

test_that("consensus iteration stops at a fixed point", {
  skip_on_cran()

  # This fixture holds a stable disagreement between the coarsest and the
  # finest resolution: pairwise ARI sticks at 0.9585 from iteration 6 and
  # never reaches 0.999, so the loop used to spend every iteration of
  # max_consensus_iter rebuilding the partition it already had.
  net <- make_ambiguous_net()
  res <- seq(0.25, 2.5, by = 0.25)
  run <- function(cap) {
    detect_modules(net,
      resolution = res,
      objective_function = "modularity", seed = 42L,
      max_consensus_iter = cap, test_k1 = FALSE
    )
  }
  r_short <- run(40L)
  r_long <- run(200L)

  expect_lt(r_short$params$n_consensus_iterations, 40L)
  expect_identical(
    r_short$params$n_consensus_iterations,
    r_long$params$n_consensus_iterations
  )
  # stopping early must not move the answer
  expect_identical(r_short$modules, r_long$modules)
})


test_that(".partition_id compares groupings, not Leiden's labels", {
  pid <- rcomplex:::.partition_id

  # Leiden hands back arbitrary module ids. Two sweeps that found the same
  # grouping under different labels must compare identical, or the `settled`
  # break in detect_modules_consensus() never fires and the loop burns every
  # iteration of max_consensus_iter rebuilding the partition it already had.
  expect_identical(pid(c(3, 3, 7, 7, 1)), pid(c(1, 1, 2, 2, 3)))
  expect_identical(pid(c(9, 4, 4)), pid(c(1, 2, 2)))

  # ...and it must stay injective on the grouping itself, or the loop would
  # stop on a partition that is still moving.
  expect_false(identical(pid(c(1, 1, 2, 2, 3)), pid(c(1, 2, 1, 2, 3))))
  expect_false(identical(pid(c(1, 1, 2)), pid(c(1, 2, 2))))

  # igraph::membership() carries vertex names; identical() is name-sensitive,
  # so the canonical form must drop them.
  a <- stats::setNames(c(2, 2, 5), c("g1", "g2", "g3"))
  b <- stats::setNames(c(1, 1, 4), c("g1", "g2", "g3"))
  expect_identical(pid(a), pid(b))
  expect_null(names(pid(a)))

  # Labels are assigned by first appearance, so the canonical form is always
  # 1..K in order of first occurrence.
  expect_identical(pid(c(8, 3, 8, 3, 5)), c(1L, 2L, 1L, 2L, 3L))
})
