# Regression test: consensus module detection must be bit-reproducible,
# and identical across core counts (see R/rng.R .task_seed()).

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
      seed = 42L, n_cores = nc, max_consensus_iter = 10L
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
})

test_that("forked consensus workers survive the parent's OpenMP threads", {
  skip_on_cran()
  skip_on_os("windows") # mclapply falls back to serial
  # Regression for a Linux-only deadlock. The parent runs the
  # co-classification scan with n_cores = 2, which starts a libgomp thread
  # pool; the next Leiden sweep forks workers, and a worker that enters any
  # OpenMP region inherits that pool without its threads and blocks for
  # good. macOS (LLVM libomp) never hung, so only Linux CI can fail this.
  # Unlike the core-count test above, n_cores = 2 stays inside the limit
  # R CMD check --as-cran imposes, so R-CMD-check runs it too.
  net <- make_ambiguous_net()
  res <- detect_modules(net,
    resolution = seq(0.25, 2.5, by = 0.25),
    objective_function = "modularity",
    seed = 42L, n_cores = 2L, max_consensus_iter = 10L
  )
  expect_gt(res$n_modules, 0L)
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
      seed = 42L, n_cores = nc
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
      objective_function = "modularity", max_consensus_iter = 1L
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
      max_consensus_iter = cap
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
