# Fork safety: .blas_fork_safe() gates the K = 1 test's fork on the BLAS,
# .check_fork_results() turns failed or crashed workers into an error at
# every fork site. The end-to-end core-count check of the K = 1 test lives
# in test-module-determinism.R.

test_that(".blas_fork_safe() refuses Accelerate unless pinned to one thread", {
  acc <- paste0(
    "/System/Library/Frameworks/Accelerate.framework/Versions/A/",
    "Frameworks/vecLib.framework/Versions/A/libBLAS.dylib"
  )
  expect_false(rcomplex:::.blas_fork_safe(acc, ""))
  expect_false(rcomplex:::.blas_fork_safe(acc, "4"))
  expect_true(rcomplex:::.blas_fork_safe(acc, "1"))
  expect_true(rcomplex:::.blas_fork_safe("/usr/lib/libopenblas.so.0", ""))
  expect_true(rcomplex:::.blas_fork_safe(
    "/Library/Frameworks/R.framework/Resources/lib/libRblas.dylib", ""
  ))
  expect_true(rcomplex:::.blas_fork_safe("", ""))
  expect_type(rcomplex:::.blas_fork_safe(), "logical")
  expect_length(rcomplex:::.blas_fork_safe(), 1L)
})

test_that("VECLIB_MAXIMUM_THREADS set after loading does not unlock forks", {
  acc <- "/System/Library/Frameworks/Accelerate.framework/libBLAS.dylib"
  old <- Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA)
  on.exit(if (is.na(old)) {
    Sys.unsetenv("VECLIB_MAXIMUM_THREADS")
  } else {
    Sys.setenv(VECLIB_MAXIMUM_THREADS = old)
  })
  at_load <- get("veclib_threads", rcomplex:::.load_env)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = if (identical(at_load, "1")) "4" else "1")
  expect_identical(
    rcomplex:::.blas_fork_safe(acc),
    rcomplex:::.blas_fork_safe(acc, at_load)
  )
})

test_that(".check_fork_results() names failed tasks and crashed workers", {
  check <- rcomplex:::.check_fork_results
  ok <- list(1, "a", list(2))
  expect_identical(check(ok, 1:3, "task"), ok)
  failed <- try(stop("boom"), silent = TRUE)
  expect_error(
    check(list(1, failed, 3), c(0.5, 1, 2), "resolution"),
    "resolution 1 failed: boom"
  )
  expect_error(
    check(list(1, NULL, NULL), 1:3, "permutation"),
    "no result for permutation 2, 3 \\(on a forked run a worker may have"
  )
  expect_error(check(list(NULL), 1L, "p"), "rerun with n_cores = 1")
  bare <- structure("Error in f() : raw failure\n", class = "try-error")
  expect_error(check(list(1, bare), 1:2, "task"), "task 2 failed: .*raw")
  classed <- try(rlang::abort("bad arg", class = "my_error"), silent = TRUE)
  expect_error(
    check(list(classed), 3L, "permutation"),
    "permutation 3 failed: bad arg",
    class = "my_error"
  )
  bullets <- try(rlang::abort(c("bad arg", i = "detail")), silent = TRUE)
  msg <- tryCatch(check(list(bullets), 1L, "t"), error = conditionMessage)
  expect_match(msg, "t 1 failed: bad arg")
  expect_identical(lengths(regmatches(msg, gregexpr("detail", msg))), 1L)
})

test_that("a condition with an empty header is wrapped, not duplicated", {
  check <- rcomplex:::.check_fork_results
  cnd <- rlang::error_cnd("custom_error", message = "", body = c(i = "detail"))
  e <- structure("Error\n", class = "try-error", condition = cnd)
  err <- tryCatch(check(list(e), 5L, "task"), error = identity)
  expect_s3_class(err, "custom_error")
  expect_match(conditionMessage(err), "task 5 failed")
})

test_that("a worker killed mid-batch stops with a named error", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(identical(Sys.getenv("R_COVR"), "true"), "covr disables forking")
  res <- suppressWarnings(parallel::mclapply(1:4, function(i) {
    if (i == 3L) tools::pskill(Sys.getpid())
    i
  }, mc.cores = 2L))
  expect_error(
    rcomplex:::.check_fork_results(res, 1:4, "permutation"),
    "no result for permutation .*3 \\(on a forked run"
  )
})

test_that("the serial fallback under Accelerate says so once", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(identical(Sys.getenv("R_COVR"), "true"), "covr disables forking")
  skip_if(
    rcomplex:::.blas_fork_safe(),
    "BLAS is fork-safe here, so the K = 1 test forks"
  )
  rlang::reset_message_verbosity("rcomplex_k1_serial")
  set.seed(1)
  x <- matrix(stats::rnorm(120 * 20), 120, 20)
  x <- x + 1.5 * matrix(stats::rnorm(2 * 20), 2, 20)[rep(1:2, 60), ]
  rownames(x) <- sprintf("g%03d", seq_len(120))
  net <- compute_network(x)
  run <- function() {
    detect_modules(net,
      resolution = c(0.5, 1), seed = 1L, n_cores = 2L,
      test_k1 = TRUE, n_perm_k1 = 10L
    )
  }
  expect_message(run(), "VECLIB_MAXIMUM_THREADS=1")
  expect_no_message(run(), message = "VECLIB")
})
