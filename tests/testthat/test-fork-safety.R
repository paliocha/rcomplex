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
    "no result \\(a worker may have crashed\\) for permutation 2, 3"
  )
  expect_error(check(list(NULL), 1L, "p"), "rerun with n_cores = 1")
})
