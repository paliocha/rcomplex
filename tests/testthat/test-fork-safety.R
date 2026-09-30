# Fork safety of the K = 1 test: .blas_fork_safe() gates forking on the
# BLAS, .collect_perm_batch() turns failed or crashed workers into an error.

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

test_that(".collect_perm_batch() passes scalars and stops on failed workers", {
  collect <- rcomplex:::.collect_perm_batch
  expect_identical(collect(list(1, 2.5, 3), "t"), c(1, 2.5, 3))
  failed <- try(stop("boom"), silent = TRUE)
  expect_error(
    collect(list(1, failed, 3), "K = 1 test"),
    "K = 1 test: 1 of 3 forked workers failed -- .*boom.*n_cores = 1"
  )
  expect_error(collect(list(1, NULL), "t"), "no result \\(it may have")
  expect_error(collect(list(1, c(2, 3)), "t"), "1 of 2 forked workers failed")
  expect_error(collect(list(), "t"), "0 of 0")
})

test_that("the K = 1 test gives the serial result on any fork path", {
  skip_on_os("windows")
  skip_on_cran()
  set.seed(1)
  x <- matrix(stats::rnorm(300 * 20), 300, 20)
  f <- matrix(stats::rnorm(4 * 20), 4, 20)
  x <- x + 1.5 * f[rep(1:4, length.out = 300), ]
  rownames(x) <- sprintf("g%03d", seq_len(300))
  net <- compute_network(x)
  run <- function(nc) {
    detect_modules(net,
      resolution = c(0.5, 1, 2), seed = 1L, n_cores = nc,
      test_k1 = TRUE, n_perm_k1 = 20L
    )$k1_test
  }
  serial <- run(1L)
  forked <- run(2L)
  expect_identical(forked$lambda_null, serial$lambda_null)
  expect_identical(forked$p_value, serial$p_value)
})
