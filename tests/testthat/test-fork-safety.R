# Fork safety: .check_fork_results() turns failed or crashed workers into
# an error at every fork site.

test_that(".check_fork_results() names failed tasks and crashed workers", {
  check <- rcomplex:::.check_fork_results
  ok <- list(1, "a", list(2))
  expect_identical(check(ok, 1:3, "task"), ok)
  err_of <- function(x) {
    e <- structure("Error\n", class = "try-error", condition = x)
    tryCatch(check(list(1, e), c(0.5, 2), "resolution"), error = identity)
  }
  count <- function(pat, msg) lengths(regmatches(msg, gregexpr(pat, msg)))
  # base error: task named, worker message kept once
  err <- err_of(simpleError("boom"))
  expect_match(conditionMessage(err), "resolution 2 failed")
  expect_identical(count("boom", conditionMessage(err)), 1L)
  # several lines survive
  expect_match(conditionMessage(err_of(simpleError(c("a", "zz")))), "zz")
  # rlang error with bullets: class kept, bullet rendered once
  cnd <- tryCatch(
    rlang::abort(c("bad arg", i = "detail"), class = "my_error"),
    error = identity
  )
  err <- err_of(cnd)
  expect_s3_class(err, "my_error")
  expect_identical(count("detail", conditionMessage(err)), 1L)
  # header from a cnd_header() method (empty message)
  err <- err_of(rlang::error_cnd("custom_error", message = "", body = "x"))
  expect_s3_class(err, "custom_error")
  # empty base message: just the task line
  expect_match(conditionMessage(err_of(simpleError(""))), "resolution 2 failed")
  # a try-error without a condition
  bare <- structure("Error in f() : raw failure\n", class = "try-error")
  expect_error(check(list(1, bare), 1:2, "task"), "task 2 failed: .*raw")
  # crashed workers
  expect_error(
    check(list(1, NULL, NULL), 1:3, "permutation"),
    "no result for permutation 2, 3 \\(on a forked run a worker may have"
  )
  expect_error(check(list(NULL), 1L, "p"), "rerun with n_cores = 1")
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
