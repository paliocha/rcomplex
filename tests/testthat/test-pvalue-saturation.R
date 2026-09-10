# Fixture reproducing the counts the eight-species Pooideae run reported:
# 511 module-directions, 172 distinct q-values, 35 tied at the floor of
# 0.00071, 20 at exactly 1.
pooideae_q <- function() {
  c(
    rep(0.00071, 35), rep(1, 20),
    rep(seq(0.002, 0.99, length.out = 170), length.out = 456)
  )
}

# The permutation p-values behind them: everything on the 1/2001 grid.
pooideae_p <- function() {
  c(
    rep(1 / 2001, 35), rep(1, 20),
    rep(round(seq(4, 1981, length.out = 170)) / 2001, length.out = 456)
  )
}


test_that("counts match the Pooideae run", {
  res <- pvalue_resolution(pooideae_q())

  expect_s3_class(res, "pvalue_resolution")
  expect_equal(res$n, 511L)
  expect_equal(res$n_missing, 0L)
  expect_equal(res$n_distinct, 172L)
  expect_equal(res$min, 0.00071)
  expect_equal(res$n_at_min, 35L)
  expect_equal(res$n_at_one, 20L)
})


test_that("without n_perm the floor fields stay unset", {
  res <- pvalue_resolution(pooideae_q())

  expect_null(res$n_perm)
  expect_true(is.na(res$floor))
  expect_true(is.na(res$floor_status))
  expect_true(is.na(res$permutation_limited))
  expect_true(is.na(res$n_off_grid))
})


test_that("the floor is detected when the minimum sits on it", {
  res <- pvalue_resolution(pooideae_p(), n_perm = 2000)

  expect_equal(res$floor, 1 / 2001)
  expect_true(res$permutation_limited)
  expect_identical(res$floor_status, "at")
  expect_equal(res$n_off_grid, 0L)
  expect_equal(res$n_distinct, 172L)
  expect_equal(res$n_at_min, 35L)
})


test_that("a minimum above the floor is called evidence-limited", {
  res <- pvalue_resolution(c(0.02, 0.11, 0.4, 0.9), n_perm = 999)

  expect_false(res$permutation_limited)
  expect_identical(res$floor_status, "above")
  expect_equal(res$n_at_min, 1L)
  # More permutations cannot separate an untied minimum, so no suggestion.
  expect_true(is.na(res$suggested_n_perm))
})


test_that("an evidence-limited tie gets no permutation suggestion", {
  # Tied, but well above the 1/1001 floor: the tie is in the data.
  res <- pvalue_resolution(c(rep(0.2, 5), 0.4, 0.8), n_perm = 1000)

  expect_false(res$permutation_limited)
  expect_equal(res$n_at_min, 5L)
  expect_true(is.na(res$suggested_n_perm))
})


test_that("dropping n_perm does not turn that tie into advice", {
  # Same vector, no n_perm: a minimum of 0.2 would mean a run of four
  # permutations, so it cannot stand in for a floor and the answer must
  # stay the NA the n_perm branch gives.
  v <- c(rep(0.2, 5), 0.4, 0.8)

  expect_true(is.na(pvalue_resolution(v)$suggested_n_perm))
  expect_true(is.na(pvalue_resolution(rep(1, 5))$suggested_n_perm))

  txt <- paste(capture.output(print(pvalue_resolution(v))), collapse = " ")
  expect_false(grepl("n_perm >=", txt, fixed = TRUE))
  expect_match(txt, "too large to stand in for a floor")

  # A minimum small enough to be a floor still gets one.
  expect_equal(
    pvalue_resolution(c(rep(0.001, 5), 0.4))$suggested_n_perm, 5000
  )
})


test_that("suggested_n_perm gives the tied block room to separate", {
  res <- pvalue_resolution(pooideae_p(), n_perm = 2000)

  # 35 values on one grid point need >= 35 grid points below it.
  expect_lte(1 / (res$suggested_n_perm + 1), res$floor / res$n_at_min)
  # ... and it is the smallest value on the 1-2-5 ladder that does: the
  # rung below leaves the block no room.
  expect_gt(res$n_at_min / (5e4 + 1), res$floor)
  expect_equal(res$suggested_n_perm, 1e5)

  # Without n_perm the observed minimum plays the floor's role.
  expect_equal(pvalue_resolution(pooideae_q())$suggested_n_perm, 5e4)
})


test_that("the suggestion counts grid points, not grid spacings", {
  # 7 values tied at the 1/143 floor need 7 grid points below it, so
  # n_perm + 1 >= 7 * 143 and n_perm >= 1000 -- exactly a rung of the
  # 1-2-5 ladder, so an off-by-one in the bound shows up as 2000.
  res <- pvalue_resolution(c(rep(1 / 143, 7), 0.5, 0.9), n_perm = 142)

  expect_equal(res$n_at_min, 7L)
  expect_equal(res$suggested_n_perm, 1000)
})


test_that("the floor test survives last-bit differences", {
  # A p-value formed as (1 + count) / (n_perm + 1) elsewhere need not be
  # bit-identical to 1/(n_perm + 1) recomputed here.
  wobbled <- (1 / 2001) * (1 + c(0, 1, -1, 2) * .Machine$double.eps)
  res <- pvalue_resolution(c(wobbled, 0.5, 0.9), n_perm = 2000)

  expect_equal(res$n_at_min, 4L)
  expect_true(res$permutation_limited)
  # The distinct count must read equality the same way, or the printout
  # says "6 distinct (100%)" and "4 tied at the minimum" at once.
  expect_equal(res$n_distinct, 3L)
  expect_lte(res$n_distinct, res$n - res$n_at_min + 1L)
})


test_that("values at 1 are counted with the same tolerance", {
  res <- pvalue_resolution(c(0.5, 1, 1 - .Machine$double.eps))

  expect_equal(res$n_at_one, 2L)
  expect_equal(res$n_distinct, 2L)
})


test_that("the tie tolerance comes from the shared helper", {
  # The constant used to be written out here as well as in .tie_tol(),
  # so the two had to be kept in step by hand. Pin the width of the
  # n_at_one window against the helper: a change to .tie_tol() must move
  # this function's counts along with preservation_matrix_test()'s.
  tol <- rcomplex:::.tie_tol()

  expect_identical(pvalue_resolution(c(0.5, 1, 1 - tol / 2))$n_at_one, 2L)
  expect_identical(pvalue_resolution(c(0.5, 1, 1 - tol * 2))$n_at_one, 1L)
})


test_that("a value above 1 is rejected, never counted at 1", {
  # The @return text for n_at_one used to promise that values above 1
  # are counted here rather than dropped. They cannot be: validation
  # refuses the vector first. Only preservation_matrix_test(), which
  # does not validate, can report one.
  expect_error(pvalue_resolution(c(0.5, 1.0000001)), "\\[0, 1\\]")
})


test_that("the floor comparison is relative, not absolute", {
  # 1.001e-6 is a thousandth above the 1/1000001 floor -- a mile away in
  # relative terms, but only 1e-9 in absolute ones, well inside
  # sqrt(.Machine$double.eps). Only a relative test calls this above.
  res <- pvalue_resolution(c(1.001e-6, 0.3, 0.9), n_perm = 1e6)

  expect_identical(res$floor_status, "above")
  expect_false(res$permutation_limited)
})


test_that("the off-grid test scales with the grid index", {
  # grid = p * (n_perm + 1), so a relative wobble in p becomes an absolute
  # error proportional to the index. k = 1000 has to be judged on the same
  # relative scale as k = 1, or arithmetic on large p reads as off grid.
  res <- pvalue_resolution(
    c(1 / 2001, (1000 + 1e-7) / 2001, 1),
    n_perm = 2000
  )

  expect_equal(res$n_off_grid, 0L)
})


test_that("off-grid values are counted", {
  # q-values are not multiples of the floor, which is how a reader learns
  # the vector is not raw permutation output.
  res <- pvalue_resolution(c(1 / 2001, 0.0333, 0.4171, 0.9), n_perm = 2000)

  expect_equal(res$n_off_grid, 3L)
  expect_true(res$permutation_limited)
})


test_that("values below the floor warn", {
  expect_warning(
    pvalue_resolution(c(1e-6, 0.3, 0.9), n_perm = 100),
    "below the permutation floor"
  )
})


test_that("a minimum below the floor is its own state, not 'above'", {
  res <- suppressWarnings(pvalue_resolution(c(0, 0, 0.5), n_perm = 1000))

  expect_identical(res$floor_status, "below")
  # Neither "the minimum IS the floor" nor "above it" is true here, so the
  # two-way flag must abstain rather than pick the wrong one.
  expect_true(is.na(res$permutation_limited))
  expect_true(is.na(res$suggested_n_perm))
})


test_that("print does not claim a below-floor minimum is above it", {
  # The object is the deliverable: the construction-time warning is gone
  # by the time a saved result is reprinted, so print has to carry it.
  res <- suppressWarnings(pvalue_resolution(c(0, 0, 0.5), n_perm = 1000))
  txt <- paste(capture.output(print(res)), collapse = " ")

  expect_match(txt, "below the floor")
  expect_false(grepl("minimum is above the floor", txt, fixed = TRUE))
  expect_false(grepl("More permutations will not help", txt, fixed = TRUE))
  expect_false(grepl("permutation-limited", txt, fixed = TRUE))
})


test_that("NAs are counted and dropped", {
  res <- pvalue_resolution(c(0.001, NA, 0.001, 0.5, NA))

  expect_equal(res$n, 3L)
  expect_equal(res$n_missing, 2L)
  expect_equal(res$n_distinct, 2L)
  expect_equal(res$n_at_min, 2L)
})


test_that("invalid input is rejected", {
  expect_error(pvalue_resolution("a"), "non-empty numeric")
  expect_error(pvalue_resolution(numeric(0)), "non-empty numeric")
  expect_error(pvalue_resolution(c(NA_real_, NA_real_)), "no non-missing")
  expect_error(pvalue_resolution(c(0.5, 1.5)), "\\[0, 1\\]")
  expect_error(pvalue_resolution(c(0.5, -0.1)), "\\[0, 1\\]")
  expect_error(pvalue_resolution(0.5, n_perm = 0), "positive whole number")
  expect_error(pvalue_resolution(0.5, n_perm = 10.5), "positive whole number")
  expect_error(
    pvalue_resolution(0.5, n_perm = c(10, 20)),
    "positive whole number"
  )
  expect_error(
    pvalue_resolution(0.5, n_perm = "1000"),
    "positive whole number"
  )
})


test_that("print names the tie and points at the effect size", {
  out <- capture.output(print(pvalue_resolution(pooideae_q())))
  txt <- paste(out, collapse = " ")

  expect_match(txt, "511")
  expect_match(txt, "172")
  expect_match(txt, "35 values tied there")
  expect_match(txt, "NO ordering")
  expect_match(txt, "Zsummary_std")
  expect_match(txt, "n_perm >= 50000")
  # No floor lines without n_perm.
  expect_false(any(grepl("floor 1/", out, fixed = TRUE)))
})


test_that("print reports the floor and its status", {
  out <- capture.output(
    print(pvalue_resolution(pooideae_p(), n_perm = 2000))
  )
  txt <- paste(out, collapse = " ")

  expect_match(txt, "floor 1/\\(n_perm \\+ 1\\)")
  expect_match(txt, "permutation-limited")
  expect_match(txt, "the minimum IS the floor")
})


test_that("print says so when nothing is tied", {
  out <- capture.output(
    print(pvalue_resolution(c(0.01, 0.2, 0.6), n_perm = 999))
  )
  txt <- paste(out, collapse = " ")

  expect_match(txt, "\\(unique\\)")
  expect_match(txt, "still order the tests")
  expect_match(txt, "evidence-limited")
  expect_false(grepl("NO ordering", txt))
})


test_that("print returns its input invisibly", {
  res <- pvalue_resolution(c(0.1, 0.2))
  expect_output(out <- withVisible(print(res)))

  expect_false(out$visible)
  expect_identical(out$value, res)
})


test_that("round-up lands on 1, 2 or 5 times a power of ten", {
  expect_equal(rcomplex:::.round_up_nice(70034), 1e5)
  expect_equal(rcomplex:::.round_up_nice(1), 1)
  expect_equal(rcomplex:::.round_up_nice(1.5), 2)
  expect_equal(rcomplex:::.round_up_nice(3), 5)
  expect_equal(rcomplex:::.round_up_nice(6), 10)
  expect_true(is.na(rcomplex:::.round_up_nice(0)))
  expect_true(is.na(rcomplex:::.round_up_nice(NaN)))
})
