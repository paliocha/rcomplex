# DiscreteQvalue::DQ as permutation_hog_test() calls it (R/summary.R):
# method = "Liang" on the known p-value support.

test_that("Liang discrete q-values are bounded, monotone, and beat BH", {
  # Upper-tail hypergeometric p-values P(X >= x) for one urn.
  x <- 0:8
  ss <- stats::phyper(x - 1, 10, 30, 8, lower.tail = FALSE)
  # Null counts follow dhyper; a few strong overlaps sit in the tail.
  counts <- round(200 * stats::dhyper(x, 10, 30, 8))
  counts[7:9] <- counts[7:9] + c(6, 8, 10)
  p <- rep(ss, counts)

  q <- DiscreteQvalue::DQ(p, ss = ss, method = "Liang")$q.values

  expect_true(all(q >= 0 & q <= 1))
  o <- order(p)
  expect_true(all(diff(q[o]) >= -1e-12))
  expect_true(all(q <= stats::p.adjust(p, "BH") + 1e-12))
})
