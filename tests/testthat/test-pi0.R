# pi0 estimation from randomized p-values (D4)
#
# Exact hypergeometric p-values are discrete and pile up at 1 (small
# expected overlaps, the x > 1 gate), so Storey's estimator returns
# pi0 = 1 even with 30% planted alternatives. The randomized p-value
# P(X > x) + U * P(X = x) is exactly U(0, 1) under H0 (Dickhaus et al.
# 2012), so pi0 estimated on it is unbiased; it is then applied to the
# exact p-values.

# Appendix B recipe: 2000 pairs, N = 30000, lognormal neighbourhood sizes,
# 30% planted alternatives (true pi0 ~ 0.69)
sim_pairs <- function(seed = 1, n = 2000, N = 30000, alt_rate = 0.3) {
  set.seed(seed)
  m <- pmin(N - 1, pmax(3, round(exp(rnorm(n, log(600), 0.9)))))
  k <- pmin(N - 1, pmax(2, round(exp(rnorm(n, log(500), 0.9)))))
  alt <- runif(n) < alt_rate
  x <- rhyper(n, m, N - m, k)
  x[alt] <- pmin(pmin(m, k)[alt],
                 x[alt] + rpois(sum(alt), 3 + 0.5 * m[alt] * k[alt] / N))
  list(
    alt = alt,
    p = phyper(x - 1, m, N - m, k, lower.tail = FALSE),
    p_gt = phyper(x, m, N - m, k, lower.tail = FALSE),
    p_eq = dhyper(x, m, N - m, k)
  )
}

rand_fn <- function(s) function() s$p_gt + runif(length(s$p_gt)) * s$p_eq


test_that("randomized pi0 recovers the true null proportion; Storey gives 1", {
  s <- sim_pairs(1)
  true_pi0 <- mean(!s$alt)
  expect_lt(abs(true_pi0 - 0.69), 0.05)

  set.seed(2)
  rnd <- compute_qvalues(s$p, rand_fn(s), pi0_method = "randomized", B = 20L)
  expect_named(rnd, c("qvalues", "pi0"))
  expect_lt(abs(rnd$pi0 - true_pi0), 0.05)
  expect_length(rnd$qvalues, length(s$p))
  # q-values are the exact p-values' BH values scaled by pi0
  expect_equal(rnd$qvalues, rnd$pi0 * p.adjust(s$p, "BH"))

  st <- compute_qvalues(s$p, pi0_method = "storey")
  expect_equal(st$pi0, 1)
  expect_equal(st$qvalues, qvalue::qvalue(s$p)$qvalues)
})


test_that("pi0_method = 'none' is Benjamini-Hochberg", {
  s <- sim_pairs(1)
  none <- compute_qvalues(s$p, pi0_method = "none")
  expect_equal(none$pi0, 1)
  expect_equal(none$qvalues, p.adjust(s$p, "BH"))
  # p_rand_fn is not needed for the non-randomized methods
  expect_silent(compute_qvalues(s$p, NULL, pi0_method = "none"))
})


test_that("null-only simulation yields no q-value below 0.05", {
  s <- sim_pairs(1, alt_rate = 0)
  expect_false(any(s$alt))
  set.seed(3)
  rnd <- compute_qvalues(s$p, rand_fn(s), pi0_method = "randomized", B = 20L)
  expect_gt(rnd$pi0, 0.9)
  expect_equal(sum(rnd$qvalues < 0.05), 0L)
})


test_that("randomized pi0 is reproducible under set.seed and needs p_rand_fn", {
  s <- sim_pairs(1)
  set.seed(11)
  a <- compute_qvalues(s$p, rand_fn(s), pi0_method = "randomized", B = 5L)
  set.seed(11)
  b <- compute_qvalues(s$p, rand_fn(s), pi0_method = "randomized", B = 5L)
  expect_identical(a, b)
  set.seed(12)
  c <- compute_qvalues(s$p, rand_fn(s), pi0_method = "randomized", B = 5L)
  expect_false(identical(a$pi0, c$pi0))

  expect_error(compute_qvalues(s$p, pi0_method = "randomized"),
               "p_rand_fn")
  expect_error(compute_qvalues(s$p, rand_fn(s), pi0_method = "bogus"))
})


test_that("degenerate inputs fall back gracefully", {
  # fewer than 2 p-values: returned unchanged, pi0 undefined
  one <- compute_qvalues(0.2, function() 0.5, pi0_method = "randomized")
  expect_equal(one$qvalues, 0.2)
  expect_true(is.na(one$pi0))

  # all randomized p-values tiny: pi0est errors (max p below the lambda
  # range) -> fallback pi0 = 1 (BH) instead of an error
  p <- c(0.001, 0.002, 0.003, 0.5, 0.6)
  set.seed(4)
  r <- compute_qvalues(p, function() rep(1e-4, 5), pi0_method = "randomized",
                       B = 3L)
  expect_equal(r$pi0, 1)
  expect_equal(r$qvalues, p.adjust(p, "BH"))
})


test_that("filter_zero does not deflate the randomized pi0 (null-only, small E[X])", {
  # Null-only pairs with small expected overlap (median E[X] ~ 0.08 at
  # m ~ 300, k ~ 5, N = 20000, lognormal spread): most rows have overlap
  # 0 and are removed by filter_zero. The randomized p-value is U(0, 1)
  # under H0 only unconditionally -- conditioning on X > 0 truncates it
  # to [0, P(X > 0)) and collapses pi0, so an all-null data set comes out
  # significant. pi0 must therefore be estimated on the unfiltered
  # comparison rows. (Other data seeds make every pi0est() draw error and
  # silently fall back to pi0 = 1 -- the failure is data-mix dependent;
  # this seed gives a stable deflated estimate under the filtered draw.)
  set.seed(1)
  n <- 2000L
  N <- 20000L
  m <- pmin(N - 1, pmax(3, round(exp(rnorm(n, log(300), 1)))))
  k <- pmin(N - 1, pmax(1, round(exp(rnorm(n, log(5), 0.8)))))
  x <- rhyper(n, m, N - m, k)
  p_con <- ifelse(x <= 1, 1, phyper(x - 1, m, N - m, k, lower.tail = FALSE))
  p_div <- phyper(x, m, N - m, k)
  p_gt <- phyper(x, m, N - m, k, lower.tail = FALSE)
  p_eq <- dhyper(x, m, N - m, k)
  cmp <- data.frame(
    Species1 = paste0("A", seq_len(n)), Species2 = paste0("B", seq_len(n)),
    hog = paste0("HOG", seq_len(n)),
    Species1.neigh.overlap = x, Species1.p.val.con = p_con,
    Species1.p.val.div = p_div, Species1.p.val.gt = p_gt,
    Species1.p.val.eq = p_eq,
    Species2.neigh.overlap = x, Species2.p.val.con = p_con,
    Species2.p.val.div = p_div, Species2.p.val.gt = p_gt,
    Species2.p.val.eq = p_eq,
    stringsAsFactors = FALSE
  )
  # the filter bites (most rows dropped) but enough survive to test on
  kept <- sum(cmp$Species1.neigh.overlap > 0)
  expect_gt(kept, 100L)
  expect_lt(kept, n / 2L)

  set.seed(32)
  s <- summarize_comparison(cmp)
  expect_equal(nrow(s$results), kept)
  expect_gt(min(s$summary$pi0), 0.9)
  expect_equal(sum(s$results$Species1.q.val.con < 0.05), 0L)
  expect_equal(sum(s$results$Species2.q.val.con < 0.05), 0L)
})
