# Tests for module_auroc_reciprocal(). Fixture: two random 200-gene
# species, 1:1 orthologs, the same 20 blocks of 10 as modules in both;
# species 2 splits the last block in two, so one half has no partner.

recip_fixture <- function() {
  withr::with_seed(1L, {
    n <- 200L
    x1 <- matrix(stats::rnorm(n * 20L), n)
    x2 <- matrix(stats::rnorm(n * 20L), n)
  })
  rownames(x1) <- paste0("A", seq_len(n))
  rownames(x2) <- paste0("B", seq_len(n))
  lab <- paste0("m", rep(1:20, each = 10L))
  lab2 <- lab
  lab2[196:200] <- "m21"
  list(
    net1 = compute_network(x1, density = 0.05),
    net2 = compute_network(x2, density = 0.05),
    ortho = data.frame(
      Species1 = rownames(x1), Species2 = rownames(x2),
      hog = paste0("H", seq_len(n)), stringsAsFactors = FALSE
    ),
    mods1 = as_modules(stats::setNames(lab, rownames(x1))),
    mods2 = as_modules(stats::setNames(lab2, rownames(x2)))
  )
}

rfx <- recip_fixture()

run_recip <- function(..., n_cores = 1L) {
  module_auroc_reciprocal(rfx$mods1, rfx$mods2, rfx$net1, rfx$net2,
    rfx$ortho, ...,
    n_null = 50L, max_draws = 50L, n_cores = n_cores, seed = 1L
  )
}


test_that("blocks pair with themselves and p.val is pmax", {
  r <- run_recip()
  expect_named(r, c(
    "module1", "module2", "jaccard", "p.val.1to2", "p.val.2to1", "p.val",
    "q.val", "z.1to2", "z.2to1"
  ))
  expect_identical(r$module1, r$module2)
  expect_identical(nrow(r), 20L)
  expect_equal(r$jaccard, ifelse(r$module1 == "m20", 0.5, 1))
  expect_identical(r$p.val, pmax(r$p.val.1to2, r$p.val.2to1))
  expect_true(all(r$q.val >= 0 & r$q.val <= 1))
  expect_s3_class(attr(r, "resolution"), "pvalue_resolution")
})


test_that("pval_combine = 'min' gives pmin with the same directions", {
  mx <- run_recip()
  mn <- run_recip(pval_combine = "min")
  expect_identical(mn$p.val, pmin(mn$p.val.1to2, mn$p.val.2to1))
  expect_identical(mn$p.val.1to2, mx$p.val.1to2)
})


test_that("a module without a reciprocal partner is reported", {
  u <- attr(run_recip(), "unmatched")
  expect_identical(u$species1, character(0))
  expect_identical(u$species2, "m21")
})


test_that("the result does not depend on n_cores", {
  expect_identical(run_recip(n_cores = 1L), run_recip(n_cores = 2L))
})
