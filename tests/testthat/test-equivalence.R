# Equivalence against canonical ComPlEx.
#
# Fixtures live in fixtures/complex_py/; see the README.md there for
# provenance (make_fixture.R, seed 20260901; expected calls from
# natstreet/ComPlEx_python commit 1838a5c run with a BH shim identical to
# R p.adjust(method = "BH"), density 0.03) and the regeneration recipe.
#
# The fixture is only valid under the canonical gene universe, i.e. genes
# restricted to those with an ortholog (ComPlEx_python builds the network on
# the ortholog-restricted expression tables).

fx <- function(...) testthat::test_path("fixtures", "complex_py", ...)

read_expr <- function(file) {
  d <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  x <- as.matrix(d[, -1])
  rownames(x) <- d$Genes
  x
}

fixture_inputs <- c("sp1_expr.tsv", "sp2_expr.tsv", "ortho_pairs.tsv")

skip_if_no_fixture <- function() {
  skip_if_not(all(file.exists(fx(c("expected_calls.tsv", "ortho_pairs.tsv",
                                   "sp1_expr.tsv", "sp2_expr.tsv",
                                   "make_fixture.R")))))
}

load_complex_py <- function() {
  skip_if_no_fixture()
  ortho <- read.delim(fx("ortho_pairs.tsv"), stringsAsFactors = FALSE)
  expected <- read.delim(fx("expected_calls.tsv"), stringsAsFactors = FALSE)
  x1 <- read_expr(fx("sp1_expr.tsv"))
  x2 <- read_expr(fx("sp2_expr.tsv"))
  x1 <- x1[rownames(x1) %in% ortho$Species1, ]
  x2 <- x2[rownames(x2) %in% ortho$Species2, ]
  list(ortho = ortho, expected = expected,
       n1 = compute_network(x1, density = 0.03),
       n2 = compute_network(x2, density = 0.03))
}

test_that("density thresholds are stable and match pure-R reference", {
  d <- load_complex_py()
  expect_equal(d$n1$threshold, reference_density_threshold(d$n1$network, 0.03),
               tolerance = 1e-9)
  expect_equal(d$n2$threshold, reference_density_threshold(d$n2$network, 0.03),
               tolerance = 1e-9)
  # Stability pin (rcomplex's own values, not written by ComPlEx_python).
  # Raw MR = sqrt(rank_ij * rank_ji), so a threshold is the square root of
  # a rank product.
  expect_equal(d$n1$threshold, sqrt(95172), tolerance = 1e-9)
  expect_equal(d$n2$threshold, sqrt(115258), tolerance = 1e-9)
})

test_that("co-expressolog calls match canonical ComPlEx", {
  d <- load_complex_py()
  cmp <- compare_neighborhoods(d$n1, d$n2, d$ortho)
  pool <- cmp$Species1.neigh.overlap > 0 & cmp$Species2.neigh.overlap > 0
  cmp <- cmp[pool, ]
  bh1 <- p.adjust(cmp$Species1.p.val.con, method = "BH")
  bh2 <- p.adjust(cmp$Species2.p.val.con, method = "BH")
  maxbh <- pmax(bh1, bh2)

  expected <- d$expected
  key <- paste(cmp$Species1, cmp$Species2)
  key_expected <- paste(expected$Species1, expected$Species2)
  expect_equal(nrow(expected), 149L)
  expect_true(setequal(key[maxbh < 0.05], key_expected))

  idx <- match(key_expected, key)
  expect_false(anyNA(idx))
  expect_lt(max(abs(maxbh[idx] - expected$Max.p.val)), 1e-12)
  expect_equal(cmp$Species1.neigh.overlap[idx], expected$Species1.neigh.overlap)
  expect_equal(cmp$Species2.neigh.overlap[idx], expected$Species2.neigh.overlap)
})

test_that("make_fixture.R reproduces the committed fixture inputs", {
  skip_if_not_installed("withr")
  skip_if_no_fixture()
  script <- normalizePath(fx("make_fixture.R"))
  committed <- normalizePath(fx(fixture_inputs))
  withr::local_preserve_seed()
  tmp <- withr::local_tempdir()
  withr::with_dir(tmp, source(script, local = new.env()))
  for (i in seq_along(fixture_inputs)) {
    expect_identical(readLines(file.path(tmp, fixture_inputs[i])),
                     readLines(committed[i]))
  }
})
