# Equivalence against canonical ComPlEx.
#
# Provenance: fixtures/complex_py/ holds synthetic 8-module expression data for
# two species (seed 20260901, 320 one-to-one orthologs + 32 sp2 paralogs, 64
# shuffled assignments) and the co-expressolog calls produced for it by
# natstreet/ComPlEx_python (commit 1838a5c) run with a BH shim for
# statsmodels.stats.multitest.multipletests (identical to R
# p.adjust(method = "BH")), density 0.03, min-expr filter disabled.
# The full generation recipe is in PLAN_sparse_network.md, section P0.
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

ortho <- read.delim(fx("ortho_pairs.tsv"), stringsAsFactors = FALSE)
expected <- read.delim(fx("expected_calls.tsv"), stringsAsFactors = FALSE)

x1 <- read_expr(fx("sp1_expr.tsv"))
x2 <- read_expr(fx("sp2_expr.tsv"))
x1 <- x1[rownames(x1) %in% ortho$Species1, ]
x2 <- x2[rownames(x2) %in% ortho$Species2, ]

n1 <- compute_network(x1, density = 0.03)
n2 <- compute_network(x2, density = 0.03)

test_that("density thresholds match canonical ComPlEx", {
  # Raw MR = sqrt(rank_ij * rank_ji), so a threshold is the square root of a
  # rank product. ComPlEx_python does not write the threshold; these guard
  # the values under which the calls below were reproduced.
  expect_equal(n1$threshold, sqrt(95172), tolerance = 1e-9)
  expect_equal(n2$threshold, sqrt(115258), tolerance = 1e-9)
  expect_equal(n1$threshold, reference_density_threshold(n1$network, 0.03),
               tolerance = 1e-9)
  expect_equal(n2$threshold, reference_density_threshold(n2$network, 0.03),
               tolerance = 1e-9)
})

test_that("co-expressolog calls match canonical ComPlEx", {
  cmp <- compare_neighborhoods(n1, n2, ortho)
  pool <- cmp$Species1.neigh.overlap > 0 & cmp$Species2.neigh.overlap > 0
  cmp <- cmp[pool, ]
  bh1 <- p.adjust(cmp$Species1.p.val.con, method = "BH")
  bh2 <- p.adjust(cmp$Species2.p.val.con, method = "BH")
  maxbh <- pmax(bh1, bh2)

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
