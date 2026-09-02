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
       n1 = compute_network(x1, density = 0.03, sparse = FALSE),
       n2 = compute_network(x2, density = 0.03, sparse = FALSE))
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

  idx <- match(key_expected, key)
  expect_false(anyNA(idx))

  # The self-excluded urn (D5, rcomplex >= 0.2.0) shifts p-values by
  # O(1/N) relative to canonical ComPlEx: the call sets are identical, or
  # differ only by pairs whose combined BH value sits within 1e-3 of alpha.
  # Measured on this fixture (N = 320): 0 flipped calls, max |diff| of the
  # combined BH value 2.7e-3, mean relative difference 7e-3.
  called <- key[maxbh < 0.05]
  flipped <- c(setdiff(called, key_expected), setdiff(key_expected, called))
  if (length(flipped) > 0L) {
    p_rc <- maxbh[match(flipped, key)]
    p_py <- expected$Max.p.val[match(flipped, key_expected)]
    expect_true(all(abs(p_rc - 0.05) < 1e-3))
    expect_true(all(abs(p_py[!is.na(p_py)] - 0.05) < 1e-3))
  }
  expect_equal(maxbh[idx], expected$Max.p.val, tolerance = 1e-2)
  # absolute guard next to the relative tolerance: measured max |diff| is
  # 2.7e-3, so 5e-3 leaves headroom without letting large shifts through
  expect_lt(max(abs(maxbh[idx] - expected$Max.p.val)), 5e-3)
  expect_equal(cmp$Species1.neigh.overlap[idx], expected$Species1.neigh.overlap)
  expect_equal(cmp$Species2.neigh.overlap[idx], expected$Species2.neigh.overlap)
})

test_that("make_fixture.R reproduces the committed fixture inputs", {
  skip_if_not_installed("withr")
  skip_if_no_fixture()
  script <- normalizePath(fx("make_fixture.R"))
  withr::local_preserve_seed()
  tmp <- withr::local_tempdir()
  withr::with_dir(tmp, source(script, local = new.env()))
  # Expression TSVs: compare parsed numbers, not bytes -- the last ULP of
  # printed doubles differs on x86/gcc (Orion); byte-identical output is
  # only guaranteed on the reference platform (arm64/clang, md5s in the
  # fixture README).
  for (f in setdiff(fixture_inputs, "ortho_pairs.tsv")) {
    regen <- read_expr(file.path(tmp, f))
    ref <- read_expr(fx(f))
    expect_equal(regen, ref, tolerance = 1e-12)
  }
  # ortho_pairs.tsv is strings only: byte-identical on every platform
  expect_identical(readLines(file.path(tmp, "ortho_pairs.tsv")),
                   readLines(fx("ortho_pairs.tsv")))
})

test_that("find_coexpressologs(pval_combine = 'max', pi0_method = 'none') reproduces the canonical calls", {
  d <- load_complex_py()
  nets <- list(sp1 = d$n1, sp2 = d$n2)
  e_max <- find_coexpressologs(nets, d$ortho, alpha = 0.05,
                               pval_combine = "max", pi0_method = "none")
  e_min <- find_coexpressologs(nets, d$ortho, alpha = 0.05,
                               pval_combine = "min", pi0_method = "none")

  key_expected <- paste(d$expected$Species1, d$expected$Species2)
  key_max <- paste(e_max$gene1, e_max$gene2)
  called <- key_max[e_max$type == "conserved"]

  # same boundary rule as the canonical-calls test above
  flipped <- c(setdiff(called, key_expected), setdiff(key_expected, called))
  if (length(flipped) > 0L) {
    p_rc <- e_max$q.value[match(flipped, key_max)]
    p_py <- d$expected$Max.p.val[match(flipped, key_expected)]
    expect_true(all(abs(p_rc - 0.05) < 1e-3))
    expect_true(all(abs(p_py[!is.na(p_py)] - 0.05) < 1e-3))
  }
  idx <- match(key_expected, key_max)
  expect_false(anyNA(idx))
  expect_equal(e_max$q.value[idx], d$expected$Max.p.val, tolerance = 1e-2)

  # "min" (permissive) calls are a superset of "max" calls
  called_min <- paste(e_min$gene1, e_min$gene2)[e_min$type == "conserved"]
  expect_true(all(called %in% called_min))
  expect_gte(length(called_min), length(called))
})


test_that("find_coexpressologs default pval_combine reproduces the canonical calls (D2)", {
  d <- load_complex_py()
  nets <- list(sp1 = d$n1, sp2 = d$n2)
  e_def <- find_coexpressologs(nets, d$ortho, alpha = 0.05,
                               pi0_method = "none")
  e_max <- find_coexpressologs(nets, d$ortho, alpha = 0.05,
                               pval_combine = "max", pi0_method = "none")
  expect_identical(e_def, e_max)

  # default = BH + pmax: the fixture's Max.p.val criterion holds by default
  key_expected <- paste(d$expected$Species1, d$expected$Species2)
  key_def <- paste(e_def$gene1, e_def$gene2)
  called <- key_def[e_def$type == "conserved"]
  flipped <- c(setdiff(called, key_expected), setdiff(key_expected, called))
  if (length(flipped) > 0L) {
    p_rc <- e_def$q.value[match(flipped, key_def)]
    p_py <- d$expected$Max.p.val[match(flipped, key_expected)]
    expect_true(all(abs(p_rc - 0.05) < 1e-3))
    expect_true(all(abs(p_py[!is.na(p_py)] - 0.05) < 1e-3))
  }
  idx <- match(key_expected, key_def)
  expect_false(anyNA(idx))
  expect_equal(e_def$q.value[idx], d$expected$Max.p.val, tolerance = 1e-2)
})
