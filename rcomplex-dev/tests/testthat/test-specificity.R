# Tests for rcomplex:::compare_specificity() against the pure-R oracle
# reference_specificity() in helper-reference.R. Fixtures: make_graded_nets(),
# make_cmp_nets(), sparse_net() from helper-reference.R.

spec_cols <- function(prefix) {
  paste0(prefix, c(".neigh", ".mapped", ".auroc", ".p_value", ".jaccard"))
}
both_cols <- c(spec_cols("species1"), spec_cols("species2"))

# nolint start: object_usage_linter. (oracle from helper-reference.R)
expect_grid_equal <- function(res, ...) {
  ref <- reference_specificity(..., grid_frac = rcomplex:::.rank_grid_frac)
  for (s in c("species1", "species2")) {
    testthat::expect_equal(unname(res[[paste0(s, ".auroc.grid")]]),
      unname(ref[[paste0(s, ".auroc.grid")]]),
      tolerance = 1e-12
    )
  }
}
# nolint end

expect_spec_equal <- function(res, ref, cols = both_cols) {
  testthat::expect_equal(as.list(res[cols]), as.list(ref[cols]),
    tolerance = 1e-12
  )
}


test_that("dense compare_specificity matches the oracle in both directions", {
  f <- make_graded_nets()
  res <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  ref <- reference_specificity(
    f$net1$network, f$net2$network, 5, 5, f$ortho
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c(
    "gene1", "gene2", "hog", spec_cols("species1"), "species1.n.cand",
    spec_cols("species2"), "species2.n.cand",
    "species1.effect_size", "species1.auroc.grid",
    "species2.effect_size", "species2.auroc.grid"
  ))
  expect_equal(res$gene1, f$ortho$gene1)
  expect_equal(res$hog, f$ortho$hog)
  expect_spec_equal(res, ref)
  expect_identical(res$species1.effect_size, res$species1.auroc)
  expect_identical(res$species2.effect_size, res$species2.auroc)
})

test_that("a store above the lowest tier matches the oracle with store", {
  f <- make_graded_nets()
  res <- rcomplex:::compare_specificity(
    sparse_net(f$net1, 5), sparse_net(f$net2, 5), f$ortho
  )
  ref <- reference_specificity(
    f$net1$network, f$net2$network, 5, 5, f$ortho,
    store1 = 5, store2 = 5
  )
  expect_spec_equal(res, ref)
  expect_grid_equal(
    res, f$net1$network, f$net2$network, 5, 5, f$ortho,
    store1 = 5, store2 = 5
  )
  # tier 4 really is unstored: collapsing it to the bottom block changes
  # other candidates' AUROCs (e.g. B16 for anchor A04), hence the p-values
  dense <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  expect_false(isTRUE(all.equal(res$species1.p_value, dense$species1.p_value)))
})

test_that("a store holding every nonzero entry equals the dense result", {
  f <- make_graded_nets()
  dense <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  sp <- rcomplex:::compare_specificity(
    sparse_net(f$net1, 4), sparse_net(f$net2, 4), f$ortho
  )
  expect_identical(sp, dense)
})

test_that("own-HOG orthologs leave the mapped set (paralog row)", {
  f <- make_cmp_nets()
  res <- rcomplex:::compare_specificity(
    sparse_net(f$net1), sparse_net(f$net2), f$ortho
  )
  ref <- reference_specificity(
    f$net1$network, f$net2$network, f$net1$threshold, f$net2$threshold,
    f$ortho,
    store1 = f$net1$threshold, store2 = f$net2$threshold
  )
  expect_equal(nrow(res), 31L)
  expect_spec_equal(res, ref)
  # several rows per anchor: the grid is written per row from one per-anchor
  # selection, so a stale or misindexed grid would show here
  expect_grid_equal(
    res, f$net1$network, f$net2$network, f$net1$threshold,
    f$net2$threshold, f$ortho,
    store1 = f$net1$threshold, store2 = f$net2$threshold
  )
  # A_001 maps to B_001 and B_031: both must be absent from its mapped set
  rows <- which(res$gene1 == "A_001")
  expect_length(rows, 2L)
  thr <- f$net1$threshold
  n1 <- setdiff(names(which(f$net1$network[, "A_001"] >= thr)), "A_001")
  mapped <- setdiff(
    unique(f$ortho$gene2[f$ortho$gene1 %in% n1]), c("B_001", "B_031")
  )
  expect_equal(res$species1.mapped[rows], rep(length(mapped), 2L))
})

test_that("p lies on the 1/n grid, AUROC in [0, 1], isolated genes NA", {
  f <- make_graded_nets()
  res <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  for (s in c("species1", "species2")) {
    p <- res[[paste0(s, ".p_value")]]
    a <- res[[paste0(s, ".auroc")]]
    expect_true(any(is.na(p)) && any(!is.na(p)))
    expect_equal(is.na(p), is.na(a))
    expect_true(all(abs(p * 30 - round(p * 30)) < 1e-9, na.rm = TRUE))
    expect_true(all(p > 0 & p <= 1, na.rm = TRUE))
    expect_true(all(a >= 0 & a <= 1, na.rm = TRUE))
  }
  iso <- res$gene1 %in% paste0("A", 28:30)
  expect_true(all(is.na(res$species1.p_value[iso])))
  expect_true(all(is.na(res$species2.p_value[iso])))
  expect_equal(res$species1.mapped[iso], rep(0L, 3L))
  expect_false(any(is.na(res$species1.p_value[res$gene1 %in% "A01"])))
})

test_that("directions restricts the columns and keeps the values", {
  f <- make_graded_nets()
  both <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  d12 <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho, directions = "1to2")
  d21 <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho, directions = "2to1")
  keys <- c("gene1", "gene2", "hog")
  expect_named(d12, c(
    keys, spec_cols("species1"), "species1.n.cand",
    "species1.effect_size", "species1.auroc.grid"
  ))
  expect_named(d21, c(
    keys, spec_cols("species2"), "species2.n.cand",
    "species2.effect_size", "species2.auroc.grid"
  ))
  expect_identical(d12, both[names(d12)])
  expect_identical(d21, both[names(d21)])
})

test_that("n_cores = 2 equals the serial result", {
  f <- make_graded_nets()
  serial <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  par <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho, n_cores = 2L)
  expect_identical(par, serial)
})

test_that("a non-symmetric sparse pattern is an error", {
  g <- c("g1", "g2", "g3")
  m <- Matrix::sparseMatrix(
    i = 1L, j = 2L, x = 1, dims = c(3L, 3L), dimnames = list(g, g)
  )
  net <- list(network = m, threshold = 1, store_threshold = 1)
  ortho <- data.frame(gene1 = g, gene2 = g, hog = 1:3)
  expect_error(rcomplex:::compare_specificity(net, net, ortho), "symmetric")
})

test_that("a membership-only network matches the oracle", {
  f <- make_graded_nets()
  binary <- function(net) {
    sp <- sparse_net(net, 5)
    sp$network@x[] <- 1
    sp$threshold <- 1
    sp$store_threshold <- 1
    sp
  }
  dense_binary <- function(net) {
    m <- (net$network >= 5) * 1
    diag(m) <- 0
    m
  }
  res <- rcomplex:::compare_specificity(binary(f$net1), binary(f$net2), f$ortho)
  ref <- reference_specificity(
    dense_binary(f$net1), dense_binary(f$net2), 1, 1, f$ortho,
    store1 = 1, store2 = 1
  )
  expect_spec_equal(res, ref)
})

test_that("the AUROC grid matches the oracle in both directions", {
  f <- make_graded_nets()
  res <- rcomplex:::compare_specificity(f$net1, f$net2, f$ortho)
  gf <- rcomplex:::.rank_grid_frac
  ref <- reference_specificity(
    f$net1$network, f$net2$network, 5, 5, f$ortho,
    grid_frac = gf
  )
  for (s in c("species1", "species2")) {
    g <- res[[paste0(s, ".auroc.grid")]]
    expect_identical(dim(g), c(nrow(res), length(gf)))
    expect_equal(unname(g), unname(ref[[paste0(s, ".auroc.grid")]]),
      tolerance = 1e-12
    )
    expect_identical(res[[paste0(s, ".n.cand")]], rep(30L, nrow(res)))
    # the grid is non-increasing along the fractions
    d <- t(apply(g, 1L, diff))
    expect_true(all(d <= 1e-12, na.rm = TRUE))
  }
})

test_that("the kernel rejects a grid that is not ascending in (0, 1]", {
  f <- make_graded_nets()
  run <- function(gf) {
    rcomplex:::.specificity_run(f$net1, f$net2, f$ortho, 1L, "both", gf)
  }
  expect_error(run(c(0.5, 0.1)), "strictly ascending")
  expect_error(run(c(0, 0.1)), "strictly ascending")
  expect_error(run(1.5), "strictly ascending")
  expect_false("species1.auroc.grid" %in% names(run(numeric(0))))
})
