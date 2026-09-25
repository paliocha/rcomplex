# Tests for compare_specificity() against the pure-R oracle
# reference_specificity() in helper-reference.R. Fixtures: make_graded_nets(),
# make_cmp_nets(), sparse_net() from helper-reference.R.

spec_cols <- function(prefix) {
  paste0(prefix, c(".neigh", ".mapped", ".auroc", ".p.val", ".jaccard"))
}
both_cols <- c(spec_cols("Species1"), spec_cols("Species2"))

expect_spec_equal <- function(res, ref, cols = both_cols) {
  testthat::expect_equal(as.list(res[cols]), as.list(ref[cols]),
    tolerance = 1e-12
  )
}


test_that("dense compare_specificity matches the oracle in both directions", {
  f <- make_graded_nets()
  res <- compare_specificity(f$net1, f$net2, f$ortho)
  ref <- reference_specificity(
    f$net1$network, f$net2$network, 5, 5, f$ortho
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c(
    "Species1", "Species2", "hog", both_cols,
    "Species1.effect.size", "Species2.effect.size"
  ))
  expect_equal(res$Species1, f$ortho$Species1)
  expect_equal(res$hog, f$ortho$hog)
  expect_spec_equal(res, ref)
  expect_identical(res$Species1.effect.size, res$Species1.auroc)
  expect_identical(res$Species2.effect.size, res$Species2.auroc)
})

test_that("a store above the lowest tier matches the oracle with store", {
  f <- make_graded_nets()
  res <- compare_specificity(
    sparse_net(f$net1, 5), sparse_net(f$net2, 5), f$ortho
  )
  ref <- reference_specificity(
    f$net1$network, f$net2$network, 5, 5, f$ortho,
    store1 = 5, store2 = 5
  )
  expect_spec_equal(res, ref)
  # tier 4 really is unstored: collapsing it to the bottom block changes
  # other candidates' AUROCs (e.g. B16 for anchor A04), hence the p-values
  dense <- compare_specificity(f$net1, f$net2, f$ortho)
  expect_false(isTRUE(all.equal(res$Species1.p.val, dense$Species1.p.val)))
})

test_that("a store holding every nonzero entry equals the dense result", {
  f <- make_graded_nets()
  dense <- compare_specificity(f$net1, f$net2, f$ortho)
  sp <- compare_specificity(
    sparse_net(f$net1, 4), sparse_net(f$net2, 4), f$ortho
  )
  expect_identical(sp, dense)
})

test_that("own-HOG orthologs leave the mapped set (paralog row)", {
  f <- make_cmp_nets()
  res <- compare_specificity(
    sparse_net(f$net1), sparse_net(f$net2), f$ortho
  )
  ref <- reference_specificity(
    f$net1$network, f$net2$network, f$net1$threshold, f$net2$threshold,
    f$ortho,
    store1 = f$net1$threshold, store2 = f$net2$threshold
  )
  expect_equal(nrow(res), 31L)
  expect_spec_equal(res, ref)
  # A_001 maps to B_001 and B_031: both must be absent from its mapped set
  rows <- which(res$Species1 == "A_001")
  expect_length(rows, 2L)
  thr <- f$net1$threshold
  n1 <- setdiff(names(which(f$net1$network[, "A_001"] >= thr)), "A_001")
  mapped <- setdiff(
    unique(f$ortho$Species2[f$ortho$Species1 %in% n1]), c("B_001", "B_031")
  )
  expect_equal(res$Species1.mapped[rows], rep(length(mapped), 2L))
})

test_that("p lies on the 1/n grid, AUROC in [0, 1], isolated genes NA", {
  f <- make_graded_nets()
  res <- compare_specificity(f$net1, f$net2, f$ortho)
  for (s in c("Species1", "Species2")) {
    p <- res[[paste0(s, ".p.val")]]
    a <- res[[paste0(s, ".auroc")]]
    expect_true(any(is.na(p)) && any(!is.na(p)))
    expect_equal(is.na(p), is.na(a))
    expect_true(all(abs(p * 30 - round(p * 30)) < 1e-9, na.rm = TRUE))
    expect_true(all(p > 0 & p <= 1, na.rm = TRUE))
    expect_true(all(a >= 0 & a <= 1, na.rm = TRUE))
  }
  iso <- res$Species1 %in% paste0("A", 28:30)
  expect_true(all(is.na(res$Species1.p.val[iso])))
  expect_true(all(is.na(res$Species2.p.val[iso])))
  expect_equal(res$Species1.mapped[iso], rep(0L, 3L))
  expect_false(any(is.na(res$Species1.p.val[res$Species1 %in% "A01"])))
})

test_that("directions restricts the columns and keeps the values", {
  f <- make_graded_nets()
  both <- compare_specificity(f$net1, f$net2, f$ortho)
  d12 <- compare_specificity(f$net1, f$net2, f$ortho, directions = "1to2")
  d21 <- compare_specificity(f$net1, f$net2, f$ortho, directions = "2to1")
  keys <- c("Species1", "Species2", "hog")
  expect_named(d12, c(keys, spec_cols("Species1"), "Species1.effect.size"))
  expect_named(d21, c(keys, spec_cols("Species2"), "Species2.effect.size"))
  expect_identical(d12, both[names(d12)])
  expect_identical(d21, both[names(d21)])
})

test_that("n_cores = 2 equals the serial result", {
  f <- make_graded_nets()
  serial <- compare_specificity(f$net1, f$net2, f$ortho)
  par <- compare_specificity(f$net1, f$net2, f$ortho, n_cores = 2L)
  expect_identical(par, serial)
})

test_that("a non-symmetric sparse pattern is an error", {
  g <- c("g1", "g2", "g3")
  m <- Matrix::sparseMatrix(
    i = 1L, j = 2L, x = 1, dims = c(3L, 3L), dimnames = list(g, g)
  )
  net <- list(network = m, threshold = 1, store_threshold = 1)
  ortho <- data.frame(Species1 = g, Species2 = g, hog = 1:3)
  expect_error(compare_specificity(net, net, ortho), "symmetric")
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
  res <- compare_specificity(binary(f$net1), binary(f$net2), f$ortho)
  ref <- reference_specificity(
    dense_binary(f$net1), dense_binary(f$net2), 1, 1, f$ortho,
    store1 = 1, store2 = 1
  )
  expect_spec_equal(res, ref)
})
