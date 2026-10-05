sub_net <- function(seed, shuffle = FALSE, n = 200L, s = 20L) {
  x <- withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * s), n)
    for (b in 1:4) {
      i <- (b - 1L) * 30L + 1:30
      x[i, ] <- x[i, ] + 1.5 * matrix(stats::rnorm(s), 30L, s, byrow = TRUE)
    }
    if (shuffle) x <- t(apply(x, 1L, sample))
    x
  })
  rownames(x) <- paste0("g", seq_len(n))
  compute_network(x, density = 0.05)
}
sub_map <- function(sp, n = 200L) {
  data.frame(
    species = rep(sp, each = n),
    gene = rep(paste0("g", seq_len(n)), length(sp)),
    hog = rep(paste0("H", seq_len(n)), length(sp))
  )
}
sp3 <- c("A", "B", "C")

test_that("a planted shared block puts S far above the null", {
  nets <- stats::setNames(lapply(1:3, sub_net), sp3)
  r <- subspace_preservation(nets, sub_map(sp3),
    K = 5L, n_null = 30L,
    seed = 1L
  )
  expect_s3_class(r, "subspace_preservation")
  expect_equal(nrow(r$pairs), 6L)
  expect_true(all(r$pairs$z > 3))
  expect_true(all(is.na(diag(r$S))))
  # a 1:1 map makes S symmetric
  expect_equal(r$S, t(r$S), tolerance = 1e-8)
  expect_output(print(r), "Subspace preservation over 3 species")
})

test_that("shuffled expression gives z about 0", {
  nets <- stats::setNames(lapply(1:3, sub_net, shuffle = TRUE), sp3)
  r <- subspace_preservation(nets, sub_map(sp3),
    K = 5L, n_null = 30L,
    seed = 1L
  )
  expect_true(all(abs(r$pairs$z) < 3))
})

test_that("identity map between a network and itself gives S = 1", {
  net <- sub_net(1L)
  r <- subspace_preservation(list(A = net, B = net), sub_map(c("A", "B")),
    K = 10L, n_null = 5L, seed = 1L
  )
  expect_equal(r$pairs$S, c(1, 1), tolerance = 1e-8)
  expect_error(
    subspace_preservation(list(A = net, B = net), sub_map(c("A", "B")),
      K = 500L
    ),
    "largest connected component"
  )
})

test_that("the result does not depend on n_cores", {
  nets <- stats::setNames(lapply(1:3, sub_net), sp3)
  a <- subspace_preservation(nets, sub_map(sp3),
    K = 5L, n_null = 10L,
    seed = 3L
  )
  b <- subspace_preservation(nets, sub_map(sp3),
    K = 5L, n_null = 10L,
    n_cores = 2L, seed = 3L
  )
  expect_identical(a, b)
})

test_that("as_preservation_matrix() feeds preservation_matrix_test()", {
  sp8 <- paste0("S", 1:8)
  nets <- stats::setNames(lapply(1:8, sub_net, n = 150L), sp8)
  r <- subspace_preservation(nets, sub_map(sp8, 150L),
    K = 4L,
    n_null = 10L, seed = 1L
  )
  pm <- as_preservation_matrix(r)
  expect_equal(nrow(pm), choose(8, 2))
  expect_equal(
    pm$Zsummary_std[1], (r$z["S1", "S2"] + r$z["S2", "S1"]) / 2
  )
  group <- stats::setNames(rep(c("annual", "perennial"), 4), sp8)
  block <- stats::setNames(rep(1:4, each = 2), sp8)
  # four blocks of two cannot reach p < 0.05; the test says so
  expect_warning(
    res <- preservation_matrix_test(pm, group, block = block),
    "smallest attainable"
  )
  expect_true(is.finite(res$observed))
  expect_true(res$p_free > 0 && res$p_free <= 1)
})
