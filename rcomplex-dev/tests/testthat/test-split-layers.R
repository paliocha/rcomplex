# split_layers(): OLS projection on a per-sample block.

split_fixture <- function() {
  x <- withr::with_seed(1L, matrix(stats::rnorm(30L * 20L), 30L))
  dimnames(x) <- list(paste0("g", 1:30), paste0("s", 1:20))
  block <- factor(rep(paste0("t", 1:5), each = 4L))
  x <- x + outer(stats::rnorm(30L), as.integer(block))
  list(x = x, block = block, sl = split_layers(x, block))
}


test_that("the layers are the OLS fit on the block", {
  f <- split_fixture()
  sl <- f$sl
  expect_s3_class(sl, "split_layers")
  expect_identical(dimnames(sl$wiring), dimnames(f$x))
  expect_identical(colnames(sl$deployment), levels(f$block))
  res_means <- t(rowsum(t(sl$wiring), f$block)) / 4
  expect_lt(max(abs(res_means)), 1e-12)
  expect_equal(sl$wiring + sl$deployment[, f$block], f$x,
    ignore_attr = TRUE
  )
  for (g in c("g1", "g7", "g30")) {
    expect_equal(
      unname(sl$r2[g]), summary(stats::lm(f$x[g, ] ~ f$block))$r.squared
    )
  }
  expect_identical(sl$df_residual, 15L)
  expect_identical(sl$block$sizes, c(
    t1 = 4L, t2 = 4L, t3 = 4L, t4 = 4L,
    t5 = 4L
  ))
})


test_that("one gene's residuals correlate -1/(n_b - 1) within a block", {
  f <- split_fixture()
  e <- f$sl$wiring["g1", ]
  r <- vapply(split(e, f$block), function(v) {
    p <- outer(v, v)
    mean(p[upper.tri(p)]) / mean(v^2)
  }, numeric(1))
  expect_equal(mean(r), -1 / 3)
})


test_that("bad blocks are refused", {
  f <- split_fixture()
  b <- f$block
  b[1] <- NA
  expect_error(split_layers(f$x, b), "NA")
  expect_error(split_layers(f$x, f$block[-1]), "one entry per sample")
  expect_error(split_layers(f$x, rep("a", 20L)), "two levels")
  expect_error(
    split_layers(f$x, c("x", rep(c("a", "b"), each = 10L)[-1])),
    "singleton: x"
  )
})


test_that("the wiring layer feeds compute_network and null_network", {
  f <- split_fixture()
  net <- compute_network(f$sl$wiring, density = 0.1)
  nn <- null_network(f$sl$wiring, net, block = f$block, seed = 1)
  expect_identical(rownames(nn$network), rownames(net$network))
  expect_output(print(f$sl), "df_residual: 15")
})


test_that("SummarizedExperiment input matches the matrix", {
  skip_if_not_installed("SummarizedExperiment")
  f <- split_fixture()
  se <- SummarizedExperiment::SummarizedExperiment(list(vst = f$x))
  expect_equal(split_layers(se, f$block), f$sl)
})
