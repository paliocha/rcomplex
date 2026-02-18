test_that("summarize_comparison returns correct structure", {
  comparison <- data.frame(
    Species1 = paste0("A_", 1:10),
    Species2 = paste0("B_", 1:10),
    OrthoGroup = rep(1:5, each = 2),
    Species1.neigh = rep(10, 10),
    Species1.ortho.neigh = rep(5, 10),
    Species1.neigh.overlap = c(3, 0, 2, 4, 1, 3, 2, 0, 1, 5),
    Species1.p.val.con = c(
      0.001, 1, 0.01, 0.0001, 0.5, 0.005, 0.05, 1, 0.3, 0.0001
    ),
    Species1.p.val.div = c(
      0.99, 0.01, 0.9, 0.999, 0.5, 0.99, 0.9, 0.01, 0.7, 0.999
    ),
    Species1.effect.size = c(5, 1, 3, 8, 1, 4, 2, 1, 1.5, 10),
    Species2.neigh = rep(8, 10),
    Species2.ortho.neigh = rep(4, 10),
    Species2.neigh.overlap = c(2, 0, 1, 3, 0, 2, 1, 0, 1, 4),
    Species2.p.val.con = c(0.01, 1, 0.1, 0.001, 1, 0.01, 0.1, 1, 0.5, 0.0001),
    Species2.p.val.div = c(
      0.9, 0.01, 0.8, 0.99, 0.01, 0.9, 0.8, 0.01, 0.5, 0.999
    ),
    Species2.effect.size = c(4, 1, 2, 6, 1, 3, 1.5, 1, 1, 8),
    stringsAsFactors = FALSE
  )

  result <- summarize_comparison(comparison)

  expect_type(result, "list")
  expect_named(result, c("results", "summary"))
  expect_s3_class(result$results, "data.frame")
  expect_type(result$summary, "list")
  expect_named(result$summary, c("gene_pairs", "genes", "orthogroups"))
})

test_that("zero-overlap rows are filtered by default", {
  comparison <- data.frame(
    Species1 = paste0("A_", 1:5),
    Species2 = paste0("B_", 1:5),
    OrthoGroup = 1:5,
    Species1.neigh = rep(10, 5),
    Species1.ortho.neigh = rep(5, 5),
    Species1.neigh.overlap = c(3, 0, 2, 0, 1),
    Species1.p.val.con = c(0.001, 1, 0.01, 1, 0.5),
    Species1.p.val.div = c(0.99, 0.01, 0.9, 0.01, 0.5),
    Species1.effect.size = c(5, 1, 3, 1, 2),
    Species2.neigh = rep(8, 5),
    Species2.ortho.neigh = rep(4, 5),
    Species2.neigh.overlap = c(2, 0, 1, 0, 1),
    Species2.p.val.con = c(0.01, 1, 0.1, 1, 0.5),
    Species2.p.val.div = c(0.9, 0.01, 0.8, 0.01, 0.5),
    Species2.effect.size = c(4, 1, 2, 1, 1),
    stringsAsFactors = FALSE
  )

  result <- summarize_comparison(comparison)
  expect_equal(nrow(result$results), 3)  # rows 2 and 4 filtered

  result_no_filter <- summarize_comparison(comparison, filter_zero = FALSE)
  expect_equal(nrow(result_no_filter$results), 5)
})

test_that("FDR correction is applied", {
  comparison <- data.frame(
    Species1 = paste0("A_", 1:5),
    Species2 = paste0("B_", 1:5),
    OrthoGroup = 1:5,
    Species1.neigh = rep(10, 5),
    Species1.ortho.neigh = rep(5, 5),
    Species1.neigh.overlap = rep(2, 5),
    Species1.p.val.con = c(0.001, 0.01, 0.02, 0.03, 0.04),
    Species1.p.val.div = c(0.9, 0.8, 0.7, 0.6, 0.5),
    Species1.effect.size = rep(3, 5),
    Species2.neigh = rep(8, 5),
    Species2.ortho.neigh = rep(4, 5),
    Species2.neigh.overlap = rep(2, 5),
    Species2.p.val.con = c(0.002, 0.02, 0.03, 0.04, 0.05),
    Species2.p.val.div = c(0.9, 0.8, 0.7, 0.6, 0.5),
    Species2.effect.size = rep(2, 5),
    stringsAsFactors = FALSE
  )

  result <- summarize_comparison(comparison)

  # FDR-adjusted p-values should be >= raw p-values
  expect_true(all(result$results$Species1.p.val.con >= c(
    0.001, 0.01, 0.02, 0.03, 0.04
  )))
  expect_true(all(result$results$Species2.p.val.con >= c(
    0.002, 0.02, 0.03, 0.04, 0.05
  )))
})

test_that("summary counts are correct", {
  comparison <- data.frame(
    Species1 = c("A_1", "A_1", "A_2"),
    Species2 = c("B_1", "B_2", "B_2"),
    OrthoGroup = c(1, 1, 2),
    Species1.neigh = rep(10, 3),
    Species1.ortho.neigh = rep(5, 3),
    Species1.neigh.overlap = rep(5, 3),
    Species1.p.val.con = c(0.001, 0.5, 0.001),
    Species1.p.val.div = c(0.99, 0.5, 0.99),
    Species1.effect.size = c(5, 1, 5),
    Species2.neigh = rep(8, 3),
    Species2.ortho.neigh = rep(4, 3),
    Species2.neigh.overlap = rep(4, 3),
    Species2.p.val.con = c(0.001, 0.5, 0.001),
    Species2.p.val.div = c(0.99, 0.5, 0.99),
    Species2.effect.size = c(4, 1, 4),
    stringsAsFactors = FALSE
  )

  result <- summarize_comparison(comparison, alpha = 0.05)

  expect_equal(result$summary$gene_pairs$total, 3)
  expect_equal(result$summary$orthogroups$total, 2)
})

test_that("empty comparison handled gracefully", {
  comparison <- data.frame(
    Species1 = character(0),
    Species2 = character(0),
    OrthoGroup = integer(0),
    Species1.neigh = integer(0),
    Species1.ortho.neigh = integer(0),
    Species1.neigh.overlap = integer(0),
    Species1.p.val.con = numeric(0),
    Species1.p.val.div = numeric(0),
    Species1.effect.size = numeric(0),
    Species2.neigh = integer(0),
    Species2.ortho.neigh = integer(0),
    Species2.neigh.overlap = integer(0),
    Species2.p.val.con = numeric(0),
    Species2.p.val.div = numeric(0),
    Species2.effect.size = numeric(0),
    stringsAsFactors = FALSE
  )

  result <- summarize_comparison(comparison)
  expect_equal(nrow(result$results), 0)
  expect_equal(result$summary$gene_pairs$total, 0L)
})

test_that("alternative='less' uses divergence p-values", {
  comparison <- data.frame(
    Species1 = paste0("A_", 1:5),
    Species2 = paste0("B_", 1:5),
    OrthoGroup = 1:5,
    Species1.neigh = rep(10, 5),
    Species1.ortho.neigh = rep(5, 5),
    Species1.neigh.overlap = c(0, 0, 0, 3, 5),
    Species1.p.val.con = c(1, 1, 1, 0.01, 0.001),
    Species1.p.val.div = c(0.001, 0.01, 0.02, 0.9, 0.99),
    Species1.effect.size = c(0, 0, 0, 3, 5),
    Species2.neigh = rep(8, 5),
    Species2.ortho.neigh = rep(4, 5),
    Species2.neigh.overlap = c(0, 0, 0, 2, 4),
    Species2.p.val.con = c(1, 1, 1, 0.01, 0.001),
    Species2.p.val.div = c(0.001, 0.01, 0.02, 0.9, 0.99),
    Species2.effect.size = c(0, 0, 0, 3, 5),
    stringsAsFactors = FALSE
  )

  # With alternative="less", should use .p.val.div for thresholding
  result <- summarize_comparison(comparison, alternative = "less", alpha = 0.05)

  # filter_zero defaults to FALSE for "less"
  expect_equal(nrow(result$results), 5)

  # FDR-adjusted divergence p-values for first three rows should be significant
  expect_true(result$results$Species1.p.val.div[1] < 0.05)
  expect_true(result$results$Species1.p.val.div[2] < 0.05)
  # Rows 4 and 5 have high div p-values, should not be significant
  expect_true(result$results$Species1.p.val.div[4] > 0.05)
})

test_that("alternative='less' disables zero-overlap filtering by default", {
  comparison <- data.frame(
    Species1 = paste0("A_", 1:3),
    Species2 = paste0("B_", 1:3),
    OrthoGroup = 1:3,
    Species1.neigh = rep(10, 3),
    Species1.ortho.neigh = rep(5, 3),
    Species1.neigh.overlap = c(0, 0, 2),
    Species1.p.val.con = c(1, 1, 0.01),
    Species1.p.val.div = c(0.001, 0.01, 0.9),
    Species1.effect.size = c(0, 0, 3),
    Species2.neigh = rep(8, 3),
    Species2.ortho.neigh = rep(4, 3),
    Species2.neigh.overlap = c(0, 0, 1),
    Species2.p.val.con = c(1, 1, 0.01),
    Species2.p.val.div = c(0.001, 0.01, 0.9),
    Species2.effect.size = c(0, 0, 2),
    stringsAsFactors = FALSE
  )

  # Zero-overlap rows are kept for divergence (the strongest signal)
  result <- summarize_comparison(comparison, alternative = "less")
  expect_equal(nrow(result$results), 3)

  # But can be overridden
  result_filtered <- summarize_comparison(comparison, alternative = "less",
                                          filter_zero = TRUE)
  expect_equal(nrow(result_filtered$results), 1)
})
