test_that("compute_network returns correct structure", {
  set.seed(42)
  expr <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)

  result <- compute_network(expr, density = 0.1)

  expect_type(result, "list")
  expect_named(result, c("network", "threshold", "n_genes", "n_removed", "params"))
  expect_equal(result$n_genes, 10)
  expect_true(is.matrix(result$network))
  expect_equal(nrow(result$network), 10)
  expect_equal(ncol(result$network), 10)
  expect_equal(rownames(result$network), paste0("gene", 1:10))
  expect_equal(colnames(result$network), paste0("gene", 1:10))
})

test_that("network is symmetric with zero diagonal", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, density = 0.05)

  expect_equal(result$network, t(result$network))
  expect_equal(unname(diag(result$network)), rep(0, 20))
})

test_that("MR log_transform produces values in [0,1]", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, norm_method = "MR",
                            mr_log_transform = TRUE, density = 0.05)
  vals <- result$network[upper.tri(result$network)]
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1))
})

test_that("MR raw mode matches R reference", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  # Compute via package (raw MR mode)
  result <- compute_network(expr, cor_method = "pearson",
                            norm_method = "MR", mr_log_transform = FALSE,
                            density = 0.05)

  # Compute via R reference
  cor_mat <- cor(t(expr), method = "pearson")
  ref_net <- reference_mr_raw(cor_mat)

  # Values should match closely
  expect_equal(result$network, ref_net, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("CLR normalization produces non-negative values", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, norm_method = "CLR", density = 0.05)
  vals <- result$network[upper.tri(result$network)]
  expect_true(all(vals >= 0))
})

test_that("CLR matches R reference", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, cor_method = "pearson",
                            norm_method = "CLR", density = 0.05)

  cor_mat <- cor(t(expr), method = "pearson")
  ref_net <- reference_clr(cor_mat)

  expect_equal(result$network, ref_net, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("density threshold is in valid range", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, density = 0.05)
  vals <- result$network[upper.tri(result$network)]

  expect_true(result$threshold >= min(vals))
  expect_true(result$threshold <= max(vals))
})

test_that("density threshold matches R reference", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, norm_method = "MR",
                            mr_log_transform = FALSE, density = 0.05)

  ref_thr <- reference_density_threshold(result$network, 0.05)

  expect_equal(result$threshold, ref_thr, tolerance = 1e-10)
})

test_that("spearman correlation method works", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, cor_method = "spearman", density = 0.05)
  expect_true(is.matrix(result$network))
  expect_equal(result$params$cor_method, "spearman")
})

test_that("abs_cor option works", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, abs_cor = TRUE, density = 0.05)
  expect_true(result$params$abs_cor)
})

test_that("input validation works", {
  expect_error(compute_network(data.frame(a = 1:5)), "must be a matrix")
  mat <- matrix(1:10, nrow = 5)
  expect_error(compute_network(mat), "must have row names")
  rownames(mat) <- paste0("g", 1:5)
  expect_error(compute_network(mat, density = 0), "between 0 and 1")
  expect_error(compute_network(mat, density = 1), "between 0 and 1")
})

test_that("min_var removes constant genes", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)
  # Make 3 genes constant (use exact binary values to avoid fp accumulation)
  expr[1, ] <- 5.0
  expr[2, ] <- 0.0
  expr[3, ] <- -3.0

  result <- compute_network(expr, density = 0.1)
  expect_equal(result$n_removed, 3L)
  expect_equal(result$n_genes, 17)
  expect_equal(nrow(result$network), 17)
  expect_false("gene1" %in% rownames(result$network))
  expect_false("gene2" %in% rownames(result$network))
  expect_false("gene3" %in% rownames(result$network))
})

test_that("min_var threshold filters near-invariant genes", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)
  # Make gene1 nearly constant (tiny variance)
  expr[1, ] <- 5.0 + rnorm(10, sd = 1e-6)

  result_strict <- compute_network(expr, density = 0.1, min_var = 1e-8)
  expect_equal(result_strict$n_removed, 1L)
  expect_false("gene1" %in% rownames(result_strict$network))

  # Default min_var=0 keeps near-invariant genes (variance > 0)
  result_default <- compute_network(expr, density = 0.1)
  expect_equal(result_default$n_removed, 0L)
  expect_true("gene1" %in% rownames(result_default$network))
})

test_that("min_var=NULL disables filtering", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)
  expr[1, ] <- 5.0  # constant gene

  result <- compute_network(expr, density = 0.1, min_var = NULL)
  expect_equal(result$n_removed, 0L)
  expect_equal(result$n_genes, 20)
})

test_that("min_var errors when too few genes remain", {
  set.seed(42)
  expr <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(expr) <- paste0("gene", 1:5)
  # Huge threshold removes all genes
  expect_error(compute_network(expr, density = 0.1, min_var = 1e6),
               "Fewer than 3 genes")
})

test_that("min_var is stored in params", {
  set.seed(42)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)

  result <- compute_network(expr, density = 0.1, min_var = 0.5)
  expect_equal(result$params$min_var, 0.5)
})
