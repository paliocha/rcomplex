# coexpressolog_null(): degree-preserving edge-swap null (P6).
#
# Uses the complex_py fixture data (sparse networks, ortholog-restricted
# gene universe) with n_perm = 19 and fixed seeds. Settings passed through
# `...` (pval_combine = "max", pi0_method = "none") are deterministic and
# reproduce the canonical fixture calls in the observed run.

null_fx <- function(...) testthat::test_path("fixtures", "complex_py", ...)

load_null_fixture <- local({
  cache <- NULL
  function() {
    skip_if_not(all(file.exists(null_fx(c("ortho_pairs.tsv", "sp1_expr.tsv",
                                          "sp2_expr.tsv")))))
    if (!is.null(cache)) return(cache)
    read_expr <- function(file) {
      d <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
      x <- as.matrix(d[, -1])
      rownames(x) <- d$Genes
      x
    }
    ortho <- read.delim(null_fx("ortho_pairs.tsv"), stringsAsFactors = FALSE)
    x1 <- read_expr(null_fx("sp1_expr.tsv"))
    x2 <- read_expr(null_fx("sp2_expr.tsv"))
    x1 <- x1[rownames(x1) %in% ortho$Species1, ]
    x2 <- x2[rownames(x2) %in% ortho$Species2, ]
    networks <- list(
      sp1 = compute_network(x1, density = 0.03, sparse = TRUE),
      sp2 = compute_network(x2, density = 0.03, sparse = TRUE)
    )
    cache <<- list(networks = networks, ortho = ortho)
    cache
  }
})

test_that("coexpressolog_null requires sparse networks", {
  d <- make_cmp_nets()
  nets_dense <- list(A = d$net1, B = d$net2)
  expect_error(
    coexpressolog_null(nets_dense, d$ortho, n_perm = 2L),
    "as_sparse_network"
  )
  # mixed dense/sparse is rejected too
  nets_mixed <- list(A = sparse_net(d$net1), B = d$net2)
  expect_error(
    coexpressolog_null(nets_mixed, d$ortho, n_perm = 2L),
    "as_sparse_network"
  )
})

test_that("observed conserved calls exceed the rewired null", {
  d <- load_null_fixture()
  res <- coexpressolog_null(
    d$networks, d$ortho, n_perm = 19L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  )

  expect_s3_class(res, "data.frame")
  # one row per species pair plus the total
  expect_identical(res$statistic, c("sp1~sp2", "total"))
  expect_identical(
    names(res),
    c("statistic", "observed", "null_mean", "null_sd", "null_max",
      "fold", "p_emp")
  )

  total <- res[res$statistic == "total", ]
  expect_gt(total$observed, total$null_max)
  expect_equal(total$p_emp, 1 / 20)
  # two-species fixture: the pair row equals the total row
  expect_equal(res$observed[1], total$observed)
  expect_equal(res$p_emp[1], 1 / 20)

  null_mat <- attr(res, "null")
  expect_identical(dim(null_mat), c(19L, 2L))
  expect_identical(colnames(null_mat), c("sp1~sp2", "total"))
  expect_equal(res$null_mean, vapply(1:2, function(j) mean(null_mat[, j]),
                                     numeric(1)))
  expect_equal(res$null_max, vapply(1:2, function(j) max(null_mat[, j]),
                                    numeric(1)))
})

test_that("shuffled orthologs give a non-significant null", {
  d <- load_null_fixture()
  ortho_shuf <- d$ortho
  set.seed(99)
  ortho_shuf$Species2 <- sample(ortho_shuf$Species2)
  res <- coexpressolog_null(
    d$networks, ortho_shuf, n_perm = 19L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  )
  expect_gt(res$p_emp[res$statistic == "total"], 0.05)
})

test_that("n_cores = 2 reproduces the serial result", {
  skip_if(.Platform$OS.type != "unix")
  d <- load_null_fixture()
  res1 <- coexpressolog_null(
    d$networks, d$ortho, n_perm = 5L, seed = 7L, n_cores = 1L,
    pval_combine = "max", pi0_method = "none"
  )
  res2 <- coexpressolog_null(
    d$networks, d$ortho, n_perm = 5L, seed = 7L, n_cores = 2L,
    pval_combine = "max", pi0_method = "none"
  )
  expect_equal(res1, res2)
})
