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


# Tiny matched-networks fixture for the missing-pair and RNG-state tests:
# a disjoint perfect matching gives every gene exactly one neighbour, so a
# rewired permutation in which no ortholog pair shares a mapped neighbour
# has zero overlap > 0 rows and the species pair is absent from the null
# edges entirely. With seed = 1 the first of 6 permutations is such a run.
make_match_nets <- function() {
  mk <- function(prefix) {
    n <- 8L
    g <- paste0(prefix, seq_len(n))
    m <- matrix(0, n, n, dimnames = list(g, g))
    for (e in list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))) {
      m[e[1], e[2]] <- m[e[2], e[1]] <- 10
    }
    list(network = m, threshold = 5)
  }
  list(networks = list(A = sparse_net(mk("A")), B = sparse_net(mk("B"))),
       ortho = data.frame(Species1 = paste0("A", 1:8),
                          Species2 = paste0("B", 1:8),
                          hog = paste0("H", 1:8),
                          stringsAsFactors = FALSE))
}


test_that("null runs missing a species pair record 0 for the built-in statistic", {
  d <- make_match_nets()
  res <- coexpressolog_null(
    d$networks, d$ortho, n_perm = 6L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  )
  expect_identical(res$statistic, c("A~B", "total"))
  null_mat <- attr(res, "null")
  expect_true(all(is.finite(null_mat)))
  # the pair-less rewiring is a null observation of 0, not an abort
  expect_equal(unname(null_mat[1L, "A~B"]), 0)
})


test_that("a user statistic missing a name still errors", {
  d <- make_match_nets()
  per_pair <- function(edges) {
    if (is.null(edges) || nrow(edges) == 0L) return(c(total = 0))
    pair <- paste(edges$species1, edges$species2, sep = "~")
    counts <- vapply(split(edges$type == "conserved", pair), sum,
                     numeric(1))
    c(counts, total = sum(counts))
  }
  expect_error(
    coexpressolog_null(d$networks, d$ortho, statistic = per_pair,
                       n_perm = 6L, seed = 1L,
                       pval_combine = "max", pi0_method = "none"),
    "statistic is missing"
  )
})


test_that("the serial path restores the caller's RNG state", {
  d <- make_match_nets()
  set.seed(11)
  before <- runif(5)
  set.seed(11)
  invisible(coexpressolog_null(
    d$networks, d$ortho, n_perm = 2L, seed = 3L,
    pval_combine = "max", pi0_method = "none"
  ))
  after <- runif(5)
  # set.seed(seed + b) inside the serial loop must not leak: the caller's
  # stream continues exactly where the observed run left it
  expect_equal(after, before)
})
