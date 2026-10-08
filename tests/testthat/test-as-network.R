# Tests for as_network(): networks built outside rcomplex.

an_expr <- function(seed, prefix, n = 60L) {
  set.seed(seed)
  x <- matrix(rnorm(n * 12L), n, 12L)
  rownames(x) <- sprintf("%s%02d", prefix, seq_len(n))
  x
}

test_that("a compute_network() store round-trips to identical edges", {
  na <- compute_network(an_expr(1, "A"))
  nb <- compute_network(an_expr(2, "B"))
  ra <- as_network(na$network)
  rb <- as_network(nb$network)
  expect_identical(ra$network, na$network)
  expect_identical(ra$threshold, na$threshold)
  expect_identical(ra$store_threshold, na$store_threshold)
  ortho <- data.frame(
    gene1 = rownames(na$network),
    gene2 = rownames(nb$network),
    hog = paste0("H", seq_len(60))
  )
  e1 <- find_coexpressologs(list(A = na, B = nb), ortho, seed = 1)
  e2 <- find_coexpressologs(list(A = ra, B = rb), ortho, seed = 1)
  expect_identical(e2, e1)
})

test_that("dense, sparse and edge-list input give one network", {
  dense <- compute_network(an_expr(3, "g"), sparse = FALSE)$network
  from_dense <- as_network(dense)
  keep <- c("network", "threshold", "store_density", "store_threshold")
  expect_identical(
    from_dense[keep],
    compute_network(an_expr(3, "g"))[keep]
  )
  m <- from_dense$network
  up <- Matrix::summary(Matrix::triu(m))
  el <- data.frame(
    gene1 = rownames(m)[up$i], gene2 = rownames(m)[up$j],
    weight = up$x
  )
  from_el <- as_network(el, genes = rownames(m))
  expect_identical(from_el$network, m)
  expect_identical(from_el$threshold, from_dense$threshold)
  sym <- methods::as(m, "symmetricMatrix")
  expect_identical(as_network(sym)$network, m)
})

test_that("genes subsets a matrix and adds isolated genes to an edge list", {
  el <- data.frame(
    gene1 = c("a", "b", "c"), gene2 = c("b", "c", "a"),
    weight = c(3, 2, 1)
  )
  net <- as_network(el, density = 0.5)
  expect_identical(rownames(net$network), c("a", "b", "c"))
  expect_identical(net$threshold, 2)
  wide <- as_network(el, density = 0.2, genes = c("a", "b", "c", "d"))
  expect_identical(dim(wide$network), c(4L, 4L))
  expect_identical(wide$threshold, 3)
  sub <- as_network(wide$network, density = 0.5, genes = c("c", "a", "b"))
  expect_identical(rownames(sub$network), c("c", "a", "b"))
  expect_identical(sub$threshold, 2)
})

test_that("as_network refuses bad input", {
  m <- matrix(c(0, 1, 2, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(as_network(m), "symmetric")
  expect_error(as_network(unname(m + t(m))), "gene names")
  expect_error(as_network(m + t(m), genes = "z"), "z")
  expect_error(as_network(m + t(m), density = 1), "density")
  expect_error(as_network(data.frame(a = 1)), "gene1, gene2, weight")
  el <- data.frame(
    gene1 = c("a", "b"), gene2 = c("b", "a"),
    weight = c(1, 2)
  )
  expect_error(as_network(el), "two weights")
  el4 <- data.frame(gene1 = "a", gene2 = "b", weight = 1)
  expect_error(
    as_network(el4, density = 0.9, genes = letters[1:4]),
    "fewer"
  )
  expect_error(as_network("x"), "matrix")
  expect_error(as_network(m + t(m)), "3 genes")
})

test_that("find_cliques and module_preservation run on as_network input", {
  fx <- pres_fixture()
  net_a <- as_network(fx$netA$network)
  net_b <- as_network(fx$netB$network)
  mods <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(mods, net_a, net_b, fx$ortho,
    n_perm = 30L, seed = 2
  )
  expect_identical(nrow(pres$preservation), fx$n_mod)
  edges <- find_coexpressologs(list(A = net_a, B = net_b), fx$ortho,
    seed = 1
  )
  cl <- find_cliques(edges, c("A", "B"), min_species = 2L)
  expect_s3_class(cl, "data.frame")
  expect_gt(nrow(cl), 0L)
})
