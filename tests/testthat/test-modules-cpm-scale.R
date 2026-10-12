# CPM resolution is read on the [0, 1] weight scale (0.3.2). Raw MR weights
# (~2e4) made resolution 1 effectively 0 and CPM returned one module.

# Three dense blocks of 15 genes at weight ~2e4, sparse cross-block edges
# just above the threshold so the graph is connected.
cpm_scale_fixture <- function() {
  set.seed(1)
  n <- 45L
  lab <- rep(1:3, each = 15L)
  m <- matrix(0, n, n)
  same <- outer(lab, lab, "==")
  m[same] <- stats::runif(sum(same), 19500, 20000)
  cross <- !same & matrix(stats::runif(n * n) < 0.1, n, n)
  m[cross] <- stats::runif(sum(cross), 18500, 19000)
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 0
  rownames(m) <- colnames(m) <- paste0("G", seq_len(n))
  list(net = list(network = m, threshold = 18000), lab = lab)
}

ari <- function(a, b) {
  igraph::compare(as.integer(a), as.integer(b), method = "adjusted.rand")
}

test_that("CPM resolution is read on the unit weight scale", {
  fx <- cpm_scale_fixture()
  res <- detect_modules(fx$net,
    objective_function = "CPM", resolution = 0.5, seed = 1
  )
  expect_gte(res$n_modules, 2L)
  expect_gt(ari(res$modules, fx$lab), 0.8)
})

test_that("modularity partitions do not depend on the weight scale", {
  fx <- cpm_scale_fixture()
  unit <- fx$net
  unit$network <- unit$network / max(unit$network)
  unit$threshold <- fx$net$threshold / max(fx$net$network)
  a <- detect_modules(fx$net,
    objective_function = "modularity", seed = 3
  )
  b <- detect_modules(unit,
    objective_function = "modularity", seed = 3
  )
  expect_equal(ari(a$modules, b$modules), 1)
})

test_that("the default is modularity at 1; CPM defaults to the density", {
  fx <- cpm_scale_fixture()
  mod <- detect_modules(fx$net, seed = 1)
  expect_identical(mod$params$objective_function, "modularity")
  expect_identical(mod$params$resolution, 1)
  res <- detect_modules(fx$net, objective_function = "CPM", seed = 1)
  g <- res$graph
  dens <- igraph::ecount(g) / choose(igraph::vcount(g), 2)
  expect_gt(dens, 1 / 3)
  expect_equal(res$params$resolution, dens)
  expect_gt(res$n_modules, 1L)
  expect_lt(res$n_modules, 45L)
  expect_gt(ari(res$modules, fx$lab), 0.8)
})
