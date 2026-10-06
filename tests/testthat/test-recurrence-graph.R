# Tests for recurrence_graph() and recurrence_modules(). Fixture: three
# species, 180 ortholog groups, 20 of them with two copies (200 genes per
# species, H161-H180). Groups H1-H7 share one expression factor in every
# species; the rest is noise. With three species a pair's smallest p is
# about d^3, so the block is single-copy and the density low (0.03).

recur_fixture <- function(shuffle = FALSE) {
  hogs <- paste0("H", 1:180)
  two <- 161:180
  hog_of <- c(hogs, hogs[two])
  block <- hog_of %in% hogs[1:7]
  sp <- c("S1", "S2", "S3")
  nets <- list()
  map <- list()
  withr::with_seed(1L, {
    for (k in seq_along(sp)) {
      x <- matrix(stats::rnorm(200 * 20), 200)
      x[block, ] <- x[block, ] +
        4 * outer(rep(1, sum(block)), stats::rnorm(20))
      if (shuffle) x <- t(apply(x, 1L, sample))
      g <- paste0(sp[k], "_g", 1:200)
      rownames(x) <- g
      nets[[sp[k]]] <- compute_network(x, density = 0.03)
      map[[k]] <- data.frame(species = sp[k], gene = g, hog = hog_of)
    }
  })
  list(nets = nets, map = do.call(rbind, map), block = hogs[1:7])
}

fx <- recur_fixture()
rg <- recurrence_graph(fx$nets, fx$map)


test_that("planted co-expressed groups recur in all species", {
  e <- rg$edges
  blk <- e[e$hog1 %in% fx$block & e$hog2 %in% fx$block, ]
  expect_identical(nrow(blk), as.integer(choose(7L, 2L)))
  expect_true(all(blk$K == 3L))
  expect_true(all(blk$profile == "111"))
  expect_true(all(blk$q.val < 0.05))
  expect_identical(rg$copies["H161", "S1"], 2L)
  expect_true(all(rg$species$n_genes == 200L))
  expect_output(print(rg), "significant")
})


test_that("shuffled expression gives no calls and the predicted K counts", {
  sh <- recur_fixture(shuffle = TRUE)
  r <- recurrence_graph(sh$nets, sh$map)
  expect_identical(sum(r$edges$q.val < 0.1), 0L)
  # Poisson-binomial expectation over every group pair
  cp <- r$copies
  pr <- utils::combn(nrow(cp), 2L)
  v <- cp[pr[1L, ], ] * cp[pr[2L, ], ]
  p <- v
  for (s in 1:3) p[, s] <- 1 - (1 - r$species$density[s])^v[, s]
  ex2 <- sum(rcomplex:::.poibin_upper(p, rep(2L, nrow(p))))
  ex3 <- sum(rcomplex:::.poibin_upper(p, rep(3L, nrow(p))))
  obs2 <- sum(r$edges$K >= 2L)
  obs3 <- sum(r$edges$K == 3L)
  expect_lt(abs(obs2 - ex2), 4 * sqrt(ex2))
  expect_lt(abs(obs3 - ex3), 4 * sqrt(ex3) + 2)
})


test_that("the Poisson-binomial tail matches brute-force enumeration", {
  p <- withr::with_seed(1L, matrix(stats::runif(20), 4))
  grid <- as.matrix(expand.grid(rep(list(0:1), 5)))
  bf <- vapply(1:4, function(r) {
    pr <- apply(grid, 1L, function(b) {
      prod(ifelse(b == 1, p[r, ], 1 - p[r, ]))
    })
    sum(pr[rowSums(grid) >= 3])
  }, numeric(1))
  expect_equal(rcomplex:::.poibin_upper(p, rep(3L, 4)), bf)
})


test_that("Leiden recovers the planted block and feeds module_auroc()", {
  m <- recurrence_modules(rg, fx$map, min_size = 5L, seed = 1L)
  hit <- vapply(m$hog_modules, function(h) all(fx$block %in% h), TRUE)
  expect_identical(sum(hit), 1L)
  expect_lte(length(m$hog_modules[[which(hit)]]), 15L)
  expect_named(m$modules, c("S1", "S2", "S3"))
  m1 <- m$modules$S1
  expect_identical(as_modules(m1), m1)
  expect_gte(length(m1$module_genes[[names(which(hit))]]), 7L)
  ortho <- merge(
    fx$map[fx$map$species == "S1", c("gene", "hog")],
    fx$map[fx$map$species == "S2", c("gene", "hog")],
    by = "hog"
  )
  names(ortho) <- c("hog", "Species1", "Species2")
  r <- module_auroc(m1, fx$nets$S1, fx$nets$S2, ortho,
    n_null = 5L, max_draws = 10L, batch = 5L, seed = 1L
  )
  expect_true(names(which(hit)) %in% r$module)
})


test_that("the anchored densest subgraph contains the planted block", {
  m <- recurrence_modules(rg, fx$map,
    method = "anchored", anchors = "H3", min_size = 5L
  )
  expect_named(m$hog_modules, "H3")
  expect_true(all(fx$block %in% m$hog_modules$H3))
  expect_message(
    recurrence_modules(rg, fx$map,
      method = "anchored", anchors = c("H3", "nope"), min_size = 5L
    ),
    "not in rg"
  )
})


test_that("the densest subgraph is exact on a K4 with a pendant anchor", {
  g <- igraph::graph_from_literal(
    A - B, B - C, B - D, B - E, C - D, C - E, D - E
  )
  igraph::E(g)$weight <- 1
  expect_setequal(rcomplex:::.anchored_densest(g, "A"), LETTERS[1:5])
})
