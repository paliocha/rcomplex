# Tests for module_auroc(). Fixture: two random 400-gene species, 1:1
# orthologs, optionally a planted module (genes 1-20 share one factor in
# both species).

auroc_fixture <- function(plant = FALSE) {
  withr::with_seed(1L, {
    n <- 400L
    s <- 20L
    x1 <- matrix(stats::rnorm(n * s), n)
    x2 <- matrix(stats::rnorm(n * s), n)
    if (plant) {
      x1[1:20, ] <- x1[1:20, ] + 2 * outer(rep(1, 20), stats::rnorm(s))
      x2[1:20, ] <- x2[1:20, ] + 2 * outer(rep(1, 20), stats::rnorm(s))
    }
  })
  rownames(x1) <- paste0("A", seq_len(n))
  rownames(x2) <- paste0("B", seq_len(n))
  list(
    net1 = compute_network(x1, density = 0.05),
    net2 = compute_network(x2, density = 0.05),
    ortho = data.frame(
      Species1 = rownames(x1), Species2 = rownames(x2),
      hog = paste0("H", seq_len(n)), stringsAsFactors = FALSE
    )
  )
}

# 40 modules of 10 genes over species-1 genes in fixed blocks
block_modules <- function(fx) {
  g <- rownames(fx$net1$network)
  as_modules(stats::setNames(rep(paste0("m", 1:40), each = 10L), g))
}

fx0 <- auroc_fixture()
mods0 <- block_modules(fx0)


test_that("random modules are calibrated under H0", {
  r <- module_auroc(mods0, fx0$net1, fx0$net2, fx0$ortho,
    max_draws = 300L, seed = 1L
  )
  expect_identical(nrow(r), 40L)
  expect_lte(mean(r$p.val < 0.05), 0.12)
  expect_lt(abs(mean(r$z)), 0.3)
  expect_true(all(abs(r$degree_auroc - 0.5) < 0.4))
})


test_that("a planted conserved module reaches the floor with z > 3", {
  fx <- auroc_fixture(plant = TRUE)
  mods <- as_modules(list(plant = paste0("A", 1:20)))
  r <- module_auroc(mods, fx$net1, fx$net2, fx$ortho,
    max_draws = 1000L, seed = 1L
  )
  expect_identical(r$n_draws, 1000L)
  expect_identical(r$n_exceed, 0L)
  expect_equal(r$p.val, 1 / 1001)
  expect_gt(r$z, 3)
})


test_that("sequential stopping stops null modules early", {
  r <- module_auroc(mods0, fx0$net1, fx0$net2, fx0$ortho,
    n_null = 50L, batch = 50L, max_draws = 2000L, h = 3L, seed = 2L
  )
  expect_true(all(r$n_draws %% 50L == 0L))
  expect_gt(mean(r$n_draws < 2000L), 0.8)
  stopped <- r$n_draws < 2000L
  expect_true(all(r$n_exceed[stopped] >= 3L))
  expect_equal(r$p.val, r$p.val.gt + r$p.val.eq)
  expect_equal(r$p.val, (1 + r$n_exceed) / (1 + r$n_draws))
})


test_that("the result does not depend on n_cores", {
  run <- function(k) {
    module_auroc(mods0, fx0$net1, fx0$net2, fx0$ortho,
      n_null = 50L, batch = 50L, max_draws = 200L, n_cores = k, seed = 3L
    )
  }
  expect_identical(run(1L), run(2L))
})


test_that("attributes are attached and q-values lie in [0, 1]", {
  r <- module_auroc(mods0, fx0$net1, fx0$net2, fx0$ortho,
    n_null = 50L, max_draws = 50L, seed = 4L
  )
  expect_named(r, c(
    "module", "n_hogs", "n_genes_2", "auroc", "degree_auroc", "null_mean",
    "null_sd", "z", "n_draws", "n_exceed", "p.val", "p.val.gt", "p.val.eq",
    "q.val", "n_exhausted"
  ))
  expect_identical(r$module, names(mods0$module_genes))
  expect_s3_class(attr(r, "resolution"), "pvalue_resolution")
  expect_identical(attr(r, "params")$max_draws, 50L)
  expect_true(all(r$q.val >= 0 & r$q.val <= 1))
})


test_that("modules with too few ortholog groups are dropped", {
  mods <- as_modules(list(small = c("A1", "A2"), big = paste0("A", 3:12)))
  expect_message(
    r <- module_auroc(mods, fx0$net1, fx0$net2, fx0$ortho,
      n_null = 20L, max_draws = 20L, seed = 1L
    ),
    "1 module"
  )
  expect_identical(r$module, "big")
})


test_that("drop_within_hog is inert 1:1 and removes paralog edges", {
  run <- function(ortho, drop) {
    module_auroc(mods0, fx0$net1, fx0$net2, ortho,
      n_null = 30L, max_draws = 30L, drop_within_hog = drop, seed = 5L
    )
  }
  on <- run(fx0$ortho, TRUE)
  off <- run(fx0$ortho, FALSE)
  attr(on, "params") <- attr(off, "params") <- NULL
  expect_identical(on, off)

  # species-2 genes paired into two-copy groups
  par <- fx0$ortho
  par$hog <- paste0("H", (seq_len(nrow(par)) + 1L) %/% 2L)
  op <- rcomplex:::.ortholog_pair_index(fx0$net1, fx0$net2, par)
  tgt <- function(drop) {
    rcomplex:::.module_auroc_target(fx0$net2, op$orthologs, op$sp2_idx, drop)
  }
  kept <- tgt(FALSE)
  dropped <- tgt(TRUE)
  col <- rep.int(seq_along(kept$deg), kept$deg)
  n_within <- sum((kept$i + 1L + 1L) %/% 2L == (col + 1L) %/% 2L)
  expect_gt(n_within, 0L)
  expect_identical(length(kept$i) - length(dropped$i), n_within)
  dcol <- rep.int(seq_along(dropped$deg), dropped$deg)
  expect_false(any((dropped$i + 2L) %/% 2L == (dcol + 1L) %/% 2L))
})
