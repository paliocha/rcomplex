# Tests for mr_block_network_cpp(): the blockwise sparse mutual-rank
# network must reproduce compute_network(sparse = TRUE) slot for slot.

block_fixture <- function(n, s, seed) {
  set.seed(seed)
  x <- matrix(rnorm(n * s), n, s)
  rownames(x) <- paste0("g", seq_len(n))
  x
}

# Standardise like Rfast::cora(): crossprod(zt) is the correlation matrix
block_zt <- function(x, cor_method) {
  rcomplex:::.standardise_for_cor(x, cor_method)
}

# No two distinct correlations in one column closer than 1e-12, so BLAS
# last-bit differences between the two paths cannot swap ranks
has_near_ties <- function(zt, abs_cor) {
  cm <- pmin(pmax(crossprod(zt), -1), 1)
  if (abs_cor) cm <- abs(cm)
  near <- apply(cm, 2, function(v) {
    d <- diff(sort(v))
    any(d > 0 & d < 1e-12)
  })
  any(near)
}

block_modes <- list(
  pearson_raw = list(cor = "pearson", log = FALSE, abs = FALSE),
  pearson_log = list(cor = "pearson", log = TRUE, abs = FALSE),
  spearman_raw = list(cor = "spearman", log = FALSE, abs = FALSE),
  pearson_abs = list(cor = "pearson", log = FALSE, abs = TRUE)
)

run_block <- function(zt, m, density = 0.03, store_density = 0.05,
                      block_size = 7L, n_cores = 1L) {
  rcomplex:::mr_block_network_cpp(
    zt, m$log, m$abs, density, store_density, block_size, n_cores
  )
}

# The slots both paths must agree on (block result or dense network)
net_slots <- function(net) {
  if (!is.null(net$network)) {
    net <- c(
      list(i = net$network@i, p = net$network@p, x = net$network@x),
      net[c("threshold", "store_threshold")]
    )
  }
  net[c("i", "p", "x", "threshold", "store_threshold")]
}

# Spearman correlations are rationals, so tied pairs land on neighbouring
# floats rather than equal ones; many samples (and a seed checked tie-free)
# keep the near-tie guard satisfiable.
block_fixtures <- list(
  pearson = list(c(60, 12, 1), c(300, 20, 2)),
  spearman = list(c(60, 1000, 1), c(300, 1000, 4))
)

for (mode in names(block_modes)) {
  m <- block_modes[[mode]]
  for (fx in block_fixtures[[m$cor]]) {
    test_that(paste("block network matches dense:", fx[1], mode), {
      x <- block_fixture(fx[1], fx[2], fx[3])
      ref <- compute_network(
        x,
        cor_method = m$cor, density = 0.03, store_density = 0.05,
        mr_log_transform = m$log, abs_cor = m$abs, sparse = TRUE
      )
      zt <- block_zt(x, m$cor)
      expect_false(has_near_ties(zt, m$abs))
      blk <- run_block(zt, m)
      expect_identical(net_slots(blk), net_slots(ref))
      for (bs in c(1L, ncol(zt))) {
        expect_identical(run_block(zt, m, block_size = bs), blk)
      }
      expect_identical(run_block(zt, m, n_cores = 2L), blk)
    })
  }
}

# The rank fraction starts at 3 * store_density: raw MR passes there on
# this fixture, log MR's looser bound forces widening (by 1.5x, up to 1).
test_that("block network widens the candidate fraction when needed", {
  x <- block_fixture(60, 12, 1)
  for (m in block_modes[c("pearson_raw", "pearson_log")]) {
    ref <- compute_network(
      x,
      cor_method = "pearson", density = 0.03, store_density = 0.1,
      mr_log_transform = m$log, abs_cor = FALSE, sparse = TRUE
    )
    blk <- run_block(block_zt(x, "pearson"), m, store_density = 0.1)
    expect_identical(net_slots(blk), net_slots(ref))
    if (m$log) {
      expect_gt(blk$fraction, 3 * 0.1)
    } else {
      expect_equal(blk$fraction, 3 * 0.1)
    }
  }
})

test_that("block network errors on NaN input", {
  zt <- block_zt(block_fixture(60, 12, 1), "pearson")
  zt[, 5] <- NaN
  expect_error(run_block(zt, block_modes$pearson_raw), "NaN")
})

# Through compute_network(): block_size must return the dense sparse object
for (mode in names(block_modes)) {
  m <- block_modes[[mode]]
  for (fx in block_fixtures[[m$cor]]) {
    test_that(paste("compute_network block_size matches:", fx[1], mode), {
      x <- block_fixture(fx[1], fx[2], fx[3])
      build <- function(bs) {
        # log MR can fall back to all pairs; the fallback test covers it
        suppressMessages(compute_network(
          x,
          cor_method = m$cor, density = 0.03, store_density = 0.05,
          mr_log_transform = m$log, abs_cor = m$abs, block_size = bs
        ))
      }
      ref <- build(NULL)
      for (bs in c(1L, 7L, nrow(x))) {
        expect_identical(build(bs), ref)
      }
    })
  }
}

test_that("compute_network block_size validates its arguments", {
  x <- block_fixture(60, 12, 1)
  expect_error(compute_network(x, block_size = 0), "positive whole")
  expect_error(compute_network(x, block_size = 2.5), "positive whole")
  expect_error(compute_network(x, block_size = "a"), "positive whole")
  expect_error(compute_network(x, block_size = c(1, 2)), "positive whole")
  expect_error(
    compute_network(x, sparse = FALSE, block_size = 7), "sparse = TRUE"
  )
  expect_error(
    compute_network(x, norm_method = "CLR", block_size = 7), "MR"
  )
  expect_error(
    compute_network(x, use_torch = TRUE, block_size = 7), "use_torch"
  )
})

test_that("compute_network block_size keeps the variance filter", {
  x <- block_fixture(60, 12, 1)
  x[3, ] <- 5
  x[10, ] <- x[10, ] * 1e-3
  for (mv in list(0, 1e-3)) {
    ref <- compute_network(x, min_var = mv)
    blk <- compute_network(x, min_var = mv, block_size = 7)
    expect_identical(blk$n_removed, ref$n_removed)
    expect_identical(rownames(blk$network), rownames(ref$network))
    expect_identical(blk, ref)
  }
  expect_identical(compute_network(x, block_size = 7)$n_removed, 1L)
})

test_that("mr_block() agrees on block and dense networks", {
  x <- block_fixture(60, 12, 1)
  genes <- paste0("g", c(2, 9, 17, 40))
  ref <- compute_network(x)
  blk <- compute_network(x, block_size = 7)
  expect_identical(mr_block(x, genes, blk), mr_block(x, genes, ref))
})

test_that("compute_network block_size reports the all-pairs fallback", {
  x <- block_fixture(60, 12, 1)
  expect_message(
    blk <- compute_network(
      x,
      store_density = 0.1, mr_log_transform = TRUE, block_size = 7
    ),
    "all pairs"
  )
  ref <- compute_network(x, store_density = 0.1, mr_log_transform = TRUE)
  expect_identical(blk, ref)
  expect_silent(compute_network(x, block_size = 7))
})
