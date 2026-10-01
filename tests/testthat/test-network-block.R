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
  # exact ties count too: they may round apart under another BLAS
  near <- vapply(seq_len(ncol(cm)), function(j) {
    any(diff(sort(cm[-j, j])) < 1e-12)
  }, logical(1))
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

# The rank fraction starts at 3 * store_density (raw) or 2 * store_density
# (log); when the store threshold is too low to prove the candidates
# complete it widens to the fraction that threshold implies. On this
# fixture raw passes at 0.1 and log at 0.05; raw at 0.02 and 0.01 and log
# at 0.02 widen.
test_that("block network widens the candidate fraction when needed", {
  x <- block_fixture(60, 12, 1)
  zt <- block_zt(x, "pearson")
  cases <- list(
    list(m = block_modes$pearson_raw, sd = 0.1, widen = FALSE),
    list(m = block_modes$pearson_raw, sd = 0.02, widen = TRUE),
    list(m = block_modes$pearson_raw, sd = 0.01, widen = TRUE),
    list(m = block_modes$pearson_log, sd = 0.05, widen = FALSE),
    list(m = block_modes$pearson_log, sd = 0.02, widen = TRUE)
  )
  for (cs in cases) {
    ref <- compute_network(
      x,
      cor_method = "pearson", density = 0.01, store_density = cs$sd,
      mr_log_transform = cs$m$log, abs_cor = FALSE, sparse = TRUE
    )
    blk <- run_block(zt, cs$m, density = 0.01, store_density = cs$sd)
    expect_identical(net_slots(blk), net_slots(ref))
    f0 <- (if (cs$m$log) 2 else 3) * cs$sd
    expect_equal(blk$start_fraction, f0)
    expect_identical(blk$fraction > f0, cs$widen)
    expect_lt(blk$fraction, 1)
    if (cs$widen) {
      expect_message(
        compute_network(
          x,
          density = 0.01, store_density = cs$sd,
          mr_log_transform = cs$m$log, block_size = 7
        ),
        "widened"
      )
    }
  }
})

# Entries in {-0.5, 0, 0.5} over 4 samples make every cross-product an
# exact multiple of 0.25 in [-1, 1], identical under any summation order,
# so the kernel meets heavy rank ties (half-integer ranks, tie groups across
# the list cut-off, tied MR at the thresholds) with nothing left to BLAS.
test_that("block network matches the dense kernels under exact ties", {
  set.seed(5)
  zt <- matrix(sample(c(-0.5, 0, 0.5), 4 * 40, replace = TRUE), 4, 40)
  for (lg in c(FALSE, TRUE)) {
    m <- crossprod(zt) + 0
    rcomplex:::mutual_rank_inplace_cpp(m, lg, FALSE, 1L)
    t_s <- rcomplex:::density_threshold_cpp(m, 0.1)
    ref <- c(
      rcomplex:::extract_sparse_cpp(m, t_s, 1L),
      list(
        threshold = rcomplex:::density_threshold_cpp(m, 0.05),
        store_threshold = t_s
      )
    )
    m_mode <- list(log = lg, abs = FALSE)
    for (bs in c(1L, 7L, 40L)) {
      blk <- run_block(zt, m_mode,
        density = 0.05, store_density = 0.1,
        block_size = bs
      )
      expect_identical(net_slots(blk), net_slots(ref))
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
        compute_network(
          x,
          cor_method = m$cor, density = 0.03, store_density = 0.05,
          mr_log_transform = m$log, abs_cor = m$abs, block_size = bs
        )
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
  # beyond the integer range is capped at the gene count, not NA
  expect_identical(
    compute_network(x, block_size = 1e12),
    compute_network(x, block_size = 60L)
  )
  expect_error(
    compute_network(x, sparse = FALSE, block_size = 7), "sparse = TRUE"
  )
  expect_error(
    compute_network(x, norm_method = "CLR", block_size = 7), "MR"
  )
  expect_error(
    compute_network(x, use_torch = TRUE, block_size = 7),
    "block_size requires use_torch"
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

# Only a fraction of 1 falls back to all pairs: store_density 0.5 starts
# log MR there and 0.4 raw MR, the defaults stay below it for both.
test_that("compute_network block_size reports the all-pairs fallback", {
  x <- block_fixture(60, 12, 1)
  expect_message(
    blk <- compute_network(
      x,
      store_density = 0.5, mr_log_transform = TRUE, block_size = 7
    ),
    "all pairs"
  )
  ref <- compute_network(x, store_density = 0.5, mr_log_transform = TRUE)
  expect_identical(blk, ref)
  # raw MR starts at all pairs from store_density 1/3
  expect_message(
    blk <- compute_network(x, store_density = 0.4, block_size = 7),
    "all pairs"
  )
  expect_identical(blk, compute_network(x, store_density = 0.4))
  expect_silent(compute_network(x, block_size = 7))
  expect_silent(compute_network(x, mr_log_transform = TRUE, block_size = 7))
})
