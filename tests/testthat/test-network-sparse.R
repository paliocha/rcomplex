# Tests for the sparse network object (P2):
# compute_network(sparse = TRUE) as the default and as_sparse_network().


# ---- (a) compute_network(sparse = TRUE) and as_sparse_network() ----

test_that(
  "compute_network(sparse = TRUE) is the default and matches as_sparse_network",
  {
    set.seed(42)
    expr <- matrix(rnorm(500), nrow = 50, ncol = 10)
    rownames(expr) <- paste0("g", sprintf("%02d", 1:50))

    dense <- compute_network(expr, density = 0.03, sparse = FALSE)
    sp_default <- compute_network(expr, density = 0.03)
    sp_explicit <- compute_network(expr, density = 0.03, sparse = TRUE)

    expect_s4_class(sp_default$network, "dgCMatrix")
    expect_equal(sp_default, sp_explicit)
    expect_equal(sp_default, rcomplex:::as_sparse_network(dense, 0.05))

    expect_named(sp_default, c(
      "network", "threshold", "n_genes", "n_removed",
      "params", "store_density", "store_threshold"
    ))
    # dense object is exactly today's: no store_* fields
    expect_named(dense, c(
      "network", "threshold", "n_genes", "n_removed",
      "params"
    ))
    expect_true(is.matrix(dense$network))

    # analysis threshold and store fields
    expect_equal(sp_default$threshold, dense$threshold)
    expect_true(sp_default$store_threshold <= sp_default$threshold)
    expect_equal(sp_default$store_density, 0.05)
    expect_equal(sp_default$params$store_density, 0.05)

    # stored entries: >= store_threshold, no diagonal, both triangles
    m <- sp_default$network
    expect_true(all(m@x >= sp_default$store_threshold))
    expect_false(any(m@i == rep.int(seq_len(ncol(m)) - 1L, diff(m@p))))
    expect_equal(m, Matrix::t(m))

    # sparse matrix equals the thresholded dense matrix
    expect_equal(m, dense_to_dgc(dense$network, sp_default$store_threshold))
  }
)


test_that("store_density defaults to max(density, 0.05) and is validated", {
  set.seed(42)
  expr <- matrix(rnorm(300), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("g", sprintf("%02d", 1:30))

  sp <- compute_network(expr, density = 0.1)
  expect_equal(sp$store_density, 0.1)
  expect_equal(sp$store_threshold, sp$threshold)

  species2 <- compute_network(expr, density = 0.1, store_density = 0.2)
  expect_equal(species2$store_density, 0.2)
  expect_true(species2$store_threshold < species2$threshold)
  expect_true(length(species2$network@x) > length(sp$network@x))

  expect_error(
    compute_network(expr, density = 0.1, store_density = 0.05),
    "store_density"
  )
  expect_error(
    compute_network(expr, density = 0.1, store_density = 1),
    "store_density"
  )
  expect_error(
    compute_network(expr, density = 0.1, sparse = FALSE, store_density = 0.2),
    "sparse = TRUE"
  )
})


test_that("as_sparse_network validates its input", {
  set.seed(42)
  expr <- matrix(rnorm(300), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("g", sprintf("%02d", 1:30))
  dense <- compute_network(expr, density = 0.1, sparse = FALSE)

  sp <- rcomplex:::as_sparse_network(dense, 0.1)
  expect_s4_class(sp$network, "dgCMatrix")
  expect_equal(sp$store_threshold, dense$threshold)

  expect_error(rcomplex:::as_sparse_network(sp, 0.2), "already sparse")
  expect_error(rcomplex:::as_sparse_network(dense, 0.05), "store_density")
  expect_error(rcomplex:::as_sparse_network(dense, 0), "store_density")
  expect_error(rcomplex:::as_sparse_network(dense, 1), "store_density")
  expect_error(
    rcomplex:::as_sparse_network(list(network = dense$network)),
    "network object"
  )

  # hand-built net without params$density: converts at any store_density
  hand <- list(network = dense$network, threshold = dense$threshold)
  hs <- rcomplex:::as_sparse_network(hand, 0.2)
  expect_s4_class(hs$network, "dgCMatrix")
  expect_equal(hs$store_density, 0.2)

  # ... unless the store threshold lands ABOVE the analysis threshold:
  # that store would drop analysis edges, so it must fail here (naming
  # store_density), not one call later in .net_check()
  expect_error(rcomplex:::as_sparse_network(hand, 0.01), "store_density 0.01")
})
