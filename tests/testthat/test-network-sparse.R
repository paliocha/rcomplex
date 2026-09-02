# Tests for the sparse network object (P2):
# compute_network(sparse = TRUE) as the default, as_sparse_network(), and
# dense-vs-sparse equality for every R-side consumer of a network object.
# Fixtures: make_cmp_nets(), make_graded_nets(), sparse_net() from
# helper-reference.R; make_clique_fixture() from helper-clique-fixtures.R.


# Local fixture: 30-gene network with three clean modules
make_module_net <- function() {
  n <- 30
  mat <- matrix(0, n, n)
  rownames(mat) <- colnames(mat) <- paste0("g", sprintf("%02d", 1:n))
  mat[1:10, 1:10] <- 0.9
  mat[11:20, 11:20] <- 0.9
  mat[21:30, 21:30] <- 0.9
  mat[2, 15] <- mat[15, 2] <- 0.6
  diag(mat) <- 1
  list(network = mat, threshold = 0.5)
}

# Undirected edge data frame of a graph, orientation- and order-normalized
edge_df <- function(g) {
  e <- igraph::as_data_frame(g)
  out <- data.frame(from = pmin(e$from, e$to), to = pmax(e$from, e$to),
                    weight = e$weight, stringsAsFactors = FALSE)
  out <- out[order(out$from, out$to), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Everything but the graph (igraph objects differ in edge order between
# dense and sparse construction; edge content is compared via edge_df())
no_graph <- function(m) m[setdiff(names(m), "graph")]


# ---- (a) compute_network(sparse = TRUE) and as_sparse_network() ----

test_that("compute_network(sparse = TRUE) is the default and matches as_sparse_network", {
  set.seed(42)
  expr <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("g", sprintf("%02d", 1:50))

  dense <- compute_network(expr, density = 0.03, sparse = FALSE)
  sp_default <- compute_network(expr, density = 0.03)
  sp_explicit <- compute_network(expr, density = 0.03, sparse = TRUE)

  expect_s4_class(sp_default$network, "dgCMatrix")
  expect_equal(sp_default, sp_explicit)
  expect_equal(sp_default, as_sparse_network(dense, 0.05))

  expect_named(sp_default, c("network", "threshold", "n_genes", "n_removed",
                             "params", "store_density", "store_threshold"))
  # dense object is exactly today's: no store_* fields
  expect_named(dense, c("network", "threshold", "n_genes", "n_removed",
                        "params"))
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
})


test_that("store_density defaults to max(density, 0.05) and is validated", {
  set.seed(42)
  expr <- matrix(rnorm(300), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("g", sprintf("%02d", 1:30))

  sp <- compute_network(expr, density = 0.1)
  expect_equal(sp$store_density, 0.1)
  expect_equal(sp$store_threshold, sp$threshold)

  sp2 <- compute_network(expr, density = 0.1, store_density = 0.2)
  expect_equal(sp2$store_density, 0.2)
  expect_true(sp2$store_threshold < sp2$threshold)
  expect_true(length(sp2$network@x) > length(sp$network@x))

  expect_error(compute_network(expr, density = 0.1, store_density = 0.05),
               "store_density")
  expect_error(compute_network(expr, density = 0.1, store_density = 1),
               "store_density")
  expect_error(
    compute_network(expr, density = 0.1, sparse = FALSE, store_density = 0.2),
    "sparse = TRUE")
})


test_that("as_sparse_network validates its input", {
  set.seed(42)
  expr <- matrix(rnorm(300), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("g", sprintf("%02d", 1:30))
  dense <- compute_network(expr, density = 0.1, sparse = FALSE)

  sp <- as_sparse_network(dense, 0.1)
  expect_s4_class(sp$network, "dgCMatrix")
  expect_equal(sp$store_threshold, dense$threshold)

  expect_error(as_sparse_network(sp, 0.2), "already sparse")
  expect_error(as_sparse_network(dense, 0.05), "store_density")
  expect_error(as_sparse_network(dense, 0), "store_density")
  expect_error(as_sparse_network(dense, 1), "store_density")
  expect_error(as_sparse_network(list(network = dense$network)),
               "network object")

  # hand-built net without params$density: converts at any store_density
  hand <- list(network = dense$network, threshold = dense$threshold)
  hs <- as_sparse_network(hand, 0.2)
  expect_s4_class(hs$network, "dgCMatrix")
  expect_equal(hs$store_density, 0.2)

  # ... unless the store threshold lands ABOVE the analysis threshold:
  # that store would drop analysis edges, so it must fail here (naming
  # store_density), not one call later in .net_check()
  expect_error(as_sparse_network(hand, 0.01), "store_density 0.01")
})


test_that(".net_check rejects sparse networks with stored diagonal entries", {
  td <- make_cmp_nets()
  net_s <- sparse_net(td$net1)
  net_s$network[1, 1] <- td$net1$threshold + 1
  expect_error(compare_neighborhoods(net_s, sparse_net(td$net2), td$ortho),
               "diagonal")
})


# ---- (b) dense vs sparse equality for every consumer ----

test_that("compare_neighborhoods and summarize_comparison are identical dense vs sparse", {
  td <- make_cmp_nets()
  net1_s <- as_sparse_network(td$net1, td$net1$params$density)
  net2_s <- as_sparse_network(td$net2, td$net2$params$density)

  cmp_d <- compare_neighborhoods(td$net1, td$net2, td$ortho)
  cmp_s <- compare_neighborhoods(net1_s, net2_s, td$ortho)
  expect_equal(cmp_s, cmp_d)

  set.seed(3)
  sum_d <- summarize_comparison(cmp_d)
  set.seed(3)
  sum_s <- summarize_comparison(cmp_s)
  expect_equal(sum_s, sum_d)
})


test_that("permutation_hog_test is identical dense vs sparse (seeded)", {
  td <- make_cmp_nets()
  net1_s <- as_sparse_network(td$net1, td$net1$params$density)
  net2_s <- as_sparse_network(td$net2, td$net2$params$density)
  cmp <- compare_neighborhoods(td$net1, td$net2, td$ortho)

  set.seed(7)
  r_d <- permutation_hog_test(td$net1, td$net2, cmp,
                              max_permutations = 200L)
  set.seed(7)
  r_s <- permutation_hog_test(net1_s, net2_s, cmp,
                              max_permutations = 200L)
  expect_equal(r_s, r_d)
})


test_that("find_coexpressologs is identical dense vs sparse", {
  td <- make_cmp_nets()
  nets_d <- list(SP_A = td$net1, SP_B = td$net2)
  nets_s <- list(SP_A = as_sparse_network(td$net1, 0.1),
                 SP_B = as_sparse_network(td$net2, 0.1))

  set.seed(5)
  e_d <- find_coexpressologs(nets_d, td$ortho)
  set.seed(5)
  e_s <- find_coexpressologs(nets_s, td$ortho)
  expect_equal(e_s, e_d)
})


test_that("density_sweep is identical dense vs sparse for multipliers within the store", {
  td <- make_cmp_nets()
  nets_d <- list(SP_A = td$net1, SP_B = td$net2)
  nets_s <- list(SP_A = as_sparse_network(td$net1, 0.3),
                 SP_B = as_sparse_network(td$net2, 0.3))

  sw_d <- suppressMessages(density_sweep(
    nets_d, td$ortho, multipliers = c(0.9, 1, 1.1),
    method = "analytical", pi0_method = "storey"))
  sw_s <- suppressMessages(density_sweep(
    nets_s, td$ortho, multipliers = c(0.9, 1, 1.1),
    method = "analytical", pi0_method = "storey"))
  expect_equal(sw_s, sw_d)
})


test_that("get_coexpressed_hogs is identical dense vs sparse", {
  td <- make_graded_nets()
  nets_d <- list(SP_A = td$net1, SP_B = td$net2)
  nets_s <- list(SP_A = sparse_net(td$net1), SP_B = sparse_net(td$net2))
  ortho <- td$ortho
  ortho$Species1 <- as.character(ortho$Species1)

  r_d <- get_coexpressed_hogs("HOG1", nets_d, ortho, min_species = 2L)
  r_s <- get_coexpressed_hogs("HOG1", nets_s, ortho, min_species = 2L)
  expect_true(nrow(r_d) > 0)

  # asymmetric HOG1 copies (see make_graded_nets()): gene 1 -> gene 16 and
  # gene 2 -> gene 17 at tier 10, cross entries at the sub-store tier 4.
  # The per-copy weight must average only edges of the copy's own
  # neighbourhood (values >= threshold, hence always stored), so HOG6 gets
  # mean_weight 10 regardless of storage; a union-of-copies mean would read
  # the tier-4 cross entry from the dense matrix but an implicit 0 from the
  # dgCMatrix.
  expect_true("HOG6" %in% r_d$partner_hog)
  expect_equal(r_d$mean_weight[r_d$partner_hog == "HOG6"], 10)
  expect_equal(r_s$mean_weight[r_s$partner_hog == "HOG6"], 10)
  expect_equal(r_s, r_d)
})


test_that("detect_modules (single) is identical dense vs sparse", {
  net_d <- make_module_net()
  net_s <- sparse_net(net_d)

  m_d <- detect_modules(net_d, method = "leiden",
                        objective_function = "modularity", seed = 42)
  m_s <- detect_modules(net_s, method = "leiden",
                        objective_function = "modularity", seed = 42)
  expect_equal(no_graph(m_s), no_graph(m_d))
  expect_equal(edge_df(m_s$graph), edge_df(m_d$graph))
  expect_equal(m_d$n_modules, 3L)

  i_d <- detect_modules(net_d, method = "infomap", seed = 42)
  i_s <- detect_modules(net_s, method = "infomap", seed = 42)
  expect_equal(no_graph(i_s), no_graph(i_d))
})


test_that("detect_modules (consensus) is identical dense vs sparse", {
  net_d <- make_module_net()
  net_s <- sparse_net(net_d)

  m_d <- detect_modules(net_d, resolution = c(0.5, 1, 2),
                        objective_function = "modularity", seed = 42,
                        test_k1 = FALSE)
  m_s <- detect_modules(net_s, resolution = c(0.5, 1, 2),
                        objective_function = "modularity", seed = 42,
                        test_k1 = FALSE)
  expect_equal(no_graph(m_s), no_graph(m_d))
  expect_equal(edge_df(m_s$graph), edge_df(m_d$graph))
})


test_that("detect_modules densifies for SBM with a warning", {
  skip_if_not_installed("sbm")
  net_s <- sparse_net(make_module_net())
  expect_warning(
    m <- detect_modules(net_s, method = "sbm", seed = 42),
    "densifying")
  expect_true(m$n_modules >= 1L)
})


test_that("identify_module_hubs is identical dense vs sparse", {
  net_d <- make_module_net()
  net_s <- sparse_net(net_d)

  m_d <- detect_modules(net_d, method = "leiden",
                        objective_function = "modularity", seed = 42)
  m_s <- detect_modules(net_s, method = "leiden",
                        objective_function = "modularity", seed = 42)

  set.seed(1)
  h_d <- identify_module_hubs(m_d, net_d)
  set.seed(1)
  h_s <- identify_module_hubs(m_s, net_s)
  expect_equal(h_s, h_d)
})


test_that("clique_persistence is identical dense vs sparse", {
  set.seed(99)
  setup <- make_clique_fixture()
  skip_if(nrow(setup$cliques) == 0, "No baseline cliques found")
  nets_s <- lapply(setup$networks, sparse_net)

  p_d <- clique_persistence(setup$cliques, setup$target_species,
                            setup$networks, setup$edges)
  p_s <- clique_persistence(setup$cliques, setup$target_species,
                            nets_s, setup$edges)
  expect_equal(p_s, p_d)
})


test_that("clique_threshold_sweep is identical dense vs sparse", {
  set.seed(99)
  setup <- make_clique_fixture()
  skip_if(nrow(setup$cliques) == 0, "No baseline cliques found")
  nets_s <- lapply(setup$networks, sparse_net)

  sw_d <- suppressMessages(clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(1.5, 2)))
  sw_s <- suppressMessages(clique_threshold_sweep(
    setup$cliques, setup$target_species, nets_s,
    setup$orthologs, multipliers = c(1.5, 2)))
  expect_equal(sw_s, sw_d)

  # the sparse store travels through the multiplier rescale (modifyList),
  # so the guard stays armed inside the sweep
  expect_true(!is.null(nets_s$SP_A$store_threshold))
})


test_that("clique_intensity_test is identical dense vs sparse (seeded)", {
  set.seed(99)
  setup <- make_clique_fixture()
  skip_if(nrow(setup$cliques) == 0, "No baseline cliques found")
  nets_s <- lapply(setup$networks, sparse_net)

  it_d <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks, setup$orthologs,
    n_perm = 3L, seed = 11L, edges = setup$edges)
  it_s <- clique_intensity_test(
    setup$cliques, setup$target_species, nets_s, setup$orthologs,
    n_perm = 3L, seed = 11L, edges = setup$edges)
  expect_equal(it_s, it_d)
})


# ---- (c) clique_perturbation_test on sparse networks ----

test_that("clique_perturbation_test works on sparse networks", {
  set.seed(99)
  setup <- make_clique_fixture()
  skip_if(nrow(setup$cliques) == 0, "No baseline cliques found")
  nets_s <- lapply(setup$networks, sparse_net)

  # zero noise: every clique survives every bootstrap
  r0 <- clique_perturbation_test(
    setup$cliques, setup$target_species, nets_s, setup$orthologs,
    n_boot = 3L, noise_sd = 0, seed = 42L)
  expect_true(all(r0$survival_rate == 1.0))
  expect_true(all(r0$mean_jaccard == 1.0))

  # seeded runs are reproducible
  r1 <- suppressWarnings(clique_perturbation_test(
    setup$cliques, setup$target_species, nets_s, setup$orthologs,
    n_boot = 3L, noise_sd = 0.5, seed = 7L))
  r2 <- suppressWarnings(clique_perturbation_test(
    setup$cliques, setup$target_species, nets_s, setup$orthologs,
    n_boot = 3L, noise_sd = 0.5, seed = 7L))
  expect_equal(r1, r2)
})


# ---- (d) store guard through density_sweep ----

test_that("density_sweep errors when a multiplier falls below the store", {
  td <- make_cmp_nets()
  # store_density = density: no headroom below multiplier 1
  nets_s <- list(SP_A = as_sparse_network(td$net1, td$net1$params$density),
                 SP_B = as_sparse_network(td$net2, td$net2$params$density))
  expect_error(
    suppressMessages(density_sweep(nets_s, td$ortho, multipliers = 0.1,
                                   method = "analytical")),
    "store_density")
})


# ---- (e) rcomplex container with sparse networks ----

test_that("rcomplex accepts sparse networks and print shows storage", {
  td <- make_cmp_nets()
  nets_s <- list(SP_A = as_sparse_network(td$net1, 0.1),
                 SP_B = as_sparse_network(td$net2, 0.1))
  ortho <- td$ortho
  ortho$hog <- as.character(ortho$hog)

  rcx <- rcomplex(
    species = c("SP_A", "SP_B"),
    traits = c(SP_A = "annual", SP_B = "perennial"),
    networks = nets_s,
    orthologs = ortho
  )
  expect_s3_class(rcx, "rcomplex")

  out <- capture.output(print(rcx))
  expect_true(any(grepl("dgCMatrix", out)))
  expect_true(any(grepl("stored", out)))
  expect_true(any(grepl("store_density", out)))

  s <- summary(rcx)
  expect_true(any(grepl("dgCMatrix", s$network_storage)))

  # other Matrix classes are rejected (dsCMatrix from Matrix::Matrix())
  bad <- nets_s
  bad$SP_A$network <- Matrix::Matrix(td$net1$network, sparse = TRUE)
  expect_error(
    rcomplex(species = c("SP_A", "SP_B"),
             traits = c(SP_A = "annual", SP_B = "perennial"),
             networks = bad, orthologs = ortho),
    "dgCMatrix")
})
