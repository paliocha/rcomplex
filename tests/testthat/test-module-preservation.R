# Tests for module_preservation() and classify_preservation()

# Fixture: two species sharing module structure. Gene loadings on each
# module's latent factor are heavy-tailed and SHARED between species, so hub
# identity is conserved and cor.degree has signal. A fixture where every gene
# in a module is exchangeable (one factor, iid noise) correctly yields
# cor.degree ~ 0 even for preserved modules, and would look like a bug.
pres_expr <- function(seed, n, prefix, loadings, per, n_samp = 40) {
  set.seed(seed)
  n_mod <- length(loadings)
  e <- matrix(stats::rnorm(n * n_samp), n, n_samp)
  for (k in seq_len(n_mod)) {
    f <- stats::rnorm(n_samp)
    rows <- ((k - 1) * per + 1):(k * per)
    for (j in seq_along(rows)) {
      lam <- loadings[[k]][j]
      e[rows[j], ] <- lam * f +
        stats::rnorm(n_samp, sd = sqrt(max(1e-6, 1 - lam^2)))
    }
  }
  rownames(e) <- paste0(prefix, sprintf("%04d", seq_len(n)))
  e
}

pres_fixture <- function() {
  n_mod <- 4L
  per <- 40L
  set.seed(77)
  loadings <- lapply(seq_len(n_mod), function(k) {
    l <- stats::rlnorm(per, 0, 0.9)
    l / max(l)
  })
  eA <- pres_expr(31, 300, "A", loadings, per)
  eB <- pres_expr(32, 500, "B", loadings, per)
  np <- n_mod * per

  # Map background genes too, not just module genes. The mappable universe is
  # what the permutation draws from, so an ortholog table covering only module
  # genes makes the null "genes from the other modules" -- and since every
  # module is dense, avg.weight then has almost no contrast. Real ortholog
  # tables include background, and so must the fixture.
  n_map <- 300L

  list(
    n_mod = n_mod, per = per,
    netA = compute_network(eA, density = 0.03, sparse = FALSE),
    netB = compute_network(eB, density = 0.03, sparse = FALSE),
    netB_sparse = compute_network(eB,
      density = 0.03, sparse = TRUE,
      store_density = 0.03
    ),
    ortho = data.frame(
      Species1 = paste0("A", sprintf("%04d", seq_len(n_map))),
      Species2 = paste0("B", sprintf("%04d", seq_len(n_map))),
      hog = paste0("H", seq_len(n_map)),
      stringsAsFactors = FALSE
    ),
    mods = lapply(
      seq_len(n_mod),
      function(k) ((k - 1) * per + 1):(k * per)
    )
  )
}

# Module labels matching the simulated block structure, in detect_modules()
# shape, so the tests do not depend on Leiden's partition.
true_modules <- function(net, mods) {
  genes <- rownames(net$network)
  membership <- stats::setNames(rep(NA_integer_, length(genes)), genes)
  for (k in seq_along(mods)) membership[mods[[k]]] <- k
  membership <- membership[!is.na(membership)]
  list(
    modules = membership,
    module_genes = split(names(membership), membership),
    n_modules = length(mods)
  )
}


# ---- C++ kernel against the pure-R reference ----

test_that("kernel statistics match the R reference implementation", {
  fx <- pres_fixture()
  net <- fx$netA
  n <- nrow(net$network)
  keep <- as.integer(seq_len(n) - 1L)
  mm <- lapply(fx$mods, function(z) as.integer(z - 1L))

  got <- module_gene_stats_dense_cpp(
    net$network, net$threshold, keep, mm,
    FALSE
  )
  adj <- reference_adjacency(net, rownames(net$network))

  for (k in seq_along(fx$mods)) {
    idx <- fx$mods[[k]]
    want <- reference_module_stats(adj, idx)
    expect_equal(got$kIM[idx], unname(want$kIM), tolerance = 1e-10)
    expect_equal(got$CC[idx], unname(want$CC), tolerance = 1e-10)
    expect_equal(got$MAR[idx], unname(want$MAR), tolerance = 1e-10)
  }
})

test_that("observed preservation statistics match the R reference", {
  fx <- pres_fixture()
  net <- fx$netA
  n <- nrow(net$network)
  keep <- as.integer(seq_len(n) - 1L)
  mm <- lapply(fx$mods, function(z) as.integer(z - 1L))
  adj <- reference_adjacency(net, rownames(net$network))
  rs <- lapply(fx$mods, function(i) reference_module_stats(adj, i))

  got <- module_preservation_dense_cpp(
    net$network, net$threshold, keep, mm,
    lapply(rs, `[[`, "kIM"), lapply(rs, `[[`, "CC"), lapply(rs, `[[`, "MAR"),
    n_perm = 0L, n_cores = 1L, binary = FALSE
  )

  for (k in seq_along(fx$mods)) {
    want <- reference_preservation_stats(adj, fx$mods[[k]], adj, fx$mods[[k]])
    expect_equal(unname(got$observed[k, ]), unname(want), tolerance = 1e-10)
  }
})

test_that("dense and sparse kernels agree exactly", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  set.seed(3)
  dense <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, seed = 11
  )
  set.seed(3)
  sparse <- module_preservation(tm, fx$netA, fx$netB_sparse, fx$ortho,
    n_perm = 100L, seed = 11
  )

  expect_equal(dense$preservation, sparse$preservation)
})

test_that("results do not depend on n_cores", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  one <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, n_cores = 1L, seed = 5
  )
  four <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, n_cores = 4L, seed = 5
  )

  expect_equal(one$preservation$Zsummary, four$preservation$Zsummary,
    tolerance = 1e-8
  )
  expect_equal(one$preservation$p.value, four$preservation$p.value)
})


# ---- calibration: the checks that catch a broken permutation scheme ----

test_that("permutation p-values are uniform under the null", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  genes_b <- rownames(fx$netB$network)

  set.seed(4)
  pvals <- c()
  for (rep in 1:8) {
    shuffled <- sample(genes_b, nrow(fx$ortho))
    ortho_rand <- fx$ortho
    ortho_rand$Species2 <- shuffled
    pr <- module_preservation(tm, fx$netA, fx$netB, ortho_rand,
      n_perm = 100L, n_cores = 2L, seed = rep
    )
    pvals <- c(
      pvals, pr$preservation$p.avg.weight,
      pr$preservation$p.cor.degree
    )
  }

  # Under a random ortholog map neither statistic should be significant more
  # often than chance.
  expect_gt(mean(pvals, na.rm = TRUE), 0.3)
  expect_lt(mean(pvals, na.rm = TRUE), 0.7)
  expect_lt(mean(pvals < 0.05, na.rm = TRUE), 0.2)
})

test_that("Zsummary centres near zero under the null", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  genes_b <- rownames(fx$netB$network)

  set.seed(6)
  zs <- c()
  for (rep in 1:8) {
    ortho_rand <- fx$ortho
    ortho_rand$Species2 <- sample(genes_b, nrow(fx$ortho))
    pr <- module_preservation(tm, fx$netA, fx$netB, ortho_rand,
      n_perm = 100L, n_cores = 2L, seed = rep
    )
    zs <- c(zs, pr$preservation$Zsummary)
  }

  expect_lt(abs(mean(zs, na.rm = TRUE)), 1.5)
})

test_that("shared module structure is detected as preserved", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 200L, n_cores = 2L, seed = 1
  )

  expect_equal(nrow(pres$preservation), fx$n_mod)
  expect_true(all(pres$preservation$p.value < 0.05))
  expect_true(all(pres$preservation$Zsummary > 2))
})

test_that("an unrelated test network gives no preservation", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  set.seed(99)
  e <- matrix(stats::rnorm(500 * 40), 500, 40)
  rownames(e) <- paste0("B", sprintf("%04d", seq_len(500)))
  net_rand <- compute_network(e, density = 0.03, sparse = FALSE)

  pres <- module_preservation(tm, fx$netA, net_rand, fx$ortho,
    n_perm = 200L, n_cores = 2L, seed = 1
  )
  cls <- classify_preservation(pres)

  expect_true(all(cls$classification == "diverged"))
})


# ---- structure, options, validation ----

test_that("module_preservation returns the documented structure", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, seed = 1
  )

  expect_named(pres, c(
    "preservation", "observed", "projection", "map",
    "params"
  ))
  expect_true(all(c(
    "module", "size", "size_mapped", "avg.weight",
    "cor.degree", "p.value", "q.value", "Zsummary",
    "medianRank"
  ) %in% names(pres$preservation)))
  # The diagnostics are reported but take no part in the call.
  expect_true(all(c("meanMAR", "meanClusterCoeff") %in% names(pres$observed)))
})

test_that("binary mode makes avg.weight the module edge density", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, binary = TRUE, seed = 1
  )

  expect_true(all(pres$preservation$avg.weight >= 0))
  expect_true(all(pres$preservation$avg.weight <= 1))
  expect_equal(pres$params$scale, 1)
})

test_that("min_module_size excludes small modules", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  big <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, min_module_size = 10L, seed = 1
  )
  expect_equal(nrow(big$preservation), fx$n_mod)

  expect_error(
    module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = 50L, min_module_size = 500L, seed = 1
    ),
    "no module has at least min_module_size"
  )
})

test_that("module_preservation validates its inputs", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  expect_error(
    module_preservation(list(a = 1), fx$netA, fx$netB, fx$ortho),
    "must be output from detect_modules"
  )
  expect_error(
    module_preservation(tm, fx$netA, fx$netB, fx$ortho, n_perm = 0L),
    "n_perm must be >= 1"
  )
  expect_error(
    module_preservation(tm, fx$netA, fx$netB, fx$ortho, min_module_size = 2L),
    "min_module_size must be >= 3"
  )
  expect_error(
    module_preservation(tm, fx$netA, fx$netB, orthologs = NULL),
    "supply either 'orthologs' or a pre-built 'map'"
  )
})


# ---- projection ----

test_that("resolved map entries win over unresolved ones", {
  map <- data.frame(
    gene1 = c("a1", "a2", "a3"),
    gene2 = c("B1", "B1", "B1"),
    module = c("1", "2", "2"),
    source = c("clique", "unresolved", "unresolved"),
    stringsAsFactors = FALSE
  )
  proj <- .pres_project(map)

  expect_equal(nrow(proj), 1L)
  expect_equal(proj$module, "1")
  expect_equal(proj$source, "clique")
})

test_that("unresolved genes take the modal label and drop on ties", {
  majority <- data.frame(
    gene1 = c("a1", "a2", "a3"),
    gene2 = "B1",
    module = c("2", "2", "3"),
    source = "unresolved",
    stringsAsFactors = FALSE
  )
  expect_equal(.pres_project(majority)$module, "2")

  tied <- data.frame(
    gene1 = c("a1", "a2"),
    gene2 = "B1",
    module = c("2", "3"),
    source = "unresolved",
    stringsAsFactors = FALSE
  )
  expect_null(.pres_project(tied))
})


# ---- classification ----

test_that("classify_preservation applies the documented criteria", {
  pres <- list(preservation = data.frame(
    module = c("1", "2", "3"),
    size = c(30L, 30L, 30L),
    size_mapped = c(20L, 20L, 20L),
    Zsummary = c(15, 5, 0.2),
    q.value = c(0.001, 0.001, 0.9),
    stringsAsFactors = FALSE
  ))

  cls <- classify_preservation(pres, alpha = 0.05, z_conserved = 10)
  expect_equal(cls$classification, c("conserved", "moderate", "diverged"))
})

test_that("classify_preservation records species and pair labels", {
  pres <- list(preservation = data.frame(
    module = "1", size = 30L, size_mapped = 20L,
    Zsummary = 15, q.value = 0.001, stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(pres, species = "SP_A", pair_name = "A.B")

  expect_equal(cls$species, "SP_A")
  expect_equal(cls$pair_name, "A.B")
})

test_that("classify_preservation validates its input", {
  expect_error(
    classify_preservation(list(a = 1)),
    "must be output from module_preservation"
  )
})
