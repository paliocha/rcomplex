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


# ---- sensitivity: the circularity guard ----

# A HOG that is multi-copy on the reference side, spanning two modules: each
# listed species-2 gene has one partner in module 1 and one in module 2, so the
# naive map ties and .pres_project() drops it, while a clique resolves it to
# module 1. Both maps still cover the same species-2 genes -- only the copy
# choice differs -- which is exactly the case sensitivity exists to measure.
ambiguous_fixture <- function(fx, n_amb = 10L) {
  amb <- seq_len(n_amb)
  list(
    ortho = rbind(fx$ortho, data.frame(
      Species1 = paste0("A", sprintf("%04d", fx$per + amb)),
      Species2 = paste0("B", sprintf("%04d", amb)),
      hog = paste0("H", amb),
      stringsAsFactors = FALSE
    )),
    cliques = data.frame(
      hog = paste0("H", amb),
      A = paste0("A", sprintf("%04d", amb)),
      B = paste0("B", sprintf("%04d", amb)),
      n_species = 2L, mean_q = 0.01,
      stringsAsFactors = FALSE
    )
  )
}

test_that("sensitivity detects a copy choice that changes the result", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)

  pres <- module_preservation(tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques, sp_ref = "A", sp_test = "B",
    n_perm = 100L, sensitivity = TRUE, seed = 1
  )

  expect_true("sensitivity" %in% names(pres))
  expect_named(pres$sensitivity, c(
    "module", "Zsummary", "Zsummary_naive",
    "q.value", "q.value_naive", "Zsummary_delta"
  ))

  # Resolution may only change which copy carries a label, never which genes
  # are mappable.
  expect_true(attr(pres$sensitivity, "same_gene_set"))

  # The clique rescues genes the naive majority vote drops on a tie, so at
  # least one module must actually move.
  expect_true(any(pres$sensitivity$Zsummary_delta != 0))
})

test_that("sensitivity warns when the two maps cover different genes", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  naive <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  # Drop a mappable gene: a resolution layer that filtered like this would be
  # selecting the tested genes on the statistic being tested.
  trimmed <- naive[naive$gene2 != naive$gene2[1], , drop = FALSE]

  expect_warning(
    pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      map = trimmed, n_perm = 100L, sensitivity = TRUE, seed = 1
    ),
    "cover different"
  )
  expect_false(attr(pres$sensitivity, "same_gene_set"))
})

test_that("sensitivity is skipped when there is nothing to resolve", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  expect_warning(
    pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = 50L, sensitivity = TRUE, seed = 1
    ),
    "nothing to compare"
  )
  expect_false("sensitivity" %in% names(pres))
})

test_that("sensitivity is absent unless requested", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, seed = 1
  )
  expect_false("sensitivity" %in% names(pres))
})

test_that("sensitivity warns and is skipped without orthologs", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )

  expect_warning(
    pres <- module_preservation(tm, fx$netA, fx$netB,
      orthologs = NULL, map = map, n_perm = 50L, sensitivity = TRUE, seed = 1
    ),
    "needs 'orthologs'"
  )
  expect_false("sensitivity" %in% names(pres))
})


# ---- module_correspondence ----

test_that("module_correspondence matches the modules that correspond", {
  fx <- pres_fixture()
  tm_a <- true_modules(fx$netA, fx$mods)
  tm_b <- true_modules(fx$netB, fx$mods)
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )

  corr <- module_correspondence(tm_a, tm_b, map)
  p <- corr$pairs

  expect_true(all(c(
    "module_sp1", "module_sp2", "overlap", "jaccard",
    "p.value", "q.value"
  ) %in% names(p)))
  expect_equal(nrow(p), tm_a$n_modules * tm_b$n_modules)

  # The fixture maps module k of species A onto module k of species B, so the
  # diagonal must be the significant part of the table.
  diag_rows <- p$module_sp1 == p$module_sp2
  expect_true(all(p$q.value[diag_rows] < 0.05))
  expect_gt(min(p$jaccard[diag_rows]), max(p$jaccard[!diag_rows]))
})

test_that("module_correspondence validates its inputs", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )

  expect_error(module_correspondence(list(a = 1), tm, map),
    "must be output from detect_modules"
  )
  expect_error(module_correspondence(tm, tm, data.frame(x = 1)),
    "must be a data frame from resolve_ortholog_map"
  )
})


# ---- preservation_paired ----

test_that("preservation_paired runs both directions per contrast", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    n_perm = 100L, seed = 1
  )

  expect_named(res, c("classification", "summary", "raw"))
  expect_setequal(names(res$raw), c("A.B", "B.A"))
  expect_setequal(unique(res$classification$reference), c("A", "B"))

  # tag_permutation() reads exactly these columns.
  expect_true(all(c("pair_name", "module", "species", "classification") %in%
    names(res$classification)))
})

test_that("preservation_paired tags trait groups", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    group = c(A = "annual", B = "perennial"), n_perm = 100L, seed = 1
  )

  expect_true("group" %in% names(res$classification))
  expect_true(all(res$classification$group %in%
    c("conserved", "annual", "perennial")))
})

test_that("preservation_paired validates its inputs", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods))
  nets <- list(A = fx$netA)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  expect_error(
    preservation_paired(mods, nets, fx$ortho, pairs, n_perm = 10L),
    "modules and networks must both cover"
  )
  expect_error(
    preservation_paired(mods, nets, fx$ortho, data.frame(x = 1)),
    "must have columns"
  )
})


# ---- an uncomputable statistic must not be scored as significant ----

test_that("a constant reference degree gives NA, not a floor p-value", {
  fx <- pres_fixture()
  net <- fx$netA
  keep <- as.integer(seq_len(nrow(net$network)) - 1L)
  mm <- lapply(fx$mods, function(z) as.integer(z - 1L))

  # cor.degree is undefined when the reference vector is constant. The
  # exceedance counter can never fire for an NA observed value, so a naive
  # (0 + 1) / (n + 1) would report the most significant p-value attainable
  # for a statistic that could not be computed.
  flat <- lapply(fx$mods, function(i) rep(1, length(i)))
  got <- module_preservation_dense_cpp(
    net$network, net$threshold, keep, mm, flat, flat, flat,
    n_perm = 100L, n_cores = 1L, binary = FALSE
  )

  expect_true(all(is.na(got$observed[, 4])))
  expect_true(all(is.na(got$p_value[, 4])))
  # All three outputs must agree, or a consumer could recompute the floor
  # p-value from the counts.
  expect_true(all(is.na(got$n_perm_used[, 4])))
  expect_true(all(is.na(got$n_exceed[, 4])))
  # The density statistic is unaffected and still computed.
  expect_false(any(is.na(got$p_value[, 1])))
  expect_false(any(is.na(got$n_perm_used[, 1])))
})

test_that("an NA combined p-value is reported as untested, not diverged", {
  pres <- list(preservation = data.frame(
    module = c("1", "2"),
    size = c(30L, 30L), size_mapped = c(20L, 20L),
    Zsummary = c(NA_real_, 15),
    q.value = c(NA_real_, 0.001),
    stringsAsFactors = FALSE
  ))
  # Nothing was measured for module 1, so calling it diverged would assert
  # something the data does not support.
  expect_warning(cls <- classify_preservation(pres), "could not be tested")
  expect_equal(cls$classification, c("untested", "conserved"))
})

test_that("q-value correction passes NA through", {
  p <- c(0.001, NA, 0.5, 0.9)
  q <- rcomplex:::.pres_qvalues(p, n_perm = 1000L, method = "bh")

  expect_true(is.na(q[2]))
  expect_false(any(is.na(q[-2])))
})


test_that("preservation_paired rejects a repeated or self species pair", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)

  # Both directions of each contrast are run, so (A, B) and (B, A) would
  # collide on the "<reference>.<test>" key and silently overwrite each other.
  expect_error(
    preservation_paired(mods, nets, fx$ortho,
      data.frame(sp1 = c("A", "B"), sp2 = c("B", "A"),
                 stringsAsFactors = FALSE),
      n_perm = 10L
    ),
    "same species pair more than once"
  )
  expect_error(
    preservation_paired(mods, nets, fx$ortho,
      data.frame(sp1 = "A", sp2 = "A", stringsAsFactors = FALSE),
      n_perm = 10L
    ),
    "not compare a species with itself"
  )
})
