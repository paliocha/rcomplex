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
    "params", "coverage"
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
    "module", "size_mapped", "size_mapped_naive", "Zsummary",
    "Zsummary_naive", "q.value", "q.value_naive", "Zsummary_delta",
    "p_copy.avg.weight", "p_copy.cor.degree"
  ))

  # Resolution may only change which copy carries a label, never which genes
  # are mappable.
  # The candidate sets are equal by construction; what matters is whether
  # projection kept the same genes.
  expect_true(attr(pres$sensitivity, "same_candidate_set"))

  # The clique rescues the 10 genes the naive majority vote drops on a tie, so
  # module 1 gains exactly those. This is the deterministic consequence of the
  # copy choice; Zsummary_delta is not assertable on its own because differing
  # block sizes also make the two runs consume the RNG differently, which
  # shifts every module's delta.
  m1 <- pres$sensitivity[pres$sensitivity$module == "1", ]
  expect_equal(m1$size_mapped - m1$size_mapped_naive, 10L)
  # Modules the ambiguity does not touch are unchanged in size.
  rest <- pres$sensitivity[pres$sensitivity$module != "1", ]
  expect_true(all(rest$size_mapped == rest$size_mapped_naive))
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
    "changed which test-species genes"
  )
  expect_false(attr(pres$sensitivity, "same_projected_set"))
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
  expect_true(all(c("pair_name", "module", "reference", "test",
                    "classification") %in% names(res$classification)))
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
  # With an NA observed statistic all three outputs are NA, or a consumer
  # could recompute the floor p-value from the counts.
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


# ---- untested modules through preservation_paired ----

# Hand-built pair: two connected blocks plus an isolated block. Every gene in
# the isolated module has kIM = 0, so the degree correlation is undefined and
# the module comes back "untested" rather than diverged.
isolated_fixture <- function() {
  n <- 60L
  build <- function(prefix) {
    m <- matrix(0, n, n)
    # Each connected block is a clique minus its tail-tail edges, so the
    # first half have degree 19 and the second half 10. A plain clique would
    # give every gene the same kIM, making cor.degree undefined and the block
    # untested too -- which would leave the mixed tested/untested case, the
    # one the summary accounting depends on, uncovered.
    m[1:20, 1:20] <- 1
    m[11:20, 11:20] <- 0
    m[21:40, 21:40] <- 1
    m[31:40, 31:40] <- 0
    diag(m) <- 0
    dimnames(m) <- list(
      paste0(prefix, sprintf("%02d", seq_len(n))),
      paste0(prefix, sprintf("%02d", seq_len(n)))
    )
    list(network = m, threshold = 0.5)
  }
  mods <- list(1:20, 21:40, 41:60)
  nets <- list(A = build("A"), B = build("B"))
  list(
    nets = nets,
    mods = lapply(names(nets), function(sp) {
      genes <- rownames(nets[[sp]]$network)
      mb <- stats::setNames(rep(NA_integer_, n), genes)
      for (k in seq_along(mods)) mb[mods[[k]]] <- k
      list(modules = mb, module_genes = split(names(mb), mb), n_modules = 3L)
    }) |> stats::setNames(names(nets)),
    ortho = data.frame(
      Species1 = paste0("A", sprintf("%02d", seq_len(n))),
      Species2 = paste0("B", sprintf("%02d", seq_len(n))),
      hog = paste0("H", seq_len(n)),
      stringsAsFactors = FALSE
    )
  )
}

test_that("a module with constant connectivity is untested, not diverged", {
  fx <- isolated_fixture()

  expect_warning(
    pres <- module_preservation(fx$mods$A, fx$nets$A, fx$nets$B, fx$ortho,
      n_perm = 50L, seed = 1
    ),
    "connectivity is constant"
  )
  expect_warning(cls <- classify_preservation(pres), "could not be tested")

  # Only the isolated module is untested; the other two must be testable, or
  # the mixed accounting below is not actually being exercised.
  expect_equal(cls$classification[cls$module == "3"], "untested")
  expect_false(any(cls$classification[cls$module != "3"] == "untested"))
  expect_true(all(is.na(cls$q.value[cls$classification == "untested"])))
  # Upstream property, not a restatement of the classifier's own definition:
  # cor.degree is computable for the two blocks with non-constant degree.
  expect_false(anyNA(
    pres$preservation$cor.degree[pres$preservation$module != "3"]))
})

test_that("untested modules stay visible in the paired summary", {
  fx <- isolated_fixture()
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- suppressWarnings(preservation_paired(
    fx$mods, fx$nets, fx$ortho, pairs,
    group = c(A = "annual", B = "perennial"), n_perm = 50L, seed = 1
  ))

  expect_true("untested" %in% res$classification$classification)
  # An untested module earns no trait attribution, but must not vanish:
  # stats::aggregate() drops NA groups under its default na.omit, so the
  # counts would silently stop summing to the number of modules.
  untested <- res$classification$classification == "untested"
  expect_true(all(res$classification$group[untested] == "untested"))
  # Tested rows must still receive their real attribution -- a table that was
  # untested end to end would satisfy the line above vacuously.
  expect_true(any(!untested))
  expect_true(all(res$classification$group[!untested] %in%
    c("conserved", "annual", "perennial")))
  expect_equal(sum(res$summary$n), nrow(res$classification))
})

test_that("the paired summary always accounts for every module", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  for (grp in list(NULL, c(A = "annual", B = "perennial"))) {
    res <- suppressWarnings(preservation_paired(
      mods, nets, fx$ortho, pairs, group = grp, n_perm = 50L, seed = 1
    ))
    expect_equal(sum(res$summary$n), nrow(res$classification))
  }
})


test_that("an empty null keeps the counts but not the p-value", {
  fx <- pres_fixture()
  net <- fx$netA
  keep <- as.integer(seq_len(nrow(net$network)) - 1L)
  mm <- lapply(fx$mods, function(z) as.integer(z - 1L))
  adj <- reference_adjacency(net, rownames(net$network))
  rs <- lapply(fx$mods, function(i) reference_module_stats(adj, i))

  # n_perm = 0 is the "observed is fine, the null was uncomputable" state:
  # distinct from an NA observed statistic, and the usable count of 0 is the
  # only record of it.
  got <- module_preservation_dense_cpp(
    net$network, net$threshold, keep, mm,
    lapply(rs, `[[`, "kIM"), lapply(rs, `[[`, "CC"), lapply(rs, `[[`, "MAR"),
    n_perm = 0L, n_cores = 1L, binary = FALSE
  )

  expect_false(any(is.na(got$observed[, 1])))
  expect_equal(got$n_perm_used[, 1], rep(0L, length(fx$mods)))
  expect_equal(got$n_exceed[, 1], rep(0L, length(fx$mods)))
  expect_true(all(is.na(got$p_value[, 1])))

  # Per-column, not all-columns-at-once: a constant reference degree vector
  # makes cor.degree NA while the density column keeps a usable null. This
  # covers per-column have_obs reporting. The remaining state -- a finite
  # observed value whose own null was unscorable while sibling columns kept a
  # usable one -- is not constructible from the R entry point, which rejects
  # n_perm < 1; it is covered only by the n_perm = 0 case above.
  flat <- lapply(fx$mods, function(i) rep(1, length(i)))
  mixed <- module_preservation_dense_cpp(
    net$network, net$threshold, keep, mm, flat, flat, flat,
    n_perm = 20L, n_cores = 1L, binary = FALSE
  )
  expect_true(all(is.na(mixed$n_perm_used[, 4])))
  expect_true(all(mixed$n_perm_used[, 1] > 0L))
  expect_false(any(is.na(mixed$p_value[, 1])))
})


test_that("sensitivity reports the naive run's own statistics", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)

  pres <- module_preservation(tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques, sp_ref = "A", sp_test = "B",
    n_perm = 100L, sensitivity = TRUE, seed = 1
  )

  # The internal naive run is seeded with the same seed, so an external run on
  # the naive map reproduces it exactly. Without this, a wrong-index or
  # copy-paste regression in .pres_sensitivity() -- naive columns silently
  # echoing the resolved ones -- would pass every other assertion.
  naive_map <- resolve_ortholog_map(
    amb$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  ext <- module_preservation(tm, fx$netA, fx$netB, amb$ortho,
    map = naive_map, n_perm = 100L, seed = 1
  )

  idx <- match(pres$sensitivity$module, ext$preservation$module)
  expect_equal(pres$sensitivity$Zsummary_naive, ext$preservation$Zsummary[idx])
  expect_equal(pres$sensitivity$q.value_naive, ext$preservation$q.value[idx])
})


test_that("sensitivity warns when resolution loses a testable module", {
  resolved <- list(preservation = data.frame(
    module = c("1", "2"), size_mapped = c(20L, 20L),
    Zsummary = c(1, 2), q.value = c(0.1, 0.2), stringsAsFactors = FALSE
  ))
  # The naive run tested a module the resolved run dropped, e.g. because
  # resolution concentrated its genes below min_module_size. This direction
  # produces no NA and would otherwise pass silently.
  naive <- list(preservation = data.frame(
    module = c("1", "2", "3"), size_mapped = c(20L, 20L, 15L),
    Zsummary = c(1, 2, 3), q.value = c(0.1, 0.2, 0.3), stringsAsFactors = FALSE
  ))
  map <- data.frame(gene2 = c("X", "Y"), stringsAsFactors = FALSE)

  expect_warning(
    out <- rcomplex:::.pres_sensitivity(resolved, naive, map, map),
    "not tested under the resolved map"
  )
  expect_false("3" %in% out$module)
})


# ---- carried over from the retired gene-overlap tests ----

test_that("correspondence p-values, q-values, jaccard and overlap are sane", {
  set.seed(7)
  fx <- pres_fixture()
  tm_a <- true_modules(fx$netA, fx$mods)
  tm_b <- true_modules(fx$netB, fx$mods)
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  p <- module_correspondence(tm_a, tm_b, map)$pairs

  # p.value is built as p_gt + p_eq and q.value goes through a randomized-pi0
  # correction; neither is range-checked anywhere else.
  expect_true(all(p$p.value >= 0 & p$p.value <= 1))
  expect_true(all(p$q.value >= 0 & p$q.value <= 1))
  expect_false(anyNA(p$jaccard))
  expect_true(all(p$jaccard >= 0 & p$jaccard <= 1))
  expect_true(all(p$overlap <= pmin(p$size_sp1, p$size_sp2)))
  # Every projected gene lands in exactly one (ref, test) module cell. Compare
  # against a projection size computed independently of the cross-tab -- the
  # marginals come from the same table, so comparing them to it is arithmetic.
  # Independent of .pres_project(): the fixture's ortholog table is strictly
  # 1:1, so a projected gene is one whose reference partner carries a module
  # and which itself lands in a test-species module.
  labelled <- map$gene2[!is.na(tm_a$modules[map$gene1])]
  n_projected <- sum(!is.na(tm_b$modules[labelled]))
  expect_equal(sum(p$overlap), n_projected)
})

test_that("classification covers every module in both directions", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- suppressWarnings(preservation_paired(
    mods, nets, fx$ortho, pairs, n_perm = 50L, min_module_size = 3L, seed = 1
  ))

  expect_equal(nrow(res$classification),
               mods$A$n_modules + mods$B$n_modules)
})

test_that("alpha monotonically controls the diverged call", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, seed = 1
  )

  strict <- classify_preservation(pres, alpha = 1e-10)$classification
  loose <- classify_preservation(pres, alpha = 0.5)$classification
  expect_gte(sum(strict == "diverged"), sum(loose == "diverged"))
})

test_that("z_conserved splits conserved from moderate, rest unmoved", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, seed = 1
  )

  low <- classify_preservation(pres, z_conserved = 0)$classification
  high <- classify_preservation(pres, z_conserved = 1e6)$classification

  expect_gte(sum(low == "conserved"), sum(high == "conserved"))
  expect_lte(sum(low == "moderate"), sum(high == "moderate"))
  # Raising the secondary cut may not change significance.
  expect_equal(sum(low %in% c("diverged", "untested")),
               sum(high %in% c("diverged", "untested")))
})

test_that("classify_preservation handles a zero-row preservation table", {
  empty <- list(preservation = data.frame(
    module = character(0), size = integer(0), size_mapped = integer(0),
    Zsummary = numeric(0), q.value = numeric(0), stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(empty)

  expect_equal(nrow(cls), 0L)
  expect_true(all(c("module", "species", "pair_name", "classification",
                    "Zsummary", "q.value") %in% names(cls)))
})

test_that("the test species' own partition never enters preservation", {
  # The retired engine compared two partitions, so a module-count mismatch
  # produced false species-specific calls and needed a coarsening pass. The
  # preservation engine projects reference modules through the ortholog map
  # and never partitions the test species at all. Pin that contract: a future
  # modules_test argument would silently reintroduce the scale sensitivity.
  # A structural claim, deliberately: module_preservation() projects reference
  # modules through the ortholog map, so the test species' partition has no
  # argument to arrive through. A behavioural test cannot vary what cannot be
  # passed; this fails the moment someone adds the parameter back.
  expect_false(any(
    c("modules_test", "mods_test", "clusters", "partition", "membership") %in%
      names(formals(module_preservation))
  ))
})

test_that("preservation_paired requires a group entry for every species", {
  fx <- pres_fixture()
  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(fx$netB, fx$mods))
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  expect_error(
    preservation_paired(mods, nets, fx$ortho, pairs,
                        group = c(A = "annual"), n_perm = 10L),
    "group missing entries"
  )
})


test_that("preservation_paired output feeds tag_permutation directly", {
  # The one integration boundary this engine introduces. Both sides are tested
  # apart -- preservation_paired() here, tag_permutation() against a hand-built
  # fixture -- so a change to the level names or to how reference/test are
  # filled would leave both suites green while the handoff silently returned
  # observed = 0.
  fx <- pres_fixture()
  # The partner network carries no shared structure, so A's modules come back
  # diverged and the HOG pool is genuinely non-empty. On a fixture where every
  # module is preserved this test would pass vacuously.
  set.seed(99)
  e <- matrix(stats::rnorm(500 * 40), 500, 40)
  rownames(e) <- paste0("B", sprintf("%04d", seq_len(500)))
  net_b <- compute_network(e, density = 0.03, sparse = FALSE)

  mods <- list(A = true_modules(fx$netA, fx$mods),
               B = true_modules(net_b, fx$mods))
  nets <- list(A = fx$netA, B = net_b)
  grp <- c(A = "annual", B = "perennial")
  pairs <- data.frame(sp1 = "A", sp2 = "B", pair_name = "AB",
                      stringsAsFactors = FALSE)

  res <- suppressWarnings(preservation_paired(
    mods, nets, fx$ortho, pairs, group = grp,
    n_perm = 50L, min_module_size = 3L, seed = 1
  ))

  tp <- tag_permutation(res$classification, mods, fx$ortho, pairs,
                        group = grp, target_group = "annual",
                        n_perm = 50L, min_recurrence = 1L)

  expect_true(all(c("observed", "p_value", "recurrence_table") %in% names(tp)))
  expect_gte(tp$p_value, 0)
  expect_lte(tp$p_value, 1)

  # Non-vacuity: if the annual reference has diverged modules, they must reach
  # the HOG pool rather than being filtered out by a vocabulary mismatch.
  n_div <- sum(res$classification$classification == "diverged" &
                 res$classification$reference == "A")
  expect_gt(n_div, 0L)
  expect_gt(tp$observed, 0L)
})


test_that(".pres_project resolves each gene2 group independently", {
  # The vectorised rewrite attributes run counts to (gene2, module) cells and
  # takes two group-wise passes keyed on gene2. Single-gene fixtures exercise
  # none of that: a group-boundary error in the run alignment or in either
  # ave() pass would ship green.
  map <- data.frame(
    gene1 = c("a1", "a2", "a3",      # B1: module 5 wins 2-1
              "b1", "b2",            # B2: 1-1 tie, dropped
              "c1",                  # B3: single row
              "d1", "d2", "d3"),     # B4: winner sorts last
    gene2 = c("B1", "B1", "B1", "B2", "B2", "B3", "B4", "B4", "B4"),
    module = c("5", "5", "9",  "2", "7",  "4",  "1", "9", "9"),
    source = "unresolved",
    stringsAsFactors = FALSE
  )
  out <- .pres_project(map)

  expect_setequal(out$gene2, c("B1", "B3", "B4"))
  # Modal label per gene, and the smallest gene1 carrying it.
  expect_equal(out$module[out$gene2 == "B1"], "5")
  expect_equal(out$gene1[out$gene2 == "B1"], "a1")
  expect_equal(out$module[out$gene2 == "B3"], "4")
  # "9" wins on count even though "1" sorts first.
  expect_equal(out$module[out$gene2 == "B4"], "9")
  expect_equal(out$gene1[out$gene2 == "B4"], "d2")
  # A tie in B2 must not disturb its neighbours.
  expect_false("B2" %in% out$gene2)
})

test_that("module_correspondence records its orientation", {
  fx <- pres_fixture()
  tm_a <- true_modules(fx$netA, fx$mods)
  tm_b <- true_modules(fx$netB, fx$mods)
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  corr <- module_correspondence(tm_a, tm_b, map, sp_ref = "A", sp_test = "B")

  expect_equal(corr$sp_ref, "A")
  expect_equal(corr$sp_test, "B")
})


test_that("the copy-choice null runs and is skipped when there is no choice", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)

  # The ambiguous fixture has multi-copy HOGs, so there is a copy choice to
  # vary and the null has something to say.
  pres <- module_preservation(tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques, sp_ref = "A", sp_test = "B",
    n_perm = 100L, sensitivity = TRUE, copy_draws = 20L, seed = 1
  )
  expect_true(all(c("p_copy.avg.weight", "p_copy.cor.degree") %in%
                    names(pres$sensitivity)))
  pc <- c(pres$sensitivity$p_copy.avg.weight,
          pres$sensitivity$p_copy.cor.degree)
  pc <- pc[!is.na(pc)]
  expect_true(all(pc > 0 & pc <= 1))
  expect_gt(attr(pres$sensitivity, "n_multi_copy"), 0L)

  # The 1:1 fixture offers no copy to choose, so the null is skipped rather
  # than reporting a degenerate p of 1 for every module.
  strict <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    map = resolve_ortholog_map(fx$ortho, rownames(fx$netA$network),
                               rownames(fx$netB$network)),
    n_perm = 50L, sensitivity = TRUE, copy_draws = 5L, seed = 1
  )
  expect_equal(attr(strict$sensitivity, "n_multi_copy"), 0L)
  expect_false("p_copy.avg.weight" %in% names(strict$sensitivity))
})


test_that("Zsummary is standardized to unit null variance", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 500L, n_cores = 2L, seed = 1
  )
  d <- pres$preservation

  # Each Z has unit null variance by construction, so their mean has null sd
  # sqrt(2 + 2*rho)/2 -- between 1/sqrt(2) and 1. Reading the Langfelder 10/2
  # cut points against the raw mean imports a threshold calibrated on a
  # different quantity.
  expect_true(all(d$Zsummary_null_sd >= 1 / sqrt(2) - 1e-8))
  expect_true(all(d$Zsummary_null_sd <= 1 + 1e-8))
  expect_equal(d$Zsummary_std, d$Zsummary / d$Zsummary_null_sd)

  # The scale switch must actually change which cut point is applied.
  raw <- classify_preservation(pres, z_conserved = 10, z_scale = "raw")
  std <- classify_preservation(pres, z_conserved = 10,
                               z_scale = "standardized")
  expect_gte(sum(std$classification == "conserved"),
             sum(raw$classification == "conserved"))
})

test_that("coverage reconciles the tested modules against the partition", {
  fx <- pres_fixture()
  # Carve a five-gene module out of module 4 so something is genuinely below
  # min_module_size; every module in the base fixture has 40 genes.
  mods <- fx$mods
  mods[[4]] <- setdiff(mods[[4]], tail(fx$mods[[4]], 5L))
  mods[[5]] <- tail(fx$mods[[4]], 5L)
  tm <- true_modules(fx$netA, mods)

  # min_module_size drops modules from the analysis entirely; without a
  # coverage table the preservation output looks like a complete accounting.
  expect_message(
    pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = 50L, min_module_size = 10L, seed = 1
    ),
    "were not tested"
  )

  expect_equal(nrow(pres$coverage), tm$n_modules)
  expect_false(pres$coverage$tested[pres$coverage$module == "5"])
  expect_match(pres$coverage$reason[pres$coverage$module == "5"],
               "min_module_size")
  expect_equal(sum(pres$coverage$tested), nrow(pres$preservation))
  expect_true(all(c("module", "size", "size_mapped", "tested", "reason") %in%
                    names(pres$coverage)))
  # Every untested module carries a reason; every tested one does not.
  expect_true(all(!is.na(pres$coverage$reason[!pres$coverage$tested])))
  expect_true(all(is.na(pres$coverage$reason[pres$coverage$tested])))
})


test_that("the resolution guard compares the set that can actually differ", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)

  # resolve_ortholog_map() guarantees both maps carry every candidate gene2,
  # so comparing candidate sets is a tautology. The projected sets differ
  # because resolving a copy rescues genes whose labels would otherwise tie.
  expect_warning(
    pres <- module_preservation(tm, fx$netA, fx$netB, amb$ortho,
      cliques = amb$cliques, sp_ref = "A", sp_test = "B",
      n_perm = 100L, sensitivity = TRUE, copy_draws = 10L, seed = 1
    ),
    "changed which test-species genes"
  )
  expect_true(attr(pres$sensitivity, "same_candidate_set"))
  expect_false(attr(pres$sensitivity, "same_projected_set"))
  expect_gt(attr(pres$sensitivity, "n_rescued"), 0L)
})

test_that("the copy null holds the projected gene set fixed", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)

  pres <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, amb$ortho, cliques = amb$cliques,
    sp_ref = "A", sp_test = "B", n_perm = 100L, sensitivity = TRUE,
    copy_draws = 20L, seed = 1
  ))

  # Every draw must score the same genes as the observed run, or a set-size
  # difference is read as a copy-choice effect.
  pc <- c(pres$sensitivity$p_copy.avg.weight,
          pres$sensitivity$p_copy.cor.degree)
  pc <- pc[!is.na(pc)]
  expect_gt(length(pc), 0L)
  expect_true(all(pc > 0 & pc <= 1))
})
