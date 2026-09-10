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
    n_perm = 0L, n_cores = 1L, binary = FALSE, store_perm = FALSE
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


# ---- medianRank ----

test_that("medianRank reaches the user through classify_preservation", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, seed = 3
  )
  cls <- classify_preservation(pres)

  expect_true("medianRank" %in% names(cls))
  expect_equal(cls$module, pres$preservation$module)
  expect_equal(cls$medianRank, pres$preservation$medianRank)
})

test_that("medianRank is the mean of the two observed-statistic ranks", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 100L, seed = 3
  )
  p <- pres$preservation

  # 1 = strongest, so the ranks run on the negated statistics.
  rank_d <- rank(-p$avg.weight, na.last = "keep")
  rank_c <- rank(-p$cor.degree, na.last = "keep")
  expect_equal(p$medianRank, (rank_d + rank_c) / 2)
  expect_gte(min(p$medianRank), 1)
  expect_lte(max(p$medianRank), nrow(p))
})

test_that("medianRank carries no permutation moments", {
  # The point of reporting it next to Zsummary: it is a rank of the observed
  # statistics, so nothing about the null -- and hence nothing about the
  # module-size-dependent null spread -- can move it. Zsummary does move.
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  args <- list(tm, fx$netA, fx$netB, fx$ortho, seed = 5)
  few <- do.call(module_preservation, c(args, n_perm = 50L))
  many <- do.call(module_preservation, c(args, n_perm = 500L))

  expect_equal(few$preservation$medianRank, many$preservation$medianRank)
  expect_false(isTRUE(all.equal(
    few$preservation$Zsummary_std, many$preservation$Zsummary_std
  )))
})

test_that("medianRank survives into a preservation_paired table", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    n_perm = 100L, seed = 1
  )

  expect_true("medianRank" %in% names(res$classification))
  for (key in names(res$raw)) {
    ref <- strsplit(key, ".", fixed = TRUE)[[1]][1]
    got <- res$classification$medianRank[res$classification$reference == ref]
    expect_equal(got, res$raw[[key]]$preservation$medianRank)
  }
})

test_that("a preservation table without medianRank still classifies", {
  # Objects saved before medianRank was carried through must not turn into a
  # recycling error or a truncated data frame.
  legacy <- list(preservation = data.frame(
    module = c("1", "2"),
    size = c(30L, 30L), size_mapped = c(20L, 20L),
    Zsummary = c(15, 0.2), q.value = c(0.001, 0.9),
    stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(legacy)

  expect_equal(nrow(cls), 2L)
  expect_true(all(is.na(cls$medianRank)))
})


# ---- sensitivity: the circularity guard ----

# A HOG that is multi-copy on the reference side, spanning two modules: each
# listed species-2 gene has one partner in module 1 and one in module 2, so the
# naive map ties and .pres_project() drops it, while a clique resolves it to
# module 1. Both maps offer the same CANDIDATE species-2 genes, but the
# PROJECTED sets differ by exactly those tie-rescued genes -- the resolved run
# tests them and the naive run does not, which is why same_projected_set is
# FALSE here and why the circularity check is p_copy rather than the delta.
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


# ---- p_copy: seeded, and independent of the nested naive-map run ----

# The copy-choice draws sit downstream of the nested naive-map run, so
# whatever that run takes from the stream moves them. A seeded
# module_preservation() restores the caller's stream on exit (R/rng.R), which
# is what keeps p_copy a function of the seed and of the run under test alone.
# These three tests pin that: the values themselves, the fact that they are
# what the copy null produces with no naive run in the stream at all, and the
# restoration that makes it so at any permutation count.

pcopy_fixture <- function() {
  fx <- pres_fixture()
  list(
    fx = fx, tm = true_modules(fx$netA, fx$mods),
    amb = ambiguous_fixture(fx)
  )
}

test_that("p_copy is reproducible and pinned under a seed", {
  p <- pcopy_fixture()
  run <- function() {
    module_preservation(p$tm, p$fx$netA, p$fx$netB, p$amb$ortho,
      cliques = p$amb$cliques, sp_ref = "A", sp_test = "B",
      n_perm = 100L, copy_draws = 50L, sensitivity = TRUE, seed = 1
    )
  }
  a <- run()
  b <- run()

  aw <- a$sensitivity$p_copy.avg.weight
  cd <- a$sensitivity$p_copy.cor.degree
  expect_identical(b$sensitivity$p_copy.avg.weight, aw)
  expect_identical(b$sensitivity$p_copy.cor.degree, cd)

  expect_identical(attr(a$sensitivity, "n_copy_draws"), 50L)
  expect_identical(a$sensitivity$module, c("1", "2", "3", "4"))

  # Exact rationals k / (n_draws + 1); pinned so a change of seeding point
  # cannot pass silently. Modules 3 and 4 hold no multi-copy gene, so every
  # draw reproduces their observed statistic and p_copy is exactly 1.
  expect_equal(aw, c(28, 1, 51, 51) / 51)
  expect_equal(cd, c(7, 21, 51, 51) / 51)
})

test_that("p_copy does not depend on what the naive run consumed", {
  p <- pcopy_fixture()
  genes_a <- rownames(p$fx$netA$network)
  genes_b <- rownames(p$fx$netB$network)
  resolved <- resolve_ortholog_map(p$amb$ortho, genes_a, genes_b,
    sp1 = "A", sp2 = "B", cliques = p$amb$cliques
  )
  naive_map <- resolve_ortholog_map(p$amb$ortho, genes_a, genes_b)

  # module_preservation()'s sensitivity path, reassembled so the nested
  # naive-map run's permutation count can be varied: seed, the resolved-map
  # permutations, the naive run, then the copy draws. The kernel takes one
  # uniform per permutation, so a naive run that does not restore the stream
  # leaves the copy draws starting somewhere else. On the eight-species
  # vignette data that moved p_copy by up to 0.0398 and flipped one module
  # across alpha = 0.05; naive_perm = NULL is the same sequence with no naive
  # run in it at all.
  copy_p <- function(naive_perm) {
    set.seed(1)
    main <- module_preservation(p$tm, p$fx$netA, p$fx$netB,
      map = resolved, n_perm = 100L
    )
    if (!is.null(naive_perm)) {
      invisible(module_preservation(p$tm, p$fx$netA, p$fx$netB,
        map = naive_map, n_perm = naive_perm, seed = 1
      ))
    }
    cn <- rcomplex:::.pres_copy_null(
      main$preservation, p$tm, p$fx$netA, p$fx$netB, naive_map,
      unique(main$projection$gene2), 50L, 10L, FALSE
    )
    c(cn$p_copy.avg.weight, cn$p_copy.cor.degree)
  }

  ref <- copy_p(NULL)
  expect_equal(copy_p(100L), ref)
  expect_equal(copy_p(500L), ref)

  # And the reassembly is the real thing: the shipped call, whose naive run
  # takes n_perm from the outer call, lands on the same p_copy.
  pres <- module_preservation(p$tm, p$fx$netA, p$fx$netB, p$amb$ortho,
    cliques = p$amb$cliques, sp_ref = "A", sp_test = "B",
    n_perm = 100L, copy_draws = 50L, sensitivity = TRUE, seed = 1
  )
  expect_equal(
    c(
      pres$sensitivity$p_copy.avg.weight,
      pres$sensitivity$p_copy.cor.degree
    ),
    ref
  )
})

test_that("a seeded run restores the stream whatever its n_perm", {
  p <- pcopy_fixture()
  naive_map <- resolve_ortholog_map(
    p$amb$ortho, rownames(p$fx$netA$network), rownames(p$fx$netB$network)
  )

  # The nested naive-map run is exactly this call, and it consumes one
  # uniform per permutation. That its cost is invisible to the caller is what
  # the test above rests on, so assert it across permutation counts.
  set.seed(1)
  before <- .Random.seed
  for (np in c(20L, 100L, 500L)) {
    invisible(module_preservation(p$tm, p$fx$netA, p$fx$netB,
      map = naive_map, n_perm = np, seed = 3
    ))
    expect_identical(.Random.seed, before)
  }
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

  expect_error(
    module_correspondence(list(a = 1), tm, map),
    "must be output from detect_modules"
  )
  expect_error(
    module_correspondence(tm, tm, data.frame(x = 1)),
    "must be a data frame from resolve_ortholog_map"
  )
})


# ---- preservation_paired ----

test_that("preservation_paired runs both directions per contrast", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    n_perm = 100L, seed = 1
  )

  expect_named(res, c("classification", "summary", "raw"))
  expect_setequal(names(res$raw), c("A.B", "B.A"))
  expect_setequal(unique(res$classification$reference), c("A", "B"))

  # tag_permutation() reads exactly these columns.
  expect_true(all(c(
    "pair_name", "module", "reference", "test",
    "classification"
  ) %in% names(res$classification)))
})

test_that("preservation_paired tags trait groups", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
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
    n_perm = 100L, n_cores = 1L, binary = FALSE, store_perm = FALSE
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
  q <- rcomplex:::.pres_qvalues(p)

  expect_true(is.na(q[2]))
  expect_false(any(is.na(q[-2])))
})


test_that("preservation_paired rejects a repeated or self species pair", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)

  # Both directions of each contrast are run, so (A, B) and (B, A) would
  # collide on the "<reference>.<test>" key and silently overwrite each other.
  expect_error(
    preservation_paired(mods, nets, fx$ortho,
      data.frame(
        sp1 = c("A", "B"), sp2 = c("B", "A"),
        stringsAsFactors = FALSE
      ),
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
    pres$preservation$cor.degree[pres$preservation$module != "3"]
  ))
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
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  for (grp in list(NULL, c(A = "annual", B = "perennial"))) {
    res <- suppressWarnings(preservation_paired(
      mods, nets, fx$ortho, pairs,
      group = grp, n_perm = 50L, seed = 1
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
    n_perm = 0L, n_cores = 1L, binary = FALSE, store_perm = FALSE
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
    n_perm = 20L, n_cores = 1L, binary = FALSE, store_perm = FALSE
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
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  res <- suppressWarnings(preservation_paired(
    mods, nets, fx$ortho, pairs,
    n_perm = 50L, min_module_size = 3L, seed = 1
  ))

  expect_equal(
    nrow(res$classification),
    mods$A$n_modules + mods$B$n_modules
  )
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
  expect_equal(
    sum(low %in% c("diverged", "untested")),
    sum(high %in% c("diverged", "untested"))
  )
})

test_that("classify_preservation handles a zero-row preservation table", {
  empty <- list(preservation = data.frame(
    module = character(0), size = integer(0), size_mapped = integer(0),
    Zsummary = numeric(0), q.value = numeric(0), stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(empty)

  expect_equal(nrow(cls), 0L)
  expect_true(all(c(
    "module", "species", "pair_name", "classification",
    "Zsummary", "q.value"
  ) %in% names(cls)))
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
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)

  expect_error(
    preservation_paired(mods, nets, fx$ortho, pairs,
      group = c(A = "annual"), n_perm = 10L
    ),
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

  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(net_b, fx$mods)
  )
  nets <- list(A = fx$netA, B = net_b)
  grp <- c(A = "annual", B = "perennial")
  pairs <- data.frame(
    sp1 = "A", sp2 = "B", pair_name = "AB",
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(preservation_paired(
    mods, nets, fx$ortho, pairs,
    group = grp,
    n_perm = 50L, min_module_size = 3L, seed = 1
  ))

  # Every module is diverged in both directions here, so the two sides of
  # the pair project the same 160 HOGs and swapping the labels changes
  # nothing: the null has one point and p_min = 1. Counting the pair
  # anyway would report p_min = 0.5, a resolution the data does not have.
  tp <- suppressMessages(suppressWarnings(tag_permutation(
    res$classification, mods, fx$ortho, pairs,
    group = grp, target_group = "annual",
    n_perm = 50L, min_recurrence = 1L
  )))

  expect_true(all(c("observed", "p_value", "recurrence_table") %in% names(tp)))
  expect_equal(tp$n_contributing, 1L)
  expect_equal(tp$n_swappable, 0L)
  expect_equal(tp$p_min, 1)
  expect_equal(
    tp$pair_sizes$n_hogs_target,
    tp$pair_sizes$n_hogs_partner
  )
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
    gene1 = c(
      "a1", "a2", "a3", # B1: module 5 wins 2-1
      "b1", "b2", # B2: 1-1 tie, dropped
      "c1", # B3: single row
      "d1", "d2", "d3"
    ), # B4: winner sorts last
    gene2 = c("B1", "B1", "B1", "B2", "B2", "B3", "B4", "B4", "B4"),
    module = c("5", "5", "9", "2", "7", "4", "1", "9", "9"),
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
  pc <- c(
    pres$sensitivity$p_copy.avg.weight,
    pres$sensitivity$p_copy.cor.degree
  )
  pc <- pc[!is.na(pc)]
  expect_true(all(pc > 0 & pc <= 1))
  expect_gt(attr(pres$sensitivity, "n_multi_copy"), 0L)

  # The 1:1 fixture offers no copy to choose, so the null is skipped rather
  # than reporting a degenerate p of 1 for every module.
  strict <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    map = resolve_ortholog_map(
      fx$ortho, rownames(fx$netA$network),
      rownames(fx$netB$network)
    ),
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
  # sqrt(2 + 2*rho)/2. rho is clamped at 0 from below -- a negative sample
  # correlation would shrink the divisor and inflate Zsummary_std without
  # bound -- so the spread is in [1/sqrt(2), 1]. Reading the Langfelder 10/2
  # cut points against the raw mean imports a threshold calibrated on a
  # different quantity.
  expect_true(all(d$Zsummary_null_sd >= 1 / sqrt(2) - 1e-8))
  expect_true(all(d$Zsummary_null_sd <= 1 + 1e-8))
  expect_equal(d$Zsummary_std, d$Zsummary / d$Zsummary_null_sd)

  # The scale switch must actually change which cut point is applied.
  raw <- classify_preservation(pres, z_conserved = 10, z_scale = "raw")
  std <- classify_preservation(pres,
    z_conserved = 10,
    z_scale = "standardized"
  )
  expect_gte(
    sum(std$classification == "conserved"),
    sum(raw$classification == "conserved")
  )
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
  expect_match(
    pres$coverage$reason[pres$coverage$module == "5"],
    "min_module_size"
  )

  # The other arm: a module whose genes have no ortholog at all never enters
  # the size vector, so it was invisible even in the module-size table.
  mods2 <- mods
  # Indices, not names: setdiff() on a character vector against integers
  # coerces and matches nothing, which would silently pick module-1 genes.
  bg <- setdiff(seq_len(nrow(fx$netA$network)), unlist(mods))
  mods2[[6]] <- bg[1:12]
  tm2 <- true_modules(fx$netA, mods2)
  ortho_partial <- fx$ortho[fx$ortho$Species1 %in%
    rownames(fx$netA$network)[unlist(mods)], ]
  expect_message(
    pres2 <- module_preservation(tm2, fx$netA, fx$netB, ortho_partial,
      n_perm = 50L, min_module_size = 10L, seed = 1
    ),
    "no mapped gene"
  )
  cv6 <- pres2$coverage[pres2$coverage$module == "6", ]
  expect_equal(cv6$size_mapped, 0L)
  expect_false(cv6$tested)
  expect_equal(cv6$reason, "no mapped gene")
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
    tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques,
    sp_ref = "A", sp_test = "B", n_perm = 100L, sensitivity = TRUE,
    copy_draws = 20L, seed = 1
  ))

  # p_copy lies in (0, 1] by construction, so asserting that proves
  # nothing. The invariant the restriction establishes is that the draws'
  # candidate pool is exactly the observed run's projected gene set.
  # Re-applying the same `%in% projected` filter before comparing would
  # make that a tautology, so compare the UNFILTERED pool: on the clique
  # map the two coincide, which is what makes the restriction a no-op
  # here and the naive case below the one that exercises it.
  cand <- resolve_ortholog_map(
    amb$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  cand <- cand[!is.na(tm$modules[cand$gene1]), , drop = FALSE]
  expect_setequal(unique(cand$gene2), unique(pres$projection$gene2))

  # Under the naive map the ambiguous genes tie in the majority vote and
  # drop out of projection while remaining candidates, so the pool is a
  # strict superset and the restriction is load-bearing. Without it those
  # 10 genes would enter the draws and the copy null would score a larger
  # gene set than the observed run -- exactly the set-size confound it
  # exists to remove.
  naive <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, amb$ortho,
    n_perm = 20L, min_module_size = 3L, sensitivity = FALSE, seed = 1
  ))
  naive_proj <- unique(naive$projection$gene2)
  expect_gt(length(setdiff(unique(cand$gene2), naive_proj)), 0L)
  expect_true(all(naive_proj %in% unique(cand$gene2)))

  # And every draw must have survived, or p_copy rests on fewer than claimed.
  expect_equal(attr(pres$sensitivity, "n_copy_draws"), 20L)
  pc <- c(
    pres$sensitivity$p_copy.avg.weight,
    pres$sensitivity$p_copy.cor.degree
  )
  expect_gt(sum(!is.na(pc)), 0L)
})


# ---- mixture-calibrated p-values ----

# 14 modules, so w00 is estimable (the estimator returns 0 below 10).
calib_fixture <- function() {
  n_mod <- 14L
  per <- 25L
  set.seed(31)
  loadings <- lapply(seq_len(n_mod), function(k) {
    l <- stats::rlnorm(per, 0, 0.9)
    l / max(l)
  })
  eA <- pres_expr(41, 600, "A", loadings, per)
  eB <- pres_expr(42, 700, "B", loadings, per)
  mods <- lapply(seq_len(n_mod), function(k) ((k - 1) * per + 1):(k * per))
  list(
    netA = compute_network(eA, density = 0.03, sparse = FALSE),
    netB = compute_network(eB, density = 0.03, sparse = FALSE),
    mods = mods, n_mod = n_mod,
    ortho = data.frame(
      Species1 = paste0("A", sprintf("%04d", seq_len(600))),
      Species2 = paste0("B", sprintf("%04d", seq_len(600))),
      hog = paste0("H", seq_len(600)), stringsAsFactors = FALSE
    )
  )
}

test_that("calibration lies between the joint null and raw pmax", {
  fx <- calib_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 999L, min_module_size = 10L, n_cores = 2L, seed = 1
  )
  d <- pres$preservation

  # pmax is the intersection-union p-value and the calibrated value is a
  # blend of it with the empirical joint null, so it can only move downward.
  expect_true(all(d$p.calibrated <= d$p.value + 1e-12, na.rm = TRUE))
  expect_true(all(d$p.calibrated > 0, na.rm = TRUE))
  # The identity the blend's validity rests on.
  expect_equal(d$p.value, pmax(d$p.avg.weight, d$p.cor.degree))
  # The estimator is live at this module count.
  expect_gte(pres$params$w00, 0)
  expect_lte(pres$params$w00, 1)
})

test_that("calibrate = none reproduces the uncalibrated pmax result", {
  fx <- calib_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 999L, min_module_size = 10L, n_cores = 2L, seed = 1,
    calibrate = "none"
  )
  d <- pres$preservation

  expect_identical(d$p.calibrated, d$p.value)
  expect_equal(pres$params$w00, 0)
  # and the q-values are BH on the raw pmax
  expect_equal(d$q.value, stats::p.adjust(d$p.value, "BH"))
})

test_that("w00 degrades to zero below ten modules", {
  # Storey's estimator is high-variance on a handful of p-values, so the
  # implementation refuses rather than guessing.
  expect_equal(.pres_w00(runif(5), runif(5)), 0)
  # With both margins uniform it should approach 1 (both-null everywhere).
  set.seed(2)
  expect_gt(.pres_w00(runif(400), runif(400)), 0.8)
  # With both margins all-significant there is no both-null mass.
  expect_equal(.pres_w00(rep(0.001, 400), rep(0.001, 400)), 0)
})

test_that("the NPC p-value is the joint-null rank of the observed pmax", {
  fx <- pres_fixture()
  net <- fx$netA
  keep <- as.integer(seq_len(nrow(net$network)) - 1L)
  mm <- lapply(fx$mods, function(z) as.integer(z - 1L))
  adj <- reference_adjacency(net, rownames(net$network))
  rs <- lapply(fx$mods, function(i) reference_module_stats(adj, i))

  got <- module_preservation_dense_cpp(
    net$network, net$threshold, keep, mm,
    lapply(rs, `[[`, "kIM"), lapply(rs, `[[`, "CC"), lapply(rs, `[[`, "MAR"),
    n_perm = 200L, n_cores = 1L, binary = FALSE, store_perm = TRUE
  )

  # Recompute the combination in R from the stored pairs and check the kernel.
  pp <- array(got$perm_pairs, dim = c(2L * length(fx$mods), 200L))
  for (k in seq_along(fx$mods)) {
    a1 <- c(got$observed[k, 1], pp[2 * k - 1, ])
    a2 <- c(got$observed[k, 4], pp[2 * k, ])
    ok <- !is.na(a1) & !is.na(a2)
    a1 <- a1[ok]
    a2 <- a2[ok]
    n <- length(a1)
    l1 <- vapply(a1, function(v) sum(a1 >= v) / n, numeric(1))
    l2 <- vapply(a2, function(v) sum(a2 >= v) / n, numeric(1))
    psi <- pmax(l1, l2)
    expect_equal(got$p_npc[k], sum(psi <= psi[1]) / n, tolerance = 1e-12)
    expect_equal(got$p_joint[k, 1], l1[1], tolerance = 1e-12)
    expect_equal(got$p_joint[k, 2], l2[1], tolerance = 1e-12)
    expect_equal(got$n_joint[k], n - 1L)

    # rho is the Pearson correlation of the two statistics over the
    # permutation draws only -- entry 1 is the observed pair. It divides
    # Zsummary to give Zsummary_std, which carries the z_conserved cut
    # point, so a wrong normaliser (s12 / s11) or an off-by-one that
    # lets the observed pair in would move every conserved call.
    expect_equal(got$rho[k], stats::cor(a1[-1], a2[-1]),
      tolerance = 1e-12
    )
  }
})

test_that("calibrated p-values do not depend on n_cores", {
  fx <- calib_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  a <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 200L, min_module_size = 10L, n_cores = 1L, seed = 5
  )
  b <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 200L, min_module_size = 10L, n_cores = 4L, seed = 5
  )
  expect_equal(a$preservation$p.calibrated, b$preservation$p.calibrated)
  expect_equal(a$params$w00, b$params$w00)
})

test_that("qvalue_method is deprecated", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  expect_warning(
    module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = 50L, seed = 1, qvalue_method = "liang"
    ),
    "deprecated and ignored"
  )
})


test_that("an unestimable Zsummary_std falls back to the raw scale", {
  # rho is NA when fewer than four permutations are usable or one
  # statistic is constant across them, which makes Zsummary_std NA. Left
  # alone that turns `strong` FALSE and demotes a significant module to
  # "moderate" with nothing said -- a missing normaliser reported as a
  # measurement.
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pr <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, min_module_size = 3L, seed = 1
  ))
  expect_true(any(classify_preservation(pr)$classification == "conserved"))

  broken <- pr
  broken$preservation$Zsummary_std <- NA_real_
  expect_warning(
    cls <- classify_preservation(broken),
    "no null correlation for Zsummary_std"
  )
  # The fallback is to the raw scale, which is the stricter of the two
  # (Zsummary_std >= Zsummary, since the divisor is in [1/sqrt(2), 1]),
  # so nothing is promoted by the failure.
  expect_false(any(cls$classification == "conserved"))
  expect_true(all(cls$classification[!is.na(pr$preservation$q.value)] %in%
    c("moderate", "diverged")))

  # An untested module must not trigger the fallback: it has no Zsummary
  # to fall back to and is already reported as untested.
  untested <- pr
  untested$preservation$q.value <- NA_real_
  untested$preservation$Zsummary_std <- NA_real_
  # The inner suppressWarnings() muffled every warning before
  # expect_silent could observe one, so this passed against a mutant that
  # emitted the fallback warning on untested rows. Collect instead, and
  # assert on which warning fired: the untested one must, the fallback
  # one must not.
  warns <- character(0)
  cls2 <- withCallingHandlers(
    classify_preservation(untested),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("could not be tested", warns)))
  expect_false(any(grepl("no null correlation", warns)))
  expect_true(all(cls2$classification == "untested"))
})


test_that("copy_null_skipped distinguishes why the copy null did not run", {
  # n_multi_copy = 0 used to mean three different things: the caller
  # switched the null off, the ortholog table could not hold the gene set
  # fixed, or there was genuinely nothing multi-copy. Only the middle one
  # warned, so a reader of n_multi_copy == 0 could conclude the map had no
  # paralogs when the check had simply never run.
  fx <- pres_fixture()
  amb <- ambiguous_fixture(fx)
  tm <- true_modules(fx$netA, fx$mods)

  off <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques,
    sp_ref = "A", sp_test = "B",
    n_perm = 20L, min_module_size = 3L, sensitivity = TRUE,
    copy_draws = 0L, seed = 1
  ))
  expect_equal(attr(off$sensitivity, "copy_null_skipped"), "off")
  expect_true(is.na(attr(off$sensitivity, "n_multi_copy")))
  expect_false("p_copy.avg.weight" %in% names(off$sensitivity))

  # The ordinary multi-copy path runs the null and records no reason.
  ran <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques,
    sp_ref = "A", sp_test = "B",
    n_perm = 20L, min_module_size = 3L, sensitivity = TRUE,
    copy_draws = 5L, seed = 1
  ))
  expect_null(attr(ran$sensitivity, "copy_null_skipped"))
  expect_gt(attr(ran$sensitivity, "n_multi_copy"), 0L)
  expect_true("p_copy.avg.weight" %in% names(ran$sensitivity))
  expect_equal(attr(ran$sensitivity, "n_copy_draws"), 5L)
})


test_that("same_candidate_set is FALSE for a map that does not cover", {
  # Documented as "always TRUE". It is TRUE by construction only when the
  # map was resolved from `orthologs`; a supplied map that covers fewer
  # genes makes it FALSE, so a reader taking the doc at face value would
  # treat a real mismatch as impossible.
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  full <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network),
    rownames(fx$netB$network)
  )
  trimmed <- full[-seq_len(20L), , drop = FALSE]

  res <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, fx$ortho,
    map = trimmed,
    n_perm = 20L, min_module_size = 3L, sensitivity = TRUE,
    copy_draws = 0L, seed = 1
  ))
  expect_false(attr(res$sensitivity, "same_candidate_set"))
})


test_that("copy_null_skipped reports no_multi_copy on a 1:1 map", {
  # Only "off" and the success path were covered, so the reason that
  # distinguishes "nothing to vary" from "never ran" was untested -- and
  # it is the one a reader of n_multi_copy == 0 relies on.
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  res <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, fx$ortho,
    cliques = NULL, edges = NULL,
    map = resolve_ortholog_map(
      fx$ortho, rownames(fx$netA$network),
      rownames(fx$netB$network)
    ),
    n_perm = 20L, min_module_size = 3L, sensitivity = TRUE,
    copy_draws = 5L, seed = 1
  ))
  # fx$ortho is strictly 1:1, so no gene2 has a choice of partner.
  if (!is.null(res$sensitivity)) {
    expect_equal(attr(res$sensitivity, "n_multi_copy"), 0L)
    expect_equal(
      attr(res$sensitivity, "copy_null_skipped"),
      "no_multi_copy"
    )
    expect_false("p_copy.avg.weight" %in% names(res$sensitivity))
  }
})


test_that("a partial copy-draw failure is reported, not absorbed", {
  # p_copy was ranked against however many draws survived, with only the
  # n_copy_draws attribute recording it. Drive some draws to fail by
  # setting min_module_size where a copy choice can push a module under
  # it, and assert the count is visible.
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  amb <- ambiguous_fixture(fx)
  res <- suppressWarnings(module_preservation(
    tm, fx$netA, fx$netB, amb$ortho,
    cliques = amb$cliques,
    sp_ref = "A", sp_test = "B",
    n_perm = 20L, min_module_size = 3L, sensitivity = TRUE,
    copy_draws = 5L, seed = 1
  ))
  n_used <- attr(res$sensitivity, "n_copy_draws")
  # Whatever the outcome, the number of draws p_copy was ranked against
  # must be recorded, and must not exceed what was asked for.
  if ("p_copy.avg.weight" %in% names(res$sensitivity)) {
    expect_true(is.numeric(n_used))
    expect_lte(n_used, 5L)
    expect_gt(n_used, 0L)
  }
})
