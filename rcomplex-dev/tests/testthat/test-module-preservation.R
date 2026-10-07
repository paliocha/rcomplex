# Tests for module_preservation() and classify_preservation(). Fixture:
# pres_fixture() / true_modules() in helper-preservation.R.

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
  expect_equal(one$preservation$p_value, four$preservation$p_value)
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
    ortho_rand$gene2 <- shuffled
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
    ortho_rand$gene2 <- sample(genes_b, nrow(fx$ortho))
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
  expect_true(all(pres$preservation$p_value < 0.05))
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
    "cor.degree", "p_value", "q_value", "Zsummary",
    "medianRank"
  ) %in% names(pres$preservation)))
  # The diagnostics are reported but take no part in the call.
  expect_true(all(c("meanMAR", "meanClusterCoeff") %in% names(pres$observed)))
})

test_that("module_preservation validates its inputs", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  expect_error(
    module_preservation(list(a = 1), fx$netA, fx$netB, fx$ortho),
    "must be a module assignment"
  )
  expect_error(
    module_preservation(tm, fx$netA, fx$netB, fx$ortho, n_perm = 0L),
    "n_perm must be >= 1"
  )
  expect_error(
    module_preservation(tm, fx$netA, fx$netB, orthologs = NULL),
    "must be a data.frame"
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
    q_value = c(0.001, 0.001, 0.9),
    stringsAsFactors = FALSE
  ))

  cls <- classify_preservation(pres)
  expect_equal(cls$classification, c("conserved", "moderate", "diverged"))
})

test_that("classify_preservation records species and pair labels", {
  pres <- list(preservation = data.frame(
    module = "1", size = 30L, size_mapped = 20L,
    Zsummary = 15, q_value = 0.001, stringsAsFactors = FALSE
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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    n_perm = 100L, seed = 1
  )

  expect_true("medianRank" %in% names(res$classification))
  for (key in names(res$raw)) {
    # raw is keyed "<reference>.<test>"; species names that would make
    # that split ambiguous are rejected by preservation_paired() itself.
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
    Zsummary = c(15, 0.2), q_value = c(0.001, 0.9),
    stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(legacy)

  expect_equal(nrow(cls), 2L)
  expect_true(all(is.na(cls$medianRank)))
})


test_that("a seeded run restores the stream whatever its n_perm", {
  fx <- pres_fixture()
  tm <- true_modules(fx$netA, fx$mods)

  set.seed(1)
  before <- .Random.seed
  for (np in c(20L, 100L, 500L)) {
    invisible(module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = np, seed = 3
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
    "module1", "module2", "overlap", "jaccard",
    "p_value", "q_value"
  ) %in% names(p)))
  expect_equal(nrow(p), tm_a$n_modules * tm_b$n_modules)

  # The fixture maps module k of species A onto module k of species B, so the
  # diagonal must be the significant part of the table.
  diag_rows <- p$module1 == p$module2
  expect_true(all(p$q_value[diag_rows] < 0.05))
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
    "must be a module assignment"
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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

  res <- preservation_paired(mods, nets, fx$ortho, pairs,
    n_perm = 100L, seed = 1
  )

  expect_named(res, c("classification", "summary", "raw"))
  expect_setequal(names(res$raw), c("A.B", "B.A"))
  expect_setequal(unique(res$classification$reference), c("A", "B"))

  # Downstream consumers read exactly these columns.
  expect_true(all(c(
    "pair_name", "module", "reference", "test",
    "classification"
  ) %in% names(res$classification)))
})

test_that("colliding '<reference>.<test>' keys are refused up front", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  # ref "A" / test "B.C" and ref "A.B" / test "C" both key as "A.B.C".
  # The names only have to collide; the run must stop before computing
  # anything, so the modules behind them are irrelevant.
  mods4 <- stats::setNames(
    list(mods$A, mods$B, mods$A, mods$B),
    c("A", "B.C", "A.B", "C")
  )
  nets4 <- stats::setNames(
    list(fx$netA, fx$netB, fx$netA, fx$netB),
    c("A", "B.C", "A.B", "C")
  )
  pairs <- data.frame(
    species1 = c("A", "A.B"), species2 = c("B.C", "C"),
    stringsAsFactors = FALSE
  )
  expect_error(
    preservation_paired(mods4, nets4, fx$ortho, pairs, n_perm = 5L),
    "colliding contrast keys"
  )
})

test_that("preservation_paired tags trait groups", {
  fx <- pres_fixture()
  mods <- list(
    A = true_modules(fx$netA, fx$mods),
    B = true_modules(fx$netB, fx$mods)
  )
  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

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
    q_value = c(NA_real_, 0.001),
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
        species1 = c("A", "B"), species2 = c("B", "A"),
        stringsAsFactors = FALSE
      ),
      n_perm = 10L
    ),
    "same species pair more than once"
  )
  expect_error(
    preservation_paired(mods, nets, fx$ortho,
      data.frame(species1 = "A", species2 = "A", stringsAsFactors = FALSE),
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
      gene1 = paste0("A", sprintf("%02d", seq_len(n))),
      gene2 = paste0("B", sprintf("%02d", seq_len(n))),
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
  expect_true(all(is.na(cls$q_value[cls$classification == "untested"])))
  # Upstream property, not a restatement of the classifier's own definition:
  # cor.degree is computable for the two blocks with non-constant degree.
  expect_false(anyNA(
    pres$preservation$cor.degree[pres$preservation$module != "3"]
  ))
})

test_that("untested modules stay visible in the paired summary", {
  fx <- isolated_fixture()
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

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

  # p_value is built as p_gt + p_eq and q_value goes through a randomized-pi0
  # correction; neither is range-checked anywhere else.
  expect_true(all(p$p_value >= 0 & p$p_value <= 1))
  expect_true(all(p$q_value >= 0 & p$q_value <= 1))
  expect_false(anyNA(p$jaccard))
  expect_true(all(p$jaccard >= 0 & p$jaccard <= 1))
  expect_true(all(p$overlap <= pmin(p$size1, p$size2)))
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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

  res <- suppressWarnings(preservation_paired(
    mods, nets, fx$ortho, pairs,
    n_perm = 50L, seed = 1
  ))

  expect_equal(
    nrow(res$classification),
    mods$A$n_modules + mods$B$n_modules
  )
})

test_that("classify_preservation handles a zero-row preservation table", {
  empty <- list(preservation = data.frame(
    module = character(0), size = integer(0), size_mapped = integer(0),
    Zsummary = numeric(0), q_value = numeric(0), stringsAsFactors = FALSE
  ))
  cls <- classify_preservation(empty)

  expect_equal(nrow(cls), 0L)
  expect_true(all(c(
    "module", "species", "pair_name", "classification",
    "Zsummary", "q_value"
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
  pairs <- data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE)

  expect_error(
    preservation_paired(mods, nets, fx$ortho, pairs,
      group = c(A = "annual"), n_perm = 10L
    ),
    "group missing entries"
  )
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
  corr <- module_correspondence(tm_a, tm_b, map,
    species_ref = "A", species_test = "B"
  )

  expect_equal(corr$species_ref, "A")
  expect_equal(corr$species_test, "B")
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
})

test_that("coverage reconciles the tested modules against the partition", {
  fx <- pres_fixture()
  # Carve a five-gene module out of module 4 so something is genuinely below
  # 10 mapped genes; every module in the base fixture has 40 genes.
  mods <- fx$mods
  mods[[4]] <- setdiff(mods[[4]], tail(fx$mods[[4]], 5L))
  mods[[5]] <- tail(fx$mods[[4]], 5L)
  tm <- true_modules(fx$netA, mods)

  # Small modules drop out of the analysis entirely; without a
  # coverage table the preservation output looks like a complete accounting.
  expect_message(
    pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
      n_perm = 50L, seed = 1
    ),
    "were not tested"
  )

  expect_equal(nrow(pres$coverage), tm$n_modules)
  expect_false(pres$coverage$tested[pres$coverage$module == "5"])
  expect_match(
    pres$coverage$reason[pres$coverage$module == "5"],
    "fewer than 10 mapped genes"
  )

  # The other arm: a module whose genes have no ortholog at all never enters
  # the size vector, so it was invisible even in the module-size table.
  mods2 <- mods
  # Indices, not names: setdiff() on a character vector against integers
  # coerces and matches nothing, which would silently pick module-1 genes.
  bg <- setdiff(seq_len(nrow(fx$netA$network)), unlist(mods))
  mods2[[6]] <- bg[1:12]
  tm2 <- true_modules(fx$netA, mods2)
  ortho_partial <- fx$ortho[fx$ortho$gene1 %in%
                              rownames(fx$netA$network)[unlist(mods)], ]
  expect_message(
    pres2 <- module_preservation(tm2, fx$netA, fx$netB, ortho_partial,
      n_perm = 50L, seed = 1
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
  eA <- pres_expr(41, 600, "A", loadings, per)  # nolint
  eB <- pres_expr(42, 700, "B", loadings, per)  # nolint
  mods <- lapply(seq_len(n_mod), function(k) ((k - 1) * per + 1):(k * per))
  list(
    netA = compute_network(eA, density = 0.03, sparse = FALSE),
    netB = compute_network(eB, density = 0.03, sparse = FALSE),
    mods = mods, n_mod = n_mod,
    ortho = data.frame(
      gene1 = paste0("A", sprintf("%04d", seq_len(600))),
      gene2 = paste0("B", sprintf("%04d", seq_len(600))),
      hog = paste0("H", seq_len(600)), stringsAsFactors = FALSE
    )
  )
}

test_that("calibration lies between the joint null and raw pmax", {
  fx <- calib_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 999L, n_cores = 2L, seed = 1
  )
  d <- pres$preservation

  # pmax is the intersection-union p-value and the calibrated value is a
  # blend of it with the empirical joint null, so it can only move downward.
  expect_true(all(d$p_calibrated <= d$p_value + 1e-12, na.rm = TRUE))
  expect_true(all(d$p_calibrated > 0, na.rm = TRUE))
  # The identity the blend's validity rests on.
  expect_equal(d$p_value, pmax(d$p.avg.weight, d$p.cor.degree))
  # The estimator is live at this module count.
  expect_gte(pres$params$w00, 0)
  expect_lte(pres$params$w00, 1)
})

test_that("calibrate = none reproduces the uncalibrated pmax result", {
  fx <- calib_fixture()
  tm <- true_modules(fx$netA, fx$mods)
  pres <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 999L, n_cores = 2L, seed = 1,
    calibrate = "none"
  )
  d <- pres$preservation

  expect_identical(d$p_calibrated, d$p_value)
  expect_equal(pres$params$w00, 0)
  # and the q-values are BH on the raw pmax
  expect_equal(d$q_value, stats::p.adjust(d$p_value, "BH"))
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
    n_perm = 200L, n_cores = 1L, seed = 5
  )
  b <- module_preservation(tm, fx$netA, fx$netB, fx$ortho,
    n_perm = 200L, n_cores = 4L, seed = 5
  )
  expect_equal(a$preservation$p_calibrated, b$preservation$p_calibrated)
  expect_equal(a$params$w00, b$params$w00)
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
    n_perm = 50L, seed = 1
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
  expect_true(all(cls$classification[!is.na(pr$preservation$q_value)] %in%
                    c("moderate", "diverged")))

  # An untested module must not trigger the fallback: it has no Zsummary
  # to fall back to and is already reported as untested.
  untested <- pr
  untested$preservation$q_value <- NA_real_
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
