# Tests for clique_intensity_test()
# Fixtures: make_clique_fixture() and make_clique_fixture_3sp()
# from helper-clique-fixtures.R


test_that("clique_intensity_test output has correct structure", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L, seed = 42L
  )

  expect_true(is.data.frame(result))
  expected_cols <- c(
    "clique_idx", "hog", "observed_intensity",
    "null_mean", "null_sd", "z_score", "p_value",
    "n_perm", "n_matched"
  )
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true(all(result$clique_idx >= 1L))
  expect_true(all(result$n_perm == 3L))
  expect_true(all(result$n_matched >= 0L & result$n_matched <= 3L))
  expect_true(is.numeric(result$observed_intensity))
  expect_true(is.numeric(result$null_mean))
  expect_true(is.numeric(result$null_sd))
})


test_that("clique_intensity_test seed produces reproducible results", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  r1 <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 5L, seed = 123L
  )
  r2 <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 5L, seed = 123L
  )

  expect_equal(r1$null_mean, r2$null_mean)
  expect_equal(r1$null_sd, r2$null_sd)
  expect_equal(r1$z_score, r2$z_score)
  expect_equal(r1$p_value, r2$p_value)
  expect_equal(r1$n_matched, r2$n_matched)
})


test_that("clique_intensity_test empty cliques returns 0-row dataframe", {
  setup <- make_clique_fixture()
  empty_cliques <- setup$cliques[0, , drop = FALSE]
  result <- clique_intensity_test(
    empty_cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L
  )

  expect_equal(nrow(result), 0)
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test n_perm = 1 gives valid output", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 1L, seed = 42L
  )

  expect_equal(nrow(result), nrow(setup$cliques))
  expect_equal(result$n_perm[1], 1L)
  # With n_perm=1, sd is NA (need >= 2 values), so z_score should be NA
  expect_true(all(is.na(result$z_score) | is.finite(result$z_score)))
})


test_that("clique_intensity_test n_perm = 0 returns empty", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 0L
  )
  expect_equal(nrow(result), 0)
})


test_that(
  "clique_intensity_test observed_intensity matches compute_clique_edge_stats",
  {
    setup <- make_clique_fixture()
    if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

    result <- clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 2L, seed = 42L
    )
    stats <- rcomplex:::compute_clique_edge_stats(
      setup$cliques, setup$edges, setup$target_species
    )
    expect_equal(result$observed_intensity, stats$intensity)
  }
)


test_that(
  "alternative = 'less' runs without error and returns valid p-values",
  {
    setup <- make_clique_fixture()
    if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

    r_less <- clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 5L, seed = 42L, alternative = "less"
    )

    expect_true(is.data.frame(r_less))
    expect_true("p_value" %in% names(r_less))
    # p-values should be in [0, 1] or NA
    valid <- !is.na(r_less$p_value)
    if (any(valid)) {
      expect_true(all(r_less$p_value[valid] >= 0 &
                        r_less$p_value[valid] <= 1))
    }
  }
)


test_that("clique_intensity_test works with 3+ species", {
  setup <- make_clique_fixture_3sp()
  if (nrow(setup$cliques) == 0) skip("No 3-species cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 3L, seed = 42L, min_species = 2L
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
  expect_true("n_matched" %in% names(result))
})


test_that("clique_intensity_test validates inputs", {
  setup <- make_clique_fixture()

  expect_error(
    clique_intensity_test(
      "not_df", setup$target_species,
      setup$networks, setup$orthologs
    ),
    "data frame"
  )
  expect_error(
    clique_intensity_test(
      data.frame(x = 1), setup$target_species,
      setup$networks, setup$orthologs
    ),
    "data frame"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, "SP_A",
      setup$networks, setup$orthologs
    ),
    "at least 2"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, setup$target_species,
      list(), setup$orthologs
    ),
    "named list"
  )
  expect_error(
    clique_intensity_test(
      setup$cliques, setup$target_species,
      setup$networks, data.frame(x = 1)
    ),
    "Species1, Species2, hog"
  )
})


test_that("max_missing_edges is forwarded to find_cliques", {
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_intensity_test(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs,
    n_perm = 2L, seed = 1L,
    max_missing_edges = 1L
  )
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(setup$cliques))
})


test_that("pval_combine/pi0_method reach the baseline and null reruns", {
  # Baseline cliques built with pval_combine = "min" (one conserved edge).
  # Observed intensity is a Jaccard-index percentile (#11), and Jaccard
  # does not depend on how directional q-values are combined, so it cannot
  # show whether pval_combine was forwarded: it matches under either value.
  # The null reruns can. Under "min" a permutation that keeps the A1-B1
  # hub mapping rebuilds the clique and counts as a match; under "max" the
  # diluted B-side direction never reaches alpha, so no permutation can
  # produce a conserved edge and nothing matches.
  setup <- make_asym_clique_fixture()
  expect_equal(nrow(setup$cliques), 1L)

  run <- function(pc) {
    clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 20L, seed = 7L,
      pval_combine = pc, pi0_method = "none"
    )
  }
  r_min <- run("min")
  r_max <- run("max")

  stats_min <- rcomplex:::compute_clique_edge_stats(
    setup$cliques, setup$edges_min, setup$target_species
  )
  expect_equal(r_min$observed_intensity, stats_min$intensity)
  expect_gt(r_min$n_matched, 0L)
  expect_identical(r_max$n_matched, 0L)
})


test_that("a pre-filtered edges argument warns", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  setup <- make_clique_fixture()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")
  conserved <- setup$edges[setup$edges$type == "conserved", , drop = FALSE]
  expect_gt(nrow(conserved), 0L)

  expect_warning(
    clique_intensity_test(
      setup$cliques, setup$target_species, setup$networks,
      setup$orthologs,
      n_perm = 1L, seed = 1L, edges = conserved
    ),
    class = "rcomplex_prefiltered_edges"
  )
})



test_that("null intensities are the permuted runs' own Jaccard weights", {
  # A calibration test (strong fixture rejects, null fixture does not) is
  # not feasible at this scale: the null reshuffles every ortholog, so
  # permuted runs rebuild a baseline clique only by chance and small
  # fixtures match none. What can be pinned exactly is the null itself.
  # With pi0_method = "none" the only draw in a permutation is the
  # ortholog shuffle, so replaying the loop under the same seed rebuilds
  # every permuted table; null_mean must equal the mean intensity of the
  # matched permuted cliques, weighted over each permuted run's own table.
  setup <- make_asym_clique_fixture()
  sp <- setup$target_species
  n_perm <- 20L
  res <- clique_intensity_test(
    setup$cliques, sp, setup$networks, setup$orthologs,
    n_perm = n_perm, seed = 7L,
    pval_combine = "min", pi0_method = "none"
  )
  expect_gt(res$n_matched, 0L)

  set.seed(7L)
  null_int <- numeric(0)
  for (b in seq_len(n_perm)) {
    sh <- setup$orthologs
    sh$Species2 <- sample(sh$Species2)
    e_p <- find_coexpressologs(setup$networks, sh,
      method = "analytical", pval_combine = "min", pi0_method = "none"
    )
    if (nrow(e_p) == 0L) next
    cl_p <- find_cliques(e_p, sp, min_species = length(sp))
    hit <- cl_p[cl_p$hog == setup$cliques$hog[1L], , drop = FALSE]
    if (nrow(hit) == 0L) next
    st <- rcomplex:::compute_clique_edge_stats(hit[1L, ], e_p, sp)
    if (!is.na(st$intensity)) null_int <- c(null_int, st$intensity)
  }
  expect_length(null_int, res$n_matched)
  expect_equal(res$null_mean, mean(null_int), tolerance = 1e-12)
})


test_that("within_hog null matches cliques where global cannot", {
  # The global shuffle moves genes out of their own HOG, so a permuted
  # run almost never rebuilds the observed clique's HOG and nothing can
  # be matched -- measured as 0 of 204 cliques on the eight-species
  # Pooideae run, which left every z undefined. The within-HOG shuffle
  # keeps the grouping, so the same cliques do match.
  #
  # make_clique_fixture() cannot show this: it puts exactly one gene in
  # every HOG, which makes the within-HOG shuffle a no-op and the null a
  # point mass. This fixture pairs genes two to a HOG instead, so the
  # shuffle has something to permute.
  n <- 20L
  ga <- paste0("A", seq_len(n))
  gb <- paste0("B", seq_len(n))
  make_net <- function(genes) {
    m <- matrix(0, n, n, dimnames = list(genes, genes))
    for (i in 2:10) m[1, i] <- m[i, 1] <- 10
    m[1, 11] <- m[11, 1] <- 3
    for (i in 13:15) m[12, i] <- m[i, 12] <- 4
    m
  }
  networks <- list(
    SP_A = list(network = make_net(ga), threshold = 2),
    SP_B = list(network = make_net(gb), threshold = 2)
  )
  sp <- c("SP_A", "SP_B")
  orthologs <- data.frame(
    Species1 = ga, Species2 = gb,
    hog = paste0("HOG", ceiling(seq_len(n) / 2)),
    stringsAsFactors = FALSE
  )
  edges <- find_coexpressologs(networks, orthologs,
    method = "analytical", pi0_method = "storey"
  )
  cliques <- find_cliques(edges, sp, min_species = 2L)
  expect_gt(nrow(cliques), 0L)

  run <- function(nm) {
    clique_intensity_test(cliques, sp, networks, orthologs,
      n_perm = 30L, seed = 11L, null_model = nm, pi0_method = "storey"
    )
  }
  global <- run("global")
  within <- run("within_hog")

  # global destroys the HOGs: nothing matches, so no z is computable
  expect_true(all(global$n_matched == 0L))
  expect_true(all(is.na(global$z_score)))
  # within_hog keeps them: every clique matches and gets a usable z
  expect_true(all(within$n_matched > 0L))
  expect_true(all(!is.na(within$z_score)))
})


test_that("within_hog shuffle stays inside one species pair", {
  # A stacked all-pairs ortholog table puts A-B, A-C and B-C rows in the
  # same HOG. Grouping the shuffle by hog alone moves a C gene into an
  # A-B row; compare_neighborhoods() drops it by network membership, so
  # the mapping vanishes silently (measured: 33% of rows on a stacked
  # 3-species table) and multi-species cliques stop rebuilding.
  #
  # make_clique_fixture_3sp() is exactly that shape, one gene per species
  # per HOG, so every (hog, species-pair) group holds a single row and a
  # correctly partitioned shuffle is the identity: every permutation
  # reproduces the observed table and matches, with no null spread. That
  # is what pins the partition -- under the hog-only shuffle rows are
  # lost and the permutations cannot all match.
  fx <- make_clique_fixture_3sp()
  n_perm <- 20L
  res <- clique_intensity_test(fx$cliques, fx$target_species, fx$networks,
    fx$orthologs,
    n_perm = n_perm, seed = 3L, null_model = "within_hog",
    pi0_method = "storey"
  )
  # nothing dropped: every permutation rebuilds the observed clique
  expect_true(all(res$n_matched == n_perm))
  # identity shuffle, so the null has no spread and z is undefined
  expect_true(all(res$null_sd == 0))
  expect_true(all(is.na(res$z_score)))

  # the global shuffle destroys the HOGs, so nothing matches at all
  glb <- clique_intensity_test(fx$cliques, fx$target_species, fx$networks,
    fx$orthologs,
    n_perm = n_perm, seed = 3L, null_model = "global",
    pi0_method = "storey"
  )
  expect_true(all(glb$n_matched == 0L))
})


test_that("matched_edges null needs no networks or orthologs", {
  # The point of this null: it never re-runs find_coexpressologs(), so a
  # caller with an edge table does not have to hold every species
  # network in memory. Passing neither must still work.
  fx <- make_clique_fixture_3sp()
  res <- clique_intensity_test(fx$cliques, fx$target_species,
    edges = fx$edges, n_perm = 200L, seed = 1L,
    null_model = "matched_edges"
  )
  expect_equal(nrow(res), nrow(fx$cliques))
  # every draw is usable -- nothing has to be rebuilt, so nothing is lost
  expect_true(all(res$n_matched == 200L))
  # this fixture has one conserved edge per species pair, so each pool is
  # a single constant and the resample reproduces the observed value: a
  # spreadless null here is correct, not a degenerate one
  expect_true(all(res$null_sd == 0))
})


test_that("matched_edges null has spread when the weights do", {
  # This null reads only `edges`, so the fixture is an edge table --
  # no networks, no find_coexpressologs(). Varying effect_size is what
  # gives each species pair's pool something to resample over; the "ns"
  # rows keep the table unfiltered and must stay out of the pools.
  sps <- c("SP_A", "SP_B", "SP_C")
  n_hog <- 30L
  hogs <- paste0("HOG", seq_len(n_hog))
  set.seed(4)
  edges <- do.call(rbind, lapply(
    utils::combn(sps, 2, simplify = FALSE),
    function(p) {
      rbind(
        data.frame(
          species1 = p[1], species2 = p[2],
          gene1 = paste0(p[1], "_", hogs), gene2 = paste0(p[2], "_", hogs),
          hog = hogs, effect_size = stats::runif(n_hog, 1.5, 12),
          type = "conserved", stringsAsFactors = FALSE
        ),
        data.frame(
          species1 = p[1], species2 = p[2],
          gene1 = paste0(p[1], "_ns", seq_len(5L)),
          gene2 = paste0(p[2], "_ns", seq_len(5L)),
          hog = paste0("NS", seq_len(5L)),
          effect_size = stats::runif(5L, 0.1, 0.9),
          type = "ns", stringsAsFactors = FALSE
        )
      )
    }
  ))
  cliques <- data.frame(hog = hogs, stringsAsFactors = FALSE)
  for (s in sps) cliques[[s]] <- paste0(s, "_", hogs)

  res <- clique_intensity_test(cliques, sps,
    edges = edges, n_perm = 500L, seed = 3L,
    null_model = "matched_edges"
  )
  expect_equal(nrow(res), n_hog)
  expect_true(all(res$n_matched == 500L))
  # a real null: spread, and a usable z for every clique
  expect_true(all(res$null_sd > 0))
  expect_true(all(is.finite(res$z_score)))
  expect_true(all(res$p_value > 0 & res$p_value <= 1))
})


test_that("matched_edges pools honour a multi-value edge_type", {
  # edge_type is vector-valued across the package (find_cliques(),
  # clique_stability()) and the documented contract is
  # `type %in% edge_type`. Filtering the pool with `==` instead recycles
  # element-wise, which discards roughly half of EVERY type and so
  # leaves the pool's composition -- and its mean -- almost unchanged.
  # That is why a mean-based assertion cannot see the bug. What it does
  # do is warn, but only when the row count is not a multiple of
  # length(edge_type), so the row count here is deliberately odd.
  sps <- c("SP_A", "SP_B")
  n_con <- 20L
  n_div <- 21L
  mkrows <- function(type, eff, tag, k) {
    data.frame(
      species1 = "SP_A", species2 = "SP_B",
      gene1 = paste0("A_", tag, seq_len(k)),
      gene2 = paste0("B_", tag, seq_len(k)),
      hog = paste0(tag, seq_len(k)),
      effect_size = eff, type = type, stringsAsFactors = FALSE
    )
  }
  set.seed(9)
  edges <- rbind(
    mkrows("conserved", stats::runif(n_con, 6, 12), "C", n_con),
    mkrows("diverged", stats::runif(n_div, 1, 3), "D", n_div)
  )
  expect_true(nrow(edges) %% 2L == 1L)
  cliques <- data.frame(
    hog = paste0("C", seq_len(n_con)),
    SP_A = paste0("A_C", seq_len(n_con)),
    SP_B = paste0("B_C", seq_len(n_con)), stringsAsFactors = FALSE
  )
  # `==` recycles over an odd row count and warns; `%in%` does not
  expect_no_warning(
    res <- clique_intensity_test(cliques, sps,
      edges = edges, n_perm = 300L, seed = 5L,
      null_model = "matched_edges",
      edge_type = c("conserved", "diverged")
    )
  )
  expect_true(all(res$n_matched == 300L))
  expect_true(all(is.finite(res$z_score)))
})


test_that("matched_edges builds edges when none are supplied", {
  # needs_networks keeps the network-driven path reachable for this
  # null: with edges = NULL the function still computes them from
  # networks and orthologs, then resamples. Every other test here hands
  # it an edge table, so this is the one covering that interaction.
  fx <- make_clique_fixture_3sp()
  res <- clique_intensity_test(fx$cliques, fx$target_species,
    fx$networks, fx$orthologs,
    n_perm = 100L, seed = 2L, null_model = "matched_edges"
  )
  expect_equal(nrow(res), nrow(fx$cliques))
  expect_true(all(res$n_matched == 100L))
})


test_that("matched_edges pools can be matched on clique size", {
  # The pair-only pool mixes edges from every clique size, and a clique
  # exists only because all its edges passed together -- so big cliques
  # are built from stronger edges and beat a pool dominated by small
  # ones. This fixture makes that concrete: for the A-B pair, 3-clique
  # edges are strong and 2-clique edges are weak, so the size-matched
  # pool for a 3-clique excludes the weak half and its null sits higher.
  sps <- c("SP_A", "SP_B", "SP_C")
  n <- 20L
  h3 <- paste0("T", seq_len(n))
  h2 <- paste0("D", seq_len(n))
  set.seed(8)
  mk <- function(s1, s2, hogs, eff) {
    data.frame(
      species1 = s1, species2 = s2,
      gene1 = paste0(s1, "_", hogs), gene2 = paste0(s2, "_", hogs),
      hog = hogs, effect_size = eff, type = "conserved",
      stringsAsFactors = FALSE
    )
  }
  edges <- rbind(
    mk("SP_A", "SP_B", h3, stats::runif(n, 9, 12)),
    mk("SP_A", "SP_C", h3, stats::runif(n, 9, 12)),
    mk("SP_B", "SP_C", h3, stats::runif(n, 9, 12)),
    mk("SP_A", "SP_B", h2, stats::runif(n, 1.2, 2.0))
  )
  mkclq <- function(hogs, third) {
    data.frame(
      hog = hogs,
      SP_A = paste0("SP_A_", hogs),
      SP_B = paste0("SP_B_", hogs),
      SP_C = third,
      stringsAsFactors = FALSE
    )
  }
  cliques <- rbind(
    mkclq(h3, paste0("SP_C_", h3)),
    mkclq(h2, NA_character_)
  )
  run <- function(match) {
    clique_intensity_test(cliques, sps,
      edges = edges, n_perm = 400L, seed = 4L,
      null_model = "matched_edges", match_clique_size = match
    )
  }
  matched <- run(TRUE)
  pooled <- run(FALSE)
  is3 <- !is.na(cliques$SP_C)

  # every clique scored either way
  expect_true(all(matched$n_matched == 400L))
  expect_true(all(pooled$n_matched == 400L))
  # the 3-cliques face a stronger null once the weak 2-clique edges are
  # excluded from their pool, so their null mean rises and their z falls
  expect_gt(mean(matched$null_mean[is3]), mean(pooled$null_mean[is3]))
  expect_lt(mean(matched$z_score[is3]), mean(pooled$z_score[is3]))
  # the size-2 cliques only ever had the A-B pair, so matching on size
  # removes the strong 3-clique edges and lowers their null instead
  expect_lt(mean(matched$null_mean[!is3]), mean(pooled$null_mean[!is3]))
})


test_that("min_pool_size rejects a non-integer threshold", {
  # as.integer() truncates, so validating the coerced value would take
  # 1.5 as 1 and quietly change the null instead of rejecting an
  # argument documented as a positive integer. Validation has to happen
  # before the coercion, and only a value that survives it may be used.
  fx <- make_clique_fixture_3sp()
  run <- function(v) {
    clique_intensity_test(fx$cliques, fx$target_species,
      edges = fx$edges, n_perm = 10L, seed = 1L,
      null_model = "matched_edges", min_pool_size = v
    )
  }
  for (bad in list(1.5, 0.9, -1, 0, NA_real_, Inf, c(1, 2), "2")) {
    expect_error(run(bad), "min_pool_size must be a single positive integer")
  }
  # whole numbers pass whether given as double or integer
  expect_s3_class(run(1), "data.frame")
  expect_s3_class(run(2L), "data.frame")
})


test_that("min_pool_size leaves a thin pool unscored", {
  fx <- make_clique_fixture_3sp()
  # this fixture has one conserved edge per species pair, so any
  # threshold above 1 refuses every draw
  thin <- clique_intensity_test(fx$cliques, fx$target_species,
    edges = fx$edges, n_perm = 50L, seed = 1L,
    null_model = "matched_edges", min_pool_size = 5L
  )
  expect_true(all(thin$n_matched == 0L))
  expect_true(all(is.na(thin$z_score)))
  # and the default of 1 scores it, with no spread, as before
  kept <- clique_intensity_test(fx$cliques, fx$target_species,
    edges = fx$edges, n_perm = 50L, seed = 1L,
    null_model = "matched_edges"
  )
  expect_true(all(kept$n_matched == 50L))
  expect_true(all(kept$null_sd == 0))
})


test_that("null_model is validated and defaults to global", {
  fx <- formals(rcomplex:::clique_intensity_test.default)
  expect_identical(eval(fx$null_model)[1], "global")
  expect_setequal(
    eval(fx$null_model),
    c("global", "within_hog", "matched_edges")
  )
})
