# Tests for clique_threshold_sweep()

# Helper: build 2-species test setup with enough genes for significant results.
# 20 genes per species, dense neighborhoods for HOG1 (A1-B1), sparse for HOG2.
make_sweep_setup <- function() {
  n <- 20
  ga <- paste0("A", seq_len(n))
  gb <- paste0("B", seq_len(n))

  # SP_A network: A1 strongly connected to A2..A10; A11 weakly
  mat_a <- matrix(0, n, n, dimnames = list(ga, ga))
  for (i in 2:10) {
    mat_a[1, i] <- mat_a[i, 1] <- 10  # strong
  }
  mat_a[1, 11] <- mat_a[11, 1] <- 3   # moderate
  # A12 connected to A13..A15
  for (i in 13:15) {
    mat_a[12, i] <- mat_a[i, 12] <- 4
  }

  # SP_B network: mirror structure
  mat_b <- matrix(0, n, n, dimnames = list(gb, gb))
  for (i in 2:10) {
    mat_b[1, i] <- mat_b[i, 1] <- 10
  }
  mat_b[1, 11] <- mat_b[11, 1] <- 3
  for (i in 13:15) {
    mat_b[12, i] <- mat_b[i, 12] <- 4
  }

  networks <- list(
    SP_A = list(network = mat_a, threshold = 2.0),
    SP_B = list(network = mat_b, threshold = 2.0)
  )

  # 1:1 orthologs
  orthologs <- data.frame(
    Species1 = ga, Species2 = gb,
    hog = paste0("HOG", seq_len(n)),
    stringsAsFactors = FALSE
  )

  # Run baseline pipeline
  comparison <- compare_neighborhoods(networks$SP_A, networks$SP_B, orthologs)
  summary <- summarize_comparison(comparison)
  edges <- comparison_to_edges(summary$results, "SP_A", "SP_B")
  baseline <- find_cliques(edges, c("SP_A", "SP_B"), min_species = 2L)

  list(networks = networks, orthologs = orthologs,
       cliques = baseline, target_species = c("SP_A", "SP_B"),
       edges = edges)
}


test_that("clique_threshold_sweep returns correct structure", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(2, 5))

  expect_true(is.list(result))
  expect_true(all(c("survival", "sweep_cliques", "sweep_edges",
                     "persistence") %in% names(result)))

  surv <- result$survival
  expect_true(all(c("clique_idx", "hog", "multiplier", "survived",
                     "jaccard", "n_species_orig", "n_species_new") %in%
                    names(surv)))
  # +1 for the injected multiplier=1.0 baseline rows
  expect_equal(nrow(surv), nrow(setup$cliques) * 3)
  # clique_idx should be 1-based
  expect_true(all(surv$clique_idx >= 1L))
})


test_that("clique_threshold_sweep survival decreases at stricter thresholds", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(1.5, 100))

  surv <- result$survival
  survive_low <- mean(surv$survived[surv$multiplier == 1.5])
  survive_high <- mean(surv$survived[surv$multiplier == 100])
  expect_true(survive_high <= survive_low)
})


test_that("clique_threshold_sweep sweep_cliques keys match multipliers", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  mults <- c(2, 5)
  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = mults)

  expect_equal(sort(names(result$sweep_cliques)), sort(as.character(mults)))
  expect_equal(sort(names(result$sweep_edges)), sort(as.character(mults)))
})


test_that("clique_threshold_sweep with empty cliques returns empty", {
  setup <- make_sweep_setup()
  empty <- setup$cliques[0, , drop = FALSE]

  result <- clique_threshold_sweep(
    empty, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(2))

  expect_equal(nrow(result$survival), 0)
  expect_equal(nrow(result$persistence), 0)
  expect_true(all(c("clique_idx", "hog", "birth", "death", "persistence") %in%
                    names(result$persistence)))
})


test_that("clique_threshold_sweep with empty multipliers returns empty", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = numeric(0))

  expect_equal(nrow(result$survival), 0)
  expect_equal(length(result$sweep_cliques), 0)
  expect_equal(nrow(result$persistence), 0)
})


test_that("clique_threshold_sweep validates inputs", {
  setup <- make_sweep_setup()

  expect_error(
    clique_threshold_sweep("bad", setup$target_species, setup$networks,
                            setup$orthologs),
    "cliques must be")
  expect_error(
    clique_threshold_sweep(setup$cliques, c("SP_A"), setup$networks,
                            setup$orthologs),
    "at least 2 species")
  expect_error(
    clique_threshold_sweep(setup$cliques, setup$target_species, list(),
                            setup$orthologs),
    "networks must be a named list")
  expect_error(
    clique_threshold_sweep(setup$cliques, setup$target_species,
                            list(SP_A = list()), setup$orthologs),
    "networks missing")
  expect_error(
    clique_threshold_sweep(setup$cliques, setup$target_species,
                            setup$networks, data.frame(x = 1)),
    "orthologs must have columns")
})


test_that("jaccard_clique_match computes per-species-slot Jaccard", {
  row1 <- data.frame(SP_A = "A1", SP_B = "B1", SP_C = NA_character_,
                     stringsAsFactors = FALSE)
  row2 <- data.frame(SP_A = "A1", SP_B = "B2", SP_C = "C1",
                     stringsAsFactors = FALSE)

  # A1 matches, B1!=B2, C: only in row2 -> intersect=1, union=3
  jac <- rcomplex:::jaccard_clique_match(row1, row2,
                                          c("SP_A", "SP_B", "SP_C"))
  expect_equal(jac, 1 / 3, tolerance = 1e-10)

  # Identical rows -> 1.0
  jac2 <- rcomplex:::jaccard_clique_match(row1, row1,
                                           c("SP_A", "SP_B", "SP_C"))
  expect_equal(jac2, 1.0)
})


test_that("persistence dataframe has correct columns and structure", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(2, 5))

  persist <- result$persistence
  expect_true(is.data.frame(persist))
  expect_true(all(c("clique_idx", "hog", "birth", "death", "persistence") %in%
                    names(persist)))
  # One row per baseline clique
  expect_equal(nrow(persist), nrow(setup$cliques))
  # clique_idx should be 1-based and cover all baseline cliques
  expect_equal(sort(persist$clique_idx), seq_len(nrow(setup$cliques)))
})


test_that("persistence birth = 1.0 for all baseline cliques", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(2, 5))

  persist <- result$persistence
  # All baseline cliques exist at multiplier=1.0, so birth = 1.0
  expect_true(all(persist$birth == 1.0))
})


test_that("persistence = death - birth where death is not NA", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(1.5, 2, 5, 100))

  persist <- result$persistence
  has_death <- !is.na(persist$death)
  if (any(has_death)) {
    expect_equal(persist$persistence[has_death],
                 persist$death[has_death] - persist$birth[has_death])
  }
  # Where death is NA, persistence should also be NA
  no_death <- is.na(persist$death)
  if (any(no_death)) {
    expect_true(all(is.na(persist$persistence[no_death])))
  }
})


test_that("survival dataframe includes multiplier=1.0 baseline rows", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(2, 5))

  surv <- result$survival
  baseline <- surv[surv$multiplier == 1.0, , drop = FALSE]
  # One baseline row per clique

  expect_equal(nrow(baseline), nrow(setup$cliques))
  # All baseline rows survived with jaccard = 1.0
  expect_true(all(baseline$survived))
  expect_true(all(baseline$jaccard == 1.0))
})


test_that("cliques surviving all multipliers have death = NA, persistence = NA", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  # Use mild multiplier to maximize survival
  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(1.5))

  persist <- result$persistence
  surv <- result$survival
  # Check cliques that survived at all tested multipliers (1.0 and 1.5)
  for (i in seq_len(nrow(persist))) {
    ci_surv <- surv[surv$clique_idx == persist$clique_idx[i], , drop = FALSE]
    if (all(ci_surv$survived)) {
      expect_true(is.na(persist$death[i]))
      expect_true(is.na(persist$persistence[i]))
    }
  }
})


test_that("cliques dying early have smaller persistence than long-lived ones", {
  setup <- make_sweep_setup()
  if (nrow(setup$cliques) == 0) skip("No baseline cliques found")

  result <- clique_threshold_sweep(
    setup$cliques, setup$target_species, setup$networks,
    setup$orthologs, multipliers = c(1.5, 2, 5, 100))

  persist <- result$persistence
  # Among cliques with finite persistence, death should be > birth
  has_death <- !is.na(persist$death)
  if (any(has_death)) {
    expect_true(all(persist$death[has_death] > persist$birth[has_death]))
    expect_true(all(persist$persistence[has_death] > 0))
  }
})
