# Tests for clique_stability()

# Helper: build 4-species binary trait edge data with known clique structure
#
# 4 species, binary trait:
#   SP_A, SP_B = "annual"
#   SP_C, SP_D = "perennial"
#
# HOG1: annual-exclusive clique (SP_A <-> SP_B, both conserved, low FDR)
#   Edges: A1-B1 (SP_A-SP_B)
#   -> forms 2-species annual-exclusive clique
#   Stability at k=1: removing SP_C or SP_D doesn't affect it
#     -> stability_score = 1.0
#
# HOG2: perennial-exclusive clique (SP_C <-> SP_D)
#   Edges: C1-D1 (SP_C-SP_D)
#   -> forms 2-species perennial-exclusive clique
#
# HOG3: mixed clique (all 4 species connected)
#   Edges: A2-B2, A2-C2, A2-D2, B2-C2, B2-D2, C2-D2
#   -> forms 4-species mixed clique
#   Should be excluded from stability analysis

make_stability_edges_binary <- function() {
  data.frame(
    gene1 = c("A1", "C1",
              "A2", "A2", "A2", "B2", "B2", "C2"),
    gene2 = c("B1", "D1",
              "B2", "C2", "D2", "C2", "D2", "D2"),
    species1 = c("SP_A", "SP_C",
                 "SP_A", "SP_A", "SP_A", "SP_B", "SP_B", "SP_C"),
    species2 = c("SP_B", "SP_D",
                 "SP_B", "SP_C", "SP_D", "SP_C", "SP_D", "SP_D"),
    hog = c("HOG1", "HOG2", rep("HOG3", 6)),
    q.value = c(0.01, 0.02, rep(0.05, 6)),
    effect_size = c(3.0, 2.5, rep(1.5, 6)),
    type = rep("conserved", 8),
    stringsAsFactors = FALSE
  )
}


# Helper: 5-species ternary trait edge data
#
# 5 species, 3 trait levels:
#   SP_A = "tropical"
#   SP_B, SP_C = "temperate"
#   SP_D, SP_E = "arctic"
#
# HOG5: temperate-exclusive (SP_B <-> SP_C)
# HOG6: arctic-exclusive (SP_D <-> SP_E)
# HOG7: mixed (SP_A <-> SP_B, SP_A <-> SP_C -- tropical + temperate)

make_stability_edges_ternary <- function() {
  data.frame(
    gene1 = c("B5", "D6", "A7", "A7"),
    gene2 = c("C5", "E6", "B7", "C7"),
    species1 = c("SP_B", "SP_D", "SP_A", "SP_A"),
    species2 = c("SP_C", "SP_E", "SP_B", "SP_C"),
    hog = c("HOG5", "HOG6", "HOG7", "HOG7"),
    q.value = c(0.01, 0.02, 0.04, 0.04),
    effect_size = c(3.0, 2.5, 1.5, 1.5),
    type = rep("conserved", 4),
    stringsAsFactors = FALSE
  )
}


# ---- Structure & basic tests ----

test_that("clique_stability returns correct structure", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 2L)

  expect_type(result, "list")
  expect_named(result, c("stability", "clique_disruption",
                          "stability_class", "novel_cliques"))

  # stability is a data frame
  expect_s3_class(result$stability, "data.frame")
  expect_true(all(c("clique_idx", "hog", "trait_value", "k",
                     "n_subsets", "n_stable", "stability_score",
                     "sole_rep") %in% names(result$stability)))

  # clique_disruption is a data frame
  expect_s3_class(result$clique_disruption, "data.frame")
  expect_true(all(c("species", "trait_value",
                     "n_cliques_disrupted") %in% names(result$clique_disruption)))

  # stability_class is an integer vector
  expect_true(is.integer(result$stability_class))

  # novel_cliques is an integer
  expect_true(is.integer(result$novel_cliques))
})


test_that("output labels are generic (trait_value, not annual)", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 1L)

  # Column names should be generic
  stab_names <- names(result$stability)
  expect_false(any(grepl("annual|perennial|life_habit", stab_names,
                         ignore.case = TRUE)))

  # trait_value column should contain actual trait labels
  expect_true(nrow(result$stability) > 0)
  expect_true(all(result$stability$trait_value %in% c("annual", "perennial")))
})


test_that("empty input returns correct empty structure", {
  edges <- data.frame(
    gene1 = character(0), gene2 = character(0),
    species1 = character(0), species2 = character(0),
    hog = character(0), q.value = numeric(0),
    effect_size = numeric(0), type = character(0),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual")
  target <- c("SP_A", "SP_B")

  result <- clique_stability(edges, target, trait, max_k = 1L)

  expect_equal(nrow(result$stability), 0)
  expect_true(all(c("clique_idx", "hog", "trait_value", "k",
                     "n_subsets", "n_stable", "stability_score",
                     "sole_rep") %in% names(result$stability)))
  expect_equal(nrow(result$clique_disruption), 0)
  expect_equal(length(result$stability_class), 0)
  expect_equal(result$novel_cliques, 0L)
})


# ---- Binary trait tests ----

test_that("stable exclusive cliques have stability_score = 1.0 at k=1", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 1L)

  # HOG1 annual-exclusive: removing SP_C or SP_D -> still exists, still annual
  # The 2 testable subsets (remove SP_C, remove SP_D) should both be stable
  hog1_k1 <- result$stability[result$stability$hog == "HOG1" &
                                 result$stability$k == 1, ]
  expect_equal(nrow(hog1_k1), 1)
  expect_equal(hog1_k1$stability_score, 1.0)
  expect_equal(hog1_k1$trait_value, "annual")

  # HOG2 perennial-exclusive: removing SP_A or SP_B -> still exists, still perennial
  hog2_k1 <- result$stability[result$stability$hog == "HOG2" &
                                 result$stability$k == 1, ]
  expect_equal(nrow(hog2_k1), 1)
  expect_equal(hog2_k1$stability_score, 1.0)
  expect_equal(hog2_k1$trait_value, "perennial")
})


test_that("mixed cliques are excluded from stability output", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 2L)

  # HOG3 mixed clique should NOT appear in stability results
  expect_false("HOG3" %in% result$stability$hog)
})


test_that("clique_disruption counts species removals correctly at k=1", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 1L)

  # Should have one row per species
  expect_equal(nrow(result$clique_disruption), 4)
  expect_setequal(result$clique_disruption$species, target)

  # HOG1 (SP_A <-> SP_B) and HOG2 (SP_C <-> SP_D) are exactly at min_species=2.
  # Removing a member drops below min_species -> clique is untestable, not disrupted.
  # Removing a non-member keeps the clique testable and stable -> no disruption.
  # So all species should have n_cliques_disrupted = 0.
  for (sp in target) {
    row <- result$clique_disruption[result$clique_disruption$species == sp, ]
    expect_equal(row$n_cliques_disrupted, 0L)
  }
})


test_that("stability_class gives highest stable k", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait, min_species = 2L,
                             max_k = 2L)

  # 2 exclusive cliques (HOG1 annual, HOG2 perennial)
  expect_equal(length(result$stability_class), 2L)
  # Both should be stable at least at k=1
  expect_true(all(result$stability_class >= 1L))
})


# ---- Ternary trait (3 levels) tests ----

test_that("three-level trait works correctly", {
  edges <- make_stability_edges_ternary()
  trait <- c(SP_A = "tropical", SP_B = "temperate", SP_C = "temperate",
             SP_D = "arctic", SP_E = "arctic")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 1L)

  # Both temperate-exclusive and arctic-exclusive cliques should be detected
  expect_true(nrow(result$stability) > 0)
  trait_vals <- unique(result$stability$trait_value)
  expect_true("temperate" %in% trait_vals)
  expect_true("arctic" %in% trait_vals)
  # trait_value should contain actual string labels
  expect_true(all(result$stability$trait_value %in%
                    c("tropical", "temperate", "arctic")))
})


test_that("removing unrelated species preserves trait exclusivity", {
  edges <- make_stability_edges_ternary()
  trait <- c(SP_A = "tropical", SP_B = "temperate", SP_C = "temperate",
             SP_D = "arctic", SP_E = "arctic")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 1L)

  # HOG5 temperate-exclusive (SP_B <-> SP_C)
  # Removing SP_A (tropical), SP_D, or SP_E should keep it stable
  hog5_k1 <- result$stability[result$stability$hog == "HOG5" &
                                 result$stability$k == 1, ]
  expect_equal(nrow(hog5_k1), 1)
  # 3 removable species (SP_A, SP_D, SP_E) that don't affect this clique
  # stability_score should be 1.0 (all testable subsets are stable)
  expect_equal(hog5_k1$stability_score, 1.0)
})


# ---- Edge cases ----

test_that("min_species edge case: clique at minimum size", {
  # 3 species: SP_A, SP_B = "annual", SP_C = "perennial"
  # HOG_MIN: SP_A <-> SP_B (2-species annual-exclusive clique)
  # min_species = 2
  # k=1: removing SP_A -> drops below min_species for this clique -> not testable
  #       removing SP_B -> drops below min_species -> not testable
  #       removing SP_C -> clique still exists (SP_A <-> SP_B) -> testable, stable
  edges <- data.frame(
    gene1 = "A1",
    gene2 = "B1",
    species1 = "SP_A",
    species2 = "SP_B",
    hog = "HOG_MIN",
    q.value = 0.01,
    effect_size = 3.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "perennial")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 1L)

  hog_min <- result$stability[result$stability$hog == "HOG_MIN" &
                                result$stability$k == 1, ]
  expect_equal(nrow(hog_min), 1)
  # Only 1 testable subset (remove SP_C), and it should be stable
  expect_equal(hog_min$n_subsets, 1L)
  expect_equal(hog_min$n_stable, 1L)
  expect_equal(hog_min$stability_score, 1.0)
})


test_that("all mixed cliques returns empty stability", {
  # All cliques span multiple trait groups -> no exclusive cliques
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    q.value = rep(0.01, 3),
    effect_size = rep(2.0, 3),
    type = rep("conserved", 3),
    stringsAsFactors = FALSE
  )
  # All species have different traits -> any multi-species clique is mixed
  trait <- c(SP_A = "annual", SP_B = "perennial", SP_C = "biennial")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 1L)

  expect_equal(nrow(result$stability), 0)
})


test_that("sole_rep flagging works", {
  # 4 species: SP_A = "rare" (only 1 species), SP_B,SP_C,SP_D = "common" (3 species)
  # HOG1: SP_B <-> SP_C (common-exclusive, 2 species) — sole_rep FALSE (3 reps)
  edges <- data.frame(
    gene1 = "B1", gene2 = "C1",
    species1 = "SP_B", species2 = "SP_C",
    hog = "HOG1", q.value = 0.01, effect_size = 2.0,
    type = "conserved", stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "rare", SP_B = "common", SP_C = "common", SP_D = "common")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")
  result <- clique_stability(edges, target, trait, min_species = 2L, max_k = 1L)

  expect_true(nrow(result$stability) > 0)
  # "common" trait has 3 species -> sole_rep = FALSE
  expect_false(any(result$stability$sole_rep))
})


test_that("jaccard_threshold affects matching", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  # With strict threshold (1.0): only exact gene matches count as stable
  result_strict <- clique_stability(edges, target, trait,
                                    min_species = 2L,
                                    max_k = 1L, jaccard_threshold = 1.0)

  # With lenient threshold (0.0): everything matches
  result_lenient <- clique_stability(edges, target, trait,
                                     min_species = 2L,
                                     max_k = 1L, jaccard_threshold = 0.0)

  # Both should produce results
  expect_true(nrow(result_strict$stability) > 0)
  expect_true(nrow(result_lenient$stability) > 0)

  # Lenient should have >= stable cliques as strict
  strict_stable <- sum(result_strict$stability$stability_score >= 1.0)
  lenient_stable <- sum(result_lenient$stability$stability_score >= 1.0)
  expect_true(lenient_stable >= strict_stable)
})


# ---- Parallelism ----

test_that("n_cores=1 and n_cores=2 produce identical results", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result1 <- clique_stability(edges, target, trait,
                               min_species = 2L, max_k = 2L, n_cores = 1L)
  result2 <- clique_stability(edges, target, trait,
                               min_species = 2L, max_k = 2L, n_cores = 2L)

  # Stability results should be identical (deterministic, no RNG)
  expect_equal(result1$stability, result2$stability)
  expect_equal(result1$clique_disruption, result2$clique_disruption)
  expect_equal(result1$stability_class, result2$stability_class)
  expect_equal(result1$novel_cliques, result2$novel_cliques)
})


# ---- Input validation ----

test_that("missing required columns in edges raises error", {
  edges <- data.frame(
    gene1 = "A1", gene2 = "B1",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", effect_size = 2.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual")
  target <- c("SP_A", "SP_B")

  # Missing q.value column
  expect_error(clique_stability(edges, target, trait))
})


test_that("species_trait missing species raises error", {
  edges <- make_stability_edges_binary()
  # Missing SP_D from trait vector
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  expect_error(clique_stability(edges, target, trait))
})


test_that("max_k validation", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  # max_k < 1 -> error
  expect_error(clique_stability(edges, target, trait, max_k = 0L))

  # max_k >= length(target_species) -> error
  expect_error(clique_stability(edges, target, trait, max_k = 4L))
})


test_that("full_cliques parameter accepts precomputed cliques", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  # First compute cliques, then pass them in
  cliques <- find_cliques(edges, target, min_species = 2L)

  result_auto <- clique_stability(edges, target, trait,
                                  min_species = 2L, max_k = 1L)
  result_pre <- clique_stability(edges, target, trait,
                                 full_cliques = cliques,
                                 min_species = 2L, max_k = 1L)

  expect_equal(result_auto$stability, result_pre$stability)
})


test_that("max_genes_per_sp parameter is passed through", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 1L,
                             max_genes_per_sp = 5L)
  expect_type(result, "list")
  expect_true(nrow(result$stability) > 0)
})


test_that("edge_type filtering works in stability analysis", {
  edges <- make_stability_edges_binary()
  # Change HOG1 edge to diverged
  edges$type[1] <- "diverged"

  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  # With default edge_type="conserved", HOG1 clique (A1-B1) is gone
  # because that edge is now "diverged"
  result_con <- clique_stability(edges, target, trait,
                                 min_species = 2L, max_k = 1L)

  # With both types, HOG1 clique should be found
  result_both <- clique_stability(edges, target, trait,
                                  min_species = 2L, max_k = 1L,
                                  edge_type = c("conserved", "diverged"))

  # Conserved-only: HOG1 gone, only HOG2 perennial remains
  expect_false("HOG1" %in% result_con$stability$hog)
  expect_true("HOG2" %in% result_con$stability$hog)

  # Both types: HOG1 and HOG2 both present
  expect_true("HOG1" %in% result_both$stability$hog)
  expect_true("HOG2" %in% result_both$stability$hog)
})


# ---- all_species != target_species tests ----

test_that("all_species parameter draws subsets from full universe", {
  # 2 target species (annuals), 4 all_species (2 annual + 2 perennial)
  edges <- data.frame(
    gene1 = "A1", gene2 = "B1",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", q.value = 0.01, effect_size = 3.0,
    type = "conserved", stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B")
  all_sp <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait,
                             all_species = all_sp,
                             min_species = 2L, max_k = 1L)

  # k=1: C(4,1) = 4 subsets
  # Remove SP_A: clique untestable (1 member)
  # Remove SP_B: clique untestable (1 member)
  # Remove SP_C: clique testable, stable (both annuals present)
  # Remove SP_D: clique testable, stable
  hog1 <- result$stability[result$stability$hog == "HOG1" &
                             result$stability$k == 1, ]
  expect_equal(nrow(hog1), 1)
  expect_equal(hog1$n_subsets, 2L)
  expect_equal(hog1$stability_score, 1.0)

  # Disruption has 4 rows (one per all_species)
  expect_equal(nrow(result$clique_disruption), 4)
  expect_setequal(result$clique_disruption$species, all_sp)
})


test_that("non-target removal trivially preserves cliques with static edges", {
  # 3 target species in 5 all_species
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    q.value = rep(0.01, 3),
    effect_size = rep(2.0, 3),
    type = rep("conserved", 3),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "annual",
             SP_D = "perennial", SP_E = "perennial")
  target <- c("SP_A", "SP_B", "SP_C")
  all_sp <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E")

  result <- clique_stability(edges, target, trait,
                             all_species = all_sp,
                             min_species = 2L, max_k = 1L)

  # k=1: C(5,1) = 5 subsets
  # Remove SP_D or SP_E: all 3 targets present, clique survives -> stable
  # Remove SP_A/SP_B/SP_C: 2 targets remain, 2-species clique survives -> stable
  hog1 <- result$stability[result$stability$hog == "HOG1" &
                             result$stability$k == 1, ]
  expect_equal(nrow(hog1), 1)
  expect_equal(hog1$n_subsets, 5L)
  expect_equal(hog1$stability_score, 1.0)
})


test_that("all_species defaults to target_species (backward compat)", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  # Explicit all_species = target_species should match default
  result_default <- clique_stability(edges, target, trait,
                                     min_species = 2L, max_k = 1L)
  result_explicit <- clique_stability(edges, target, trait,
                                      all_species = target,
                                      min_species = 2L, max_k = 1L)

  expect_equal(result_default$stability, result_explicit$stability)
  expect_equal(result_default$stability_class, result_explicit$stability_class)
})


test_that("stability_class is non-zero for strong cliques", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_stability(edges, target, trait,
                             min_species = 2L, max_k = 2L)

  # Both exclusive cliques (HOG1 annual, HOG2 perennial) should have
  # stability_class > 0 (not universally zero)
  expect_true(all(result$stability_class > 0L))
})


test_that("target_species must be subset of all_species", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual")
  target <- c("SP_A", "SP_B")

  expect_error(
    clique_stability(edges, target, trait,
                     all_species = c("SP_A"),  # missing SP_B
                     min_species = 2L, max_k = 1L),
    "target_species must be a subset of all_species"
  )
})


test_that("max_k validated against all_species length", {
  edges <- make_stability_edges_binary()
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial",
             SP_E = "perennial", SP_F = "perennial")
  target <- c("SP_A", "SP_B")
  all_sp <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E", "SP_F")

  # max_k = 6 >= length(all_species) = 6 -> error
  expect_error(
    clique_stability(edges, target, trait,
                     all_species = all_sp,
                     min_species = 2L, max_k = 6L),
    "max_k must be"
  )

  # max_k = 5 < 6 -> should work
  expect_no_error(
    clique_stability(edges, target, trait,
                     all_species = all_sp,
                     min_species = 2L, max_k = 5L)
  )
})
