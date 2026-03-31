# Tests for classify_cliques()

# Helper: build edges for 4-species, binary trait (annual/perennial)
make_classify_edges <- function() {
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")

  # HOG1: ALL 6 edges conserved (complete)
  hog1 <- data.frame(
    gene1 = c("A1", "A1", "A1", "B1", "B1", "C1"),
    gene2 = c("B1", "C1", "D1", "C1", "D1", "D1"),
    species1 = c("SP_A", "SP_A", "SP_A", "SP_B", "SP_B", "SP_C"),
    species2 = c("SP_B", "SP_C", "SP_D", "SP_C", "SP_D", "SP_D"),
    hog = "HOG1",
    q.value = rep(0.01, 6),
    effect_size = rep(3.0, 6),
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # HOG2: 3-species clique A,B,C (partial — D missing)
  hog2 <- data.frame(
    gene1 = c("A2", "A2", "B2"),
    gene2 = c("B2", "C2", "C2"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "HOG2",
    q.value = rep(0.02, 3),
    effect_size = rep(2.5, 3),
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # HOG3: annual clique (A3-B3) + perennial clique (C3-D3), no cross-group
  hog3 <- data.frame(
    gene1 = c("A3", "C3", "A3"),
    gene2 = c("B3", "D3", "C3"),
    species1 = c("SP_A", "SP_C", "SP_A"),
    species2 = c("SP_B", "SP_D", "SP_C"),
    hog = "HOG3",
    q.value = c(0.01, 0.01, 0.5),
    effect_size = c(3.0, 3.0, 0.8),
    type = c("conserved", "conserved", "ns"),
    stringsAsFactors = FALSE
  )

  # HOG4: annual-only clique (A4-B4), no perennial
  hog4 <- data.frame(
    gene1 = "A4", gene2 = "B4",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG4",
    q.value = 0.01, effect_size = 2.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # HOG5: no conserved edges (unclassified)
  hog5 <- data.frame(
    gene1 = "A5", gene2 = "B5",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG5",
    q.value = 0.8, effect_size = 0.5,
    type = "ns",
    stringsAsFactors = FALSE
  )

  edges <- rbind(hog1, hog2, hog3, hog4, hog5)
  list(edges = edges, target = target, trait = trait)
}


test_that("classify_cliques returns correct structure", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  expect_true(is.data.frame(result))
  expected_cols <- c("hog", "classification", "n_species", "best_mean_q",
                     "trait_groups", "stability_class", "persistence", "robust")
  expect_true(all(expected_cols %in% names(result)))
})


test_that("classify_cliques identifies complete cliques", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog1 <- result[result$hog == "HOG1", ]
  expect_equal(hog1$classification, "complete")
  expect_equal(hog1$n_species, 4L)
})


test_that("classify_cliques identifies partial cliques", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog2 <- result[result$hog == "HOG2", ]
  expect_equal(hog2$classification, "partial")
  expect_equal(hog2$n_species, 3L)
})


test_that("classify_cliques identifies differentiated HOGs", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog3 <- result[result$hog == "HOG3", ]
  expect_equal(hog3$classification, "differentiated")
  expect_true(grepl("annual", hog3$trait_groups))
  expect_true(grepl("perennial", hog3$trait_groups))
})


test_that("classify_cliques identifies trait-specific HOGs", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog4 <- result[result$hog == "HOG4", ]
  expect_equal(hog4$classification, "trait_specific")
  expect_equal(hog4$trait_groups, "annual")
})


test_that("classify_cliques identifies unclassified HOGs", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog5 <- result[result$hog == "HOG5", ]
  expect_equal(hog5$classification, "unclassified")
  expect_true(is.na(hog5$n_species))
})


test_that("every HOG appears exactly once", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  all_hogs <- unique(setup$edges$hog)
  expect_equal(sort(result$hog), sort(all_hogs))
  expect_equal(nrow(result), length(all_hogs))
})


test_that("cross-group conserved edge blocks differentiated", {
  setup <- make_classify_edges()
  # Modify HOG3: make the A3-C3 edge conserved (cross-group)
  setup$edges$type[setup$edges$hog == "HOG3" &
                     setup$edges$gene1 == "A3" &
                     setup$edges$gene2 == "C3"] <- "conserved"

  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  hog3 <- result[result$hog == "HOG3", ]
  # With a cross-group conserved edge, HOG3 should NOT be differentiated
  # It would be partial (A3-B3-C3 form a 3-species clique) or similar
  expect_true(hog3$classification != "differentiated")
})


test_that("ternary trait works for differentiated", {
  target <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E", "SP_F")
  trait <- c(SP_A = "warm", SP_B = "warm",
             SP_C = "cold", SP_D = "cold",
             SP_E = "arid", SP_F = "arid")

  # HOG1: warm clique + cold clique + arid clique, no cross-group
  edges <- data.frame(
    gene1 = c("A1", "C1", "E1"),
    gene2 = c("B1", "D1", "F1"),
    species1 = c("SP_A", "SP_C", "SP_E"),
    species2 = c("SP_B", "SP_D", "SP_F"),
    hog = "HOG1",
    q.value = rep(0.01, 3),
    effect_size = rep(3.0, 3),
    type = "conserved",
    stringsAsFactors = FALSE
  )

  result <- classify_cliques(edges, target, trait)
  expect_equal(result$classification, "differentiated")
  # All 3 groups should be listed
  groups <- strsplit(result$trait_groups, ",")[[1]]
  expect_equal(length(groups), 3)
})


test_that("classify_cliques validates inputs", {
  setup <- make_classify_edges()

  expect_error(
    classify_cliques(data.frame(x = 1), setup$target, setup$trait),
    "edges missing required columns")

  expect_error(
    classify_cliques(setup$edges, c("SP_A"), setup$trait),
    "at least 2 species")

  expect_error(
    classify_cliques(setup$edges, setup$target, c("a", "b")),
    "species_trait must be a named")

  expect_error(
    classify_cliques(setup$edges, setup$target,
                     c(SP_A = "x", SP_B = "y")),
    "species_trait missing entries")
})


test_that("classify_cliques with empty edges returns empty", {
  setup <- make_classify_edges()
  empty <- setup$edges[0, , drop = FALSE]

  result <- classify_cliques(empty, setup$target, setup$trait)
  expect_equal(nrow(result), 0)
})


test_that("stability annotation populates stability_class", {
  setup <- make_classify_edges()

  # Mock stability output
  stab <- list(
    stability = data.frame(
      clique_idx = 0L, hog = "HOG1",
      trait_value = "annual", k = 1L,
      n_subsets = 2L, n_stable = 2L,
      stability_score = 1.0, sole_rep = FALSE,
      stringsAsFactors = FALSE
    ),
    stability_class = 1L
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
                              stability = stab)

  hog1 <- result[result$hog == "HOG1", ]
  expect_false(is.na(hog1$stability_class))
})


test_that("stability_class uses max across multi-clique HOGs", {
  setup <- make_classify_edges()

  # Mock: HOG3 has two exclusive cliques with different stability_class
  stab <- list(
    stability = data.frame(
      clique_idx = c(2L, 3L),
      hog = c("HOG3", "HOG3"),
      trait_value = c("annual", "perennial"),
      k = c(1L, 1L),
      n_subsets = c(2L, 2L),
      n_stable = c(2L, 0L),
      stability_score = c(1.0, 0.0),
      sole_rep = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    stability_class = c(2L, 0L)  # first clique very stable, second not
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
                              stability = stab)

  hog3 <- result[result$hog == "HOG3", ]
  # max(2, 0) = 2 — best clique wins
  expect_equal(hog3$stability_class, 2L)
})


test_that("robust flag without annotations is NA", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  expect_true(all(is.na(result$robust)))
})


test_that("robust flag is TRUE when both stability and sweep pass thresholds", {
  setup <- make_classify_edges()

  stab <- list(
    stability = data.frame(
      clique_idx = 0L, hog = "HOG1",
      trait_value = "annual", k = 1L,
      n_subsets = 2L, n_stable = 2L,
      stability_score = 1.0, sole_rep = FALSE,
      stringsAsFactors = FALSE
    ),
    stability_class = 1L
  )

  # Mock sweep: HOG1 survived at multiplier 2
  sweep <- list(
    survival = data.frame(
      clique_idx = 0L, hog = "HOG1",
      multiplier = 2, survived = TRUE,
      jaccard = 1.0, n_species_orig = 4L,
      n_species_new = 4L,
      stringsAsFactors = FALSE
    )
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
                              stability = stab, sweep = sweep,
                              min_stability_class = 1L,
                              min_persistence = 1.5)

  hog1 <- result[result$hog == "HOG1", ]
  expect_false(is.na(hog1$stability_class))
  expect_false(is.na(hog1$persistence))
  expect_equal(hog1$persistence, 2)
  expect_true(hog1$robust)

  # HOGs without stability/sweep data should not be robust
  hog5 <- result[result$hog == "HOG5", ]
  expect_false(isTRUE(hog5$robust))
})


test_that("robust is FALSE when stability passes but sweep fails", {
  setup <- make_classify_edges()

  stab <- list(
    stability = data.frame(
      clique_idx = 0L, hog = "HOG1",
      trait_value = "annual", k = 1L,
      n_subsets = 2L, n_stable = 2L,
      stability_score = 1.0, sole_rep = FALSE,
      stringsAsFactors = FALSE
    ),
    stability_class = 1L
  )

  # Sweep: HOG1 did NOT survive at any multiplier
  sweep <- list(
    survival = data.frame(
      clique_idx = 0L, hog = "HOG1",
      multiplier = 2, survived = FALSE,
      jaccard = 0.0, n_species_orig = 4L,
      n_species_new = NA_integer_,
      stringsAsFactors = FALSE
    )
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
                              stability = stab, sweep = sweep,
                              min_stability_class = 1L,
                              min_persistence = 1.5)

  hog1 <- result[result$hog == "HOG1", ]
  # Stability passes but persistence is NA (never survived) -> robust FALSE
  expect_false(hog1$robust)
})


test_that("stability=list() is rejected by input validation", {
  setup <- make_classify_edges()
  expect_error(
    classify_cliques(setup$edges, setup$target, setup$trait,
                     stability = list()),
    "stability must be output of clique_stability")
})


test_that("single-species trait groups are handled gracefully", {
  # Only 3 species: SP_A alone in its group
  target <- c("SP_A", "SP_B", "SP_C")
  trait <- c(SP_A = "rare", SP_B = "common", SP_C = "common")

  edges <- data.frame(
    gene1 = c("B1"), gene2 = c("C1"),
    species1 = c("SP_B"), species2 = c("SP_C"),
    hog = "HOG1",
    q.value = 0.01, effect_size = 3.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # SP_A has no genes -> "rare" group can't form a clique
  # "common" group (B,C) has a clique -> trait_specific
  result <- classify_cliques(edges, target, trait)

  hog1 <- result[result$hog == "HOG1", ]
  expect_equal(hog1$classification, "trait_specific")
  expect_equal(hog1$trait_groups, "common")
})


test_that("end-to-end: real clique_stability output feeds classify_cliques", {
  setup <- make_classify_edges()
  target <- setup$target
  trait <- setup$trait

  # Find cliques — need trait-exclusive ones for stability to track
  cliques <- find_cliques(setup$edges, target, min_species = 2L,
                          edge_type = "conserved")
  if (nrow(cliques) == 0) skip("No cliques found for stability test")

  # clique_stability only tracks trait-exclusive cliques.
  # HOG3 has annual-exclusive (A3-B3) and perennial-exclusive (C3-D3) cliques.
  # HOG4 has annual-exclusive (A4-B4).
  stab <- clique_stability(setup$edges, target, trait,
                            all_species = target,
                            full_cliques = cliques,
                            max_k = 1L)

  result <- classify_cliques(setup$edges, target, trait,
                              stability = stab)

  # HOG3 is differentiated — its annual clique should have stability data
  hog3 <- result[result$hog == "HOG3", ]
  expect_equal(hog3$classification, "differentiated")

  # HOG4 is trait_specific — its annual clique should have stability data
  hog4 <- result[result$hog == "HOG4", ]
  expect_equal(hog4$classification, "trait_specific")

  # At least one of HOG3/HOG4 should have non-NA stability_class
  # (they are trait-exclusive, so stability tracks them)
  has_sc <- !is.na(hog3$stability_class) || !is.na(hog4$stability_class)
  expect_true(has_sc)

  # HOG1 (complete, mixed-trait) should have NA stability_class
  # (not trait-exclusive, so stability doesn't track it)
  hog1 <- result[result$hog == "HOG1", ]
  expect_true(is.na(hog1$stability_class))
})
