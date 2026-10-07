# Tests for classify_cliques()

# Helper: build edges for 4-species, binary trait (annual/perennial)
make_classify_edges <- function() {
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")
  trait <- c(
    SP_A = "annual", SP_B = "annual",
    SP_C = "perennial", SP_D = "perennial"
  )

  # HOG1: ALL 6 edges conserved (complete)
  hog1 <- data.frame(
    gene1 = c("A1", "A1", "A1", "B1", "B1", "C1"),
    gene2 = c("B1", "C1", "D1", "C1", "D1", "D1"),
    species1 = c("SP_A", "SP_A", "SP_A", "SP_B", "SP_B", "SP_C"),
    species2 = c("SP_B", "SP_C", "SP_D", "SP_C", "SP_D", "SP_D"),
    hog = "HOG1",
    q_value = rep(0.01, 6),
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
    q_value = rep(0.02, 3),
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
    q_value = c(0.01, 0.01, 0.5),
    effect_size = c(3.0, 3.0, 0.8),
    type = c("conserved", "conserved", "ns"),
    stringsAsFactors = FALSE
  )

  # HOG4: annual-only clique (A4-B4), no perennial
  hog4 <- data.frame(
    gene1 = "A4", gene2 = "B4",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG4",
    q_value = 0.01, effect_size = 2.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # HOG5: no conserved edges (unclassified)
  hog5 <- data.frame(
    gene1 = "A5", gene2 = "B5",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG5",
    q_value = 0.8, effect_size = 0.5,
    type = "ns",
    stringsAsFactors = FALSE
  )

  edges <- rbind(hog1, hog2, hog3, hog4, hog5)
  list(
    edges = edges, target = target,
    trait = split(names(trait), trait)
  )
}


test_that("classify_cliques returns correct structure", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  expect_true(is.data.frame(result))
  expected_cols <- c(
    "hog", "classification", "n_species", "best_mean_q",
    "trait_groups", "clade", "stability_class", "robust"
  )
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
  trait <- c(
    SP_A = "warm", SP_B = "warm",
    SP_C = "cold", SP_D = "cold",
    SP_E = "arid", SP_F = "arid"
  )

  # HOG1: warm clique + cold clique + arid clique, no cross-group
  edges <- data.frame(
    gene1 = c("A1", "C1", "E1"),
    gene2 = c("B1", "D1", "F1"),
    species1 = c("SP_A", "SP_C", "SP_E"),
    species2 = c("SP_B", "SP_D", "SP_F"),
    hog = "HOG1",
    q_value = rep(0.01, 3),
    effect_size = rep(3.0, 3),
    type = "conserved",
    stringsAsFactors = FALSE
  )

  clades <- split(names(trait), trait)
  result <- classify_cliques(edges, target, clades)
  expect_equal(result$classification, "differentiated")
  # All 3 groups should be listed
  groups <- strsplit(result$trait_groups, ",")[[1]]
  expect_equal(length(groups), 3)
})


test_that("classify_cliques validates inputs", {
  setup <- make_classify_edges()

  expect_error(
    classify_cliques(data.frame(x = 1), setup$target, setup$trait),
    "edges missing required columns"
  )

  expect_error(
    classify_cliques(setup$edges, c("SP_A"), setup$trait),
    "at least 2 species"
  )

  expect_error(
    classify_cliques(setup$edges, setup$target, c("a", "b")),
    "clades must be a named list"
  )
})


test_that("classify_cliques with empty edges returns empty", {
  setup <- make_classify_edges()
  empty <- setup$edges[0, , drop = FALSE]

  result <- classify_cliques(empty, setup$target, setup$trait)
  expect_equal(nrow(result), 0)
})


test_that("stability annotation populates stability_class", {
  setup <- make_classify_edges()

  # Mock stability output (structural, no trait_value)
  stab <- list(
    stability = data.frame(
      clique_idx = 1L, hog = "HOG1",
      k = 1L,
      n_subsets = 2L, n_stable = 2L,
      stability_score = 1.0,
      species_present = "SP_A,SP_B",
      traits = "annual",
      sole_rep = FALSE,
      stringsAsFactors = FALSE
    ),
    stability_class = 1L
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
    stability = stab
  )

  hog1 <- result[result$hog == "HOG1", ]
  expect_false(is.na(hog1$stability_class))
})


test_that("stability_class uses max across multi-clique HOGs", {
  setup <- make_classify_edges()

  # Mock: HOG3 has two cliques with different stability_class
  stab <- list(
    stability = data.frame(
      clique_idx = c(2L, 3L),
      hog = c("HOG3", "HOG3"),
      k = c(1L, 1L),
      n_subsets = c(2L, 2L),
      n_stable = c(2L, 0L),
      stability_score = c(1.0, 0.0),
      species_present = c("SP_A,SP_B", "SP_C,SP_D"),
      traits = c("annual", "perennial"),
      sole_rep = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    # stability_class is positional (one per clique in full_cliques):
    # clique 1 = not tested, clique 2 = very stable, clique 3 = not stable
    stability_class = c(0L, 2L, 0L)
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
    stability = stab
  )

  hog3 <- result[result$hog == "HOG3", ]
  # max(2, 0) = 2 — best clique wins
  expect_equal(hog3$stability_class, 2L)
})


test_that("robust flag without annotations is NA", {
  setup <- make_classify_edges()
  result <- classify_cliques(setup$edges, setup$target, setup$trait)

  expect_true(all(is.na(result$robust)))
})


test_that("robust flag is TRUE when stability results are provided", {
  setup <- make_classify_edges()

  stab <- list(
    stability = data.frame(
      clique_idx = 1L, hog = "HOG1",
      k = 1L,
      n_subsets = 2L, n_stable = 2L,
      stability_score = 1.0,
      species_present = "SP_A,SP_B",
      traits = "annual",
      sole_rep = FALSE,
      stringsAsFactors = FALSE
    ),
    stability_class = 1L
  )

  result <- classify_cliques(setup$edges, setup$target, setup$trait,
    stability = stab
  )

  hog1 <- result[result$hog == "HOG1", ]
  expect_false(is.na(hog1$stability_class))
  expect_true(hog1$robust)

  # HOGs without stability data should not be robust
  hog5 <- result[result$hog == "HOG5", ]
  expect_false(isTRUE(hog5$robust))
})


test_that("stability=list() is rejected by input validation", {
  setup <- make_classify_edges()
  expect_error(
    classify_cliques(setup$edges, setup$target, setup$trait,
      stability = list()
    ),
    "stability must be output of clique_stability"
  )
})


test_that("single-species trait groups are handled gracefully", {
  # Only 3 species: SP_A alone in its group
  target <- c("SP_A", "SP_B", "SP_C")
  trait <- c(SP_A = "rare", SP_B = "common", SP_C = "common")

  edges <- data.frame(
    gene1 = c("B1"), gene2 = c("C1"),
    species1 = c("SP_B"), species2 = c("SP_C"),
    hog = "HOG1",
    q_value = 0.01, effect_size = 3.0,
    type = "conserved",
    stringsAsFactors = FALSE
  )

  # SP_A has no genes -> "rare" group can't form a clique
  # "common" group (B,C) has a clique -> trait_specific
  result <- classify_cliques(edges, target, split(names(trait), trait))

  hog1 <- result[result$hog == "HOG1", ]
  expect_equal(hog1$classification, "trait_specific")
  expect_equal(hog1$trait_groups, "common")
})


test_that("end-to-end: real clique_stability output feeds classify_cliques", {
  setup <- make_classify_edges()
  target <- setup$target
  trait <- setup$trait

  # Find cliques — need trait-exclusive ones for stability to track
  cliques <- find_cliques(setup$edges, target, min_species = 2L)
  if (nrow(cliques) == 0) skip("No cliques found for stability test")

  # clique_stability now tests ALL cliques structurally.
  stab <- clique_stability(setup$edges, target, trait,
    all_species = target,
    full_cliques = cliques,
    max_k = 1L
  )

  result <- classify_cliques(setup$edges, target, trait,
    stability = stab
  )

  # HOG3 is differentiated — should have stability data
  hog3 <- result[result$hog == "HOG3", ]
  expect_equal(hog3$classification, "differentiated")

  # HOG4 is trait_specific — should have stability data
  hog4 <- result[result$hog == "HOG4", ]
  expect_equal(hog4$classification, "trait_specific")

  # All HOGs now have stability_class (structural, trait-agnostic)
  expect_false(is.na(hog3$stability_class))
  expect_false(is.na(hog4$stability_class))

  # HOG1 (complete, mixed-trait) also gets stability_class now
  hog1 <- result[result$hog == "HOG1", ]
  expect_false(is.na(hog1$stability_class))
})


# --- Underpowered specificity and divergence (#12) ---

# One extra non-conserved edge row for HOG4 (annual-only clique A4-B4).
up_row <- function(gene1, species1, gene2, species2, power) {
  data.frame(
    gene1 = gene1, gene2 = gene2, species1 = species1,
    species2 = species2, hog = "HOG4", q_value = 0.7,
    effect_size = 0.9, type = "ns", power = power,
    stringsAsFactors = FALSE
  )
}

classify_with_power <- function(extra = NULL, hog3_power = 0.99, ...) {
  res <- .cwp_result(extra, hog3_power, ...)
  stats::setNames(res$classification, res$hog)
}

# The flag, for the same inputs. underpowered no longer replaces the
# classification, so a test that only reads the label cannot tell a
# qualified call from a clean one.
underpowered_with_power <- function(extra = NULL, hog3_power = 0.99, ...) {
  res <- .cwp_result(extra, hog3_power, ...)
  stats::setNames(res$underpowered, res$hog)
}

.cwp_result <- function(extra = NULL, hog3_power = 0.99, ...) {
  setup <- make_classify_edges()
  e <- setup$edges
  e$power <- 0.99
  e$power[e$hog == "HOG3" & e$type == "ns"] <- hog3_power
  if (!is.null(extra)) e <- rbind(e, extra)
  classify_cliques(e, setup$target, setup$trait, ...)
}


test_that("powered edges leave the classification unchanged", {
  setup <- make_classify_edges()
  base <- classify_cliques(setup$edges, setup$target, setup$trait)
  cls <- classify_with_power()
  expect_equal(cls[base$hog], stats::setNames(base$classification, base$hog))
  expect_equal(
    classify_with_power(hog3_power = NA_real_)[["HOG3"]],
    "differentiated"
  )
})


test_that("an underpowered cross edge blocks differentiated", {
  setup <- make_classify_edges()
  e <- setup$edges
  e$power <- 0.99
  e$power[e$hog == "HOG3" & e$type == "ns"] <- 0.1
  res <- classify_cliques(e, setup$target, setup$trait)
  hog3 <- res[res$hog == "HOG3", ]
  expect_equal(hog3$classification, "differentiated")
  expect_true(hog3$underpowered)
  expect_equal(hog3$trait_groups, "annual,perennial")
  expect_equal(
    classify_with_power(hog3_power = 0.1, min_power = 0.05)[["HOG3"]],
    "differentiated"
  )
})


test_that("an underpowered cross edge blocks trait_specific", {
  arg <- up_row("A4", "SP_A", "C4", "SP_C", 0.1)
  low <- classify_with_power(arg)
  expect_equal(low[["HOG4"]], "trait_specific")
  expect_equal(low[["HOG3"]], "differentiated")
  # the call is kept; the flag is what separates it from a clean one
  expect_true(underpowered_with_power(arg)[["HOG4"]])

  hi_arg <- up_row("A4", "SP_A", "C4", "SP_C", 0.9)
  high <- classify_with_power(hi_arg)
  expect_equal(high[["HOG4"]], "trait_specific")
  expect_false(underpowered_with_power(hi_arg)[["HOG4"]])

  # Not deciding: the endpoint is no clique member, or the other endpoint
  # shares the clique's trait group.
  st_arg <- up_row("A4x", "SP_A", "C4", "SP_C", 0.1)
  expect_equal(classify_with_power(st_arg)[["HOG4"]], "trait_specific")
  expect_false(underpowered_with_power(st_arg)[["HOG4"]])
  sm_arg <- up_row("A4", "SP_A", "B4x", "SP_B", 0.1)
  expect_equal(classify_with_power(sm_arg)[["HOG4"]], "trait_specific")
  expect_false(underpowered_with_power(sm_arg)[["HOG4"]])

  setup <- make_classify_edges()
  e <- setup$edges
  e$power <- 0.99
  e <- rbind(e, up_row("C4", "SP_C", "B4", "SP_B", 0.1))
  res <- classify_cliques(e, setup$target, setup$trait)
  expect_equal(res$classification[res$hog == "HOG4"], "trait_specific")
  expect_true(res$underpowered[res$hog == "HOG4"])
  expect_equal(res$trait_groups[res$hog == "HOG4"], "annual")
})


test_that("classify_cliques validates min_power", {
  setup <- make_classify_edges()
  for (bad in list(-1, 2, NA_real_, c(0.1, 0.2))) {
    expect_error(
      classify_cliques(setup$edges, setup$target, setup$trait,
        min_power = bad
      ),
      "min_power must be a single number"
    )
  }
})


test_that("a clique specific to an inner clade reports that clade", {
  sp <- paste0("SP_", LETTERS[1:6])
  clades <- list(
    outer = sp[1:4], mid = sp[1:3], inner = sp[1:2], other = sp[5:6]
  )
  edge <- function(hog, s1, s2, q, power = 0.99) {
    data.frame(
      gene1 = paste0(s1, hog), gene2 = paste0(s2, hog),
      species1 = s1, species2 = s2, hog = hog, q_value = q,
      effect_size = 1, type = ifelse(q < 0.1, "conserved", "ns"),
      power = power, stringsAsFactors = FALSE
    )
  }
  e <- rbind(
    # HOG1: A-B conserved; C and E tested and rejected.
    edge("HOG1", "SP_A", "SP_B", 0.01),
    edge("HOG1", "SP_A", "SP_C", 0.7),
    edge("HOG1", "SP_B", "SP_E", 0.7),
    # HOG2: A-B-C conserved, D rejected.
    edge(
      "HOG2", c("SP_A", "SP_A", "SP_B"), c("SP_B", "SP_C", "SP_C"), 0.01
    ),
    edge("HOG2", "SP_C", "SP_D", 0.7)
  )
  res <- classify_cliques(e, sp, clades)
  cl <- stats::setNames(res$clade, res$hog)
  expect_equal(res$classification, rep("trait_specific", 2))
  expect_equal(res$trait_groups, rep("outer", 2))
  expect_equal(cl[["HOG1"]], "inner")
  expect_equal(cl[["HOG2"]], "mid")
  expect_false(any(res$underpowered))

  # A-C crosses the inner boundary, though C shares the outer group.
  e$power[2] <- 0.1
  res <- classify_cliques(e, sp, clades)
  up <- stats::setNames(res$underpowered, res$hog)
  expect_true(up[["HOG1"]])
  expect_false(up[["HOG2"]])
})


test_that("species-graph tiers read home clades, not top-level ones", {
  sp <- paste0("SP_", LETTERS[1:6])
  clades <- list(
    outer = sp[1:5], mid = sp[1:4], inner = sp[1:2], sib = sp[3:4],
    other = sp[6]
  )
  edge <- function(hog, s1, s2, q, power = 0.99) {
    data.frame(
      gene1 = paste0(s1, hog), gene2 = paste0(s2, hog),
      species1 = s1, species2 = s2, hog = hog, q_value = q,
      effect_size = 1, type = ifelse(q < 0.1, "conserved", "ns"),
      power = power, stringsAsFactors = FALSE
    )
  }
  e <- rbind(
    # HOG1: inner conserved, C rejected.
    edge("HOG1", "SP_A", "SP_B", 0.01), edge("HOG1", "SP_A", "SP_C", 0.7),
    # HOG2: inner and its sibling conserved, rejected across.
    edge("HOG2", c("SP_A", "SP_C"), c("SP_B", "SP_D"), 0.01),
    edge("HOG2", "SP_A", "SP_C", 0.7),
    # HOG3: A-B (home inner) and B-C (home mid) are nested homes.
    edge("HOG3", c("SP_A", "SP_B"), c("SP_B", "SP_C"), 0.01),
    edge("HOG3", "SP_A", "SP_C", 0.7)
  )
  res <- classify_cliques(e, sp, clades)
  got <- stats::setNames(res$classification, res$hog)
  expect_equal(
    got[c("HOG1", "HOG2", "HOG3")],
    c(
      HOG1 = "trait_specific", HOG2 = "differentiated",
      HOG3 = "trait_specific"
    )
  )
  cl <- stats::setNames(res$clade, res$hog)
  expect_equal(
    cl[c("HOG1", "HOG2", "HOG3")],
    c(HOG1 = "inner", HOG2 = "mid", HOG3 = "mid")
  )
  expect_equal(res$trait_groups, rep("outer", 3))
  expect_false(any(res$underpowered))

  # A-C crosses from inner to its sibling, inside one top-level clade.
  e$power[e$hog == "HOG2" & e$type == "ns"] <- 0.1
  res <- classify_cliques(e, sp, clades)
  expect_true(res$underpowered[res$hog == "HOG2"])
})
