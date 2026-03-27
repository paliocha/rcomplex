# Tests for clique_hubs()

make_hub_cliques <- function() {
  # Gene A1 appears in 3 HOGs, B1 in 2, C1 in 1
  data.frame(
    hog = c("HOG1", "HOG2", "HOG3"),
    SP_A = c("A1", "A1", "A1"),
    SP_B = c("B1", "B1", "B2"),
    SP_C = c("C1", "C2", "C3"),
    n_species = c(3L, 3L, 3L),
    mean_q = c(0.01, 0.02, 0.03),
    max_q = c(0.01, 0.02, 0.03),
    mean_effect_size = c(2.0, 2.0, 2.0),
    n_edges = c(3L, 3L, 3L),
    stringsAsFactors = FALSE
  )
}


test_that("hub genes are ranked by clique count", {
  cliques <- make_hub_cliques()
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, min_cliques = 1L)

  expect_equal(result$gene[1], "A1")
  expect_equal(result$n_cliques[result$gene == "A1"], 3L)
  expect_equal(result$n_cliques[result$gene == "B1"], 2L)
  # n_hogs equals n_cliques here (1 clique per HOG)
  expect_equal(result$n_hogs[result$gene == "A1"], 3L)
})


test_that("min_cliques filters genes by clique count", {
  cliques <- make_hub_cliques()
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, min_cliques = 3L)

  expect_equal(nrow(result), 1)
  expect_equal(result$gene, "A1")
  expect_equal(result$n_cliques, 3L)
})


test_that("trait annotation produces per-trait columns", {
  cliques <- make_hub_cliques()
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "annual")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, species_trait = trait, min_cliques = 1L)

  expect_true("n_exclusive" %in% names(result))
  expect_true("n_annual" %in% names(result))
  # All cliques are annual-exclusive (all species are annual)
  expect_equal(result$n_exclusive[result$gene == "A1"], 3L)
  expect_equal(result$n_annual[result$gene == "A1"], 3L)
})


test_that("mixed cliques are not counted as exclusive", {
  # HOG1: SP_A + SP_B (annual-exclusive), HOG2: SP_A + SP_C (mixed)
  cliques <- data.frame(
    hog = c("HOG1", "HOG2"),
    SP_A = c("A1", "A1"),
    SP_B = c("B1", NA),
    SP_C = c(NA, "C1"),
    n_species = c(2L, 2L),
    mean_q = c(0.01, 0.01),
    max_q = c(0.01, 0.01),
    mean_effect_size = c(2.0, 2.0),
    n_edges = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "perennial")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, species_trait = trait, min_cliques = 1L)

  a1 <- result[result$gene == "A1", ]
  expect_equal(a1$n_cliques, 2L)
  expect_equal(a1$n_exclusive, 1L)
  expect_equal(a1$n_annual, 1L)
})


test_that("multiple trait levels produce correct columns", {
  cliques <- data.frame(
    hog = c("HOG1", "HOG2", "HOG3"),
    SP_A = c("A1", NA, "A1"),
    SP_B = c("B1", NA, NA),
    SP_C = c(NA, "C1", "C1"),
    SP_D = c(NA, "D1", NA),
    n_species = c(2L, 2L, 2L),
    mean_q = c(0.01, 0.01, 0.01),
    max_q = c(0.01, 0.01, 0.01),
    mean_effect_size = c(2.0, 2.0, 2.0),
    n_edges = c(1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "warm", SP_B = "warm", SP_C = "cold", SP_D = "cold")
  target <- c("SP_A", "SP_B", "SP_C", "SP_D")

  result <- clique_hubs(cliques, target, species_trait = trait, min_cliques = 1L)

  expect_true("n_warm" %in% names(result))
  expect_true("n_cold" %in% names(result))
  # A1: HOG1 warm-excl, HOG3 mixed (warm + cold) -> n_warm = 1
  a1 <- result[result$gene == "A1", ]
  expect_equal(a1$n_warm, 1L)
  expect_equal(a1$n_cold, 0L)
})


test_that("sorted by n_exclusive when trait provided", {
  cliques <- data.frame(
    hog = c("HOG1", "HOG2", "HOG3"),
    SP_A = c("A1", "A2", "A2"),
    SP_B = c("B1", "B1", "B1"),
    n_species = c(2L, 2L, 2L),
    mean_q = c(0.01, 0.01, 0.01),
    max_q = c(0.01, 0.01, 0.01),
    mean_effect_size = c(2.0, 2.0, 2.0),
    n_edges = c(1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "x", SP_B = "x")
  target <- c("SP_A", "SP_B")

  result <- clique_hubs(cliques, target, species_trait = trait, min_cliques = 1L)

  # B1 appears in 3 exclusive cliques, A2 in 2, A1 in 1
  expect_equal(result$gene[1], "B1")
})


test_that("empty cliques returns empty result", {
  cliques <- data.frame(
    hog = character(0), SP_A = character(0), SP_B = character(0),
    n_species = integer(0), mean_q = numeric(0), max_q = numeric(0),
    mean_effect_size = numeric(0), n_edges = integer(0),
    stringsAsFactors = FALSE
  )
  result <- clique_hubs(cliques, c("SP_A", "SP_B"), min_cliques = 1L)
  expect_equal(nrow(result), 0)
})


test_that("no trait produces no exclusive columns", {
  cliques <- make_hub_cliques()
  result <- clique_hubs(cliques, c("SP_A", "SP_B", "SP_C"), min_cliques = 1L)

  expect_true("n_cliques" %in% names(result))
  expect_true("n_hogs" %in% names(result))
  expect_false("n_exclusive" %in% names(result))
})


test_that("duplicate HOG with different species compositions annotates correctly", {
  # HOG1 row1: A1+B1 (annual+annual = annual-exclusive)
  # HOG1 row2: B1+C1 (annual+perennial = mixed)
  cliques <- data.frame(
    hog = c("HOG1", "HOG1"),
    SP_A = c("A1", NA),
    SP_B = c("B1", "B1"),
    SP_C = c(NA, "C1"),
    n_species = c(2L, 2L),
    mean_q = rep(0.01, 2),
    max_q = rep(0.01, 2),
    mean_effect_size = rep(2.0, 2),
    n_edges = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "perennial")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, species_trait = trait, min_cliques = 1L)

  # A1 only in row1 (annual-excl) -> n_exclusive = 1
  expect_equal(result$n_exclusive[result$gene == "A1"], 1L)
  # C1 only in row2 (mixed) -> n_exclusive = 0
  expect_equal(result$n_exclusive[result$gene == "C1"], 0L)
  # B1 in both rows: row1 annual-excl, row2 mixed -> n_cliques = 2, n_hogs = 1
  expect_equal(result$n_exclusive[result$gene == "B1"], 1L)
  expect_equal(result$n_cliques[result$gene == "B1"], 2L)
  expect_equal(result$n_hogs[result$gene == "B1"], 1L)
})


# ---- Stability-weighted hub tests ----

test_that("stability weighting adds mean_stability and n_stable", {
  cliques <- make_hub_cliques()
  trait <- c(SP_A = "x", SP_B = "x", SP_C = "x")
  target <- c("SP_A", "SP_B", "SP_C")

  # Mock stability output: HOG1 and HOG2 stable, HOG3 not
  stab <- list(
    stability = data.frame(
      clique_idx = c(0L, 1L, 2L),
      hog = c("HOG1", "HOG2", "HOG3"),
      trait_value = rep("x", 3),
      k = c(1L, 1L, 1L),
      n_subsets = c(2L, 2L, 2L),
      n_stable = c(2L, 2L, 0L),
      stability_score = c(1.0, 1.0, 0.0),
      sole_rep = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    stability_class = c(1L, 1L, 0L)
  )

  result <- clique_hubs(cliques, target, species_trait = trait,
                        stability = stab, min_cliques = 1L)

  expect_true("mean_stability" %in% names(result))
  expect_true("n_stable" %in% names(result))

  # A1 in HOG1 (stable), HOG2 (stable), HOG3 (unstable) -> n_stable=2, mean=2/3
  a1 <- result[result$gene == "A1", ]
  expect_equal(a1$n_stable, 2L)
  expect_equal(a1$mean_stability, 2 / 3, tolerance = 1e-10)

  # B1 in HOG1 (stable), HOG2 (stable) -> n_stable=2, mean=1.0
  b1 <- result[result$gene == "B1", ]
  expect_equal(b1$n_stable, 2L)
  expect_equal(b1$mean_stability, 1.0)
})


test_that("stability sorts by n_stable first", {
  cliques <- data.frame(
    hog = c("HOG1", "HOG2", "HOG3"),
    SP_A = c("A1", "A2", "A2"),
    SP_B = c("B1", "B1", "B1"),
    n_species = c(2L, 2L, 2L),
    mean_q = rep(0.01, 3), max_q = rep(0.01, 3),
    mean_effect_size = rep(2.0, 3), n_edges = c(1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  trait <- c(SP_A = "x", SP_B = "x")
  target <- c("SP_A", "SP_B")

  # B1 in all 3 cliques but only 1 stable; A2 in 2 cliques both stable
  stab <- list(
    stability = data.frame(
      clique_idx = c(0L, 1L, 2L),
      hog = c("HOG1", "HOG2", "HOG3"),
      trait_value = rep("x", 3),
      k = c(1L, 1L, 1L),
      n_subsets = c(1L, 1L, 1L),
      n_stable = c(0L, 1L, 1L),
      stability_score = c(0.0, 1.0, 1.0),
      sole_rep = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    stability_class = c(0L, 1L, 1L)
  )

  result <- clique_hubs(cliques, target, species_trait = trait,
                        stability = stab, min_cliques = 1L)

  # A2: n_stable=2, B1: n_stable=2 (HOG2+HOG3 stable), A1: n_stable=0
  # B1 has 3 cliques total (n_exclusive=3) but only 2 stable
  # Sort: n_stable desc, then mean_stability desc
  expect_equal(result$gene[nrow(result)], "A1")  # A1 last (0 stable)
})


test_that("stability without trait raises error", {
  cliques <- make_hub_cliques()
  stab <- list(stability = data.frame(clique_idx = integer(0)))
  expect_error(clique_hubs(cliques, c("SP_A", "SP_B", "SP_C"),
                           stability = stab),
               "species_trait is required")
})


test_that("no stability columns when stability is NULL", {
  cliques <- make_hub_cliques()
  trait <- c(SP_A = "x", SP_B = "x", SP_C = "x")
  result <- clique_hubs(cliques, c("SP_A", "SP_B", "SP_C"),
                        species_trait = trait, min_cliques = 1L)
  expect_false("mean_stability" %in% names(result))
  expect_false("n_stable" %in% names(result))
})


test_that("input validation works", {
  expect_error(clique_hubs("not_a_df", c("SP_A")),
               "data frame")
  cliques <- make_hub_cliques()
  expect_error(clique_hubs(cliques, c("SP_A", "SP_MISSING")),
               "missing columns")
})
