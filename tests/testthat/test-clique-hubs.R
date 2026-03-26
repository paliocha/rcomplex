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


test_that("hub genes are ranked by HOG count", {
  cliques <- make_hub_cliques()
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, min_hogs = 1L)

  expect_equal(result$gene[1], "A1")
  expect_equal(result$n_hogs[result$gene == "A1"], 3L)
  expect_equal(result$n_hogs[result$gene == "B1"], 2L)
})


test_that("min_hogs filters genes", {
  cliques <- make_hub_cliques()
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, min_hogs = 3L)

  expect_equal(nrow(result), 1)
  expect_equal(result$gene, "A1")
})


test_that("trait annotation produces per-trait columns", {
  cliques <- make_hub_cliques()
  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "annual")
  target <- c("SP_A", "SP_B", "SP_C")

  result <- clique_hubs(cliques, target, species_trait = trait, min_hogs = 1L)

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

  result <- clique_hubs(cliques, target, species_trait = trait, min_hogs = 1L)

  a1 <- result[result$gene == "A1", ]
  expect_equal(a1$n_hogs, 2L)
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

  result <- clique_hubs(cliques, target, species_trait = trait, min_hogs = 1L)

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

  result <- clique_hubs(cliques, target, species_trait = trait, min_hogs = 1L)

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
  result <- clique_hubs(cliques, c("SP_A", "SP_B"), min_hogs = 1L)
  expect_equal(nrow(result), 0)
})


test_that("no trait produces no exclusive columns", {
  cliques <- make_hub_cliques()
  result <- clique_hubs(cliques, c("SP_A", "SP_B", "SP_C"), min_hogs = 1L)

  expect_true("n_hogs" %in% names(result))
  expect_false("n_exclusive" %in% names(result))
})


test_that("input validation works", {
  expect_error(clique_hubs("not_a_df", c("SP_A")),
               "data frame")
  cliques <- make_hub_cliques()
  expect_error(clique_hubs(cliques, c("SP_A", "SP_MISSING")),
               "missing columns")
})
