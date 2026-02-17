# --- Tests for find_coexpression_cliques() ---

test_that("simple 3-species clique found when all edges present", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 1)
  expect_equal(result$hog, "HOG1")
  expect_equal(result$SP_A, "A1")
  expect_equal(result$SP_B, "B1")
  expect_equal(result$SP_C, "C1")
  expect_equal(result$min_effect_size, 2.0)
  expect_equal(result$mean_effect_size, 3.0)
})

test_that("incomplete graph returns 0 cliques", {
  # Missing the B1-C1 edge
  edges <- data.frame(
    gene1 = c("A1", "A1"),
    gene2 = c("B1", "C1"),
    species1 = c("SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C"),
    hog = rep("HOG1", 2),
    type = rep("conserved", 2),
    effect_size = c(2.0, 3.0),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 0)
  expect_true("SP_A" %in% names(result))
  expect_true("SP_B" %in% names(result))
  expect_true("SP_C" %in% names(result))
})

test_that("multi-paralog HOG returns correct subset of cliques", {
  # SP_A has 2 genes: A1 and A2. A1 is fully connected, A2 is missing edge to C1
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A2", "A2"),
    gene2 = c("B1", "C1", "C1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
    hog = rep("HOG1", 5),
    type = c("conserved", "conserved", "conserved", "conserved", "n.s."),
    effect_size = c(2.0, 3.0, 4.0, 1.5, 0.5),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  # Only A1-B1-C1 clique should be found (A2's edge to C1 is "n.s.", filtered out)
  expect_equal(nrow(result), 1)
  expect_equal(result$SP_A, "A1")
})

test_that("multi-paralog HOG: both paralogs form cliques when fully connected", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A2", "A2"),
    gene2 = c("B1", "C1", "C1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
    hog = rep("HOG1", 5),
    type = rep("conserved", 5),
    effect_size = c(2.0, 3.0, 4.0, 1.5, 2.5),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 2)
  expect_setequal(result$SP_A, c("A1", "A2"))
})

test_that("multiple HOGs processed correctly", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "X1", "X1", "Y1"),
    gene2 = c("B1", "C1", "C1", "Y1", "Z1", "Z1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C", "SP_C"),
    hog = c(rep("HOG1", 3), rep("HOG2", 3)),
    type = rep("conserved", 6),
    effect_size = c(2, 3, 4, 5, 6, 7),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 2)
  expect_setequal(result$hog, c("HOG1", "HOG2"))
})

test_that("empty input returns empty result with correct columns", {
  edges <- data.frame(
    gene1 = character(0), gene2 = character(0),
    species1 = character(0), species2 = character(0),
    hog = character(0), type = character(0),
    effect_size = numeric(0), stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B"))

  expect_equal(nrow(result), 0)
  expect_true(all(c("hog", "SP_A", "SP_B", "min_effect_size",
                     "mean_effect_size") %in% names(result)))
})

test_that("edge order independence: gene1/gene2 swap does not affect result", {
  # Normal order
  edges1 <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  # Swap gene1/gene2 for second edge
  edges2 <- edges1
  edges2$gene1[2] <- "C1"
  edges2$gene2[2] <- "A1"
  edges2$species1[2] <- "SP_C"
  edges2$species2[2] <- "SP_A"

  r1 <- find_coexpression_cliques(edges1, c("SP_A", "SP_B", "SP_C"))
  r2 <- find_coexpression_cliques(edges2, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(r1), 1)
  expect_equal(nrow(r2), 1)
  expect_equal(r1$SP_A, r2$SP_A)
  expect_equal(r1$SP_B, r2$SP_B)
  expect_equal(r1$SP_C, r2$SP_C)
})

test_that("effect size stats are correct", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    effect_size = c(2.0, 6.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(result$min_effect_size, 2.0)
  expect_equal(result$mean_effect_size, 4.0)
})

test_that("only specified edge_type is used", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = c("conserved", "conserved", "diverged"),
    effect_size = c(2.0, 3.0, 0.3),
    stringsAsFactors = FALSE
  )

  # Default edge_type = "conserved" should miss B1-C1 (diverged)
  result_con <- find_coexpression_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  expect_equal(nrow(result_con), 0)

  # Using both types should find the clique
  result_both <- find_coexpression_cliques(
    edges, c("SP_A", "SP_B", "SP_C"),
    edge_type = c("conserved", "diverged")
  )
  expect_equal(nrow(result_both), 1)
})

test_that("2-species cliques work", {
  edges <- data.frame(
    gene1 = c("A1", "A2"),
    gene2 = c("B1", "B1"),
    species1 = c("SP_A", "SP_A"),
    species2 = c("SP_B", "SP_B"),
    hog = rep("HOG1", 2),
    type = rep("conserved", 2),
    effect_size = c(3.0, 5.0),
    stringsAsFactors = FALSE
  )

  result <- find_coexpression_cliques(edges, c("SP_A", "SP_B"))

  # A1-B1 and A2-B1 are both valid 2-species cliques
  expect_equal(nrow(result), 2)
  expect_setequal(result$SP_A, c("A1", "A2"))
  expect_true(all(result$SP_B == "B1"))
})

test_that("missing required columns raise error", {
  edges <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)
  expect_error(
    find_coexpression_cliques(edges, c("SP_A", "SP_B")),
    "edges missing required columns"
  )
})

test_that("fewer than 2 target species raises error", {
  edges <- data.frame(
    gene1 = "A1", gene2 = "B1",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", type = "conserved", effect_size = 2.0,
    stringsAsFactors = FALSE
  )
  expect_error(
    find_coexpression_cliques(edges, "SP_A"),
    "target_species must have at least 2 species"
  )
})
