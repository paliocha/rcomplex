# --- Tests for find_cliques() ---

test_that("find_cliques simple 3-species clique", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 1)
  expect_equal(result$hog, "HOG1")
  expect_equal(result$SP_A, "A1")
  expect_equal(result$SP_B, "B1")
  expect_equal(result$SP_C, "C1")
  expect_equal(result$n_species, 3L)
  expect_equal(result$n_edges, 3L)
})


test_that("find_cliques FDR optimization picks best assignment", {
  # Two paralogs: A1 has lower FDR edges, A2 has higher
  # C++ Bron-Kerbosch finds one maximal species clique, then picks the

  # single best gene assignment (lowest mean FDR) -> returns 1 row (A1)
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A2", "A2"),
    gene2 = c("B1", "C1", "C1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
    hog = rep("HOG1", 5),
    type = rep("conserved", 5),
    q.value = c(0.01, 0.01, 0.01, 0.5, 0.5),
    effect_size = c(3.0, 3.0, 3.0, 1.0, 1.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  # C++ picks the single best assignment per species clique: A1 (lower FDR)
  expect_equal(nrow(result), 1)
  expect_equal(result$SP_A, "A1")
  expect_equal(result$mean_q, mean(c(0.01, 0.01, 0.01)), tolerance = 1e-10)
})


test_that("find_cliques handles multiple HOGs", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "X1", "X1", "Y1"),
    gene2 = c("B1", "C1", "C1", "Y1", "Z1", "Z1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C", "SP_C"),
    hog = c(rep("HOG1", 3), rep("HOG2", 3)),
    type = rep("conserved", 6),
    q.value = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
    effect_size = c(2, 3, 4, 5, 6, 7),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 2)
  expect_setequal(result$hog, c("HOG1", "HOG2"))
})


test_that("find_cliques returns correct FDR and effect size summary stats", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.04, 0.07),
    effect_size = c(2.0, 6.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(result$mean_q, mean(c(0.01, 0.04, 0.07)), tolerance = 1e-10)
  expect_equal(result$max_q, 0.07)
  expect_equal(result$mean_effect_size, 4.0)
  expect_equal(result$n_edges, 3L)
})


test_that("find_cliques empty input returns correct empty structure", {
  edges <- data.frame(
    gene1 = character(0), gene2 = character(0),
    species1 = character(0), species2 = character(0),
    hog = character(0), q.value = numeric(0),
    effect_size = numeric(0), stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B"))

  expect_equal(nrow(result), 0)
  expect_true(all(c("hog", "SP_A", "SP_B", "n_species",
                     "mean_q") %in% names(result)))
})


test_that("find_cliques validates required columns", {
  edges <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)
  expect_error(find_cliques(edges, c("SP_A", "SP_B")),
               "missing required columns")
})


test_that("find_cliques edge_type filtering works", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = c("conserved", "conserved", "diverged"),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 0.3),
    stringsAsFactors = FALSE
  )

  # conserved only: missing B1-C1 edge -> no clique
  result_con <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  expect_equal(nrow(result_con), 0)

  # Both types: complete clique
  result_both <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
                              edge_type = c("conserved", "diverged"))
  expect_equal(nrow(result_both), 1)
})


test_that("find_cliques min_species allows partial cliques", {
  # 3 target species but allow 2-species cliques
  edges <- data.frame(
    gene1 = "A1",
    gene2 = "B1",
    species1 = "SP_A",
    species2 = "SP_B",
    hog = "HOG1",
    type = "conserved",
    q.value = 0.01,
    effect_size = 3.0,
    stringsAsFactors = FALSE
  )

  # min_species = 3 (default): no clique (only 2 species present)
  result_strict <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
                                min_species = 3L)
  expect_equal(nrow(result_strict), 0)

  # min_species = 2: 2-species clique found
  result_lenient <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
                                 min_species = 2L)
  expect_equal(nrow(result_lenient), 1)
  expect_equal(result_lenient$n_species, 2L)
})


test_that("find_cliques max_genes_per_sp limits combinatorial explosion", {
  # Many paralogs in SP_A, should be capped
  genes_a <- paste0("A", 1:15)
  edges_list <- lapply(genes_a, function(a) {
    data.frame(
      gene1 = a, gene2 = "B1",
      species1 = "SP_A", species2 = "SP_B",
      hog = "HOG1", type = "conserved",
      q.value = 0.01, effect_size = 2.0,
      stringsAsFactors = FALSE
    )
  })
  edges <- do.call(rbind, edges_list)

  # With default max_genes_per_sp = 10, should cap at 10 cliques
  result <- find_cliques(edges, c("SP_A", "SP_B"), max_genes_per_sp = 10L)
  expect_true(nrow(result) <= 10)
})


test_that("find_cliques fewer than 2 target species raises error", {
  edges <- data.frame(
    gene1 = "A1", gene2 = "B1",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", type = "conserved",
    q.value = 0.01, effect_size = 2.0,
    stringsAsFactors = FALSE
  )
  expect_error(find_cliques(edges, "SP_A"))
})
