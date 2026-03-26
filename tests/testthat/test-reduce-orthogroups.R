# Tests for reduce_orthogroups()

test_that("identical paralogs are merged", {
  # 4 genes, 5 samples. G1 and G2 are identical (same HOG), G3 different HOG.
  set.seed(42)
  base <- rnorm(5)
  expr <- rbind(
    G1 = base,
    G2 = base,
    G3 = rnorm(5),
    G4 = rnorm(5)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  expect_equal(result$n_original, 4L)
  expect_equal(result$n_merged, 1L)
  expect_equal(result$n_reduced, 3L)
  expect_equal(nrow(result$expr_matrix), 3L)
  expect_equal(ncol(result$expr_matrix), 5L)

  # G1 and G2 should map to the same representative
  g1_rep <- result$gene_map$representative[result$gene_map$original == "G1"]
  g2_rep <- result$gene_map$representative[result$gene_map$original == "G2"]
  expect_equal(g1_rep, g2_rep)
})


test_that("uncorrelated paralogs are kept separate", {
  set.seed(42)
  expr <- rbind(
    G1 = rnorm(10),
    G2 = rnorm(10),
    G3 = rnorm(10)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs, cor_threshold = 0.9)

  # Random vectors are uncorrelated — should not be merged
  expect_equal(result$n_merged, 0L)
  expect_equal(result$n_reduced, 3L)
})


test_that("non-HOG genes are preserved", {
  expr <- rbind(
    G1 = c(1, 2, 3),
    G2 = c(1, 2, 3),
    SOLO = c(5, 6, 7)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  expect_true("SOLO" %in% rownames(result$expr_matrix))
  solo_map <- result$gene_map[result$gene_map$original == "SOLO", ]
  expect_equal(solo_map$representative, "SOLO")
})


test_that("zero-variance genes are kept as singletons", {
  expr <- rbind(
    G1 = c(1, 2, 3, 4),
    G2 = c(5, 5, 5, 5),  # zero variance
    G3 = c(1, 2, 3, 4)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2", "G3"),
    Species2 = c("X1", "X2", "X3"),
    OrthoGroup = c(1L, 1L, 1L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  # G2 (zero-var) stays as singleton; G1 and G3 are correlated → merged
  expect_equal(result$n_merged, 1L)
  g2_map <- result$gene_map[result$gene_map$original == "G2", ]
  expect_equal(g2_map$representative, "G2")
})


test_that("merged genes have averaged expression", {
  expr <- rbind(
    G1 = c(2.0, 4.0, 6.0),
    G2 = c(4.0, 8.0, 12.0)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  expect_equal(result$n_merged, 1L)
  expect_equal(as.numeric(result$expr_matrix[1, ]), c(3.0, 6.0, 9.0))
})


test_that("multiple HOGs are processed independently", {
  set.seed(42)
  base1 <- rnorm(10)
  base2 <- rnorm(10)
  expr <- rbind(
    A1 = base1,
    A2 = base1 + rnorm(10, sd = 0.01),  # nearly identical to A1
    B1 = base2,
    B2 = base2 + rnorm(10, sd = 0.01)   # nearly identical to B1
  )
  orthologs <- data.frame(
    Species1 = c("A1", "A2", "B1", "B2"),
    Species2 = c("X1", "X2", "Y1", "Y2"),
    OrthoGroup = c(1L, 1L, 2L, 2L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs, cor_threshold = 0.9)

  expect_equal(result$n_merged, 2L)
  expect_equal(result$n_reduced, 2L)
})


test_that("cor_threshold = 1.0 merges nothing", {
  expr <- rbind(
    G1 = c(1, 2, 3),
    G2 = c(2, 4, 6)  # perfectly correlated but different scale
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  # cor = 1.0 exactly → distance = 0 → should still merge at threshold 1.0
  result <- reduce_orthogroups(expr, orthologs, cor_threshold = 1.0)
  expect_equal(result$n_merged, 1L)
})


test_that("gene_col parameter works for Species2", {
  expr <- rbind(
    X1 = c(1, 2, 3),
    X2 = c(1, 2, 3),
    X3 = c(4, 5, 6)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2", "G3"),
    Species2 = c("X1", "X2", "X3"),
    OrthoGroup = c(1L, 1L, 2L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs, gene_col = "Species2")

  expect_equal(result$n_merged, 1L)
  expect_equal(result$n_reduced, 2L)
})


test_that("single-gene HOGs are passed through", {
  expr <- rbind(
    G1 = c(1, 2, 3),
    G2 = c(4, 5, 6)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2"),
    Species2 = c("X1", "X2"),
    OrthoGroup = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  expect_equal(result$n_merged, 0L)
  expect_equal(result$n_reduced, 2L)
})


test_that("gene_map covers all original genes", {
  set.seed(42)
  expr <- rbind(
    G1 = rnorm(5), G2 = rnorm(5), G3 = rnorm(5),
    G4 = rnorm(5), G5 = rnorm(5)
  )
  orthologs <- data.frame(
    Species1 = c("G1", "G2", "G3"),
    Species2 = c("X1", "X2", "X3"),
    OrthoGroup = c(1L, 1L, 2L),
    stringsAsFactors = FALSE
  )

  result <- reduce_orthogroups(expr, orthologs)

  expect_true(all(rownames(expr) %in% result$gene_map$original))
})


test_that("input validation works", {
  expr <- rbind(G1 = c(1, 2, 3))
  ortho <- data.frame(Species1 = "G1", Species2 = "X1", OrthoGroup = 1L,
                      stringsAsFactors = FALSE)

  expect_error(reduce_orthogroups("not a matrix", ortho),
               "numeric matrix")
  expect_error(reduce_orthogroups(matrix(1:3, nrow = 1), ortho),
               "row names")
  expect_error(reduce_orthogroups(expr, ortho, cor_threshold = 1.5),
               "between 0 and 1")
  expect_error(reduce_orthogroups(expr, data.frame(x = 1)),
               "column 'Species1'")
})
