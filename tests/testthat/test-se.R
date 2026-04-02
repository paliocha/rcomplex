# Tests for SummarizedExperiment integration

skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")

# Helper: minimal long-format data for build_se
make_long_data <- function(species = "SP_A", genes = paste0("G", 1:3),
                           samples = paste0("S", 1:4),
                           hogs = c("HOG1", "HOG1", "HOG2")) {
  n <- length(genes) * length(samples)
  data.frame(
    abbrev = rep(species, n),
    gene_id = rep(genes, each = length(samples)),
    sample_id = rep(samples, length(genes)),
    vst.count = rnorm(n, mean = 10),
    HOG = rep(hogs, each = length(samples)),
    tissue = rep("leaf", n)
  )
}


test_that("build_se creates correct SE with named assay", {
  se <- build_se(make_long_data(), "SP_A", sample_metadata = "tissue")

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(nrow(se), 3)
  expect_equal(ncol(se), 4)
  expect_equal(rownames(se), paste0("G", 1:3))
  expect_equal(colnames(se), paste0("S", 1:4))
  expect_equal(SummarizedExperiment::assayNames(se), "vst.count")
  expect_true("hog" %in% names(SummarizedExperiment::rowData(se)))
  expect_true("tissue" %in% names(SummarizedExperiment::colData(se)))
  expect_equal(se@metadata$species, "SP_A")
})


test_that("build_se validates inputs", {
  expect_error(build_se(data.frame(x = 1), "SP_A"), "missing required columns")
  data <- data.frame(abbrev = "SP_B", gene_id = "G1",
                     sample_id = "S1", vst.count = 1.0)
  expect_error(build_se(data, "SP_A"), "No rows found")
})


test_that("extract_orthologs derives correct pairs from shared HOGs", {
  se1 <- build_se(make_long_data("SP_A", c("A1", "A2"), hogs = c("HOG1", "HOG2")),
                  "SP_A")
  se2 <- build_se(make_long_data("SP_B", c("B1", "B2", "B3"),
                                 hogs = c("HOG1", "HOG2", "HOG3")),
                  "SP_B")

  ortho <- extract_orthologs(se1, se2)

  expect_equal(names(ortho), c("Species1", "Species2", "hog"))
  # HOG1: A1 x B1 = 1 pair. HOG2: A2 x B2 = 1 pair. HOG3: no match.
  expect_equal(nrow(ortho), 2)
  expect_setequal(ortho$hog, c("HOG1", "HOG2"))
})


test_that("extract_orthologs handles paralogs (multi-gene HOGs)", {
  se1 <- build_se(make_long_data("SP_A", c("A1", "A2", "A3"),
                                 hogs = c("HOG1", "HOG1", "HOG2")),
                  "SP_A")
  se2 <- build_se(make_long_data("SP_B", c("B1", "B2"),
                                 hogs = c("HOG1", "HOG2")),
                  "SP_B")

  ortho <- extract_orthologs(se1, se2)

  hog1 <- ortho[ortho$hog == "HOG1", ]
  expect_equal(nrow(hog1), 2)  # A1 x B1, A2 x B1
  expect_setequal(hog1$Species1, c("A1", "A2"))
})


test_that("extract_orthologs returns empty for no shared HOGs", {
  se1 <- build_se(make_long_data("SP_A", c("A1"), hogs = c("HOG1")), "SP_A")
  se2 <- build_se(make_long_data("SP_B", c("B1"), hogs = c("HOG99")), "SP_B")

  ortho <- extract_orthologs(se1, se2)
  expect_equal(nrow(ortho), 0)
  expect_equal(names(ortho), c("Species1", "Species2", "hog"))
})


test_that("extract_orthologs validates missing hog column", {
  mat <- matrix(1:8, nrow = 2, dimnames = list(c("G1", "G2"), paste0("S", 1:4)))
  se_no_hog <- SummarizedExperiment::SummarizedExperiment(assays = list(mat))
  se_with_hog <- build_se(make_long_data("SP_A", c("A1", "A2"),
                                         hogs = c("HOG1", "HOG2")), "SP_A")

  expect_error(extract_orthologs(se_no_hog, se_with_hog), "missing 'hog' column")
  expect_error(extract_orthologs(se_with_hog, se_no_hog), "missing 'hog' column")
})


test_that("compute_network accepts SummarizedExperiment", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(mat) <- paste0("G", 1:20)
  colnames(mat) <- paste0("S", 1:10)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(vst = mat)
  )

  net_mat <- compute_network(mat)
  net_se <- compute_network(se, assay = "vst")

  expect_equal(net_mat$threshold, net_se$threshold)
  expect_equal(net_mat$network, net_se$network)
})


# --- prepare_orthologs tests ---

# Helper: build SE + fake reduction output for testing
make_se_and_reduction <- function(species, genes, hogs, samples = paste0("S", 1:4),
                                  merge_map = NULL) {
  se <- build_se(make_long_data(species, genes, samples, hogs), species)
  # Build a gene_map: by default identity mapping (no merging)
  if (is.null(merge_map)) {
    gm <- data.frame(original = genes, representative = genes)
  } else {
    gm <- merge_map
  }
  reduction <- list(
    expr_matrix = SummarizedExperiment::assay(se),
    gene_map = gm,
    n_original = length(genes),
    n_reduced = length(unique(gm$representative)),
    n_merged = length(genes) - length(unique(gm$representative))
  )
  list(se = se, reduction = reduction)
}


test_that("prepare_orthologs returns correct structure", {
  r1 <- make_se_and_reduction("SP_A", c("A1", "A2"), c("HOG1", "HOG2"))
  r2 <- make_se_and_reduction("SP_B", c("B1", "B2"), c("HOG1", "HOG2"))

  se_list <- list(SP_A = r1$se, SP_B = r2$se)
  reductions <- list(SP_A = r1$reduction, SP_B = r2$reduction)

  ortho <- prepare_orthologs(se_list, reductions)

  expect_true(is.data.frame(ortho))
  expect_equal(names(ortho), c("Species1", "Species2", "hog"))
  expect_equal(nrow(ortho), 2)
  expect_setequal(ortho$hog, c("HOG1", "HOG2"))
})


test_that("prepare_orthologs maps gene names through reductions", {
  # SP_A has paralogs A1, A2 in HOG1 that get merged to A1
  merge_map_a <- data.frame(
    original = c("A1", "A2", "A3"),
    representative = c("A1", "A1", "A3")
  )
  r1 <- make_se_and_reduction("SP_A", c("A1", "A2", "A3"),
                               c("HOG1", "HOG1", "HOG2"),
                               merge_map = merge_map_a)
  r2 <- make_se_and_reduction("SP_B", c("B1", "B2"), c("HOG1", "HOG2"))

  se_list <- list(SP_A = r1$se, SP_B = r2$se)
  reductions <- list(SP_A = r1$reduction, SP_B = r2$reduction)

  ortho <- prepare_orthologs(se_list, reductions)

  # Before reduction, HOG1 would have A1-B1 and A2-B1.

  # After reduction, both A1 and A2 map to A1, so we get A1-B1 (deduplicated).
  hog1_rows <- ortho[ortho$hog == "HOG1", ]
  expect_equal(nrow(hog1_rows), 1)
  expect_equal(hog1_rows$Species1, "A1")
  expect_equal(hog1_rows$Species2, "B1")

  # HOG2 is unchanged
  hog2_rows <- ortho[ortho$hog == "HOG2", ]
  expect_equal(nrow(hog2_rows), 1)
  expect_equal(hog2_rows$Species1, "A3")
})


test_that("prepare_orthologs handles genes not in gene_map", {
  # gene_map only has A1; A2 is not in any orthogroup so not in gene_map
  partial_map <- data.frame(
    original = c("A1"),
    representative = c("A1")
  )
  r1 <- make_se_and_reduction("SP_A", c("A1", "A2"), c("HOG1", "HOG2"),
                               merge_map = partial_map)
  r2 <- make_se_and_reduction("SP_B", c("B1", "B2"), c("HOG1", "HOG2"))

  se_list <- list(SP_A = r1$se, SP_B = r2$se)
  reductions <- list(SP_A = r1$reduction, SP_B = r2$reduction)

  ortho <- prepare_orthologs(se_list, reductions)

  # A2 not in gene_map, should stay as A2
  hog2_rows <- ortho[ortho$hog == "HOG2", ]
  expect_equal(hog2_rows$Species1, "A2")
})


test_that("prepare_orthologs validates inputs", {
  r1 <- make_se_and_reduction("SP_A", c("A1"), c("HOG1"))
  r2 <- make_se_and_reduction("SP_B", c("B1"), c("HOG1"))

  # se_list not named
  expect_error(
    prepare_orthologs(list(r1$se, r2$se),
                              list(SP_A = r1$reduction, SP_B = r2$reduction)),
    "se_list must be a named list"
  )

  # reductions not named
  expect_error(
    prepare_orthologs(list(SP_A = r1$se, SP_B = r2$se),
                              list(r1$reduction, r2$reduction)),
    "reductions must be a named list"
  )

  # reductions missing a species
  expect_error(
    prepare_orthologs(list(SP_A = r1$se, SP_B = r2$se),
                              list(SP_A = r1$reduction)),
    "reductions missing species"
  )

  # reduction without gene_map
  bad_reduction <- list(SP_A = list(expr_matrix = matrix(1)),
                        SP_B = r2$reduction)
  expect_error(
    prepare_orthologs(list(SP_A = r1$se, SP_B = r2$se), bad_reduction),
    "\\$gene_map"
  )

  # Only one species
  expect_error(
    prepare_orthologs(list(SP_A = r1$se),
                              list(SP_A = r1$reduction)),
    "at least two species"
  )
})


test_that("prepare_orthologs works with three species", {
  r1 <- make_se_and_reduction("SP_A", c("A1", "A2"), c("HOG1", "HOG2"))
  r2 <- make_se_and_reduction("SP_B", c("B1", "B2"), c("HOG1", "HOG2"))
  r3 <- make_se_and_reduction("SP_C", c("C1", "C2"), c("HOG1", "HOG3"))

  se_list <- list(SP_A = r1$se, SP_B = r2$se, SP_C = r3$se)
  reductions <- list(SP_A = r1$reduction, SP_B = r2$reduction,
                     SP_C = r3$reduction)

  ortho <- prepare_orthologs(se_list, reductions)

  # A-B: HOG1 + HOG2 = 2 rows, A-C: HOG1 = 1 row, B-C: HOG1 = 1 row
  expect_equal(nrow(ortho), 4)
  expect_setequal(ortho$hog, c("HOG1", "HOG1", "HOG1", "HOG2"))
})
