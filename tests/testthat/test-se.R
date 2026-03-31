# Tests for SummarizedExperiment integration

test_that("build_se creates correct SE from long-format data", {
  data <- data.frame(
    abbrev = rep("SP_A", 12),
    gene_id = rep(paste0("G", 1:3), each = 4),
    sample_id = rep(paste0("S", 1:4), 3),
    vst.count = rnorm(12, mean = 10),
    HOG = rep(c("HOG1", "HOG1", "HOG2"), each = 4),
    tissue = rep("leaf", 12)
  )

  se <- build_se(data, "SP_A", sample_metadata = "tissue")

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(nrow(se), 3)   # 3 genes
  expect_equal(ncol(se), 4)   # 4 samples
  expect_equal(rownames(se), paste0("G", 1:3))
  expect_equal(colnames(se), paste0("S", 1:4))
  expect_true("hog" %in% names(SummarizedExperiment::rowData(se)))
  expect_true("tissue" %in% names(SummarizedExperiment::colData(se)))
  expect_equal(se@metadata$species, "SP_A")
})


test_that("build_se validates inputs", {
  data <- data.frame(x = 1)
  expect_error(build_se(data, "SP_A"), "missing required columns")
  data <- data.frame(abbrev = "SP_B", gene_id = "G1",
                     sample_id = "S1", vst.count = 1.0)
  expect_error(build_se(data, "SP_A"), "No rows found")
})


test_that("extract_orthologs derives correct pairs from shared HOGs", {
  data1 <- data.frame(
    abbrev = rep("SP_A", 8),
    gene_id = rep(c("A1", "A2"), each = 4),
    sample_id = rep(paste0("S", 1:4), 2),
    vst.count = rnorm(8), HOG = rep(c("HOG1", "HOG2"), each = 4)
  )
  data2 <- data.frame(
    abbrev = rep("SP_B", 12),
    gene_id = rep(c("B1", "B2", "B3"), each = 4),
    sample_id = rep(paste0("S", 1:4), 3),
    vst.count = rnorm(12), HOG = rep(c("HOG1", "HOG2", "HOG3"), each = 4)
  )

  se1 <- build_se(data1, "SP_A")
  se2 <- build_se(data2, "SP_B")

  ortho <- extract_orthologs(se1, se2)

  expect_true(all(c("Species1", "Species2", "hog") %in% names(ortho)))
  # HOG1: A1 x B1 = 1 pair. HOG2: A2 x B2 = 1 pair. HOG3: no match.
  expect_equal(nrow(ortho), 2)
  expect_setequal(ortho$hog, c("HOG1", "HOG2"))
})


test_that("extract_orthologs handles paralogs (multi-gene HOGs)", {
  data1 <- data.frame(
    abbrev = rep("SP_A", 12),
    gene_id = rep(c("A1", "A2", "A3"), each = 4),
    sample_id = rep(paste0("S", 1:4), 3),
    vst.count = rnorm(12),
    HOG = rep(c("HOG1", "HOG1", "HOG2"), each = 4)  # A1, A2 in HOG1
  )
  data2 <- data.frame(
    abbrev = rep("SP_B", 8),
    gene_id = rep(c("B1", "B2"), each = 4),
    sample_id = rep(paste0("S", 1:4), 2),
    vst.count = rnorm(8),
    HOG = rep(c("HOG1", "HOG2"), each = 4)
  )

  se1 <- build_se(data1, "SP_A")
  se2 <- build_se(data2, "SP_B")

  ortho <- extract_orthologs(se1, se2)

  # HOG1: A1 x B1, A2 x B1 = 2 pairs. HOG2: A3 x B2 = 1 pair.
  hog1_pairs <- ortho[ortho$hog == "HOG1", ]
  expect_equal(nrow(hog1_pairs), 2)
  expect_setequal(hog1_pairs$Species1, c("A1", "A2"))
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
  expect_equal(dim(net_mat$network), dim(net_se$network))
  expect_equal(net_mat$network, net_se$network)
})
