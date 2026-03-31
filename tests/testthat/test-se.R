# Tests for SummarizedExperiment integration

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
