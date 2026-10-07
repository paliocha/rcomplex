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
  data <- data.frame(
    abbrev = "SP_B", gene_id = "G1",
    sample_id = "S1", vst.count = 1.0
  )
  expect_error(build_se(data, "SP_A"), "No rows found")
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


# --- SE rowData to the long table for prepare_orthologs ---

test_that("SE rowData as the long table feeds find_coexpressologs", {
  set.seed(7)
  hogs <- paste0("HOG", 1:30)
  samples <- paste0("S", 1:10)
  se_a <- build_se(
    make_long_data("SP_A", paste0("A", 1:30), samples, hogs), "SP_A"
  )
  se_b <- build_se(
    make_long_data("SP_B", paste0("B", 1:30), samples, hogs), "SP_B"
  )
  long <- do.call(rbind, lapply(list(SP_A = se_a, SP_B = se_b), function(se) {
    data.frame(
      species = se@metadata$species, gene = rownames(se),
      hog = SummarizedExperiment::rowData(se)$hog
    )
  }))
  ortho <- prepare_orthologs(long)
  expect_setequal(paste(ortho$gene1, ortho$gene2, ortho$hog),
                  paste0("A", 1:30, " B", 1:30, " HOG", 1:30))
  nets <- list(
    SP_A = compute_network(se_a, density = 0.1),
    SP_B = compute_network(se_b, density = 0.1)
  )

  edges <- find_coexpressologs(nets, ortho)

  expect_s3_class(edges, "data.frame")
  expect_true(all(
    c("gene1", "gene2", "hog", "q_value", "effect_size") %in% names(edges)
  ))
})
