# Tests for detect_modules(), compare_modules(), classify_modules()

# Helper: build two test networks with clear module structure and orthologs
make_module_test_data <- function() {
  n <- 30
  mat1 <- matrix(0, n, n)
  mat2 <- matrix(0, n, n)
  rownames(mat1) <- colnames(mat1) <- paste0("A", 1:n)
  rownames(mat2) <- colnames(mat2) <- paste0("B", 1:n)

  # Module 1 (genes 1-10): conserved in both species
  mat1[1:10, 1:10] <- mat2[1:10, 1:10] <- 0.9

  # Module 2 (genes 11-20): conserved in both species
  mat1[11:20, 11:20] <- mat2[11:20, 11:20] <- 0.9

  # Module 3 in sp1 only (genes 21-25)
  mat1[21:25, 21:25] <- 0.9

  # Module 3 in sp2 is different genes (26-30)
  mat2[26:30, 26:30] <- 0.9

  diag(mat1) <- diag(mat2) <- 1

  net1 <- list(network = mat1, threshold = 0.5)
  net2 <- list(network = mat2, threshold = 0.5)

  # 1:1 orthologs for genes 1-20; no orthologs for 21-30
  orthologs <- data.frame(
    Species1 = paste0("A", 1:20),
    Species2 = paste0("B", 1:20),
    OrthoGroup = paste0("HOG", 1:20),
    stringsAsFactors = FALSE
  )

  list(net1 = net1, net2 = net2, orthologs = orthologs)
}


# Higher-level helper: also detects modules
make_module_test_modules <- function() {
  td <- make_module_test_data()
  m1 <- detect_modules(td$net1, method = "leiden",
                        objective_function = "modularity", seed = 42)
  m2 <- detect_modules(td$net2, method = "leiden",
                        objective_function = "modularity", seed = 42)
  list(net1 = td$net1, net2 = td$net2, orthologs = td$orthologs,
       m1 = m1, m2 = m2)
}


# ---- detect_modules() tests ----

test_that("detect_modules returns correct structure", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           objective_function = "modularity", seed = 42)

  expect_type(result, "list")
  expect_named(result, c("modules", "module_genes", "n_modules",
                         "modularity", "graph", "method", "params"))
  expect_true(igraph::is_igraph(result$graph))
  expect_equal(result$method, "leiden")
  expect_true(result$n_modules >= 2)
})


test_that("detect_modules assigns all genes", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           objective_function = "modularity", seed = 42)

  expect_equal(length(result$modules), nrow(td$net1$network))
  expect_equal(sort(names(result$modules)),
               sort(rownames(td$net1$network)))
})


test_that("detect_modules module_genes is consistent with modules", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           objective_function = "modularity", seed = 42)

  all_genes <- unlist(result$module_genes, use.names = FALSE)
  expect_equal(sort(all_genes), sort(names(result$modules)))

  for (mod_id in names(result$module_genes)) {
    genes <- result$module_genes[[mod_id]]
    expect_true(all(result$modules[genes] == as.integer(mod_id)))
  }
})


test_that("detect_modules finds expected modules", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           objective_function = "modularity", seed = 42)

  mods_1_10 <- result$modules[paste0("A", 1:10)]
  expect_equal(length(unique(mods_1_10)), 1)

  mods_11_20 <- result$modules[paste0("A", 11:20)]
  expect_equal(length(unique(mods_11_20)), 1)

  expect_true(unique(mods_1_10) != unique(mods_11_20))
})


test_that("detect_modules works with infomap", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "infomap", seed = 42)

  expect_equal(result$method, "infomap")
  expect_true(result$n_modules >= 2)
  expect_equal(length(result$modules), nrow(td$net1$network))
})


test_that("detect_modules infomap finds expected modules", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "infomap", seed = 42)

  mods_1_10 <- result$modules[paste0("A", 1:10)]
  expect_equal(length(unique(mods_1_10)), 1)

  mods_11_20 <- result$modules[paste0("A", 11:20)]
  expect_equal(length(unique(mods_11_20)), 1)

  expect_true(unique(mods_1_10) != unique(mods_11_20))
})


test_that("detect_modules works with sbm", {
  skip_if_not_installed("sbm")
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "sbm", seed = 42)

  expect_equal(result$method, "sbm")
  expect_true(result$n_modules >= 2)
  expect_equal(length(result$modules), nrow(td$net1$network))
  expect_true(igraph::is_igraph(result$graph))
  expect_true(!is.null(result$params$ICL))
})


test_that("detect_modules sbm finds expected modules", {
  skip_if_not_installed("sbm")
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "sbm", seed = 42)

  mods_1_10 <- result$modules[paste0("A", 1:10)]
  expect_equal(length(unique(mods_1_10)), 1)

  mods_11_20 <- result$modules[paste0("A", 11:20)]
  expect_equal(length(unique(mods_11_20)), 1)

  expect_true(unique(mods_1_10) != unique(mods_11_20))
})


test_that("detect_modules sbm requires sbm package", {
  skip_if(requireNamespace("sbm", quietly = TRUE),
          "sbm is installed; cannot test missing-package error")
  td <- make_module_test_data()
  expect_error(detect_modules(td$net1, method = "sbm"),
               "Package 'sbm' is required")
})


test_that("detect_modules validates input", {
  expect_error(detect_modules("not_a_net"),
               "net must be a network object")
  expect_error(detect_modules(list(something = 1)),
               "net must be a network object")
})


test_that("detect_modules errors on no edges", {
  mat <- matrix(0, 5, 5)
  rownames(mat) <- colnames(mat) <- paste0("G", 1:5)
  diag(mat) <- 1
  net <- list(network = mat, threshold = 0.5)

  expect_error(detect_modules(net), "No edges above threshold")
})


test_that("detect_modules modularity is numeric", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           objective_function = "modularity", seed = 42)
  expect_true(is.numeric(result$modularity))
})


test_that("detect_modules seed produces reproducible results", {
  td <- make_module_test_data()
  r1 <- detect_modules(td$net1, method = "leiden",
                        objective_function = "modularity", seed = 123)
  r2 <- detect_modules(td$net1, method = "leiden",
                        objective_function = "modularity", seed = 123)
  expect_identical(r1$modules, r2$modules)
})


# ---- compare_modules() hypergeometric tests ----

test_that("compare_modules hypergeometric returns correct structure", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")

  expect_type(result, "list")
  expect_named(result, c("pairs", "best_matches"))
  expect_s3_class(result$pairs, "data.frame")
  expect_true(all(c("module_sp1", "module_sp2", "size_sp1", "size_sp2",
                     "overlap", "jaccard", "p.value", "q.value") %in%
                    names(result$pairs)))
})


test_that("compare_modules produces correct number of pairs", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")
  expect_equal(nrow(result$pairs), td$m1$n_modules * td$m2$n_modules)
})


test_that("compare_modules finds conserved module pairs", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")

  mod_a1 <- as.character(td$m1$modules["A1"])
  mod_b1 <- as.character(td$m2$modules["B1"])
  pair <- result$pairs[result$pairs$module_sp1 == mod_a1 &
                         result$pairs$module_sp2 == mod_b1, ]

  expect_true(pair$overlap > 0)
  expect_true(pair$p.value < 0.05)
  expect_true(pair$jaccard > 0)
})


test_that("compare_modules p-values are valid", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")

  expect_true(all(result$pairs$p.value >= 0 & result$pairs$p.value <= 1))
  expect_true(all(result$pairs$q.value >= 0 & result$pairs$q.value <= 1))
})


test_that("compare_modules best_matches has both species", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")

  expect_true("sp1" %in% result$best_matches$species)
  expect_true("sp2" %in% result$best_matches$species)
  expect_equal(sum(result$best_matches$species == "sp1"), td$m1$n_modules)
  expect_equal(sum(result$best_matches$species == "sp2"), td$m2$n_modules)
})


test_that("compare_modules validates inputs", {
  td <- make_module_test_modules()

  expect_error(compare_modules("not_modules", td$m1, td$orthologs),
               "modules1 must be output from detect_modules")
  expect_error(compare_modules(td$m1, "not_modules", td$orthologs),
               "modules2 must be output from detect_modules")
  expect_error(compare_modules(td$m1, td$m1, data.frame(x = 1)),
               "orthologs must have columns")
})


test_that("compare_modules jaccard values are in [0, 1]", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")
  expect_true(all(result$pairs$jaccard >= 0 & result$pairs$jaccard <= 1))
})


test_that("compare_modules overlap <= min(size_sp1, size_sp2)", {
  td <- make_module_test_modules()
  result <- compare_modules(td$m1, td$m2, td$orthologs,
                            method = "hypergeometric")

  for (i in seq_len(nrow(result$pairs))) {
    expect_true(result$pairs$overlap[i] <=
                  min(result$pairs$size_sp1[i], result$pairs$size_sp2[i]))
  }
})


# ---- compare_modules() jaccard + permutation tests ----

test_that("compare_modules jaccard returns correct structure", {
  td <- make_module_test_modules()

  set.seed(42)
  result <- compare_modules(td$m1, td$m2, td$orthologs, method = "jaccard",
                            max_permutations = 200L, min_exceedances = 10L)

  expect_type(result, "list")
  expect_named(result, c("pairs", "best_matches"))
  expect_true(all(c("module_sp1", "module_sp2", "overlap", "jaccard",
                     "p.value", "q.value") %in% names(result$pairs)))
})


test_that("compare_modules jaccard finds conserved modules", {
  td <- make_module_test_modules()

  set.seed(42)
  result <- compare_modules(td$m1, td$m2, td$orthologs, method = "jaccard",
                            max_permutations = 500L, min_exceedances = 20L)

  mod_a1 <- as.character(td$m1$modules["A1"])
  mod_b1 <- as.character(td$m2$modules["B1"])
  pair <- result$pairs[result$pairs$module_sp1 == mod_a1 &
                         result$pairs$module_sp2 == mod_b1, ]

  expect_equal(nrow(pair), 1L)
  expect_true(pair$p.value < 0.5)
})


test_that("compare_modules jaccard p-values are valid", {
  td <- make_module_test_modules()

  set.seed(42)
  result <- compare_modules(td$m1, td$m2, td$orthologs, method = "jaccard",
                            max_permutations = 200L, min_exceedances = 10L)

  expect_true(all(result$pairs$p.value > 0 & result$pairs$p.value <= 1))
})


# ---- classify_modules() tests ----

test_that("classify_modules returns correct structure", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")
  result <- classify_modules(comp)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("module", "species", "classification", "best_match",
                     "best_jaccard", "best_q", "n_significant") %in%
                    names(result)))
})


test_that("classify_modules has correct number of rows", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")
  result <- classify_modules(comp)
  expect_equal(nrow(result), td$m1$n_modules + td$m2$n_modules)
})


test_that("classify_modules uses three categories", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")
  result <- classify_modules(comp)

  valid_classes <- c("conserved", "partially_conserved", "species_specific")
  expect_true(all(result$classification %in% valid_classes))
})


test_that("classify_modules identifies conserved modules", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")
  result <- classify_modules(comp)

  mod_a1 <- as.character(td$m1$modules["A1"])
  sp1_result <- result[result$species == "sp1" & result$module == mod_a1, ]
  expect_equal(sp1_result$classification, "conserved")
})


test_that("classify_modules respects alpha parameter", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")

  strict <- classify_modules(comp, alpha = 1e-10)
  lenient <- classify_modules(comp, alpha = 0.5)

  expect_true(sum(strict$classification == "species_specific") >=
                sum(lenient$classification == "species_specific"))
})


test_that("classify_modules respects jaccard_threshold", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")

  high_j <- classify_modules(comp, jaccard_threshold = 0.99)
  low_j <- classify_modules(comp, jaccard_threshold = 0.01)

  expect_true(sum(high_j$classification == "conserved") <=
                sum(low_j$classification == "conserved"))
})


test_that("classify_modules validates input", {
  expect_error(classify_modules("not_a_comparison"),
               "module_comparison must be output from compare_modules")
  expect_error(classify_modules(list(x = 1)),
               "module_comparison must be output from compare_modules")
})


test_that("classify_modules handles empty comparison", {
  empty_comp <- list(
    pairs = data.frame(
      module_sp1 = character(0), module_sp2 = character(0),
      size_sp1 = integer(0), size_sp2 = integer(0),
      overlap = integer(0), jaccard = numeric(0),
      p.value = numeric(0), q.value = numeric(0),
      stringsAsFactors = FALSE
    ),
    best_matches = data.frame(
      module = character(0), species = character(0),
      best_match = character(0), overlap = integer(0),
      jaccard = numeric(0), p.value = numeric(0),
      q.value = numeric(0), n_significant = integer(0),
      stringsAsFactors = FALSE
    )
  )
  result <- classify_modules(empty_comp)
  expect_equal(nrow(result), 0)
  expect_true(all(c("module", "species", "classification") %in% names(result)))
})


test_that("classify_modules best_jaccard matches expectation", {
  td <- make_module_test_modules()
  comp <- compare_modules(td$m1, td$m2, td$orthologs,
                          method = "hypergeometric")
  result <- classify_modules(comp)

  expect_true(all(result$best_jaccard >= 0 & result$best_jaccard <= 1))
  expect_true(all(result$best_q >= 0 & result$best_q <= 1))
})


# ---- Consensus module detection tests ----

test_that("consensus detect_modules returns correct structure with 8 elements", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42)

  expect_type(result, "list")
  expected_names <- c("modules", "module_genes", "n_modules", "modularity",
                      "graph", "method", "params", "resolution_scan")
  expect_equal(length(result), 8L)
  expect_named(result, expected_names)
})


test_that("consensus detect_modules finds correct partition", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42)

  # Genes 1-10 should be in the same module
  mods_1_10 <- result$modules[paste0("A", 1:10)]
  expect_equal(length(unique(mods_1_10)), 1)

  # Genes 11-20 should be in the same module
  mods_11_20 <- result$modules[paste0("A", 11:20)]
  expect_equal(length(unique(mods_11_20)), 1)

  # The two groups should be in different modules
  expect_true(unique(mods_1_10) != unique(mods_11_20))
})


test_that("consensus resolution_scan has correct dimensions and columns", {
  td <- make_module_test_data()
  resolutions <- seq(0.5, 2.0, by = 0.5)
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = resolutions,
                           objective_function = "modularity", seed = 42)

  expect_s3_class(result$resolution_scan, "data.frame")
  expect_equal(nrow(result$resolution_scan), length(resolutions))
  expect_equal(colnames(result$resolution_scan),
               c("resolution", "n_modules", "modularity", "ari_next"))
})


test_that("consensus resolution_scan ari_next values are valid", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42)

  ari_vals <- result$resolution_scan$ari_next
  non_na <- ari_vals[!is.na(ari_vals)]
  expect_true(all(non_na >= -1 & non_na <= 1))
  expect_true(is.na(ari_vals[length(ari_vals)]))
})


test_that("scalar resolution still returns 7-element list without resolution_scan", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = 1.0,
                           objective_function = "modularity", seed = 42)

  expect_equal(length(result), 7L)
  expect_null(result$resolution_scan)
  expect_equal(result$method, "leiden")
})


test_that("consensus detect_modules returns method leiden_consensus", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42)

  expect_equal(result$method, "leiden_consensus")
})


test_that("consensus mode errors for non-leiden methods", {
  td <- make_module_test_data()
  expect_error(
    detect_modules(td$net1, method = "infomap", resolution = c(0.5, 1.0)),
    "Consensus mode"
  )
})


test_that("consensus detect_modules is reproducible with same seed", {
  td <- make_module_test_data()
  r1 <- detect_modules(td$net1, method = "leiden",
                        resolution = seq(0.5, 2.0, by = 0.5),
                        objective_function = "modularity", seed = 42)
  r2 <- detect_modules(td$net1, method = "leiden",
                        resolution = seq(0.5, 2.0, by = 0.5),
                        objective_function = "modularity", seed = 42)

  expect_identical(r1$modules, r2$modules)
})


test_that("single-element vector resolution behaves like scalar", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = c(1.0),
                           objective_function = "modularity", seed = 42)

  expect_equal(length(result), 7L)
  expect_null(result$resolution_scan)
  expect_equal(result$method, "leiden")
})


test_that("consensus_threshold validation rejects 0 and 1", {
  td <- make_module_test_data()
  expect_error(
    detect_modules(td$net1, resolution = c(0.5, 1.0),
                   consensus_threshold = 0),
    "(0, 1)", fixed = TRUE
  )
  expect_error(
    detect_modules(td$net1, resolution = c(0.5, 1.0),
                   consensus_threshold = 1),
    "(0, 1)", fixed = TRUE
  )
})


test_that("consensus graph is original network not co-classification graph", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42)

  expect_true(igraph::is_igraph(result$graph))
  expect_equal(igraph::vcount(result$graph), nrow(td$net1$network))
})
