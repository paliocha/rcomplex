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
    hog = paste0("HOG", 1:20),
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
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  expect_type(result, "list")
  expected_names <- c("modules", "module_genes", "n_modules", "modularity",
                      "graph", "method", "params", "resolution_scan")
  expect_equal(length(result), 8L)
  expect_named(result, expected_names)

  # Adaptive threshold output: expected_coclassification in resolution_scan
  expect_true("expected_coclassification" %in%
                names(result$resolution_scan))

  # Consensus iteration count in params
  expect_true("n_consensus_iterations" %in% names(result$params))
  expect_true(is.integer(result$params$n_consensus_iterations))
})


test_that("consensus detect_modules finds correct partition", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

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
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  expect_s3_class(result$resolution_scan, "data.frame")
  expect_equal(nrow(result$resolution_scan), length(resolutions))
  expect_equal(colnames(result$resolution_scan),
               c("resolution", "n_modules", "modularity", "ari_next",
                 "expected_coclassification"))
})


test_that("consensus resolution_scan ari_next values are valid", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

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
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

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
                        objective_function = "modularity", seed = 42,
                        test_k1 = FALSE)
  r2 <- detect_modules(td$net1, method = "leiden",
                        resolution = seq(0.5, 2.0, by = 0.5),
                        objective_function = "modularity", seed = 42,
                        test_k1 = FALSE)

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
  # NULL (adaptive) should NOT error
  expect_no_error(
    detect_modules(td$net1, resolution = c(0.5, 1.0),
                   consensus_threshold = NULL,
                   objective_function = "modularity", seed = 42,
                   test_k1 = FALSE)
  )
})


test_that("consensus graph is original network not co-classification graph", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  expect_true(igraph::is_igraph(result$graph))
  expect_equal(igraph::vcount(result$graph), nrow(td$net1$network))
})


# ---- Adaptive threshold tests ----

test_that("consensus diagnostics in resolution_scan and params", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  # resolution_scan has expected_coclassification column
  expect_true("expected_coclassification" %in%
                names(result$resolution_scan))

  # params has n_consensus_iterations as integer
  expect_true("n_consensus_iterations" %in% names(result$params))
  expect_true(is.integer(result$params$n_consensus_iterations))
  expect_true(result$params$n_consensus_iterations >= 1L)

  # expected_coclassification values are in [0, 1]
  ecc <- result$resolution_scan$expected_coclassification
  expect_true(all(ecc >= 0 & ecc <= 1))

  # More modules -> lower expected co-classification (monotone or at least
  # negatively correlated with n_modules)
  n_mod <- result$resolution_scan$n_modules
  # Only check where n_modules actually varies
  if (length(unique(n_mod)) > 1) {
    expect_true(cor(n_mod, ecc) < 0)
  }
})


test_that("adaptive threshold avoids single-module collapse", {
  td <- make_module_test_data()
  # Broad resolution range that would collapse to 1 module with fixed 0.5

  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.1, 5, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  # Adaptive threshold should preserve module structure
  expect_true(result$n_modules >= 2)
})


test_that("fixed consensus_threshold still works", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           consensus_threshold = 0.3,
                           objective_function = "modularity", seed = 42)

  # Returns a valid partition
  expect_true(result$n_modules >= 1)
  expect_equal(length(result$modules), nrow(td$net1$network))
  expect_equal(result$method, "leiden_consensus")

  # Fixed threshold: n_consensus_iterations should be 0
  expect_equal(result$params$n_consensus_iterations, 0L)
})


test_that("adaptive and fixed threshold can produce different results", {
  td <- make_module_test_data()
  resolutions <- seq(0.1, 5, by = 0.5)

  # Adaptive (default NULL)
  adaptive <- detect_modules(td$net1, method = "leiden",
                              resolution = resolutions,
                              objective_function = "modularity", seed = 42,
                              test_k1 = FALSE)

  # Fixed 0.5
  fixed <- detect_modules(td$net1, method = "leiden",
                           resolution = resolutions,
                           consensus_threshold = 0.5,
                           objective_function = "modularity", seed = 42)

  # Both return valid structures
  expect_equal(length(adaptive$modules), nrow(td$net1$network))
  expect_equal(length(fixed$modules), nrow(td$net1$network))
  expect_true(adaptive$n_modules >= 1)
  expect_true(fixed$n_modules >= 1)

  # Adaptive should generally find more modules than fixed 0.5 on broad range
  expect_true(adaptive$n_modules >= fixed$n_modules)
})


test_that("expected_coclassification relates to module granularity", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 3.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  scan <- result$resolution_scan
  n_mod <- scan$n_modules
  ecc <- scan$expected_coclassification

  # Higher resolution -> more modules -> lower expected co-classification
  # Check negative correlation between n_modules and expected_coclassification
  if (length(unique(n_mod)) > 1) {
    expect_true(cor(n_mod, ecc) < 0)
  }

  # expected_coclassification should be in [0, 1] range
  expect_true(all(ecc >= 0 & ecc <= 1))
})


# ---- Iterative consensus convergence tests ----

test_that("iterative consensus converges in few iterations", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  # Should converge (n_consensus_iterations >= 1)
  expect_true(result$params$n_consensus_iterations >= 1L)
  # Should converge quickly (well below max_consensus_iter default of 10)
  expect_true(result$params$n_consensus_iterations <= 10L)
})


test_that("max_consensus_iter parameter is respected", {
  td <- make_module_test_data()
  # With max_consensus_iter = 1, only one iteration
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           max_consensus_iter = 1L, test_k1 = FALSE)

  expect_equal(result$params$n_consensus_iterations, 1L)
  # Should still produce valid modules
  expect_true(result$n_modules >= 1)
  expect_equal(length(result$modules), nrow(td$net1$network))
})


test_that("consensus iteration recovers planted partition", {
  # Build a network with 5 planted modules of 20 genes each
  set.seed(123)
  n <- 100
  n_mod <- 5
  mod_size <- n / n_mod  # 20
  mat <- matrix(0.05, n, n)  # low background
  for (m in seq_len(n_mod)) {
    idx <- ((m - 1) * mod_size + 1):(m * mod_size)
    mat[idx, idx] <- 0.9
  }
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", seq_len(n))
  net <- list(network = mat, threshold = 0.1)

  # Run consensus with broad CPM resolution range
  result <- detect_modules(net, method = "leiden",
                           resolution = seq(0.1, 3.0, by = 0.3),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  # Should find approximately 5 modules
  expect_true(result$n_modules >= 4 && result$n_modules <= 6)

  # Compute ARI against ground truth
  truth <- rep(seq_len(n_mod), each = mod_size)
  names(truth) <- paste0("G", seq_len(n))
  ari <- igraph::compare(result$modules[names(truth)], truth,
                         method = "adjusted.rand")
  expect_true(ari > 0.9)
})


# ---- Sparse co-classification tests ----

test_that("build_sparse_coclassification_cpp matches dense on original edges", {
  set.seed(42)
  n <- 20
  mat <- matrix(0, n, n)
  rownames(mat) <- colnames(mat) <- paste0("G", seq_len(n))

  mat[1:7, 1:7] <- 0.8
  mat[8:14, 8:14] <- 0.8
  mat[15:20, 15:20] <- 0.8
  diag(mat) <- 1
  net <- list(network = mat, threshold = 0.3)

  adj <- mat
  adj[adj < 0.3] <- 0
  g <- igraph::graph_from_adjacency_matrix(
    adj, mode = "upper", weighted = TRUE, diag = FALSE
  )

  resolutions <- seq(0.5, 2.0, by = 0.5)
  memberships <- lapply(resolutions, function(res) {
    comm <- igraph::cluster_leiden(g, resolution = res,
                                    objective_function = "modularity",
                                    n_iterations = 2L)
    mem <- igraph::membership(comm)
    names(mem) <- igraph::V(g)$name
    mem
  })

  dense_result <- build_coclassification_cpp(memberships, n, TRUE)

  el_0 <- igraph::as_edgelist(g, names = FALSE) - 1L
  storage.mode(el_0) <- "integer"
  sparse_result <- build_sparse_coclassification_cpp(memberships, n, el_0)

  el_1 <- igraph::as_edgelist(g, names = FALSE)
  dense_at_edges <- dense_result$excess_coclassification[el_1]

  expect_equal(sparse_result$excess, dense_at_edges, tolerance = 1e-12)
  expect_equal(sparse_result$expected, dense_result$expected, tolerance = 1e-12)
})


test_that("sparse consensus produces modules", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  expect_true(result$n_modules > 1)
  expect_equal(result$method, "leiden_consensus")
  expect_equal(length(result$modules), nrow(td$net1$network))
  expect_true(all(names(td$net1$network[, 1]) %in% names(result$modules)))
})


test_that("sparse_excess_spectral_norm_cpp returns non-negative scalar", {
  set.seed(42)
  n <- 20
  mat <- matrix(0, n, n)
  rownames(mat) <- colnames(mat) <- paste0("G", seq_len(n))
  mat[1:10, 1:10] <- 0.8
  mat[11:20, 11:20] <- 0.8
  diag(mat) <- 1

  adj <- mat
  adj[adj < 0.3] <- 0
  g <- igraph::graph_from_adjacency_matrix(
    adj, mode = "upper", weighted = TRUE, diag = FALSE
  )

  memberships <- lapply(c(0.5, 1.0, 1.5), function(res) {
    comm <- igraph::cluster_leiden(g, resolution = res,
                                    objective_function = "modularity",
                                    n_iterations = 2L)
    mem <- igraph::membership(comm)
    names(mem) <- igraph::V(g)$name
    mem
  })

  el_0 <- igraph::as_edgelist(g, names = FALSE) - 1L
  storage.mode(el_0) <- "integer"
  lambda <- sparse_excess_spectral_norm_cpp(memberships, n, el_0)

  expect_true(is.numeric(lambda))
  expect_equal(length(lambda), 1L)
  expect_true(lambda >= 0)
})


# ---- K = 1 community structure tests ----

test_that("K = 1 test rejects structure on Erdos-Renyi graph", {
  skip_on_cran()

  set.seed(123)
  n <- 50
  g_er <- igraph::sample_gnp(n, p = 0.1, directed = FALSE)
  igraph::V(g_er)$name <- paste0("G", seq_len(n))
  igraph::E(g_er)$weight <- stats::runif(igraph::ecount(g_er), 0.3, 1.0)

  adj <- igraph::as_adjacency_matrix(g_er, attr = "weight", sparse = FALSE)
  net <- list(network = adj, threshold = 0.01)

  result <- detect_modules(net, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = TRUE, n_perm_k1 = 100L)

  expect_equal(result$n_modules, 1L)
  expect_true(!is.null(result$k1_test))
  expect_true(result$k1_test$p_value > 0.05)
  # No-structure early stop: should stop before 100 perms
  expect_true(result$k1_test$n_perm_completed < 100L)
  expect_equal(length(result$k1_test$lambda_null),
               result$k1_test$n_perm_completed)
})


test_that("K = 1 test accepts structure on planted partition", {
  skip_on_cran()

  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = TRUE, n_perm_k1 = 100L)

  expect_true(result$n_modules > 1)
  expect_true(!is.null(result$k1_test))
  expect_true(result$k1_test$p_value < 0.05)
  # Early stopping: clear structure should stop at ceil(1/0.05) = 20 perms
  expect_true(result$k1_test$n_perm_completed <= 20L)
  expect_equal(length(result$k1_test$lambda_null),
               result$k1_test$n_perm_completed)
})


test_that("test_k1 = FALSE skips the K = 1 test", {
  td <- make_module_test_data()
  result <- detect_modules(td$net1, method = "leiden",
                           resolution = seq(0.5, 2.0, by = 0.5),
                           objective_function = "modularity", seed = 42,
                           test_k1 = FALSE)

  expect_null(result$k1_test)
  expect_equal(length(result), 8L)
  expect_true(result$n_modules >= 1)
  expect_equal(result$method, "leiden_consensus")
})


# ---- coarsen_modules() tests ----

# Helper: asymmetric module structure (sp1: 2 large, sp2: 4 small)
make_asymmetric_module_data <- function() {
  n <- 40
  mat1 <- matrix(0, n, n)
  mat2 <- matrix(0, n, n)
  rownames(mat1) <- colnames(mat1) <- paste0("A", 1:n)
  rownames(mat2) <- colnames(mat2) <- paste0("B", 1:n)

  mat1[1:20, 1:20] <- 0.9
  mat1[21:40, 21:40] <- 0.9

  mat2[1:10, 1:10] <- 0.9
  mat2[11:20, 11:20] <- 0.9
  mat2[21:30, 21:30] <- 0.9
  mat2[31:40, 31:40] <- 0.9
  # Weak inter-module edges within parent groups (below intra but above threshold)
  mat2[1:10, 11:20] <- mat2[11:20, 1:10] <- 0.35
  mat2[21:30, 31:40] <- mat2[31:40, 21:30] <- 0.35

  diag(mat1) <- diag(mat2) <- 1

  net1 <- list(network = mat1, threshold = 0.5)
  net2 <- list(network = mat2, threshold = 0.3)

  orthologs <- data.frame(
    Species1 = paste0("A", 1:40),
    Species2 = paste0("B", 1:40),
    hog = paste0("HOG", 1:40),
    stringsAsFactors = FALSE
  )

  list(net1 = net1, net2 = net2, orthologs = orthologs)
}


test_that("coarsen_modules returns correct structure", {
  td <- make_module_test_data()
  mods <- detect_modules(td$net1, method = "leiden",
                         objective_function = "modularity", seed = 42)
  coarse <- coarsen_modules(mods, target_n_modules = 2L)

  expect_type(coarse, "list")
  expect_named(coarse, c("modules", "module_genes", "n_modules",
                          "modularity", "graph", "method", "params",
                          "merge_map", "merge_dendrogram"),
               ignore.order = TRUE)
  expect_true(igraph::is_igraph(coarse$graph))
  expect_s3_class(coarse$merge_dendrogram, "hclust")
  expect_true(is.data.frame(coarse$merge_map))
  expect_true(all(c("original_module", "coarsened_module") %in%
                    names(coarse$merge_map)))
})


test_that("coarsen_modules produces target number of modules", {
  td <- make_module_test_data()
  mods <- detect_modules(td$net1, method = "leiden",
                         objective_function = "modularity", seed = 42)
  for (target in c(2L, 1L)) {
    coarse <- coarsen_modules(mods, target_n_modules = target)
    expect_equal(coarse$n_modules, target)
    expect_equal(length(coarse$module_genes), target)
    expect_equal(length(unique(coarse$modules)), target)
  }
})


test_that("coarsen_modules merges related modules first", {
  td <- make_asymmetric_module_data()
  # resolution = 2 needed to split the 4 sub-modules at this density
  mods <- detect_modules(td$net2, method = "leiden", resolution = 2.0,
                         objective_function = "modularity", seed = 42)
  expect_true(mods$n_modules >= 3)

  coarse <- coarsen_modules(mods, target_n_modules = 2L)

  # Genes 1-20 should be in one module, 21-40 in the other
  mods_1_20 <- coarse$modules[paste0("B", 1:20)]
  mods_21_40 <- coarse$modules[paste0("B", 21:40)]
  expect_equal(length(unique(mods_1_20)), 1)
  expect_equal(length(unique(mods_21_40)), 1)
  expect_true(unique(mods_1_20) != unique(mods_21_40))
})


test_that("coarsen_modules preserves all genes", {
  td <- make_module_test_data()
  mods <- detect_modules(td$net1, method = "leiden",
                         objective_function = "modularity", seed = 42)
  coarse <- coarsen_modules(mods, target_n_modules = 2L)

  expect_equal(sort(names(coarse$modules)), sort(names(mods$modules)))
  all_genes <- unlist(coarse$module_genes, use.names = FALSE)
  expect_equal(sort(all_genes), sort(names(mods$modules)))
})


test_that("coarsen_modules output works with compare_modules", {
  td <- make_asymmetric_module_data()
  m1 <- detect_modules(td$net1, method = "leiden",
                        objective_function = "modularity", seed = 42)
  m2 <- detect_modules(td$net2, method = "leiden", resolution = 2.0,
                        objective_function = "modularity", seed = 42)
  coarse2 <- coarsen_modules(m2, target_n_modules = m1$n_modules)

  comp <- compare_modules(m1, coarse2, td$orthologs,
                           method = "hypergeometric")
  expect_true(is.data.frame(comp$pairs))
  expect_true(nrow(comp$pairs) > 0)
  expect_true(nrow(comp$best_matches) > 0)
})


test_that("coarsen_modules validates inputs", {
  td <- make_module_test_data()
  mods <- detect_modules(td$net1, method = "leiden",
                         objective_function = "modularity", seed = 42)

  expect_error(coarsen_modules(mods, target_n_modules = 0L),
               "target_n_modules must be >= 1")
  expect_error(coarsen_modules(list(x = 1), target_n_modules = 2L),
               "modules must be output from detect_modules")
})


test_that("coarsen_modules is no-op when target equals current", {
  td <- make_module_test_data()
  mods <- detect_modules(td$net1, method = "leiden",
                         objective_function = "modularity", seed = 42)
  result <- coarsen_modules(mods, target_n_modules = mods$n_modules)
  expect_identical(result, mods)
})


test_that("coarsening reduces false species-specific calls", {
  td <- make_asymmetric_module_data()
  m1 <- detect_modules(td$net1, method = "leiden",
                        objective_function = "modularity", seed = 42)
  m2 <- detect_modules(td$net2, method = "leiden", resolution = 2.0,
                        objective_function = "modularity", seed = 42)

  # Natural scale: some sp2 modules appear species-specific
  comp_natural <- compare_modules(m1, m2, td$orthologs,
                                   method = "hypergeometric")
  class_natural <- classify_modules(comp_natural)

  # Matched scale: coarsen sp2 to match sp1
  coarse2 <- coarsen_modules(m2, target_n_modules = m1$n_modules)
  comp_matched <- compare_modules(m1, coarse2, td$orthologs,
                                   method = "hypergeometric")
  class_matched <- classify_modules(comp_matched)

  # Matched scale should have fewer species-specific modules
  n_ss_natural <- sum(class_natural$classification == "species_specific")
  n_ss_matched <- sum(class_matched$classification == "species_specific")
  expect_true(n_ss_matched <= n_ss_natural)
})
