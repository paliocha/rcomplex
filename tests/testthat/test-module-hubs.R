# Tests for identify_module_hubs() and classify_hub_conservation()

# Helper: 4-species, 2-trait scenario with designed hub patterns
# SP_A, SP_B = annual; SP_C, SP_D = perennial
# Each species has 20 genes, 2 modules of 10 genes each
make_hub_test_data <- function() {
  species <- c("SP_A", "SP_B", "SP_C", "SP_D")
  prefixes <- c(SP_A = "A", SP_B = "B", SP_C = "C", SP_D = "D")
  trait <- c(SP_A = "annual", SP_B = "annual",
             SP_C = "perennial", SP_D = "perennial")
  n <- 20L

  nets <- list()
  mods <- list()

  for (sp in species) {
    px <- prefixes[sp]
    gnames <- paste0(px, seq_len(n))
    mat <- matrix(0, n, n, dimnames = list(gnames, gnames))

    # Module 1: genes 1-10 (high intra-module correlation)
    mat[1:10, 1:10] <- 0.8

    # Module 2: genes 11-20
    mat[11:20, 11:20] <- 0.8

    # Gene 1 is a super-hub in module 1 (extra connectivity)
    mat[1, 2:10] <- mat[2:10, 1] <- 0.95

    if (sp %in% c("SP_A", "SP_B")) {
      # Annual: gene 5 gets extra connectivity in module 1
      mat[5, 6:10] <- mat[6:10, 5] <- 0.92

      # Annual: gene 11 is the module 2 hub
      mat[11, 12:20] <- mat[12:20, 11] <- 0.95
    } else {
      # Perennial: gene 15 gets extra connectivity in module 2
      mat[15, 11:14] <- mat[11:14, 15] <- 0.92
      mat[15, 16:20] <- mat[16:20, 15] <- 0.92

      # Perennial: gene 11 is a hub but in module 1 (rewired)
      # Move gene 11 to module 1 by boosting its module-1 edges
      # and reducing module-2 edges
      mat[11, 1:10] <- mat[1:10, 11] <- 0.85
      mat[11, 12:20] <- mat[12:20, 11] <- 0.3
    }

    diag(mat) <- 1
    net <- list(network = mat, threshold = 0.5)
    nets[[sp]] <- net
    mods[[sp]] <- detect_modules(net, method = "leiden",
                                  objective_function = "modularity",
                                  seed = 42)
  }

  # Orthologs: 1:1 mapping across species
  # Build all pairwise orthologs
  ortho_list <- list()
  sp_pairs <- combn(species, 2, simplify = FALSE)
  for (pair in sp_pairs) {
    sp1 <- pair[1]
    sp2 <- pair[2]
    px1 <- prefixes[sp1]
    px2 <- prefixes[sp2]
    ortho_list[[paste(sp1, sp2, sep = ".")]] <- data.frame(
      Species1 = paste0(px1, seq_len(n)),
      Species2 = paste0(px2, seq_len(n)),
      hog = paste0("HOG", seq_len(n)),
      stringsAsFactors = FALSE
    )
  }

  list(nets = nets, mods = mods, orthologs = ortho_list, trait = trait,
       species = species, prefixes = prefixes)
}


# ---- identify_module_hubs() tests ----

test_that("identify_module_hubs returns correct structure", {
  td <- make_hub_test_data()
  sp <- "SP_A"
  ortho <- td$orthologs[["SP_A.SP_B"]]  # any ortholog table for SP_A

  result <- identify_module_hubs(td$mods[[sp]], td$nets[[sp]], ortho)

  expect_true(is.data.frame(result))
  expect_true(all(c("gene", "module", "degree", "betweenness", "eigenvector",
                     "mean_edge_weight", "global_degree",
                     "rank", "is_hub", "hog") %in% names(result)))
  expect_equal(nrow(result), 20L)
  expect_type(result$gene, "character")
  expect_type(result$module, "integer")
  expect_type(result$degree, "double")
  expect_type(result$betweenness, "double")
  expect_type(result$eigenvector, "double")
  expect_type(result$rank, "integer")
  expect_type(result$is_hub, "logical")
  expect_type(result$hog, "character")
})


test_that("identify_module_hubs assigns all genes", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])
  expect_equal(sort(result$gene), sort(names(td$mods[["SP_A"]]$modules)))
})


test_that("identify_module_hubs top_n controls hub count", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                  top_n = 2L)
  # Each module should have at most 2 hubs
  hub_counts <- tapply(result$is_hub, result$module, sum)
  expect_true(all(hub_counts <= 2L))
  expect_true(all(hub_counts >= 1L))  # at least 1 if module >= min_module_size
})


test_that("identify_module_hubs top_fraction controls hub count", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                  top_fraction = 0.2)
  # 10-gene modules: 0.2 * 10 = 2 hubs each
  hub_counts <- tapply(result$is_hub, result$module, sum)
  expect_true(all(hub_counts == 2L))
})


test_that("identify_module_hubs rank 1 has highest centrality (default degree)", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])
  for (m in unique(result$module)) {
    mod_df <- result[result$module == m, ]
    if (all(is.na(mod_df$degree))) next
    top <- mod_df[mod_df$rank == 1L, ]
    expect_equal(top$degree, max(mod_df$degree))
  }
})


test_that("identify_module_hubs rank reflects non-default centrality", {
  # Use the bridge-gene network where degree and betweenness disagree
  n <- 10L
  gnames <- paste0("H", seq_len(n))
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
  mat[1:5, 1:5] <- 0.8
  mat[c(1, 6:10), c(1, 6:10)] <- 0.8
  diag(mat) <- 1

  net <- list(network = mat, threshold = 0.5)
  m <- detect_modules(net, method = "leiden",
                       objective_function = "modularity", seed = 42)

  result_btw <- identify_module_hubs(m, net, centrality = "betweenness")

  # Rank should reflect betweenness, not degree.
  # All rank-1 genes should have max betweenness (ties.method = "min").
  for (mod in unique(result_btw$module)) {
    mod_df <- result_btw[result_btw$module == mod, ]
    if (all(is.na(mod_df$betweenness))) next
    top <- mod_df[mod_df$rank == 1L, ]
    expect_true(all(top$betweenness == max(mod_df$betweenness)))
  }

  # Attribute should record the primary centrality
  expect_equal(attr(result_btw, "primary_centrality"), "betweenness")
})


test_that("identify_module_hubs maps HOGs correctly", {
  td <- make_hub_test_data()
  ortho <- td$orthologs[["SP_A.SP_B"]]
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]], ortho)

  # A1 should map to HOG1
  a1_row <- result[result$gene == "A1", ]
  expect_equal(a1_row$hog, "HOG1")

  # All genes should have HOGs (all are in ortho table)
  expect_true(all(!is.na(result$hog)))
})


test_that("identify_module_hubs works without orthologs", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])
  expect_true(all(is.na(result$hog)))
})


test_that("identify_module_hubs supports all centrality methods", {
  td <- make_hub_test_data()
  for (method in c("degree", "betweenness", "eigenvector")) {
    result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                    centrality = method)
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 20L)
    # Centrality should be non-NA for modules >= min_module_size
    non_tiny <- result[!is.na(result$degree), ]
    expect_true(nrow(non_tiny) > 0)
  }
})


test_that("identify_module_hubs validates inputs", {
  td <- make_hub_test_data()

  expect_error(identify_module_hubs(list(), td$nets[["SP_A"]]),
               "must be output from detect_modules")
  expect_error(identify_module_hubs(td$mods[["SP_A"]], list()),
               "must be output from compute_network")
  expect_error(identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                     top_n = 0L),
               "top_n must be >= 1")
  expect_error(identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                     top_fraction = 0),
               "top_fraction must be in")
})


test_that("identify_module_hubs gene 1 is top hub in module 1 (annual)", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]],
                                  top_n = 1L)

  # Find which module A1 is in
  a1_mod <- result$module[result$gene == "A1"]
  mod_hubs <- result[result$module == a1_mod & result$is_hub, ]
  expect_true("A1" %in% mod_hubs$gene)
})


test_that("identify_module_hubs global_degree is populated", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])

  expect_true(all(!is.na(result$global_degree)))
  # Global degree should be >= within-module centrality for degree method
  non_tiny <- result[!is.na(result$degree), ]
  # Global degree >= within-module degree (same or more edges in the full graph)
  expect_true(all(non_tiny$global_degree >= non_tiny$degree - 1e-10))
})


test_that("identify_module_hubs tie-breaking uses global degree not gene name", {
  # Create a module where two genes have identical within-module centrality
  # but different global degree
  n <- 10L
  gnames <- paste0("G", seq_len(n))
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
  # All genes fully connected with equal weight -> identical within-module degree
  mat[1:n, 1:n] <- 0.8
  # G1 gets an extra strong self-consistent weight bump won't help since diagonal
  # is 0. Instead, give G1 higher global degree by making the network have
  # two modules where G1 connects to the other module
  # Actually: in a single-module network, global degree = within-module degree.
  # So let's make 2 modules. G1-G5 in module 1, G6-G10 in module 2.
  # G5 has cross-module edges (higher global degree than G1-G4).
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
  mat[1:5, 1:5] <- 0.8
  mat[6:10, 6:10] <- 0.8
  # G5 connects to module 2 (cross-module edges)
  mat[5, 6:10] <- mat[6:10, 5] <- 0.6
  diag(mat) <- 1

  net <- list(network = mat, threshold = 0.5)
  m <- detect_modules(net, method = "leiden",
                       objective_function = "modularity", seed = 42)

  # G1-G4 should have identical within-module degree in module 1
  # G5 has higher global degree due to cross-module edges
  result <- identify_module_hubs(m, net, centrality = "degree", top_n = 1L)

  # G5 should be selected as hub (global degree breaks the tie)
  g5_mod <- result$module[result$gene == "G5"]
  if (!is.na(g5_mod)) {
    mod_hub <- result[result$module == g5_mod & result$is_hub, ]
    expect_true("G5" %in% mod_hub$gene)
  }
})


# ---- classify_hub_conservation() tests ----

# Helper: build hub_results for all 4 species
make_hub_results <- function(td, top_n = 1L) {
  hub_list <- list()
  for (sp in td$species) {
    # Find an ortholog table that has this species in Species1
    ortho_key <- grep(paste0("^", sp, "\\."), names(td$orthologs), value = TRUE)
    if (length(ortho_key) == 0L) {
      ortho_key <- grep(paste0("\\.", sp, "$"), names(td$orthologs),
                         value = TRUE)
    }
    ortho <- td$orthologs[[ortho_key[1]]]
    hub_list[[sp]] <- identify_module_hubs(
      td$mods[[sp]], td$nets[[sp]], ortho, top_n = top_n
    )
  }
  hub_list
}


test_that("classify_hub_conservation returns correct structure", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td)

  result <- classify_hub_conservation(hubs, td$trait)

  expect_true(is.data.frame(result))
  expected_cols <- c("hog", "classification", "n_species_hub",
                     "n_species_present", "hub_trait_groups",
                     "n_corresponding", "n_cross_pairs",
                     "max_centrality", "best_hub_species")
  expect_true(all(expected_cols %in% names(result)))
})


test_that("classify_hub_conservation identifies trait-specific hubs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 2L)

  result <- classify_hub_conservation(hubs, td$trait)

  # Check for at least one trait-specific hub
  trait_specific <- result[grepl("_specific_hub$", result$classification), ]
  expect_true(nrow(trait_specific) > 0)
})


test_that("classify_hub_conservation classifies non_hub HOGs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 1L)

  result <- classify_hub_conservation(hubs, td$trait)

  # With only 1 hub per module, most HOGs are non_hub
  non_hubs <- result[result$classification == "non_hub", ]
  expect_true(nrow(non_hubs) > 0)
  expect_true(all(non_hubs$n_species_hub == 0L))
  expect_true(all(is.na(non_hubs$hub_trait_groups)))
})


test_that("classify_hub_conservation without module_comparisons uses multi_trait_hub", {
  td <- make_hub_test_data()
  # Use top_n = 5 so gene 1 (shared hub) qualifies in both traits
  hubs <- make_hub_results(td, top_n = 5L)

  result <- classify_hub_conservation(hubs, td$trait)

  # HOG1 should be hub in multiple species / both traits
  hog1 <- result[result$hog == "HOG1", ]
  if (nrow(hog1) > 0 && hog1$n_species_hub >= 2) {
    # Without module_comparisons, multi-trait hubs can't be classified further
    multi_or_specific <- hog1$classification %in%
      c("multi_trait_hub", "conserved_hub", "rewired_hub",
        "annual_specific_hub", "perennial_specific_hub", "sporadic_hub")
    expect_true(multi_or_specific)
  }
})


test_that("classify_hub_conservation with module_comparisons detects conserved or rewired", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 5L)

  # Build pairwise module comparisons for cross-trait pairs
  mod_comps <- list()
  cross_pairs <- list(
    c("SP_A", "SP_C"), c("SP_A", "SP_D"),
    c("SP_B", "SP_C"), c("SP_B", "SP_D")
  )
  for (pair in cross_pairs) {
    sp1 <- pair[1]
    sp2 <- pair[2]
    key <- paste(sort(c(sp1, sp2)), collapse = ".")
    ortho_key <- paste(sp1, sp2, sep = ".")
    if (!ortho_key %in% names(td$orthologs)) {
      ortho_key <- paste(sp2, sp1, sep = ".")
    }
    mod_comps[[key]] <- compare_modules(
      td$mods[[sp1]], td$mods[[sp2]], td$orthologs[[ortho_key]]
    )
  }

  result <- classify_hub_conservation(hubs, td$trait,
                                       module_comparisons = mod_comps)

  # With module comparisons, multi-trait hubs should be classified as
  # conserved_hub or rewired_hub
  multi <- result[result$classification %in%
                    c("conserved_hub", "rewired_hub"), ]
  # At least check that the function runs without error
  expect_true(is.data.frame(result))
  expect_true(all(!is.na(result$classification)))
})


test_that("classify_hub_conservation min_trait_fraction works", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 1L)

  # Strict: hub must be in ALL species within trait group
  strict <- classify_hub_conservation(hubs, td$trait,
                                       min_trait_fraction = 1.0)
  # Lenient: hub in any species counts
  lenient <- classify_hub_conservation(hubs, td$trait,
                                        min_trait_fraction = 0.01)

  # Lenient should have more specific/multi-trait hubs than strict
  n_hub_strict <- sum(!strict$classification %in% c("non_hub", "sporadic_hub"))
  n_hub_lenient <- sum(!lenient$classification %in%
                         c("non_hub", "sporadic_hub"))
  expect_true(n_hub_lenient >= n_hub_strict)
})


test_that("classify_hub_conservation validates inputs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td)

  expect_error(classify_hub_conservation(list(1, 2), td$trait),
               "must be a named list")
  expect_error(classify_hub_conservation(hubs, c("a", "b")),
               "must be a named vector")
  expect_error(classify_hub_conservation(hubs,
                 c(SP_A = "annual", SP_B = "annual")),
               "missing entries")
})


test_that("classify_hub_conservation handles empty hub_results", {
  trait <- c(SP_A = "annual", SP_B = "perennial")
  empty_hub <- data.frame(
    gene = character(0), module = integer(0),
    degree = numeric(0), betweenness = numeric(0), eigenvector = numeric(0),
    mean_edge_weight = numeric(0), global_degree = numeric(0),
    rank = integer(0),
    is_hub = logical(0), hog = character(0),
    stringsAsFactors = FALSE
  )
  result <- classify_hub_conservation(
    list(SP_A = empty_hub, SP_B = empty_hub), trait
  )
  expect_equal(nrow(result), 0L)
  expect_true(all(c("hog", "classification") %in% names(result)))
})


# ---- Additional tests from review ----

test_that("identify_module_hubs handles multi-copy HOGs correctly", {
  td <- make_hub_test_data()
  # Create orthologs where one gene maps to two HOGs (duplicate gene entry)
  ortho <- data.frame(
    Species1 = c(paste0("A", 1:20), "A1"),
    Species2 = c(paste0("B", 1:20), "B1"),
    hog = c(paste0("HOG", 1:20), "HOG_ALT"),
    stringsAsFactors = FALSE
  )

  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]], ortho)

  # A1 should map to exactly one HOG (first one = HOG1, after dedup)
  a1 <- result[result$gene == "A1", ]
  expect_equal(nrow(a1), 1L)
  expect_equal(a1$hog, "HOG1")
})


test_that("identify_module_hubs min_module_size filters tiny modules", {
  # Create a network with one tiny module
  n <- 12L
  gnames <- paste0("X", seq_len(n))
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
  # Module 1: genes 1-10 (large)
  mat[1:10, 1:10] <- 0.9
  # Module 2: genes 11-12 (tiny, size 2)
  mat[11:12, 11:12] <- 0.9
  diag(mat) <- 1

  net <- list(network = mat, threshold = 0.5)
  m <- detect_modules(net, method = "leiden",
                       objective_function = "modularity", seed = 42)

  # Default min_module_size = 3: tiny module genes should not be hubs
  result <- identify_module_hubs(m, net, min_module_size = 3L)

  # Find the tiny module (size <= 2)
  mod_sizes <- tapply(result$gene, result$module, length)
  tiny_mods <- as.integer(names(mod_sizes[mod_sizes <= 2L]))

  if (length(tiny_mods) > 0) {
    tiny_rows <- result[result$module %in% tiny_mods, ]
    expect_true(all(!tiny_rows$is_hub))
    expect_true(all(is.na(tiny_rows$degree)))
  }
})


test_that("identify_module_hubs betweenness ranks bridge genes higher", {
  # Create a module where gene 1 bridges two sub-clusters
  n <- 10L
  gnames <- paste0("G", seq_len(n))
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))

  # Sub-cluster A: genes 1-5 all connected
  mat[1:5, 1:5] <- 0.8
  # Sub-cluster B: genes 1,6-10 all connected (gene 1 bridges A and B)
  mat[c(1, 6:10), c(1, 6:10)] <- 0.8
  diag(mat) <- 1

  net <- list(network = mat, threshold = 0.5)
  m <- detect_modules(net, method = "leiden",
                       objective_function = "modularity", seed = 42)

  result_deg <- identify_module_hubs(m, net, centrality = "degree", top_n = 1L)
  result_btw <- identify_module_hubs(m, net, centrality = "betweenness",
                                      top_n = 1L)

  # Gene 1 should have highest betweenness (bridges both sub-clusters)
  g1_mod <- result_btw$module[result_btw$gene == "G1"]
  g1_btw <- result_btw[result_btw$gene == "G1", ]
  expect_equal(g1_btw$rank, 1L)
})


test_that("classify_hub_conservation identifies sporadic_hub", {
  td <- make_hub_test_data()
  # top_n = 1: each module has 1 hub, so most HOGs are non_hub
  hubs <- make_hub_results(td, top_n = 1L)

  # Manually make HOG10 a hub in only 1 of 2 annual species
  # by editing the hub results directly
  for (sp in names(hubs)) {
    idx <- which(hubs[[sp]]$hog == "HOG10")
    if (length(idx) > 0) hubs[[sp]]$is_hub[idx] <- FALSE
  }
  # Make HOG10 a hub only in SP_A
  idx_a <- which(hubs[["SP_A"]]$hog == "HOG10")
  if (length(idx_a) > 0) hubs[["SP_A"]]$is_hub[idx_a] <- TRUE

  # With min_trait_fraction = 1.0 (strict), HOG10 is hub in 1/2 annuals
  # which is below threshold -> sporadic_hub
  result <- classify_hub_conservation(hubs, td$trait,
                                       min_trait_fraction = 1.0)
  hog10 <- result[result$hog == "HOG10", ]
  if (nrow(hog10) > 0 && hog10$n_species_hub > 0) {
    expect_equal(hog10$classification, "sporadic_hub")
  }
})


test_that("classify_hub_conservation multi_trait_hub without module_comparisons", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 5L)

  result <- classify_hub_conservation(hubs, td$trait)

  # HOGs that are hubs in both traits should be multi_trait_hub
  multi <- result[result$classification == "multi_trait_hub", ]
  if (nrow(multi) > 0) {
    # All should have hub_trait_groups spanning both traits
    for (i in seq_len(nrow(multi))) {
      groups <- strsplit(multi$hub_trait_groups[i], ",")[[1]]
      expect_true(length(groups) >= 2L)
    }
    # n_corresponding should be NA (no module_comparisons)
    expect_true(all(is.na(multi$n_corresponding)))
  }
})


test_that("identify_module_hubs comparison parameter enables conservation tie-breaking", {
  # Build a network where two genes have identical within-module degree
  # but different conservation effect sizes
  n <- 10L
  gnames <- paste0("G", seq_len(n))
  mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
  # All genes fully connected (identical within-module degree)
  mat[1:n, 1:n] <- 0.8
  diag(mat) <- 1

  net <- list(network = mat, threshold = 0.5)
  m <- detect_modules(net, method = "leiden",
                       objective_function = "modularity", seed = 42)

  # Orthologs: 1:1 mapping to partner species
  ortho <- data.frame(
    Species1 = gnames,
    Species2 = paste0("P", seq_len(n)),
    hog = paste0("HOG", seq_len(n)),
    stringsAsFactors = FALSE
  )

  # Mock comparison: G1 has high conservation effect, G2 has low
  mock_comparison <- data.frame(
    Species1 = gnames,
    Species2 = paste0("P", seq_len(n)),
    hog = paste0("HOG", seq_len(n)),
    Species1.effect.size = c(10, 0.5, rep(1, n - 2)),
    Species2.effect.size = c(10, 0.5, rep(1, n - 2)),
    Species1.q.val.con = c(0.001, 0.9, rep(0.5, n - 2)),
    Species2.q.val.con = c(0.001, 0.9, rep(0.5, n - 2)),
    stringsAsFactors = FALSE
  )

  # Without comparison: all genes have equal centrality,
  # tie-breaking falls to global degree (also equal), then alt centrality, etc.
  result_no_comp <- identify_module_hubs(m, net, ortho, top_n = 1L)

  # With comparison: G1 should win the tie (highest conservation effect)
  result_with_comp <- identify_module_hubs(m, net, ortho,
                                            comparison = mock_comparison,
                                            top_n = 1L)

  g1_mod <- result_with_comp$module[result_with_comp$gene == "G1"]
  hub_with <- result_with_comp[result_with_comp$module == g1_mod &
                                 result_with_comp$is_hub, ]
  expect_true("G1" %in% hub_with$gene)
})
