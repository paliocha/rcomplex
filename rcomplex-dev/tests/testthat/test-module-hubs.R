# Tests for identify_module_hubs() and classify_hub_conservation()

# Helper: 4-species, 2-trait scenario with designed hub patterns
# SP_A, SP_B = annual; SP_C, SP_D = perennial
# Each species has 20 genes, 2 modules of 10 genes each
make_hub_test_data <- function() {
  species <- c("SP_A", "SP_B", "SP_C", "SP_D")
  prefixes <- c(SP_A = "A", SP_B = "B", SP_C = "C", SP_D = "D")
  trait <- c(
    SP_A = "annual", SP_B = "annual",
    SP_C = "perennial", SP_D = "perennial"
  )
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
    mods[[sp]] <- detect_modules(net,
      objective_function = "modularity",
      seed = 42
    )
  }

  # Orthologs: 1:1 mapping across species
  # Build all pairwise orthologs
  ortho_list <- list()
  sp_pairs <- combn(species, 2, simplify = FALSE)
  for (pair in sp_pairs) {
    species1 <- pair[1]
    species2 <- pair[2]
    px1 <- prefixes[species1]
    px2 <- prefixes[species2]
    ortho_list[[paste(species1, species2, sep = ".")]] <- data.frame(
      gene1 = paste0(px1, seq_len(n)),
      gene2 = paste0(px2, seq_len(n)),
      hog = paste0("HOG", seq_len(n)),
      stringsAsFactors = FALSE
    )
  }

  list(
    nets = nets, mods = mods, orthologs = ortho_list, trait = trait,
    species = species, prefixes = prefixes
  )
}


# ---- identify_module_hubs() tests ----

test_that("identify_module_hubs returns correct structure", {
  td <- make_hub_test_data()
  sp <- "SP_A"
  ortho <- td$orthologs[["SP_A.SP_B"]] # any ortholog table for SP_A

  result <- identify_module_hubs(td$mods[[sp]], td$nets[[sp]], ortho)

  expect_true(is.data.frame(result))
  expect_true(all(c(
    "gene", "module", "degree", "betweenness", "eigenvector",
    "mean_edge_weight", "global_degree",
    "rank", "is_hub", "hog"
  ) %in% names(result)))
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


test_that(
  "identify_module_hubs rank 1 has highest centrality (default degree)",
  {
    td <- make_hub_test_data()
    result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])
    for (m in unique(result$module)) {
      mod_df <- result[result$module == m, ]
      if (all(is.na(mod_df$degree))) next
      top <- mod_df[mod_df$rank == 1L, ]
      expect_equal(top$degree, max(mod_df$degree))
    }
  }
)


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


test_that("identify_module_hubs validates inputs", {
  td <- make_hub_test_data()

  expect_error(
    identify_module_hubs(list(), td$nets[["SP_A"]]),
    "must be output from detect_modules"
  )
  expect_error(
    identify_module_hubs(td$mods[["SP_A"]], list()),
    "must be output from compute_network"
  )
})


test_that("identify_module_hubs gene 1 is top hub in module 1 (annual)", {
  td <- make_hub_test_data()
  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]])

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


test_that(
  "identify_module_hubs tie-breaking uses global degree not gene name",
  {
    # Create a module where two genes have identical within-module centrality
    # but different global degree
    n <- 10L
    gnames <- paste0("G", seq_len(n))
    mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
    # All genes fully connected with equal weight -> identical
    # within-module degree
    mat[1:n, 1:n] <- 0.8
    # G1 gets an extra strong self-consistent weight bump won't help since
    # diagonal is 0. Instead, give G1 higher global degree by making the
    # network have two modules where G1 connects to the other module
    # Actually: in a single-module network, global degree = within-module
    # degree.
    # So let's make 2 modules. G1-G5 in module 1, G6-G10 in module 2.
    # G5 has cross-module edges (higher global degree than G1-G4).
    mat <- matrix(0, n, n, dimnames = list(gnames, gnames))
    mat[1:5, 1:5] <- 0.8
    mat[6:10, 6:10] <- 0.8
    # G5 connects to module 2 (cross-module edges)
    mat[5, 6:10] <- mat[6:10, 5] <- 0.6
    diag(mat) <- 1

    net <- list(network = mat, threshold = 0.5)
    m <- detect_modules(net,
      objective_function = "modularity", seed = 42
    )

    # G1-G4 should have identical within-module degree in module 1
    # G5 has higher global degree due to cross-module edges
    result <- identify_module_hubs(m, net)

    # G5 should be selected as hub (global degree breaks the tie)
    g5_mod <- result$module[result$gene == "G5"]
    if (!is.na(g5_mod)) {
      mod_hub <- result[result$module == g5_mod & result$is_hub, ]
      expect_true("G5" %in% mod_hub$gene)
    }
  }
)


# ---- classify_hub_conservation() tests ----

# Helper: build hub_results for all 4 species
make_hub_results <- function(td, top_n = 1L) {
  hub_list <- list()
  for (sp in td$species) {
    # Find an ortholog table that has this species in gene1
    ortho_key <- grep(paste0("^", sp, "\\."), names(td$orthologs), value = TRUE)
    if (length(ortho_key) == 0L) {
      ortho_key <- grep(paste0("\\.", sp, "$"), names(td$orthologs),
        value = TRUE
      )
    }
    ortho <- td$orthologs[[ortho_key[1]]]
    hub_list[[sp]] <- identify_module_hubs(td$mods[[sp]], td$nets[[sp]], ortho)
    # designed hub pattern: the top_n ranked genes of each module
    hub_list[[sp]]$is_hub <- !is.na(hub_list[[sp]]$rank) &
      hub_list[[sp]]$rank <= top_n
  }
  hub_list
}


test_that("classify_hub_conservation returns correct structure", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td)

  result <- classify_hub_conservation(hubs, as_clades(td$trait))

  expect_true(is.data.frame(result))
  expected_cols <- c(
    "hog", "classification", "n_species_hub",
    "n_species_present", "hub_trait_groups",
    "n_corresponding", "n_cross_pairs",
    "max_centrality", "best_hub_species"
  )
  expect_true(all(expected_cols %in% names(result)))
})


test_that("classify_hub_conservation identifies trait-specific hubs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 2L)

  result <- classify_hub_conservation(hubs, as_clades(td$trait))

  # Check for at least one trait-specific hub
  trait_specific <- result[grepl("_specific_hub$", result$classification), ]
  expect_true(nrow(trait_specific) > 0)
})


test_that("classify_hub_conservation classifies non_hub HOGs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 1L)

  result <- classify_hub_conservation(hubs, as_clades(td$trait))

  # With only 1 hub per module, most HOGs are non_hub
  non_hubs <- result[result$classification == "non_hub", ]
  expect_true(nrow(non_hubs) > 0)
  expect_true(all(non_hubs$n_species_hub == 0L))
  expect_true(all(is.na(non_hubs$hub_trait_groups)))
})


test_that(
  "classify_hub_conservation without module_comparisons uses multi_trait_hub",
  {
    td <- make_hub_test_data()
    # Use top_n = 5 so gene 1 (shared hub) qualifies in both traits
    hubs <- make_hub_results(td, top_n = 5L)

    result <- classify_hub_conservation(hubs, as_clades(td$trait))

    # HOG1 should be hub in multiple species / both traits
    hog1 <- result[result$hog == "HOG1", ]
    if (nrow(hog1) > 0 && hog1$n_species_hub >= 2) {
      # Without module_comparisons, multi-trait hubs can't be classified further
      multi_or_specific <- hog1$classification %in%
        c(
          "multi_trait_hub", "conserved_hub", "rewired_hub",
          "annual_specific_hub", "perennial_specific_hub", "sporadic_hub"
        )
      expect_true(multi_or_specific)
    }
  }
)


test_that(
  paste(
    "classify_hub_conservation with module_comparisons detects conserved",
    "or rewired"
  ),
  {
    # module_correspondence() q-values draw from the global RNG.
    set.seed(42)
    td <- make_hub_test_data()
    hubs <- make_hub_results(td, top_n = 5L)

    # Build pairwise module comparisons for cross-trait pairs
    mod_comps <- list()
    cross_pairs <- list(
      c("SP_A", "SP_C"), c("SP_A", "SP_D"),
      c("SP_B", "SP_C"), c("SP_B", "SP_D")
    )
    for (pair in cross_pairs) {
      species1 <- pair[1]
      species2 <- pair[2]
      key <- paste(sort(c(species1, species2)), collapse = ".")
      ortho_key <- paste(species1, species2, sep = ".")
      if (!ortho_key %in% names(td$orthologs)) {
        ortho_key <- paste(species2, species1, sep = ".")
      }
      map <- resolve_ortholog_map(
        td$orthologs[[ortho_key]],
        rownames(td$nets[[species1]]$network),
        rownames(td$nets[[species2]]$network)
      )
      mod_comps[[key]] <- module_correspondence(
        td$mods[[species1]], td$mods[[species2]], map
      )
    }

    result <- classify_hub_conservation(hubs, as_clades(td$trait),
      module_comparisons = mod_comps
    )

    # With module comparisons, multi-trait hubs should be classified as
    # conserved_hub or rewired_hub
    multi <- result[result$classification %in%
                      c("conserved_hub", "rewired_hub"), ]
    # At least check that the function runs without error
    expect_true(is.data.frame(result))
    expect_true(all(!is.na(result$classification)))
    # A live correspondence must actually reach the lookup: n_corresponding is
    # reset to NA when the key misses, which is how the old unsorted keying
    # failed. n_cross_pairs is set regardless, so it proves nothing here.
    expect_false(all(is.na(result$n_corresponding)))
    expect_true(any(result$classification %in%
                      c("conserved_hub", "rewired_hub")))
  }
)


test_that("a species with no HOG-mapped genes is called out", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 1L)
  hubs[["SP_B"]]$hog <- NA_character_

  expect_warning(
    classify_hub_conservation(hubs, as_clades(td$trait)),
    "no hog-mapped gene"
  )
})


test_that("absent species count against the trait group", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td, top_n = 1L)

  # Hub in the only annual that carries the HOG: 1 of 2 annuals, which
  # clears the 0.5 cut.
  for (sp in names(hubs)) hubs[[sp]]$is_hub <- FALSE
  hubs[["SP_A"]]$is_hub[hubs[["SP_A"]]$hog == "HOG10"] <- TRUE
  hubs[["SP_B"]] <- hubs[["SP_B"]][hubs[["SP_B"]]$hog != "HOG10", ]

  res <- classify_hub_conservation(hubs, as_clades(td$trait))
  expect_equal(
    res$classification[res$hog == "HOG10"],
    "annual_specific_hub"
  )
})


test_that("classify_hub_conservation validates inputs", {
  td <- make_hub_test_data()
  hubs <- make_hub_results(td)

  expect_error(
    classify_hub_conservation(list(1, 2), as_clades(td$trait)),
    "must be a named list"
  )
  expect_error(
    classify_hub_conservation(hubs, c("a", "b")),
    "clades must be a named list"
  )
  expect_error(
    classify_hub_conservation(hubs, as_clades(td$trait),
      module_comparisons = list(SP_A.SP_C = list(raw = 1))
    ),
    "must be a module_correspondence\\(\\) result"
  )
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
    list(SP_A = empty_hub, SP_B = empty_hub), as_clades(trait)
  )
  expect_equal(nrow(result), 0L)
  expect_true(all(c("hog", "classification") %in% names(result)))
})


# ---- Additional tests from review ----

test_that("identify_module_hubs handles multi-copy HOGs correctly", {
  td <- make_hub_test_data()
  # Create orthologs where one gene maps to two HOGs (duplicate gene entry)
  ortho <- data.frame(
    gene1 = c(paste0("A", 1:20), "A1"),
    gene2 = c(paste0("B", 1:20), "B1"),
    hog = c(paste0("HOG", 1:20), "HOG_ALT"),
    stringsAsFactors = FALSE
  )

  result <- identify_module_hubs(td$mods[["SP_A"]], td$nets[["SP_A"]], ortho)

  # A1 should map to exactly one HOG (first one = HOG1, after dedup)
  a1 <- result[result$gene == "A1", ]
  expect_equal(nrow(a1), 1L)
  expect_equal(a1$hog, "HOG1")
})


test_that(
  "classify_hub_conservation multi_trait_hub without module_comparisons",
  {
    td <- make_hub_test_data()
    hubs <- make_hub_results(td, top_n = 5L)

    result <- classify_hub_conservation(hubs, as_clades(td$trait))

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
  }
)


test_that("orientation check survives species names containing a dot", {
  # Splitting the key on "." would make "A.thaliana.O.sativa" look transposed
  # and reject a correctly oriented table.
  trait <- c(A.thaliana = "annual", O.sativa = "perennial")
  corr <- list(pairs = data.frame(
    module1 = "1", module2 = "1", jaccard = 0.5, q_value = 0.01,
    stringsAsFactors = FALSE
  ), species_ref = "A.thaliana", species_test = "O.sativa")
  hubs <- list(
    A.thaliana = data.frame(
      gene = "a1", module = 1L, is_hub = TRUE,
      hog = "H1", degree = 1, stringsAsFactors = FALSE
    ),
    O.sativa = data.frame(
      gene = "o1", module = 1L, is_hub = TRUE,
      hog = "H1", degree = 1, stringsAsFactors = FALSE
    )
  )
  expect_no_error(
    classify_hub_conservation(hubs, as_clades(trait),
      module_comparisons = list("A.thaliana.O.sativa" = corr)
    )
  )
})
