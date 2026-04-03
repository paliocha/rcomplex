# --- Synthetic fixtures for tag_permutation ---

make_tag_perm_fixtures <- function() {
  # 3 pairs, 6 species, 2 traits (3 annual, 3 perennial)
  group <- c(A1 = "annual", P1 = "perennial",
             A2 = "annual", P2 = "perennial",
             A3 = "annual", P3 = "perennial")

  pairs <- data.frame(
    sp1 = c("A1", "A2", "A3"),
    sp2 = c("P1", "P2", "P3"),
    pair_name = c("pair1", "pair2", "pair3"),
    stringsAsFactors = FALSE
  )

  # Modules: each species gets 2 modules (1, 2)
  make_modules <- function(sp, genes) {
    membership <- stats::setNames(rep(c(1L, 2L), each = length(genes) / 2),
                                  genes)
    list(
      modules = membership,
      module_genes = split(names(membership), membership),
      n_modules = 2L,
      modularity = 0.3,
      graph = NULL,
      method = "leiden",
      params = list()
    )
  }

  modules <- list(
    A1 = make_modules("A1", paste0("A1_g", 1:6)),
    P1 = make_modules("P1", paste0("P1_g", 1:6)),
    A2 = make_modules("A2", paste0("A2_g", 1:6)),
    P2 = make_modules("P2", paste0("P2_g", 1:6)),
    A3 = make_modules("A3", paste0("A3_g", 1:6)),
    P3 = make_modules("P3", paste0("P3_g", 1:6))
  )

  # Orthologs: map genes to HOGs (consistent across species)
  orthologs <- data.frame(
    Species1 = c(paste0("A1_g", 1:6), paste0("A2_g", 1:6), paste0("A3_g", 1:6)),
    Species2 = c(paste0("P1_g", 1:6), paste0("P2_g", 1:6), paste0("P3_g", 1:6)),
    hog = rep(paste0("HOG", 1:6), 3),
    stringsAsFactors = FALSE
  )

  # Classification: all annual modules are species_specific (sp1 side)
  # Module 1 in annuals is species-specific, module 2 in perennials is species-specific
  classification <- data.frame(
    pair_name = rep(c("pair1", "pair2", "pair3"), each = 4),
    module = rep(c(1L, 2L, 1L, 2L), 3),
    species = rep(c("sp1", "sp1", "sp2", "sp2"), 3),
    classification = rep(c("species_specific", "conserved",
                           "conserved", "species_specific"), 3),
    stringsAsFactors = FALSE
  )

  list(classification = classification, modules = modules,
       orthologs = orthologs, pairs = pairs, group = group)
}


test_that("tag_permutation returns correct structure", {
  fix <- make_tag_perm_fixtures()

  result <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  )

  expect_type(result, "list")
  expect_true(all(c("observed", "null_distribution", "p_value",
                     "recurrence_table", "target_group",
                     "min_recurrence", "n_perm") %in% names(result)))
  expect_type(result$observed, "integer")
  expect_length(result$null_distribution, 100L)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_s3_class(result$recurrence_table, "data.frame")
  expect_equal(result$target_group, "annual")
  expect_equal(result$min_recurrence, 2L)
})


test_that("tag_permutation detects known parallel signal", {
  fix <- make_tag_perm_fixtures()

  # All 3 annual species have module 1 as species-specific.
  # Module 1 genes are A*_g1:g3, mapping to HOG1:HOG3.
  # So HOG1, HOG2, HOG3 each recur in 3 pairs.
  result <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  )

  expect_equal(result$observed, 3L)  # HOG1, HOG2, HOG3

  # Recurrence table should have 3 HOGs with n_pairs = 3
  recurring <- result$recurrence_table[result$recurrence_table$n_pairs >= 2, ]
  expect_equal(nrow(recurring), 3L)
})


test_that("tag_permutation handles min_recurrence thresholds", {
  fix <- make_tag_perm_fixtures()

  # min_recurrence = 3: all 3 pairs must have the HOG
  r3 <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 50L, min_recurrence = 3L
  )
  expect_equal(r3$observed, 3L)

  # min_recurrence = 4: impossible with 3 pairs
  r4 <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 50L, min_recurrence = 4L
  )
  expect_equal(r4$observed, 0L)
})


test_that("tag_permutation pair exclusion: both sides same trait", {
  fix <- make_tag_perm_fixtures()

  # If we make all species "annual", no pair has exactly one annual
  all_annual <- stats::setNames(rep("annual", 6),
                                names(fix$group))

  result <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    all_annual, target_group = "annual",
    n_perm = 50L, min_recurrence = 2L
  )

  expect_equal(result$observed, 0L)
})


test_that("tag_permutation works with > 2 trait values", {
  fix <- make_tag_perm_fixtures()

  # 3 trait values, unequal frequencies
  group3 <- c(A1 = "annual", P1 = "perennial",
              A2 = "annual", P2 = "biennial",
              A3 = "biennial", P3 = "perennial")

  result <- tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    group3, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  )

  expect_type(result, "list")
  # Only pair1 and pair2 have exactly one "annual" species
  # pair3 has biennial + perennial, so it doesn't contribute
  # HOG1:3 from pair1, HOG1:3 from pair2 -> 3 HOGs in 2 pairs
  expect_equal(result$observed, 3L)
})


test_that("tag_permutation validates inputs", {
  fix <- make_tag_perm_fixtures()

  expect_error(
    tag_permutation(fix$classification[, -1], fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "annual"),
    "classification missing columns"
  )

  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "nonexistent"),
    "target_group.*not found"
  )

  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, c("annual", "perennial"), "annual"),
    "named character"
  )

  bad_group <- fix$group[-1]
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, bad_group, "annual"),
    "group missing entries"
  )
})
