# Tests for get_coexpressed_hogs()

# Helper: build test data with 3 species, 5 HOGs, designed co-expression
# HOG1: candidate (1 gene per species)
# HOG2: co-expresses with HOG1 in all 3 species
# HOG3: co-expresses with HOG1 in SP_A and SP_B only
# HOG4: co-expresses with HOG1 in SP_A only
# HOG5: does NOT co-express with HOG1 in any species
# HOG6: extra HOG for multi-copy test (gene6a, gene6b in SP_A)
make_coexpr_test_data <- function() {
  genes_a <- paste0("A", 1:6)
  genes_b <- paste0("B", 1:6)
  genes_c <- paste0("C", 1:6)

  thr <- 0.5

  # Build symmetric 6x6 matrices per species
  # Diagonal = 1, off-diag = 0 by default, then set designed edges
  make_net <- function(genes, g1s, g2s, wts) {
    n <- length(genes)
    mat <- matrix(0, n, n, dimnames = list(genes, genes))
    diag(mat) <- 1
    for (i in seq_along(g1s)) {
      mat[g1s[i], g2s[i]] <- wts[i]
      mat[g2s[i], g1s[i]] <- wts[i]
    }
    list(network = mat, threshold = thr)
  }

  # SP_A: HOG1(A1) co-expresses with HOG2(A2), HOG3(A3), HOG4(A4)
  net_a <- make_net(genes_a,
                    c("A1", "A1", "A1"),
                    c("A2", "A3", "A4"),
                    c(0.8, 0.7, 0.6))

  # SP_B: HOG1(B1) co-expresses with HOG2(B2), HOG3(B3)
  net_b <- make_net(genes_b,
                    c("B1", "B1"),
                    c("B2", "B3"),
                    c(0.9, 0.75))

  # SP_C: HOG1(C1) co-expresses with HOG2(C2) only
  net_c <- make_net(genes_c, "C1", "C2", 0.85)

  networks <- list(SP_A = net_a, SP_B = net_b, SP_C = net_c)

  # Orthologs: all pairwise species, one gene per HOG per species
  orthologs <- data.frame(
    Species1 = c("A1", "A2", "A3", "A4", "A5", "A6",
                 "A1", "A2", "A3", "A4", "A5", "A6",
                 "B1", "B2", "B3", "B4", "B5", "B6"),
    Species2 = c("B1", "B2", "B3", "B4", "B5", "B6",
                 "C1", "C2", "C3", "C4", "C5", "C6",
                 "C1", "C2", "C3", "C4", "C5", "C6"),
    hog = c(paste0("HOG", 1:6),
            paste0("HOG", 1:6),
            paste0("HOG", 1:6)),
    stringsAsFactors = FALSE
  )

  trait <- c(SP_A = "annual", SP_B = "annual", SP_C = "perennial")

  # Edges for testing conservation join
  edges <- data.frame(
    gene1 = c("A2", "A2", "B2", "A3", "A5"),
    gene2 = c("B2", "C2", "C2", "B3", "B5"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_B"),
    hog = c("HOG2", "HOG2", "HOG2", "HOG3", "HOG5"),
    q.value = c(0.01, 0.02, 0.03, 0.04, 0.80),
    effect_size = c(3.0, 2.5, 2.0, 1.5, 0.5),
    type = c("conserved", "conserved", "conserved", "conserved", "ns"),
    stringsAsFactors = FALSE
  )

  list(networks = networks, orthologs = orthologs, trait = trait, edges = edges)
}


test_that("returns correct structure with coexpressed_species and coexpressed_traits", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  species_trait = d$trait)
  expect_true(is.data.frame(result))
  expected_cols <- c("partner_hog", "n_species", "coexpressed_species",
                     "coexpressed_traits", "mean_weight")
  expect_true(all(expected_cols %in% names(result)))
  expect_type(result$partner_hog, "character")
  expect_type(result$n_species, "integer")
  expect_type(result$coexpressed_species, "character")
  expect_type(result$coexpressed_traits, "character")
  expect_type(result$mean_weight, "double")
})


test_that("HOG2 appears with n_species = 3 and all species listed", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 1L)
  hog2 <- result[result$partner_hog == "HOG2", ]
  expect_equal(nrow(hog2), 1L)
  expect_equal(hog2$n_species, 3L)
  sp_list <- sort(strsplit(hog2$coexpressed_species, ",")[[1]])
  expect_equal(sp_list, c("SP_A", "SP_B", "SP_C"))
})


test_that("HOG4 excluded when min_species = 2", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 2L)
  expect_false("HOG4" %in% result$partner_hog)
  # But HOG2 (3 species) and HOG3 (2 species) should remain
  expect_true("HOG2" %in% result$partner_hog)
  expect_true("HOG3" %in% result$partner_hog)
})


test_that("HOG5 never appears (not co-expressed)", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 1L)
  expect_false("HOG5" %in% result$partner_hog)
})


test_that("species parameter filters correctly", {
  d <- make_coexpr_test_data()
  # Only look at SP_A and SP_B
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  species = c("SP_A", "SP_B"),
                                  min_species = 1L)
  # HOG2: 2 species, HOG3: 2, HOG4: 1
  hog2 <- result[result$partner_hog == "HOG2", ]
  expect_equal(hog2$n_species, 2L)
  # No SP_C in coexpressed_species
  expect_false(grepl("SP_C", hog2$coexpressed_species))
})


test_that("coexpressed_traits reflects trait groups correctly", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  species_trait = d$trait,
                                  min_species = 1L)

  # HOG2: all 3 species -> annual + perennial
  hog2 <- result[result$partner_hog == "HOG2", ]
  expect_equal(hog2$coexpressed_traits, "annual,perennial")

  # HOG3: SP_A + SP_B -> both annual
  hog3 <- result[result$partner_hog == "HOG3", ]
  expect_equal(hog3$coexpressed_traits, "annual")

  # HOG4: SP_A only -> annual
  hog4 <- result[result$partner_hog == "HOG4", ]
  expect_equal(hog4$coexpressed_traits, "annual")
})


test_that("species_trait = NULL produces NA coexpressed_traits", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 1L)
  expect_true(all(is.na(result$coexpressed_traits)))
})


test_that("edges parameter adds conservation columns", {
  d <- make_coexpr_test_data()
  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  edges = d$edges, min_species = 1L)
  expect_true("partner_conserved" %in% names(result))
  expect_true("partner_mean_effect" %in% names(result))
  expect_true("partner_min_q" %in% names(result))

  # HOG2 is conserved in edges
  hog2 <- result[result$partner_hog == "HOG2", ]
  expect_true(hog2$partner_conserved)
  expect_equal(hog2$partner_min_q, 0.01)
  expect_equal(hog2$partner_mean_effect, mean(c(3.0, 2.5, 2.0)))

  # HOG4 has no edges -> NA
  hog4 <- result[result$partner_hog == "HOG4", ]
  expect_false(hog4$partner_conserved)
  expect_true(is.na(hog4$partner_mean_effect))
})


test_that("multi-copy HOG: union of neighbors", {
  d <- make_coexpr_test_data()
  # Add a second gene for HOG1 in SP_A: A1b
  # A1b co-expresses with A6 (HOG6) but A1 does not
  net_a <- d$networks$SP_A
  genes_new <- c(rownames(net_a$network), "A1b")
  n <- length(genes_new)
  mat_new <- matrix(0, n, n, dimnames = list(genes_new, genes_new))
  mat_new[seq_len(n - 1L), seq_len(n - 1L)] <- net_a$network
  mat_new[n, n] <- 1
  # A1b co-expresses with A6 (HOG6)
  mat_new["A1b", "A6"] <- 0.7
  mat_new["A6", "A1b"] <- 0.7
  d$networks$SP_A <- list(network = mat_new, threshold = net_a$threshold)

  # Add A1b to orthologs as HOG1
  extra_ortho <- data.frame(
    Species1 = c("A1b", "A1b"),
    Species2 = c("B1", "C1"),
    hog = "HOG1",
    stringsAsFactors = FALSE
  )
  d$orthologs <- rbind(d$orthologs, extra_ortho)

  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 1L)
  # HOG6 should appear from the A1b neighbor
  expect_true("HOG6" %in% result$partner_hog)
  # HOG2 should still appear (from A1's neighbors)
  expect_true("HOG2" %in% result$partner_hog)
})


test_that("candidate HOG not in a species is skipped gracefully", {
  d <- make_coexpr_test_data()
  # Remove HOG1 genes from SP_C by removing C1 from orthologs
  d$orthologs <- d$orthologs[!(d$orthologs$Species1 == "A1" &
                                 d$orthologs$Species2 == "C1"), ]
  d$orthologs <- d$orthologs[!(d$orthologs$Species1 == "B1" &
                                 d$orthologs$Species2 == "C1"), ]

  result <- get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                  min_species = 1L)
  # Should still work, just fewer species for HOG2
  expect_true(is.data.frame(result))
  hog2 <- result[result$partner_hog == "HOG2", ]
  expect_equal(hog2$n_species, 2L)  # only SP_A and SP_B now
})


test_that("input validation errors", {
  d <- make_coexpr_test_data()

  expect_error(get_coexpressed_hogs(42, d$networks, d$orthologs),
               "single character string")
  expect_error(get_coexpressed_hogs(c("HOG1", "HOG2"), d$networks, d$orthologs),
               "single character string")
  expect_error(get_coexpressed_hogs("HOG1", list(d$networks$SP_A), d$orthologs),
               "named list")
  expect_error(get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                     species = "NOPE"),
               "not found in networks")
  expect_error(get_coexpressed_hogs("HOG1", d$networks, d$orthologs,
                                     species_trait = c(SP_A = "annual")),
               "missing entries")
})


test_that("empty result when no partners meet threshold", {
  d <- make_coexpr_test_data()
  # Use a HOG that has no co-expression partners (HOG5)
  result <- get_coexpressed_hogs("HOG5", d$networks, d$orthologs,
                                  min_species = 1L)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0L)
  expect_true("partner_hog" %in% names(result))
  expect_true("coexpressed_traits" %in% names(result))
})
