# --- Tests for find_cliques() ---

test_that("find_cliques simple 3-species clique", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 1)
  expect_equal(result$hog, "HOG1")
  expect_equal(result$SP_A, "A1")
  expect_equal(result$SP_B, "B1")
  expect_equal(result$SP_C, "C1")
  expect_equal(result$n_species, 3L)
  expect_equal(result$n_edges, 3L)
})


test_that("find_cliques FDR optimization picks best assignment", {
  # Two paralogs: A1 has lower FDR edges, A2 has higher
  # C++ Bron-Kerbosch finds one maximal species clique, then picks the

  # single best gene assignment (lowest mean FDR) -> returns 1 row (A1)
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A2", "A2"),
    gene2 = c("B1", "C1", "C1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
    hog = rep("HOG1", 5),
    type = rep("conserved", 5),
    q.value = c(0.01, 0.01, 0.01, 0.5, 0.5),
    effect_size = c(3.0, 3.0, 3.0, 1.0, 1.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  # C++ picks the single best assignment per species clique: A1 (lower FDR)
  expect_equal(nrow(result), 1)
  expect_equal(result$SP_A, "A1")
  expect_equal(result$mean_q, mean(c(0.01, 0.01, 0.01)), tolerance = 1e-10)
})


test_that("find_cliques handles multiple HOGs", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "X1", "X1", "Y1"),
    gene2 = c("B1", "C1", "C1", "Y1", "Z1", "Z1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C", "SP_C"),
    hog = c(rep("HOG1", 3), rep("HOG2", 3)),
    type = rep("conserved", 6),
    q.value = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
    effect_size = c(2, 3, 4, 5, 6, 7),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 2)
  expect_setequal(result$hog, c("HOG1", "HOG2"))
})


test_that("find_cliques returns correct FDR and effect size summary stats", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.04, 0.07),
    effect_size = c(2.0, 6.0, 4.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(result$mean_q, mean(c(0.01, 0.04, 0.07)), tolerance = 1e-10)
  expect_equal(result$max_q, 0.07)
  expect_equal(result$mean_effect_size, 4.0)
  expect_equal(result$n_edges, 3L)
})


test_that("find_cliques empty input returns correct empty structure", {
  edges <- data.frame(
    gene1 = character(0), gene2 = character(0),
    species1 = character(0), species2 = character(0),
    hog = character(0), q.value = numeric(0),
    effect_size = numeric(0), stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B"))

  expect_equal(nrow(result), 0)
  expect_true(all(c(
    "hog", "SP_A", "SP_B", "n_species",
    "mean_q"
  ) %in% names(result)))
})


test_that("find_cliques validates required columns", {
  edges <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)
  expect_error(
    find_cliques(edges, c("SP_A", "SP_B")),
    "missing required columns"
  )
})


test_that("find_cliques edge_type filtering works", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = c("conserved", "conserved", "diverged"),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 0.3),
    stringsAsFactors = FALSE
  )

  # conserved only: missing B1-C1 edge -> no clique
  result_con <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  expect_equal(nrow(result_con), 0)

  # Both types: complete clique
  result_both <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    edge_type = c("conserved", "diverged")
  )
  expect_equal(nrow(result_both), 1)
})


test_that("find_cliques min_species allows partial cliques", {
  # 3 target species but allow 2-species cliques
  edges <- data.frame(
    gene1 = "A1",
    gene2 = "B1",
    species1 = "SP_A",
    species2 = "SP_B",
    hog = "HOG1",
    type = "conserved",
    q.value = 0.01,
    effect_size = 3.0,
    stringsAsFactors = FALSE
  )

  # min_species = 3 (default): no clique (only 2 species present)
  result_strict <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    min_species = 3L
  )
  expect_equal(nrow(result_strict), 0)

  # min_species = 2: 2-species clique found
  result_lenient <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    min_species = 2L
  )
  expect_equal(nrow(result_lenient), 1)
  expect_equal(result_lenient$n_species, 2L)
})


test_that("find_cliques max_genes_per_sp limits combinatorial explosion", {
  # Many paralogs in SP_A, should be capped
  genes_a <- paste0("A", 1:15)
  edges_list <- lapply(genes_a, function(a) {
    data.frame(
      gene1 = a, gene2 = "B1",
      species1 = "SP_A", species2 = "SP_B",
      hog = "HOG1", type = "conserved",
      q.value = 0.01, effect_size = 2.0,
      stringsAsFactors = FALSE
    )
  })
  edges <- do.call(rbind, edges_list)

  # With default max_genes_per_sp = 10, should cap at 10 cliques
  result <- find_cliques(edges, c("SP_A", "SP_B"), max_genes_per_sp = 10L)
  expect_true(nrow(result) <= 10)
})


test_that("find_cliques fewer than 2 target species raises error", {
  edges <- data.frame(
    gene1 = "A1", gene2 = "B1",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", type = "conserved",
    q.value = 0.01, effect_size = 2.0,
    stringsAsFactors = FALSE
  )
  expect_error(find_cliques(edges, "SP_A"))
})


# --- Tests for edge cases ---

test_that("find_cliques handles duplicate edges (keeps min q-value)", {
  # Same gene pair twice with different q-values
  edges <- data.frame(
    gene1 = c("A1", "A1"),
    gene2 = c("B1", "B1"),
    species1 = c("SP_A", "SP_A"),
    species2 = c("SP_B", "SP_B"),
    hog = rep("HOG1", 2),
    type = rep("conserved", 2),
    q.value = c(0.5, 0.01), # second is better
    effect_size = c(1.0, 3.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B"))
  expect_equal(nrow(result), 1)
  expect_equal(result$mean_q, 0.01, tolerance = 1e-10)
  expect_equal(result$mean_effect_size, 3.0, tolerance = 1e-10)
})


test_that("find_cliques handles many paralogs without error", {
  # 6 paralogs per species in one HOG — exercises deep backtracking
  n_par <- 6
  pairs <- expand.grid(
    a = paste0("A", seq_len(n_par)),
    b = paste0("B", seq_len(n_par)),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    gene1 = pairs$a, gene2 = pairs$b,
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", type = "conserved",
    q.value = seq(0.001, 0.05, length.out = nrow(pairs)),
    effect_size = seq(1.5, 4.0, length.out = nrow(pairs)),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B"))
  expect_equal(nrow(result), 1)
  expect_equal(result$n_species, 2L)
  # Should pick the pair with the lowest q-value
  expect_true(result$mean_q <= max(edges$q.value))
})


test_that("max_missing_edges with paralogs picks best assignment", {
  # 3 species, 2 paralogs in SP_A, missing A-C edges
  # A1 connects to B1, A2 connects to B1
  # B1 connects to C1
  # No A-C edges -> need max_missing = 1 for 3-species clique
  edges <- data.frame(
    gene1 = c("A1", "A2", "B1"),
    gene2 = c("B1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_B", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.05, 0.02),
    effect_size = c(3.0, 1.5, 2.5),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    min_species = 3L, max_missing_edges = 1L
  )

  expect_equal(nrow(result), 1)
  # Should pick A1 (lower q with B1) over A2
  expect_equal(result$SP_A, "A1")
  expect_equal(result$n_missing, 1L)
})


# --- Tests for max_missing_edges ---

test_that("max_missing_edges=0 requires all edges (default)", {
  # 3-species HOG, missing SP_A-SP_C edge
  edges <- data.frame(
    gene1 = c("A1", "B1"),
    gene2 = c("B1", "C1"),
    species1 = c("SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C"),
    hog = rep("HOG1", 2),
    type = rep("conserved", 2),
    q.value = c(0.01, 0.02),
    effect_size = c(2.0, 3.0),
    stringsAsFactors = FALSE
  )

  # With default max_missing_edges=0, no 3-species clique (A-C missing)
  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    min_species = 3L
  )
  expect_equal(nrow(result), 0)
})


test_that("max_missing_edges=1 finds clique with one missing edge", {
  # 3-species HOG, missing SP_A-SP_C edge
  edges <- data.frame(
    gene1 = c("A1", "B1"),
    gene2 = c("B1", "C1"),
    species1 = c("SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C"),
    hog = rep("HOG1", 2),
    type = rep("conserved", 2),
    q.value = c(0.01, 0.02),
    effect_size = c(2.0, 3.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    min_species = 3L, max_missing_edges = 1L
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$n_species, 3L)
  expect_equal(result$n_edges, 2L) # 2 of 3 present
  expect_equal(result$n_missing, 1L) # 1 missing
  expect_equal(result$SP_A, "A1")
  expect_equal(result$SP_B, "B1")
  expect_equal(result$SP_C, "C1")
})


test_that("max_missing_edges prefers fewer missing edges", {
  # 4-species HOG: complete 3-species clique (A,B,C) + partial 4-species
  # with D connected to A only (missing B-D and C-D)
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A1"),
    gene2 = c("B1", "C1", "C1", "D1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_D"),
    hog = rep("HOG1", 4),
    type = rep("conserved", 4),
    q.value = c(0.01, 0.01, 0.01, 0.01),
    effect_size = rep(2.0, 4),
    stringsAsFactors = FALSE
  )

  # With max_missing=1, should find 4-species subset only if it has <= 1 missing
  # {A,B,C,D} has 6 possible edges, 4 present, 2 missing -> over budget
  # So only 3-species complete cliques survive
  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C", "SP_D"),
    min_species = 3L, max_missing_edges = 1L
  )

  # Should find the complete {A,B,C} clique with 0 missing
  complete <- result[result$n_missing == 0, ]
  expect_true(nrow(complete) >= 1)
  expect_equal(complete$n_species[1], 3L)
  expect_equal(complete$n_edges[1], 3L)
})


test_that("max_missing_edges n_missing output is 0 when all edges present", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  # Default mode: n_missing should always be 0
  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  expect_equal(result$n_missing, 0L)

  # With max_missing_edges=1, complete clique still has 0 missing
  result2 <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    max_missing_edges = 1L
  )
  # The complete {A,B,C} subset should be found with 0 missing
  full <- result2[result2$n_species == 3L & result2$n_missing == 0L, ]
  expect_true(nrow(full) >= 1)
})


# --- Tests for intensity and min_effect_size ---
#
# Onnela weights are each edge's Jaccard-index percentile among every
# tested pair of its species pair. 1 - q.value pinned every weight above
# 1 - alpha, because cliques only contain edges that already passed alpha
# (#11); effect size would have ranked genes by inverse degree, because
# fold enrichment falls like 1 / degree at a fixed conserved fraction (#12).

# One conserved triangle (HOG1: A1, B1, C1) plus background rows of type
# "ns" in other HOGs. find_cliques() never builds cliques from the ns rows,
# but they belong to the tested population the percentiles are taken over.
make_percentile_edges <- function(jac_ab = 0.2, jac_ac = 0.6, jac_bc = 0.4,
                                  background = TRUE) {
  tri <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "HOG1", type = "conserved",
    q.value = 0.01,
    effect_size = c(2, 6, 4),
    jaccard = c(jac_ab, jac_ac, jac_bc),
    stringsAsFactors = FALSE
  )
  if (!background) {
    return(tri)
  }
  bg <- data.frame(
    gene1 = c("A2", "A3", "A4", "A2", "A3", "A4", "B2", "B3", "B4"),
    gene2 = c("B2", "B3", "B4", "C2", "C3", "C4", "C2", "C3", "C4"),
    species1 = rep(c("SP_A", "SP_A", "SP_B"), each = 3),
    species2 = rep(c("SP_B", "SP_C", "SP_C"), each = 3),
    hog = paste0("HOG_bg", 1:9), type = "ns",
    q.value = 0.5,
    effect_size = 1,
    jaccard = c(0.1, 0.3, 0.5, 0.1, 0.2, 0.7, 0.8, 0.9, 0.95),
    stringsAsFactors = FALSE
  )
  rbind(tri, bg)
}


test_that("intensity weights are within-pair Jaccard percentiles", {
  edges <- make_percentile_edges()
  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 1L)
  # Intensity is the geometric mean of the edges' ensemble probabilities.
  w <- rcomplex:::.onnela_weight(edges)[1:3]
  expect_equal(result$intensity, exp(mean(log(w))), tolerance = 1e-12)
  expect_equal(result$min_effect_size, 2.0, tolerance = 1e-12)
  # Coherence (geometric / arithmetic mean) was retired: on real data it
  # had no range within a clique size and scaled with sqrt(edges) across
  # sizes (#20). It must not come back as a column.
  expect_false("coherence" %in% names(result))
})


test_that("cliques with equal q-values but different overlap separate", {
  # Two conserved triangles with identical q = 0.01 on every edge. Under
  # 1 - q both had intensity 0.99; association strength tells them apart.
  strong <- make_percentile_edges(background = FALSE)
  strong$effect_size <- c(8, 9, 10)
  weak <- make_percentile_edges(background = FALSE)
  weak$effect_size <- c(1.1, 1.2, 1.3)
  weak$hog <- "HOG2"
  weak$gene1 <- sub("1$", "9", weak$gene1)
  weak$gene2 <- sub("1$", "9", weak$gene2)
  result <- find_cliques(rbind(strong, weak), c("SP_A", "SP_B", "SP_C"))

  expect_equal(nrow(result), 2L)
  expect_gt(
    result$intensity[result$hog == "HOG1"],
    result$intensity[result$hog == "HOG2"]
  )
})


test_that("the weight ignores the Jaccard index", {
  # Two genes with the same fold enrichment score alike however their
  # Jaccard indices differ: E[Jaccard] grows with neighbourhood size, so a
  # Jaccard weight ranked hub genes above equally conserved small ones
  # (#15). Association strength is observed over expected, so it does not.
  edges <- data.frame(
    gene1 = c("A1", "A2"), gene2 = c("B1", "B2"),
    species1 = "SP_A", species2 = "SP_B",
    hog = c("HOG1", "HOG2"), type = "conserved",
    q.value = 0.01, effect_size = 5, jaccard = c(0.9, 0.1),
    stringsAsFactors = FALSE
  )
  result <- find_cliques(edges, c("SP_A", "SP_B"))

  expect_equal(nrow(result), 2L)
  expect_equal(result$intensity[1], result$intensity[2], tolerance = 1e-12)
})


test_that("the weight scale is fitted within each species pair", {
  # Rescaling one pair's association strengths is absorbed by the fitted
  # scale (z -> z/lambda), so intensity is comparable across pairs whose
  # enrichments live on different scales (Garlaschelli et al. 2013, eq. 11).
  # The background rows matter here: with one edge per pair the scale is
  # degenerate and every weight collapses to 0.5, which would make this
  # assertion vacuous.
  base <- make_percentile_edges()
  scaled <- base
  ac <- scaled$species1 == "SP_A" & scaled$species2 == "SP_C"
  scaled$effect_size[ac] <- scaled$effect_size[ac] * 100
  sp <- c("SP_A", "SP_B", "SP_C")

  out <- find_cliques(base, sp)$intensity
  expect_false(is.na(out))
  expect_false(isTRUE(all.equal(out, 0.5)))
  expect_equal(find_cliques(scaled, sp)$intensity, out, tolerance = 1e-6)
})


test_that("rows removed by edge_type still shape the weight scale", {
  sp <- c("SP_A", "SP_B", "SP_C")
  with_bg <- find_cliques(make_percentile_edges(), sp)
  without_bg <- find_cliques(make_percentile_edges(background = FALSE), sp)

  # clique membership is decided by conserved edges only ...
  expect_identical(with_bg[, c("hog", sp)], without_bg[, c("hog", sp)])
  # ... but the scale is fitted on every tested row, so dropping the ns
  # background moves the weights and hence the intensity.
  expect_false(isTRUE(all.equal(with_bg$intensity, without_bg$intensity)))
})


test_that("a clique with an edge lacking effect_size gets NA intensity", {
  # Dropping the incomplete edge would score the clique on a smaller edge
  # set, a different quantity, so the clique reports NA instead.
  sp <- c("SP_A", "SP_B", "SP_C")
  one_na <- make_percentile_edges()
  one_na$effect_size[2] <- NA_real_
  res <- find_cliques(one_na, sp)
  expect_equal(nrow(res), 1L)
  expect_true(is.na(res$intensity))

  two_na <- make_percentile_edges()
  two_na$effect_size[2:3] <- NA_real_
  expect_true(is.na(find_cliques(two_na, sp)$intensity))
})


test_that("effect_size is required, and an all-NA column gives NA", {
  sp <- c("SP_A", "SP_B", "SP_C")

  # Unlike jaccard, effect_size is a required column: find_cliques()
  # refuses the table outright rather than reaching the weight.
  gone <- make_percentile_edges()
  gone$effect_size <- NULL
  expect_error(find_cliques(gone, sp), "missing required columns")

  # Present but unusable: the weight is NA, so the clique scores NA
  # rather than being silently dropped from the intensity.
  rlang::local_options(rlib_warning_verbosity = "verbose")
  empty <- make_percentile_edges()
  empty$effect_size <- NA_real_
  expect_warning(result <- find_cliques(empty, sp),
                 class = "rcomplex_missing_effect_size")
  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$intensity))
})


test_that("a pre-filtered edge table warns", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  sp <- c("SP_A", "SP_B", "SP_C")
  expect_warning(
    find_cliques(make_percentile_edges(background = FALSE), sp),
    class = "rcomplex_prefiltered_edges"
  )
  # rows outside edge_type present: the table is the tested population
  expect_no_warning(
    find_cliques(make_percentile_edges(), sp),
    class = "rcomplex_prefiltered_edges"
  )
  # no type column: nothing to judge
  no_type <- make_percentile_edges(background = FALSE)
  no_type$type <- NULL
  expect_no_warning(
    find_cliques(no_type, sp),
    class = "rcomplex_prefiltered_edges"
  )
})


test_that("clique_stability on an unfiltered table does not warn", {
  # clique_stability() filters to edge_type for its own engine; handing
  # that filtered copy to find_cliques() used to trip the pre-filter
  # warning on a caller who passed the recommended unfiltered table.
  rlang::local_options(rlib_warning_verbosity = "verbose")
  fx <- make_clique_fixture_3sp()
  expect_true(any(fx$edges$type != "conserved"))
  expect_no_warning(
    clique_stability(fx$edges, fx$target_species, min_species = 2L),
    class = "rcomplex_prefiltered_edges"
  )
})


test_that("a single-edge clique's intensity is its edge weight", {
  edges <- data.frame(
    gene1 = c("A1", "A2"), gene2 = c("B1", "B2"),
    species1 = "SP_A", species2 = "SP_B",
    hog = c("HOG1", "HOG_bg"), type = c("conserved", "ns"),
    q.value = c(0.01, 0.5), effect_size = c(2, 4), jaccard = c(0.2, 0.4),
    stringsAsFactors = FALSE
  )
  result <- find_cliques(edges, c("SP_A", "SP_B"))

  expect_equal(nrow(result), 1L)
  # The clique's one edge has effect_size 2 among {2, 4}; intensity is its
  # ensemble probability p = z*w / (1 + z*w), not a rank.
  expect_equal(
    result$intensity,
    rcomplex:::.onnela_weight(edges)[1],
    tolerance = 1e-12
  )
})


test_that("min_effect_size returns minimum of effect sizes", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(5.0, 1.5, 3.0),
    stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))

  expect_equal(result$min_effect_size, 1.5, tolerance = 1e-10)
})


test_that("empty cliques have intensity/min_effect_size columns", {
  edges <- data.frame(
    gene1 = character(0), gene2 = character(0),
    species1 = character(0), species2 = character(0),
    hog = character(0), q.value = numeric(0),
    effect_size = numeric(0), stringsAsFactors = FALSE
  )

  result <- find_cliques(edges, c("SP_A", "SP_B"))

  expect_true("intensity" %in% names(result))
  expect_true("min_effect_size" %in% names(result))
  expect_equal(nrow(result), 0)
  expect_true(is.numeric(result$intensity))
  expect_true(is.numeric(result$min_effect_size))
})


# --- Tests for clique_persistence() ---

test_that("clique_persistence uses co-expressologs not full neighbourhood", {
  # A1 has neighbours A2 (strong) and A3 (marginal)
  # B1 has neighbours B2 (strong) and B3 (marginal)
  # Only A2<->B2 are orthologs, so persistence comes from that pair alone
  net_a <- matrix(c(0, 10, 2.1, 10, 0, 0.5, 2.1, 0.5, 0),
    nrow = 3,
    dimnames = list(c("A1", "A2", "A3"), c("A1", "A2", "A3"))
  )
  net_b <- matrix(c(0, 8, 2.6, 8, 0, 0.5, 2.6, 0.5, 0),
    nrow = 3,
    dimnames = list(c("B1", "B2", "B3"), c("B1", "B2", "B3"))
  )

  networks <- list(
    SP_A = list(network = net_a, threshold = 2.0),
    SP_B = list(network = net_b, threshold = 2.5)
  )

  cliques <- data.frame(
    hog = "HOG1", SP_A = "A1", SP_B = "B1",
    n_species = 2L, mean_q = 0.01, max_q = 0.01,
    mean_effect_size = 3.0, n_edges = 1L,
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    gene1 = c("A1", "A2"), gene2 = c("B1", "B2"),
    species1 = c("SP_A", "SP_A"), species2 = c("SP_B", "SP_B"),
    hog = c("HOG1", "HOG_X"),
    q.value = c(0.01, 0.5), effect_size = c(3.0, 1.0),
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(cliques, c("SP_A", "SP_B"), networks, edges)

  # Only co-expressolog: (A2, B2). ratio = min(10/2, 8/2.5) = 3.2
  # Neighbourhood-wide would give 1.04 (dominated by marginal A3, B3)
  expect_equal(result$persistence, 3.2, tolerance = 1e-10)
  expect_equal(result$mean_persistence, 3.2, tolerance = 1e-10)
})


test_that("clique_persistence weakest co-expressolog determines score", {
  net_a <- matrix(c(0, 10, 3, 10, 0, 0.5, 3, 0.5, 0),
    nrow = 3,
    dimnames = list(c("A1", "A2", "A3"), c("A1", "A2", "A3"))
  )
  net_b <- matrix(c(0, 8, 2.6, 8, 0, 0.5, 2.6, 0.5, 0),
    nrow = 3,
    dimnames = list(c("B1", "B2", "B3"), c("B1", "B2", "B3"))
  )

  networks <- list(
    SP_A = list(network = net_a, threshold = 2.0),
    SP_B = list(network = net_b, threshold = 2.5)
  )

  cliques <- data.frame(
    hog = "HOG1", SP_A = "A1", SP_B = "B1",
    n_species = 2L, mean_q = 0.01, max_q = 0.01,
    mean_effect_size = 3.0, n_edges = 1L,
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    gene1 = c("A1", "A2", "A3"),
    gene2 = c("B1", "B2", "B3"),
    species1 = rep("SP_A", 3), species2 = rep("SP_B", 3),
    hog = c("HOG1", "HOG_X", "HOG_Y"),
    q.value = rep(0.01, 3), effect_size = rep(3.0, 3),
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(cliques, c("SP_A", "SP_B"), networks, edges)

  # (A2,B2): min(10/2, 8/2.5) = min(5, 3.2) = 3.2
  # (A3,B3): min(3/2, 2.6/2.5) = min(1.5, 1.04) = 1.04
  expect_equal(result$persistence, 1.04, tolerance = 1e-10)
  expect_equal(result$mean_persistence, mean(c(3.2, 1.04)), tolerance = 1e-10)
})


test_that("clique_persistence returns NA with no co-expressologs", {
  net_a <- matrix(c(0, 5, 5, 0),
    nrow = 2,
    dimnames = list(c("A1", "A2"), c("A1", "A2"))
  )
  net_b <- matrix(c(0, 3, 3, 0),
    nrow = 2,
    dimnames = list(c("B1", "B2"), c("B1", "B2"))
  )

  networks <- list(
    SP_A = list(network = net_a, threshold = 2.0),
    SP_B = list(network = net_b, threshold = 1.0)
  )

  cliques <- data.frame(
    hog = "HOG1", SP_A = "A1", SP_B = "B1",
    n_species = 2L, mean_q = 0.01, max_q = 0.01,
    mean_effect_size = 2.0, n_edges = 1L,
    stringsAsFactors = FALSE
  )

  # A2's ortholog B3 is not in SP_B's network -> no co-expressologs
  edges <- data.frame(
    gene1 = c("A1", "A2"), gene2 = c("B1", "B3"),
    species1 = c("SP_A", "SP_A"), species2 = c("SP_B", "SP_B"),
    hog = c("HOG1", "HOG_X"),
    q.value = c(0.01, 0.5), effect_size = c(3.0, 1.0),
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(cliques, c("SP_A", "SP_B"), networks, edges)

  expect_true(is.na(result$persistence))
  expect_true(is.na(result$mean_persistence))
})


test_that("clique_persistence handles multiple cliques", {
  net_a <- matrix(c(0, 10, 3, 10, 0, 0.5, 3, 0.5, 0),
    nrow = 3,
    dimnames = list(c("A1", "A2", "A3"), c("A1", "A2", "A3"))
  )
  net_b <- matrix(c(0, 8, 4, 8, 0, 0.5, 4, 0.5, 0),
    nrow = 3,
    dimnames = list(c("B1", "B2", "B3"), c("B1", "B2", "B3"))
  )

  networks <- list(
    SP_A = list(network = net_a, threshold = 2.0),
    SP_B = list(network = net_b, threshold = 2.5)
  )

  cliques <- data.frame(
    hog = c("HOG1", "HOG2"),
    SP_A = c("A1", "A2"), SP_B = c("B1", "B2"),
    n_species = c(2L, 2L), mean_q = c(0.01, 0.02),
    max_q = c(0.01, 0.02), mean_effect_size = c(2.0, 3.0),
    n_edges = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    gene1 = c("A1", "A2", "A3"),
    gene2 = c("B1", "B2", "B3"),
    species1 = rep("SP_A", 3), species2 = rep("SP_B", 3),
    hog = c("HOG1", "HOG2", "HOG3"),
    q.value = rep(0.01, 3), effect_size = rep(3.0, 3),
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(cliques, c("SP_A", "SP_B"), networks, edges)

  # HOG1 (A1, B1): co-expressologs (A2,B2) min(10/2,8/2.5)=3.2,
  #   (A3,B3) min(3/2,4/2.5)=min(1.5,1.6)=1.5 -> persistence=1.5
  # HOG2 (A2, B2): co-expressolog (A1,B1) min(10/2,8/2.5)=3.2
  #   (A3 not neighbour of A2) -> persistence=3.2
  expect_equal(result$persistence[1], 1.5, tolerance = 1e-10)
  expect_equal(result$persistence[2], 3.2, tolerance = 1e-10)
})


test_that("clique_persistence aggregates across 3 species and reversed edges", {
  net_a <- matrix(c(0, 10, 3, 10, 0, 0.5, 3, 0.5, 0),
    nrow = 3,
    dimnames = list(c("A1", "A2", "A3"), c("A1", "A2", "A3"))
  )
  net_b <- matrix(c(0, 8, 4, 8, 0, 0.5, 4, 0.5, 0),
    nrow = 3,
    dimnames = list(c("B1", "B2", "B3"), c("B1", "B2", "B3"))
  )
  net_c <- matrix(c(0, 6, 5, 6, 0, 0.5, 5, 0.5, 0),
    nrow = 3,
    dimnames = list(c("C1", "C2", "C3"), c("C1", "C2", "C3"))
  )

  networks <- list(
    SP_A = list(network = net_a, threshold = 2.0),
    SP_B = list(network = net_b, threshold = 2.5),
    SP_C = list(network = net_c, threshold = 3.0)
  )

  cliques <- data.frame(
    hog = "HOG1", SP_A = "A1", SP_B = "B1", SP_C = "C1",
    n_species = 3L, mean_q = 0.01, max_q = 0.01,
    mean_effect_size = 3.0, n_edges = 3L,
    stringsAsFactors = FALSE
  )

  # SP_C->SP_A edge is reversed (species1=SP_C, species2=SP_A)
  edges <- data.frame(
    gene1 = c("A2", "C2", "B3"),
    gene2 = c("B2", "A2", "C3"),
    species1 = c("SP_A", "SP_C", "SP_B"),
    species2 = c("SP_B", "SP_A", "SP_C"),
    hog = c("HOG_X", "HOG_X", "HOG_Y"),
    q.value = c(0.5, 0.5, 0.5), effect_size = c(1.0, 1.0, 1.0),
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(
    cliques, c("SP_A", "SP_B", "SP_C"),
    networks, edges
  )

  # (SP_A,SP_B): A2->B2, min(10/2, 8/2.5) = 3.2
  # (SP_A,SP_C): A2->C2 via reversed edge, min(10/2, 6/3) = 2.0
  # (SP_B,SP_C): B3->C3, min(4/2.5, 5/3) = 1.6
  expect_equal(result$persistence, 1.6, tolerance = 1e-10)
  expect_equal(result$mean_persistence, mean(c(3.2, 2.0, 1.6)),
    tolerance = 1e-10
  )
})


test_that("clique_persistence excludes self from neighbours", {
  net <- matrix(c(99, 1, 1, 99),
    nrow = 2,
    dimnames = list(c("A1", "A2"), c("A1", "A2"))
  )

  networks <- list(
    SP_A = list(network = net, threshold = 2.0),
    SP_B = list(network = net, threshold = 2.0)
  )

  cliques <- data.frame(
    hog = "HOG1", SP_A = "A1", SP_B = "A2",
    n_species = 2L, mean_q = 0.01, max_q = 0.01,
    mean_effect_size = 2.0, n_edges = 1L,
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    gene1 = "A1", gene2 = "A2",
    species1 = "SP_A", species2 = "SP_B",
    hog = "HOG1", q.value = 0.01, effect_size = 3.0,
    stringsAsFactors = FALSE
  )

  result <- clique_persistence(cliques, c("SP_A", "SP_B"), networks, edges)

  # With self-exclusion: A1's only neighbour is A2 (MR=1 < 2.0) -> no neighbours
  expect_true(is.na(result$persistence))
})


test_that("clique_persistence validates inputs", {
  dummy_edges <- data.frame(
    gene1 = "A", gene2 = "B", species1 = "SP_A", species2 = "SP_B",
    hog = "H", q.value = 0.01, effect_size = 1.0,
    stringsAsFactors = FALSE
  )

  expect_error(
    clique_persistence(data.frame(x = 1), c("A", "B"), list(), dummy_edges),
    "cliques must be a data frame from find_cliques"
  )
  expect_error(
    clique_persistence(
      data.frame(hog = "H", A = "g", B = "g"),
      c("A", "B"), list(), dummy_edges
    ),
    "networks must be a named list"
  )
  expect_error(
    clique_persistence(
      data.frame(hog = "H", A = "g"),
      c("A"), list(A = list()), dummy_edges
    ),
    "target_species must have at least 2 species"
  )
  expect_error(
    clique_persistence(
      data.frame(hog = "H", A = "g"),
      c("A", "B"), list(A = list(), B = list()), dummy_edges
    ),
    "cliques missing columns for species"
  )
  expect_error(
    clique_persistence(
      data.frame(hog = "H", A = "g", B = "g"),
      c("A", "B"), list(A = list()), dummy_edges
    ),
    "networks missing entries for species"
  )
  expect_error(
    clique_persistence(
      data.frame(hog = "H", A = "g", B = "g"),
      c("A", "B"), list(A = list(), B = list()),
      data.frame(x = 1)
    ),
    "edges missing required columns"
  )
})


# --- Tests for cost_weights (composite backtracking cost) ---

test_that(
  "cost_weights default produces identical output to current behavior",
  {
    # Regression test: default cost_weights = c(q=1, effect=0) must give the

    # same result as the original mean-q-only ranking
    edges <- data.frame(
      gene1 = c("A1", "A1", "B1", "A2", "A2"),
      gene2 = c("B1", "C1", "C1", "B1", "C1"),
      species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
      species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
      hog = rep("HOG1", 5),
      type = rep("conserved", 5),
      q.value = c(0.01, 0.01, 0.01, 0.5, 0.5),
      effect_size = c(3.0, 3.0, 3.0, 1.0, 1.0),
      stringsAsFactors = FALSE
    )

    result_default <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
    result_explicit <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
      cost_weights = c(q = 1.0, effect = 0.0)
    )

    expect_identical(result_default, result_explicit)
    expect_equal(result_default$SP_A, "A1") # lower q wins
  }
)


test_that("cost_weights effect-only selects paralog with higher effect", {
  # Two paralogs in SP_A:
  #   A1 has lower q but lower effect
  #   A2 has higher q but higher effect
  # With cost_weights = c(q=0, effect=1), A2 should win (higher effect)
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1", "A2", "A2"),
    gene2 = c("B1", "C1", "C1", "B1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C"),
    hog = rep("HOG1", 5),
    type = rep("conserved", 5),
    q.value = c(0.01, 0.01, 0.01, 0.05, 0.05),
    effect_size = c(1.0, 1.0, 1.0, 5.0, 5.0),
    stringsAsFactors = FALSE
  )

  # Default: A1 wins (lower q)
  result_q <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  expect_equal(result_q$SP_A, "A1")

  # Effect-only: A2 wins (higher effect => lower cost = 0 - 1*mean_eff)
  result_eff <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    cost_weights = c(q = 0, effect = 1)
  )
  expect_equal(result_eff$SP_A, "A2")
})


test_that("cost_weights does not affect single-copy HOGs", {
  edges <- data.frame(
    gene1 = c("A1", "A1", "B1"),
    gene2 = c("B1", "C1", "C1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = rep("HOG1", 3),
    type = rep("conserved", 3),
    q.value = c(0.01, 0.02, 0.03),
    effect_size = c(2.0, 3.0, 4.0),
    stringsAsFactors = FALSE
  )

  result_default <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
  result_weighted <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"),
    cost_weights = c(q = 0.5, effect = 2.0)
  )

  expect_equal(result_default$SP_A, result_weighted$SP_A)
  expect_equal(result_default$SP_B, result_weighted$SP_B)
  expect_equal(result_default$SP_C, result_weighted$SP_C)
  expect_equal(result_default$mean_q, result_weighted$mean_q)
  expect_equal(
    result_default$mean_effect_size, result_weighted$mean_effect_size
  )
})


test_that("cost_weights validation rejects bad input", {
  edges <- data.frame(
    gene1 = c("A1"), gene2 = c("B1"),
    species1 = c("SP_A"), species2 = c("SP_B"),
    hog = "HOG1", type = "conserved",
    q.value = 0.01, effect_size = 2.0,
    stringsAsFactors = FALSE
  )
  sp <- c("SP_A", "SP_B")

  # Not numeric
  expect_error(
    find_cliques(edges, sp, cost_weights = c("a", "b")),
    "cost_weights must be a named numeric"
  )
  # Wrong length
  expect_error(
    find_cliques(edges, sp, cost_weights = c(q = 1)),
    "cost_weights must be a named numeric vector of length 2"
  )
  # Missing names
  expect_error(
    find_cliques(edges, sp, cost_weights = c(1, 0)),
    "cost_weights must have names"
  )
  # Wrong names
  expect_error(
    find_cliques(edges, sp, cost_weights = c(x = 1, y = 0)),
    "cost_weights must have names 'q' and 'effect'"
  )
  # Negative value
  expect_error(
    find_cliques(edges, sp, cost_weights = c(q = -1, effect = 0)),
    "cost_weights values must be >= 0"
  )
})


test_that("cost_weights affects edge deduplication for duplicate gene pairs", {
  # HOG with 2 species. Two edges between the same gene pair A1-B1:
  # Edge 1: q=0.01, effect=1.0 (low q, low effect)
  # Edge 2: q=0.02, effect=10.0 (higher q, much higher effect)
  edges <- data.frame(
    gene1 = c("A1", "A1"),
    gene2 = c("B1", "B1"),
    species1 = c("SP_A", "SP_A"),
    species2 = c("SP_B", "SP_B"),
    hog = c("HOG1", "HOG1"),
    type = c("conserved", "conserved"),
    q.value = c(0.01, 0.02),
    effect_size = c(1.0, 10.0),
    stringsAsFactors = FALSE
  )
  sp <- c("SP_A", "SP_B")

  # Default (q-only): should keep q=0.01 edge
  r_default <- find_cliques(edges, sp, min_species = 2L)
  expect_equal(nrow(r_default), 1)
  expect_equal(r_default$mean_q, 0.01, tolerance = 1e-10)

  # Effect-only: should keep effect=10.0 edge (q=0.02)
  r_eff <- find_cliques(edges, sp,
    min_species = 2L,
    cost_weights = c(q = 0, effect = 1)
  )
  expect_equal(nrow(r_eff), 1)
  expect_equal(r_eff$mean_effect_size, 10.0, tolerance = 1e-10)

  # Mixed weights: cost = 1*q - 1*effect
  # Edge 1: 1*0.01 - 1*1.0  = -0.99
  # Edge 2: 1*0.02 - 1*10.0 = -9.98  <- lower cost, wins
  r_mix <- find_cliques(edges, sp,
    min_species = 2L,
    cost_weights = c(q = 1, effect = 1)
  )
  expect_equal(nrow(r_mix), 1)
  expect_equal(r_mix$mean_effect_size, 10.0, tolerance = 1e-10)
  expect_equal(r_mix$mean_q, 0.02, tolerance = 1e-10)
})


test_that(".maxent_scale finds the entropy maximum and its bracket", {
  set.seed(1)
  w <- rlnorm(200, meanlog = 0.5, sdlog = 1)
  z <- rcomplex:::.maxent_scale(w)

  # Garlaschelli et al. (2013): the root is bracketed by 1/max(w), 1/min(w).
  expect_gt(z, 1 / max(w))
  expect_lt(z, 1 / min(w))

  # It maximises the ensemble entropy, so a grid must not beat it.
  entropy <- function(zz) {
    p <- zz * w / (1 + zz * w)
    -sum(p * log(p) + (1 - p) * log1p(-p))
  }
  grid <- exp(seq(log(1 / max(w)), log(1 / min(w)), length.out = 400))
  expect_gte(entropy(z), max(vapply(grid, entropy, numeric(1))) - 1e-8)

  # A collapsed bracket (one weight, or all weights equal) has nothing to
  # separate: maximum entropy for a single probability is 1/2, which
  # z = 1/mean(w) delivers. Returning NA there would void the clique.
  expect_equal(rcomplex:::.maxent_scale(5), 1 / 5)
  expect_equal(rcomplex:::.maxent_scale(c(5, 5)), 1 / 5)
  z1 <- rcomplex:::.maxent_scale(5)
  expect_equal(5 * z1 / (1 + 5 * z1), 0.5)
  expect_true(is.na(rcomplex:::.maxent_scale(numeric(0))))
})


test_that(".onnela_weight is a probability and ignores the Jaccard column", {
  e <- data.frame(
    gene1 = paste0("a", 1:4), gene2 = paste0("b", 1:4),
    species1 = "SP_A", species2 = "SP_B", hog = "HOG1",
    q.value = 0.01, effect_size = c(1, 2, 4, 8),
    jaccard = c(0.9, 0.1, 0.5, 0.2), stringsAsFactors = FALSE
  )
  w <- rcomplex:::.onnela_weight(e)
  expect_true(all(w > 0 & w < 1))
  # Monotone in association strength, and unmoved by jaccard.
  expect_equal(order(w), order(e$effect_size))
  e2 <- e
  e2$jaccard <- rev(e$jaccard)
  expect_equal(rcomplex:::.onnela_weight(e2), w)

  # No usable effect_size: NA plus a classed warning, not a silent zero.
  # The warning is .frequency = "once", so it must be re-armed to be seen.
  e3 <- e
  e3$effect_size <- NA_real_
  rlang::local_options(rlib_warning_verbosity = "verbose")
  expect_warning(
    out <- rcomplex:::.onnela_weight(e3),
    class = "rcomplex_missing_effect_size"
  )
  expect_true(all(is.na(out)))
})
