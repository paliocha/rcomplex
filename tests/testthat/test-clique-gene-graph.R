# Tests for gene_clique_graph() and classify_gene_cliques()

# All pairs of a species/gene set as a co-expressolog edge table.
gcg_pairs <- function(sp, genes, hog, q, effect = 1) {
  cmb <- utils::combn(length(sp), 2L)
  data.frame(
    gene1 = genes[cmb[1L, ]], gene2 = genes[cmb[2L, ]],
    species1 = sp[cmb[1L, ]], species2 = sp[cmb[2L, ]],
    hog = hog, q.value = q, effect_size = effect,
    stringsAsFactors = FALSE
  )
}

gcg_six <- c("SP_A", "SP_B", "SP_C", "SP_D", "SP_E", "SP_F")
gcg_lin <- c(
  SP_A = "L1", SP_B = "L1", SP_C = "L1",
  SP_D = "L2", SP_E = "L2", SP_F = "L2"
)

# Six-species fixture, one HOG per tier. Lineages L1 = A,B,C and
# L2 = D,E,F, so S = 6, P = 15, W = 6, X = 9 and the generalised
# cross_max is choose(5,2) - 6 = 4.
make_gcg_fixture <- function() {
  five <- gcg_six[1:5]

  # HOG1: every one of the 15 pairs significant.
  h1 <- gcg_pairs(
    gcg_six, paste0(c("a", "b", "c", "d", "e", "f"), 1),
    "HOG1", 0.01
  )

  # HOG2: all 15 pairs in the 0.9 graph, 11 significant at 0.1.
  # The four dropped pairs are SP_A's edges to B, C, D and E, so SP_A
  # still keeps one significant edge -- the S - 2 tolerance at work.
  q2 <- rep(0.01, 15)
  q2[1:4] <- 0.5
  h2 <- gcg_pairs(
    gcg_six, paste0(c("a", "b", "c", "d", "e", "f"), 2),
    "HOG2", q2
  )

  # HOG3: five species, all 10 pairs significant, SP_F not in the HOG.
  h3 <- gcg_pairs(
    five, paste0(c("a", "b", "c", "d", "e"), 3),
    "HOG3", 0.01
  )

  # HOG4: same shape, but SP_F WAS tested against every member and
  # failed. This must not be read as an annotation gap.
  h4 <- gcg_pairs(
    five, paste0(c("a", "b", "c", "d", "e"), 4),
    "HOG4", 0.01
  )
  h4_f <- data.frame(
    gene1 = paste0(c("a", "b", "c", "d", "e"), 4), gene2 = "f4",
    species1 = five, species2 = "SP_F", hog = "HOG4",
    q.value = 0.95, effect_size = 1, stringsAsFactors = FALSE
  )

  # HOG5: both lineages fully significant within, all nine cross pairs
  # present in the 0.9 graph but none significant.
  q5 <- rep(0.5, 15)
  cmb <- utils::combn(6L, 2L)
  same <- gcg_lin[gcg_six[cmb[1L, ]]] == gcg_lin[gcg_six[cmb[2L, ]]]
  q5[same] <- 0.01
  h5 <- gcg_pairs(
    gcg_six, paste0(c("a", "b", "c", "d", "e", "f"), 5),
    "HOG5", q5
  )

  # HOG6: a complete L1 clique; L2 has no genes in the HOG at all.
  h6 <- gcg_pairs(
    gcg_six[1:3], paste0(c("a", "b", "c"), 6),
    "HOG6", 0.01
  )

  # HOG7: like HOG3, but SP_F is in the HOG through a gene pair that
  # never touches the clique -- present yet untested.
  h7 <- gcg_pairs(
    five, paste0(c("a", "b", "c", "d", "e"), 7),
    "HOG7", 0.01
  )
  h7_f <- data.frame(
    gene1 = "a7b", gene2 = "f7", species1 = "SP_A",
    species2 = "SP_F", hog = "HOG7", q.value = 0.02,
    effect_size = 1, stringsAsFactors = FALSE
  )

  rbind(h1, h2, h3, h4, h4_f, h5, h6, h7, h7_f)
}


test_that("gene_clique_graph returns one row per clique member", {
  e <- data.frame(
    gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "HOG1", q.value = c(0.01, 0.02, 0.03),
    effect_size = c(3, 2, 4), stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e)
  expect_equal(nrow(cl), 3L)
  expect_equal(length(unique(cl$clique_id)), 1L)
  expect_setequal(cl$species, c("SP_A", "SP_B", "SP_C"))
  expect_setequal(cl$gene, c("a1", "b1", "c1"))
  expect_true(all(cl$n_members == 3L))
  expect_true(all(cl$n_species == 3L))
  expect_true(all(cl$n_edges == 3L))
  expect_equal(unique(cl$mean_q), mean(c(0.01, 0.02, 0.03)))
  expect_equal(unique(cl$max_q), 0.03)
  expect_equal(unique(cl$mean_effect_size), 3)
})


test_that("gene_clique_graph reports every paralog combination", {
  # Two SP_A copies both clique with b1 and c1: find_cliques() would
  # pick one assignment, the gene graph keeps both.
  e <- data.frame(
    gene1 = c("a1", "a1", "a2", "a2", "b1"),
    gene2 = c("b1", "c1", "b1", "c1", "c1"),
    species1 = c("SP_A", "SP_A", "SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_B", "SP_C", "SP_C"),
    hog = "HOG1", q.value = 0.01, stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e)
  expect_equal(length(unique(cl$clique_id)), 2L)
  by_id <- split(cl$gene, cl$clique_id)
  expect_setequal(vapply(by_id, function(g) {
    intersect(g, c("a1", "a2"))
  }, character(1)), c("a1", "a2"))
  expect_true(all(cl$n_members == 3L))
})


test_that("gene_clique_graph honours min_size and alpha_graph", {
  e <- data.frame(
    gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "HOG1", q.value = c(0.01, 0.02, 0.5),
    stringsAsFactors = FALSE
  )
  # b1-c1 drops out at 0.1, leaving a star with no 3-clique.
  expect_equal(nrow(gene_clique_graph(e, alpha_graph = 0.1)), 0L)
  expect_equal(nrow(gene_clique_graph(e, alpha_graph = 0.9)), 3L)
  # At min_size 2 the star yields two 2-cliques.
  cl2 <- gene_clique_graph(e, min_size = 2L, alpha_graph = 0.1)
  expect_equal(length(unique(cl2$clique_id)), 2L)
})


test_that("gene_clique_graph collapses reversed duplicate rows", {
  e <- data.frame(
    gene1 = c("a1", "a1", "b1", "b1"),
    gene2 = c("b1", "c1", "c1", "a1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_A"),
    hog = "HOG1", q.value = c(0.04, 0.02, 0.03, 0.01),
    stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e)
  expect_true(all(cl$n_edges == 3L))
  # The duplicated A-B pair keeps its most significant copy.
  expect_equal(unique(cl$max_q), 0.03)
})


test_that("gene_clique_graph validates its inputs", {
  e <- data.frame(
    gene1 = "a", gene2 = "b", species1 = "SP_A",
    species2 = "SP_B", hog = "H", q.value = 0.01
  )
  expect_error(gene_clique_graph(e[, 1:3]), "missing required columns")
  expect_error(gene_clique_graph(e, min_size = 1L), "min_size")
  expect_error(gene_clique_graph(e, alpha_graph = NA), "alpha_graph")
  bad <- e
  bad$gene1 <- "a\x01b"
  expect_error(gene_clique_graph(bad), "must not contain")
})


test_that("gene_clique_graph returns a typed empty frame", {
  e <- data.frame(
    gene1 = "a", gene2 = "b", species1 = "SP_A",
    species2 = "SP_B", hog = "H", q.value = 0.5,
    effect_size = 1
  )
  out <- gene_clique_graph(e, alpha_graph = 0.1)
  expect_equal(nrow(out), 0L)
  expect_true(all(c(
    "clique_id", "hog", "species", "gene",
    "n_members", "n_edges", "mean_effect_size"
  ) %in%
    names(out)))
})


test_that("id_prefix keeps clique ids distinct across graph runs", {
  e <- make_gcg_fixture()
  a <- gene_clique_graph(e, alpha_graph = 0.1, id_prefix = "call_")
  b <- gene_clique_graph(e, alpha_graph = 0.9, id_prefix = "loose_")
  expect_length(intersect(a$clique_id, b$clique_id), 0L)
  expect_true(all(grepl("^call_", a$clique_id)))
})


test_that("floor diagnostics report the available q resolution", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_equal(attr(cl, "q_floor"), 0.01)
  expect_equal(attr(cl, "mean_q_floor"), 0.01)
  # Several HOGs sit at mean_q 0.01: the floor is not a unique rank.
  expect_gt(attr(cl, "n_cliques_at_q_floor"), 1L)
})


test_that("the five tiers are recovered on the six-species fixture", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  cls <- stats::setNames(res$classification, res$hog)
  expect_equal(cls[["HOG1"]], "complete_conserved")
  expect_equal(cls[["HOG2"]], "partial_significant")
  expect_equal(cls[["HOG3"]], "partial_present")
  expect_equal(cls[["HOG5"]], "differentiated")
  expect_equal(cls[["HOG6"]], "lineage_specific")
  expect_equal(cls[["HOG7"]], "partial_present")
})


test_that("the gap tier refuses a species that was tested and failed", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  h4 <- res[res$hog == "HOG4", ]
  expect_equal(h4$classification, "unclassified")
  expect_equal(h4$missing_species, "SP_F")
  expect_equal(h4$missing_reason, "tested_ns")
  # HOG3 and HOG4 differ ONLY in whether SP_F has rows, and the tier
  # call follows that difference.
  h3 <- res[res$hog == "HOG3", ]
  expect_equal(h3$classification, "partial_present")
  expect_equal(h3$missing_reason, "absent")
  expect_equal(h3$n_sig, h4$n_sig)
})


test_that("a species present but never compared reads as untested", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  h7 <- res[res$hog == "HOG7", ]
  expect_equal(h7$missing_reason, "untested")
  expect_equal(h7$classification, "partial_present")
})


test_that("a joinable species is flagged extendable, not a gap", {
  # Build the graph at 0.005 so the 0.01 edges to SP_C are invisible;
  # classifying at alpha_call 0.1 then sees a species that could join.
  e <- data.frame(
    gene1 = c("a1", "a1", "a1", "b1", "b1", "c1"),
    gene2 = c("b1", "c1", "d1", "c1", "d1", "d1"),
    species1 = c("SP_A", "SP_A", "SP_A", "SP_B", "SP_B", "SP_C"),
    species2 = c("SP_B", "SP_C", "SP_D", "SP_C", "SP_D", "SP_D"),
    hog = "HOG1",
    q.value = c(0.001, 0.01, 0.001, 0.01, 0.001, 0.01),
    stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.005)
  res <- classify_gene_cliques(
    cl, e,
    c("SP_A", "SP_B", "SP_C", "SP_D")
  )
  expect_equal(res$missing_species, "SP_C")
  expect_equal(res$missing_reason, "extendable")
  expect_equal(res$classification, "unclassified")
})


test_that("lineage_specific outranks partial_present", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  # With max_gap 3 the HOG6 clique satisfies partial_present as well;
  # the waterfall must still call it lineage_specific.
  res <- classify_gene_cliques(cl, e, gcg_six,
    lineage = gcg_lin,
    max_gap = 3L
  )
  expect_equal(
    res$classification[res$hog == "HOG6"],
    "lineage_specific"
  )
  # Without lineages the same clique falls through to the gap tier.
  res0 <- classify_gene_cliques(cl, e, gcg_six, max_gap = 3L)
  expect_equal(
    res0$classification[res0$hog == "HOG6"],
    "partial_present"
  )
})


test_that("lineage_specific refuses a tested, diverged outside species", {
  # HOG6 (L2 absent from the orthogroup) is lineage_specific; the same
  # L1 triangle with L2 compared against every member and rejected is
  # not. The published workflow requires Cross == 0 in a matrix whose
  # 1s mark a *tested* pair, so a rejected comparison disqualifies
  # there too, and its differentiated set is the one that scores it.
  l1 <- gcg_six[1:3]
  h <- gcg_pairs(l1, c("a1", "b1", "c1"), "HOG1", 0.01)
  cross <- do.call(rbind, lapply(gcg_six[4:6], function(s) {
    data.frame(
      gene1 = c("a1", "b1", "c1"), gene2 = paste0(s, "_g"),
      species1 = l1, species2 = s, hog = "HOG1", q.value = 0.95,
      effect_size = 1, stringsAsFactors = FALSE
    )
  }))
  e <- rbind(h, cross)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  expect_equal(res$missing_reason, "tested_ns,tested_ns,tested_ns")
  expect_equal(res$classification, "unclassified")

  # And the tier that does own this pattern picks it up, once the
  # outside species are testable among themselves and the cliques come
  # from the unfiltered graph.
  cmb <- utils::combn(6L, 2L)
  same <- gcg_lin[gcg_six[cmb[1L, ]]] == gcg_lin[gcg_six[cmb[2L, ]]]
  e2 <- gcg_pairs(
    gcg_six, paste0(c("a", "b", "c", "d", "e", "f"), 1),
    "HOG1", ifelse(same, 0.01, 0.95)
  )
  cl2 <- gene_clique_graph(e2, alpha_graph = Inf)
  res2 <- classify_gene_cliques(cl2, e2, gcg_six, lineage = gcg_lin)
  expect_equal(res2$n_members, 6L)
  expect_equal(res2$classification, "differentiated")
})


test_that("lineage tiers are skipped when no lineage is supplied", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six)
  expect_true(all(is.na(res$n_sig_within)))
  expect_true(all(is.na(res$n_sig_cross)))
  expect_false("differentiated" %in% res$classification)
  expect_equal(res$classification[res$hog == "HOG5"], "unclassified")
})


test_that("cross_max bounds the differentiated tier", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six,
    lineage = gcg_lin,
    cross_max = 0
  )
  expect_equal(res$classification[res$hog == "HOG5"], "differentiated")
  expect_equal(attr(res, "cross_max"), 0)

  # One significant cross-lineage pair now exceeds cross_max = 0.
  e2 <- e
  hit <- e2$hog == "HOG5" & e2$species1 == "SP_A" &
    e2$species2 == "SP_D"
  expect_equal(sum(hit), 1L)
  e2$q.value[hit] <- 0.01
  cl2 <- gene_clique_graph(e2, alpha_graph = 0.9)
  res2 <- classify_gene_cliques(cl2, e2, gcg_six,
    lineage = gcg_lin,
    cross_max = 0
  )
  expect_equal(
    res2$classification[res2$hog == "HOG5"],
    "unclassified"
  )
  expect_equal(res2$n_sig_cross[res2$hog == "HOG5"], 1L)
})


test_that("the derived constants match the published six-species run", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  expect_equal(attr(res, "n_pairs_total"), 15)
  expect_equal(attr(res, "n_within_pairs"), 6)
  expect_equal(attr(res, "n_cross_pairs"), 9)
  expect_equal(attr(res, "cross_max"), 4)
  # complete needs 15, partial_significant 11, partial_present 10.
  expect_equal(res$n_sig[res$hog == "HOG1"], 15L)
  expect_equal(res$n_sig[res$hog == "HOG2"], 11L)
  expect_equal(res$n_sig[res$hog == "HOG3"], 10L)
  expect_equal(choose(6 - 1, 2) + 1, choose(6, 2) - (6 - 2))
})


test_that("no threshold is tied to six species", {
  # The same code on four species: complete needs choose(4,2) = 6 and
  # partial_significant choose(3,2) + 1 = 4.
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  lin <- c(SP_A = "L1", SP_B = "L1", SP_C = "L2", SP_D = "L2")
  g <- paste0(c("a", "b", "c", "d"), 1)
  q <- rep(0.01, 6)
  q[1:2] <- 0.5
  e <- rbind(
    gcg_pairs(sp, g, "HOG1", 0.01),
    gcg_pairs(sp, paste0(c("a", "b", "c", "d"), 2), "HOG2", q)
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp, lineage = lin)
  cls <- stats::setNames(res$classification, res$hog)
  expect_equal(cls[["HOG1"]], "complete_conserved")
  expect_equal(cls[["HOG2"]], "partial_significant")
  expect_equal(res$n_sig[res$hog == "HOG2"], 4L)
  expect_equal(attr(res, "n_pairs_total"), 6)
  expect_equal(attr(res, "cross_max"), max(0, choose(3, 2) - 2))
})


test_that("hog_class applies the published HOG-level precedence", {
  # Two cliques of one HOG: a complete one and a paralog copy that is
  # only partially significant. The HOG takes the earlier tier.
  sp <- c("SP_A", "SP_B", "SP_C")
  e <- rbind(
    gcg_pairs(sp, c("a1", "b1", "c1"), "HOG1", 0.01),
    data.frame(
      gene1 = c("a2", "a2", "b1"), gene2 = c("b1", "c1", "c1"),
      species1 = c("SP_A", "SP_A", "SP_B"),
      species2 = c("SP_B", "SP_C", "SP_C"), hog = "HOG1",
      q.value = c(0.5, 0.5, 0.01), effect_size = 1,
      stringsAsFactors = FALSE
    )
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(length(unique(res$clique_id)), 2L)
  expect_setequal(
    res$classification,
    c("complete_conserved", "unclassified")
  )
  expect_true(all(res$hog_class == "complete_conserved"))
})


test_that("classify_gene_cliques recomputes q from the edge table", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  h2 <- res[res$hog == "HOG2", ]
  expect_equal(h2$n_present, 15L)
  expect_equal(h2$mean_q, mean(c(rep(0.5, 4), rep(0.01, 11))))
  expect_equal(h2$max_q, 0.5)
  expect_true("mean_effect_size" %in% names(res))
  # The floor is reported alongside the q-values it saturates.
  expect_equal(attr(res, "mean_q_floor"), 0.01)
  expect_gt(attr(res, "n_cliques_at_q_floor"), 1L)
})


test_that("classify_gene_cliques validates its inputs", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_error(
    classify_gene_cliques(cl[, 1:2], e, gcg_six),
    "cliques missing required columns"
  )
  expect_error(
    classify_gene_cliques(cl, e[, 1:3], gcg_six),
    "edges missing required columns"
  )
  expect_error(classify_gene_cliques(cl, e, "SP_A"), "at least 2")
  expect_error(
    classify_gene_cliques(cl, e, c("SP_A", "SP_A")),
    "at least 2"
  )
  expect_error(
    classify_gene_cliques(cl, e, gcg_six, max_gap = -1L),
    "max_gap"
  )
  expect_error(
    classify_gene_cliques(cl, e, gcg_six,
      lineage = gcg_lin,
      cross_max = -1
    ),
    "cross_max"
  )
  expect_error(
    classify_gene_cliques(cl, e, gcg_six,
      lineage = unname(gcg_lin)
    ),
    "named vector"
  )
  expect_error(
    classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin[1:3]),
    "lineage missing entries"
  )
})


test_that("classify_gene_cliques handles an empty clique table", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  out <- classify_gene_cliques(cl[0, ], e, gcg_six)
  expect_equal(nrow(out), 0L)
  expect_true(all(c(
    "classification", "missing_reason",
    "hog_class"
  ) %in% names(out)))
})


test_that("gene cliques carry at most one gene per species", {
  # Co-expressolog edges are cross-species, so adjacency inside a
  # clique cannot pair two genes of the same species.
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_true(all(cl$n_members == cl$n_species))
  per <- tapply(
    cl$species, cl$clique_id,
    function(s) anyDuplicated(s) == 0L
  )
  expect_true(all(per))
})


test_that("a one-member clique is never scored as conserved", {
  # choose(1, 2) == 0, so every tier count is vacuously satisfied
  # unless singletons are refused outright.
  cl <- data.frame(
    clique_id = "X_1", hog = "HOG9", species = "SP_A",
    gene = "a9", stringsAsFactors = FALSE
  )
  e <- data.frame(
    gene1 = "a9", gene2 = "b9", species1 = "SP_A",
    species2 = "SP_B", hog = "HOG9", q.value = 0.5,
    stringsAsFactors = FALSE
  )
  res <- classify_gene_cliques(cl, e, c("SP_A", "SP_B"),
    lineage = c(SP_A = "L1", SP_B = "L2")
  )
  expect_equal(res$classification, "unclassified")
  expect_equal(res$n_pairs, 0L)
})


test_that("a within-species edge refuses species-pair scoring", {
  # Two genes of one species inside a clique break the one-gene-per-
  # species invariant, so the pair arithmetic must not be applied.
  cl <- data.frame(
    clique_id = "X_1", hog = "HOG9",
    species = c("SP_A", "SP_A", "SP_B"),
    gene = c("a1", "a2", "b1"), stringsAsFactors = FALSE
  )
  e <- data.frame(
    gene1 = c("a1", "a1", "a2"), gene2 = c("a2", "b1", "b1"),
    species1 = c("SP_A", "SP_A", "SP_A"),
    species2 = c("SP_A", "SP_B", "SP_B"), hog = "HOG9",
    q.value = 0.01, stringsAsFactors = FALSE
  )
  res <- classify_gene_cliques(cl, e, c("SP_A", "SP_B"))
  expect_equal(res$n_members, 3L)
  expect_equal(res$n_species, 2L)
  expect_equal(res$classification, "unclassified")
})


test_that("a clique species outside `species` is refused", {
  # SP_Z was counted into the clique's own species total but never
  # into choose(S, 2), so the clique scored complete_conserved while
  # the same row reported SP_C missing.
  e <- gcg_pairs(
    c("SP_A", "SP_B", "SP_Z"), c("a1", "b1", "z1"), "HOG1", 0.01
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_error(
    classify_gene_cliques(cl, e, c("SP_A", "SP_B", "SP_C")),
    "absent from"
  )
})


test_that("repeated clique ids are refused, not silently merged", {
  e <- make_gcg_fixture()
  cl <- rbind(
    gene_clique_graph(e, alpha_graph = 0.1),
    gene_clique_graph(e, alpha_graph = 0.9)
  )
  expect_error(
    classify_gene_cliques(cl, e, gcg_six),
    "distinct id_prefix"
  )
})


test_that("colliding cliques with no shared member are still refused", {
  # The (clique_id, species, gene) triple only collides when the two
  # cliques share a member. Two runs whose same-numbered cliques are
  # disjoint used to merge silently into one six-member "clique" with
  # n_pairs 15, six of them scored and nine invented.
  sp <- c("SP_A", "SP_B", "SP_C")
  tri <- function(genes, q) {
    data.frame(
      gene1 = genes[c(1L, 1L, 2L)], gene2 = genes[c(2L, 3L, 3L)],
      species1 = c("SP_A", "SP_A", "SP_B"),
      species2 = c("SP_B", "SP_C", "SP_C"),
      hog = "HOG1", q.value = q, stringsAsFactors = FALSE
    )
  }
  e <- rbind(tri(c("a", "b", "c"), 0.5), tri(c("d", "ee", "f"), 0.01))
  tight <- gene_clique_graph(e, alpha_graph = 0.1)
  loose <- gene_clique_graph(e, alpha_graph = 0.9)
  # Pick the loose clique that shares no gene with the tight one rather
  # than assuming it is numbered HOG1_1. That numbering is
  # igraph::max_cliques()'s enumeration order, which is not part of its
  # contract, and the test would silently stop exercising the collision
  # if it flipped.
  tight_genes <- tight$gene
  share <- vapply(
    split(loose$gene, loose$clique_id),
    function(g) any(g %in% tight_genes), logical(1)
  )
  disjoint_id <- names(share)[!share][1L]
  expect_false(is.na(disjoint_id))
  # Overwrite the id rather than selecting on it. Selecting still needed
  # the disjoint loose clique to happen to carry the same number as the
  # tight one, which is the enumeration-order assumption in a different
  # spelling; if the order flipped, the fixture would build no collision
  # and the test would fail while exercising nothing.
  dis <- loose[loose$clique_id == disjoint_id, , drop = FALSE]
  dis$clique_id <- tight$clique_id[1L]
  both <- rbind(tight, dis)
  expect_equal(length(unique(both$clique_id)), 1L)
  expect_equal(anyDuplicated(paste(
    both$clique_id, both$species, both$gene
  )), 0L)
  expect_error(
    classify_gene_cliques(both, e, sp),
    "more member rows than n_members declares"
  )
})


test_that("one clique id may not span two hogs", {
  # The row-count signature needs n_members; a hand-built table without
  # it is caught by the id-to-hog check instead.
  sp <- c("SP_A", "SP_B", "SP_C")
  e <- rbind(
    gcg_pairs(sp, c("a1", "b1", "c1"), "HOG1", 0.01),
    gcg_pairs(sp, c("a2", "b2", "c2"), "HOG2", 0.01)
  )
  cl <- data.frame(
    clique_id = "C1", hog = rep(c("HOG1", "HOG2"), each = 3L),
    species = rep(sp, 2L),
    gene = c("a1", "b1", "c1", "a2", "b2", "c2"),
    stringsAsFactors = FALSE
  )
  expect_error(
    classify_gene_cliques(cl, e, sp),
    "spanning several hogs"
  )
})


test_that("a clique table without n_members is still accepted", {
  sp <- c("SP_A", "SP_B", "SP_C")
  e <- gcg_pairs(sp, c("a1", "b1", "c1"), "HOG1", 0.01)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  bare <- cl[, c("clique_id", "hog", "species", "gene")]
  res <- classify_gene_cliques(bare, e, sp)
  expect_equal(res$classification, "complete_conserved")
})


test_that("all-singleton lineages cannot make differentiated vacuous", {
  # all(logical(0)) is TRUE, so with no multi-species lineage the tier
  # collapsed to "few cross pairs significant" and passed a clique with
  # zero significant pairs.
  sp <- paste0("SP_", 1:6)
  lin <- stats::setNames(paste0("L", 1:6), sp)
  e <- gcg_pairs(sp, paste0("g", 1:6), "HOG1", 0.5)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp, lineage = lin)
  expect_equal(res$n_sig, 0L)
  expect_equal(res$n_sig_within, 0L)
  expect_equal(res$classification, "unclassified")
})


test_that("differentiated needs every member pair tested", {
  # Only the two within-lineage pairs have a row; the four cross pairs
  # were never tested. Absent evidence is not evidence of divergence.
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  lin <- c(SP_A = "L1", SP_B = "L1", SP_C = "L2", SP_D = "L2")
  e <- rbind(
    gcg_pairs(c("SP_A", "SP_B"), c("a1", "b1"), "HOG1", 0.01),
    gcg_pairs(c("SP_C", "SP_D"), c("c1", "d1"), "HOG1", 0.01)
  )
  cl <- data.frame(
    clique_id = "HOG1_1", hog = "HOG1", species = sp,
    gene = c("a1", "b1", "c1", "d1"), stringsAsFactors = FALSE
  )
  res <- classify_gene_cliques(cl, e, sp, lineage = lin)
  expect_equal(res$n_pairs, 6L)
  expect_equal(res$n_present, 2L)
  expect_equal(res$n_sig_cross, 0L)
  expect_equal(res$classification, "unclassified")
})


test_that("duplicate rows break their q tie on effect size", {
  # Both copies sit at the permutation floor, so a stable sort would
  # let input row order pick the effect carried into the ranking key.
  base <- gcg_pairs(
    c("SP_A", "SP_B", "SP_C"), c("a1", "b1", "c1"), "HOG1", 0.00071
  )
  hi <- base[1, ]
  hi$effect_size <- 50
  lo <- base[1, ]
  lo$effect_size <- 2
  e_hi <- rbind(hi, lo, base[2:3, ])
  e_lo <- rbind(lo, hi, base[2:3, ])
  m_hi <- gene_clique_graph(e_hi, alpha_graph = 0.9)
  m_lo <- gene_clique_graph(e_lo, alpha_graph = 0.9)
  expect_equal(unique(m_hi$mean_effect_size), mean(c(50, 1, 1)))
  expect_equal(unique(m_lo$mean_effect_size), mean(c(50, 1, 1)))
  # classify() recomputes from `edges` and must resolve the tie the
  # same way.
  res <- classify_gene_cliques(m_lo, e_lo, c("SP_A", "SP_B", "SP_C"))
  expect_equal(res$mean_effect_size, mean(c(50, 1, 1)))
})


test_that("duplicate rows cannot fake an extendable species", {
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  core <- gcg_pairs(sp[1:3], c("a1", "b1", "c1"), "HOG1", 0.001)
  # Three rows for ONE SP_D--a1 pair. Counting significant rows rather
  # than distinct members read that as SP_D joined to all three.
  dup <- data.frame(
    gene1 = c("d1", "a1", "d1"), gene2 = c("a1", "d1", "a1"),
    species1 = c("SP_D", "SP_A", "SP_D"),
    species2 = c("SP_A", "SP_D", "SP_A"), hog = "HOG1",
    q.value = 0.01, effect_size = 1, stringsAsFactors = FALSE
  )
  e <- rbind(core, dup)
  cl <- gene_clique_graph(e, alpha_graph = 0.005)
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(res$missing_species, "SP_D")
  expect_equal(res$missing_reason, "tested_ns")
})


test_that("floor diagnostics survive subsetting and empty results", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  # subset() drops attributes; only the columns survive it.
  keep <- subset(res, res$n_sig > 0)
  expect_equal(unique(keep$mean_q_floor), attr(res, "mean_q_floor"))
  expect_equal(
    unique(keep$n_cliques_at_q_floor),
    attr(res, "n_cliques_at_q_floor")
  )
  expect_equal(unique(cl$mean_q_floor), attr(cl, "mean_q_floor"))
  # An empty result still carries the documented attributes.
  empty_cl <- gene_clique_graph(e, alpha_graph = 0.001)
  expect_equal(nrow(empty_cl), 0L)
  expect_equal(attr(empty_cl, "alpha_graph"), 0.001)
  expect_true(is.na(attr(empty_cl, "mean_q_floor")))
  expect_equal(attr(empty_cl, "n_cliques_at_q_floor"), 0L)
  empty_res <- classify_gene_cliques(cl[0, ], e, gcg_six)
  expect_equal(attr(empty_res, "n_pairs_total"), 15)
  expect_true(is.na(attr(empty_res, "mean_q_floor")))
  expect_equal(attr(empty_res, "n_cliques_at_q_floor"), 0L)
})


test_that("a within-species paralog edge is kept, not read as a loop", {
  # Only a row whose two endpoints are the SAME node is a self-loop.
  # A within-species edge between distinct paralogs is a real edge and
  # surfaces as n_species < n_members.
  e <- data.frame(
    gene1 = c("a1", "a1", "a2"), gene2 = c("a2", "b1", "b1"),
    species1 = "SP_A", species2 = c("SP_A", "SP_B", "SP_B"),
    hog = "HOG1", q.value = 0.01, stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_equal(unique(cl$n_members), 3L)
  expect_equal(unique(cl$n_species), 2L)
  expect_equal(unique(cl$n_edges), 3L)
  loop <- rbind(e, data.frame(
    gene1 = "a1", gene2 = "a1", species1 = "SP_A", species2 = "SP_A",
    hog = "HOG1", q.value = 0.01, stringsAsFactors = FALSE
  ))
  expect_equal(unique(gene_clique_graph(loop, alpha_graph = 0.9)$n_edges), 3L)
})


test_that("complete_conserved needs every pair, not all but one", {
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  q <- rep(0.01, 6)
  q[1] <- 0.5
  e <- gcg_pairs(sp, paste0(c("a", "b", "c", "d"), 1), "HOG1", q)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(res$n_sig, choose(4, 2) - 1)
  expect_equal(res$classification, "partial_significant")
})


test_that("partial_present needs every present pair significant", {
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  e <- gcg_pairs(sp[1:3], c("a1", "b1", "c1"), "HOG1",
    q = c(0.01, 0.01, 0.5)
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(res$n_species, 3L)
  expect_equal(res$n_sig, choose(3, 2) - 1)
  expect_equal(res$missing_reason, "absent")
  expect_equal(res$classification, "unclassified")
})


test_that("partial_significant needs every pair tested", {
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  # AB, AC, BC and AD are significant; BD and CD have no row at all.
  e <- rbind(
    gcg_pairs(sp[1:3], c("a1", "b1", "c1"), "HOG1", 0.01),
    data.frame(
      gene1 = "a1", gene2 = "d1", species1 = "SP_A",
      species2 = "SP_D", hog = "HOG1", q.value = 0.01,
      effect_size = 1, stringsAsFactors = FALSE
    )
  )
  cl <- data.frame(
    clique_id = "HOG1_1", hog = "HOG1", species = sp,
    gene = c("a1", "b1", "c1", "d1"), stringsAsFactors = FALSE
  )
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(res$n_sig, choose(3, 2) + 1)
  expect_equal(res$n_present, 4L)
  expect_equal(res$classification, "unclassified")
})


test_that("partial_significant needs every edge under alpha_graph", {
  sp <- c("SP_A", "SP_B", "SP_C", "SP_D")
  q <- rep(0.01, 6)
  q[5:6] <- 0.95
  e <- gcg_pairs(sp, paste0(c("a", "b", "c", "d"), 1), "HOG1", q)
  # Built on the unfiltered graph, so the 0.95 edges are in the clique.
  cl <- gene_clique_graph(e, alpha_graph = 1)
  res <- classify_gene_cliques(cl, e, sp, alpha_graph = 0.9)
  expect_equal(res$n_present, 6L)
  expect_equal(res$n_sig, choose(3, 2) + 1)
  expect_equal(res$max_q, 0.95)
  expect_equal(res$classification, "unclassified")
})


test_that("max_gap bounds how many species may be absent", {
  sp <- paste0("SP_", LETTERS[1:5])
  e <- gcg_pairs(sp[1:3], c("a1", "b1", "c1"), "HOG1", 0.01)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res1 <- classify_gene_cliques(cl, e, sp, max_gap = 1L)
  expect_equal(res1$n_missing, 2L)
  expect_equal(res1$classification, "unclassified")
  res2 <- classify_gene_cliques(cl, e, sp, max_gap = 2L)
  expect_equal(res2$classification, "partial_present")
})


test_that("a paralog inside a clique blocks the pair arithmetic", {
  # Four members over three species: the three cross-species pairs are
  # significant, so without the one-gene-per-species guard n_sig would
  # equal choose(3, 2) and the clique would score complete_conserved.
  sp <- c("SP_A", "SP_B", "SP_C")
  e <- data.frame(
    gene1 = c("a1", "a1", "b1", "a1", "a2", "a2"),
    gene2 = c("b1", "c1", "c1", "a2", "b1", "c1"),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A", "SP_A"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_A", "SP_B", "SP_C"),
    hog = "HOG1", q.value = c(0.01, 0.01, 0.01, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_equal(unique(cl$n_members), 4L)
  res <- classify_gene_cliques(cl, e, sp)
  expect_equal(res$n_species, 3L)
  expect_equal(res$n_pairs, 6L)
  expect_equal(res$n_sig, 3L)
  expect_equal(res$classification, "unclassified")
})


test_that("clique summaries use only the edges inside the clique", {
  # a1..d1 form a K4; e1 hangs off a1 alone. n_edges must be the six
  # internal edges and mean_q must average those six only.
  q4 <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.15)
  e <- rbind(
    gcg_pairs(
      c("SP_A", "SP_B", "SP_C", "SP_D"),
      c("a1", "b1", "c1", "d1"), "HOG1", q4
    ),
    data.frame(
      gene1 = "a1", gene2 = "e1", species1 = "SP_A",
      species2 = "SP_E", hog = "HOG1", q.value = 0.5,
      effect_size = 1, stringsAsFactors = FALSE
    )
  )
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  expect_equal(unique(cl$n_members), 4L)
  expect_equal(unique(cl$n_edges), 6L)
  expect_equal(unique(cl$mean_q), mean(q4))
  expect_equal(unique(cl$max_q), max(q4))
  # mean, not median: the two differ on this fixture.
  expect_false(isTRUE(all.equal(mean(q4), stats::median(q4))))
})


test_that("within and cross lineage counts are kept apart", {
  e <- make_gcg_fixture()
  hit <- e$hog == "HOG5" & e$species1 == "SP_A" & e$species2 == "SP_D"
  e$q.value[hit] <- 0.01
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  h5 <- res[res$hog == "HOG5", ]
  expect_equal(h5$n_species, 6L)
  expect_equal(h5$n_present, 15L)
  expect_equal(h5$n_sig, 7L)
  expect_equal(h5$n_sig_within, 6L)
  expect_equal(h5$n_sig_cross, 1L)
})


test_that("both thresholds are strict at q == alpha", {
  sp <- c("SP_A", "SP_B", "SP_C")
  e <- gcg_pairs(sp, c("a1", "b1", "c1"), "HOG1", c(0.01, 0.02, 0.1))
  # An edge exactly at alpha_graph is out of the graph.
  expect_equal(nrow(gene_clique_graph(e, alpha_graph = 0.1)), 0L)
  # A pair exactly at alpha_call is not a significant pair.
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, sp, alpha_call = 0.1)
  expect_equal(res$n_sig, 2L)
  expect_equal(res$classification, "partial_significant")
})


test_that("classification cost is not multiplied by unrelated edges", {
  # Matching each clique against the whole edge table separately made
  # this O(n_cliques x n_edges) and cost 18.7x on a real run. The
  # assertion counts elements scanned by match(), not seconds: a
  # wall-clock ratio on a shared CI runner fails for reasons that have
  # nothing to do with the code, and would have to be skipped, leaving
  # the regression unguarded.
  sp <- paste0("SP_", LETTERS[1:4])
  n_hog <- 1000L
  cmb <- utils::combn(4L, 2L)
  h <- rep(seq_len(n_hog), each = 6L)
  i1 <- rep(cmb[1L, ], times = n_hog)
  i2 <- rep(cmb[2L, ], times = n_hog)
  core <- data.frame(
    gene1 = paste0("g", h, "_", i1), gene2 = paste0("g", h, "_", i2),
    species1 = sp[i1], species2 = sp[i2], hog = paste0("HOG", h),
    q.value = 0.01, stringsAsFactors = FALSE
  )
  one <- gene_clique_graph(core, alpha_graph = 0.9)
  cl <- do.call(rbind, lapply(1:5, function(k) {
    x <- one
    x$clique_id <- paste0("r", k, "_", x$clique_id)
    x
  }))
  pad <- rbind(core, do.call(rbind, lapply(1:7, function(k) {
    x <- core
    x$hog <- paste0(x$hog, "_p", k)
    x$gene1 <- paste0(x$gene1, "_p", k)
    x$gene2 <- paste0(x$gene2, "_p", k)
    x
  })))
  # Shadow match() in a child of the package namespace and total the
  # lookup-table lengths it is handed. body(f) <- body(f) drops the
  # byte-compiled body, which would otherwise resolve match() straight
  # to base and never reach the shadow.
  scanned <- function(cliques, edges, species) {
    n <- 0
    f <- classify_gene_cliques
    env <- new.env(parent = environment(f))
    env$match <- function(x, table, ...) {
      n <<- n + length(table)
      base::match(x, table, ...)
    }
    environment(f) <- env
    body(f) <- body(f)
    res <- f(cliques, edges, species)
    list(res = res, scanned = n)
  }
  n_cl <- length(unique(cl$clique_id))
  a <- scanned(cl, core, sp)
  b <- scanned(cl, pad, sp)
  expect_equal(a$res$classification, b$res$classification)
  # One pass over the pair lookup, not one per clique. The regressed
  # form scanned n_cliques x n_edges, three orders of magnitude more.
  # SCOPE: the shadowed `match` only sees calls that resolve through
  # classify_gene_cliques()'s own body. `%in%` and any lookup inside
  # .gcg_classify_one() are invisible to it, so this pins the vectorised
  # pass where the fix lives and would NOT catch a regression that moved
  # a per-clique lookup down into the helper.
  expect_lt(a$scanned, nrow(core) + 100)
  expect_lt(b$scanned, nrow(pad) + 100)
  expect_lt(b$scanned, n_cl * nrow(pad) / 100)
})


test_that("alpha_call and alpha_graph are validated", {
  # Unchecked, an NA alpha makes every q < alpha comparison NA, so the
  # tier predicates evaluate to NA and the call dies inside a helper with
  # base R's "missing value where TRUE/FALSE needed", naming neither
  # argument. A length-2 vector silently used only its first element.
  e <- data.frame(
    gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "H1", q.value = 0.01, stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e, min_size = 3L, alpha_graph = 0.9)
  sp <- c("SP_A", "SP_B", "SP_C")
  for (bad in list(NA_real_, c(0.1, 0.2), "x", numeric(0))) {
    expect_error(
      classify_gene_cliques(cl, e, sp, alpha_call = bad),
      "alpha_call must be a single non-missing number"
    )
    expect_error(
      classify_gene_cliques(cl, e, sp, alpha_graph = bad),
      "alpha_graph must be a single non-missing number"
    )
  }
})


test_that("mean_q ties are counted to tolerance, not bitwise", {
  # mean_q is a mean over a different edge subset per clique, so two
  # mathematically equal values need not be bitwise equal. An exact ==
  # under-reports the very tie the count exists to report, and would
  # disagree with pvalue_resolution() on the same vector.
  eps <- .Machine$double.eps
  v <- 0.25 * (1 + c(0, 1, -1, 2) * eps)
  ties <- rcomplex:::.tol_min_ties(c(v, 0.5, 0.9))
  expect_equal(ties$n_at_min, 4L)
  expect_equal(ties$n_distinct, 3L)
  # and the exact comparison the fix replaced would have said 1 and 6
  expect_lt(ties$n_distinct, length(unique(c(v, 0.5, 0.9))))
})


test_that("near-tied mean_q is counted as tied at both call sites", {
  # The tolerant comparison exists for the two attributes below, but the
  # only test for it called .tol_min_ties() directly, so an exact ==
  # could have come back at either call site unnoticed. mean_q is a mean
  # over a different edge subset per clique, so build two cliques whose
  # mean_q is mathematically equal and bitwise apart.
  eps <- .Machine$double.eps
  q1 <- c(0.2, 0.4, 0.6)
  q2 <- q1 * (1 + c(1, -1, 2) * eps)
  mk <- function(tag, qs) {
    data.frame(
      gene1 = c("a", "a", "b"), gene2 = c("b", "c", "c"),
      species1 = c("SP_A", "SP_A", "SP_B"),
      species2 = c("SP_B", "SP_C", "SP_C"),
      hog = tag, q.value = qs, stringsAsFactors = FALSE
    )
  }
  e <- rbind(mk("H1", q1), mk("H2", q2))
  e$gene1 <- paste0(e$hog, "_", e$gene1)
  e$gene2 <- paste0(e$hog, "_", e$gene2)

  cl <- gene_clique_graph(e, min_size = 3L, alpha_graph = 0.9)
  expect_equal(length(unique(cl$clique_id)), 2L)
  # The test only discriminates if the two mean_q really are a near-tie:
  # distinct as doubles, but within the tolerance. Assert that rather
  # than trusting the last-bit arithmetic to have come out that way.
  mq <- sort(unique(cl$mean_q))
  expect_equal(length(mq), 2L)
  expect_gt(diff(mq), 0)
  expect_lt(diff(mq), rcomplex:::.tie_tol() * abs(mq[1L]))
  # Exact equality would report 1; the two mean_q differ only in the last
  # bits, so both cliques sit at the floor.
  expect_equal(attr(cl, "n_cliques_at_q_floor"), 2L)

  res <- classify_gene_cliques(cl, e, c("SP_A", "SP_B", "SP_C"),
    alpha_call = 0.9, alpha_graph = 0.9
  )
  expect_equal(unique(res$n_cliques_at_q_floor), 2L)
})


test_that("a row-filtered clique table is accepted, not called a merge", {
  # The duplicate-id guard checked the row count in both directions, but
  # a merge can only ADD rows. Refusing a table with FEWER rows than
  # n_members declares rejected the legitimate case -- a caller who has
  # subset the member rows -- and told them to set id_prefix, which is
  # not the problem.
  e <- data.frame(
    gene1 = c("a", "a", "b"), gene2 = c("b", "c", "c"),
    species1 = c("SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C"),
    hog = "H1", q.value = 0.01, stringsAsFactors = FALSE
  )
  cl <- gene_clique_graph(e, min_size = 3L, alpha_graph = 0.9)
  expect_equal(nrow(cl), 3L)

  trimmed <- cl[cl$species != "SP_C", , drop = FALSE]
  expect_lt(nrow(trimmed), unique(trimmed$n_members))
  expect_silent(rcomplex:::.gcg_check_clique_ids(
    trimmed$clique_id, paste(trimmed$species, trimmed$gene),
    trimmed$hog, trimmed$n_members
  ))

  # The regression was reported through classify_gene_cliques(), so pin
  # that path too and not only the internal check. n_members must be
  # recomputed from the surviving rows (2), not carried over from the
  # stale declared value on the untrimmed table (3).
  e_sub <- e[e$species1 != "SP_C" & e$species2 != "SP_C", , drop = FALSE]
  res <- suppressMessages(classify_gene_cliques(
    trimmed, e_sub, c("SP_A", "SP_B"),
    alpha_call = 0.9,
    alpha_graph = 0.9
  ))
  expect_identical(unique(res$n_members), 2L)

  # More rows than declared is still refused: that is a real merge.
  doubled <- rbind(cl, cl)
  doubled$gene <- paste0(doubled$gene, rep(c("", "_x"), each = nrow(cl)))
  expect_error(
    rcomplex:::.gcg_check_clique_ids(
      doubled$clique_id, paste(doubled$species, doubled$gene),
      doubled$hog, doubled$n_members
    ),
    "more member rows than n_members declares"
  )
})


# --- Underpowered tests (#12) ---

# An L1 triangle whose L2 species were each compared against every
# member and came back non-significant, at the given power.
gcg_up_lineage <- function(power) {
  l1 <- gcg_six[1:3]
  h <- gcg_pairs(l1, c("a1", "b1", "c1"), "HOG1", 0.01)
  h$power <- 0.99
  cross <- do.call(rbind, lapply(gcg_six[4:6], function(s) {
    data.frame(
      gene1 = c("a1", "b1", "c1"), gene2 = paste0(s, "_g"),
      species1 = l1, species2 = s, hog = "HOG1", q.value = 0.95,
      effect_size = 1, power = power, stringsAsFactors = FALSE
    )
  }))
  rbind(h, cross)
}

# Six species, both lineages fully significant within, every cross pair
# tested at q = 0.95 with the power given (recycled over the 9 pairs).
gcg_up_diff <- function(cross_power, cross_q = 0.95) {
  cmb <- utils::combn(6L, 2L)
  same <- gcg_lin[gcg_six[cmb[1L, ]]] == gcg_lin[gcg_six[cmb[2L, ]]]
  e <- gcg_pairs(
    gcg_six, paste0(c("a", "b", "c", "d", "e", "f"), 1),
    "HOG1", ifelse(same, 0.01, cross_q)
  )
  e$power <- 0.99
  e$power[!same] <- cross_power
  e
}


test_that("an outside species tested without power reads underpowered", {
  e <- gcg_up_lineage(0.1)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  expect_equal(
    res$missing_reason, "underpowered,underpowered,underpowered"
  )
  expect_equal(res$classification, "underpowered")
  expect_equal(res$hog_class, "underpowered")

  # One species underpowered, the others absent from the orthogroup.
  e1 <- e[e$species2 %in% gcg_six[1:4], ]
  cl1 <- gene_clique_graph(e1, alpha_graph = 0.9)
  res1 <- classify_gene_cliques(cl1, e1, gcg_six, lineage = gcg_lin)
  expect_equal(res1$missing_reason, "underpowered,absent,absent")
  expect_equal(res1$classification, "underpowered")
})


test_that("a powered or unmeasured failed test still blocks the tier", {
  for (pw in list(0.95, NA_real_)) {
    e <- gcg_up_lineage(pw)
    cl <- gene_clique_graph(e, alpha_graph = 0.9)
    res <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
    expect_equal(res$missing_reason, "tested_ns,tested_ns,tested_ns")
    expect_equal(res$classification, "unclassified")
  }

  # Power at or above min_power is powered.
  e <- gcg_up_lineage(0.1)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res <- classify_gene_cliques(cl, e, gcg_six,
    lineage = gcg_lin,
    min_power = 0.1
  )
  expect_equal(res$classification, "unclassified")

  # One powered failure among a species' tests keeps it tested_ns, and a
  # single powered species blocks the call however many others are not.
  e2 <- gcg_up_lineage(0.1)
  e2$power[e2$species2 == "SP_D"][1] <- 0.95
  cl2 <- gene_clique_graph(e2, alpha_graph = 0.9)
  res2 <- classify_gene_cliques(cl2, e2, gcg_six, lineage = gcg_lin)
  expect_equal(
    res2$missing_reason, "tested_ns,underpowered,underpowered"
  )
  expect_equal(res2$classification, "unclassified")
})


test_that("power moves only the calls it is meant to", {
  e <- make_gcg_fixture()
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  res0 <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  e$power <- 0.1
  res1 <- classify_gene_cliques(cl, e, gcg_six, lineage = gcg_lin)
  want <- res0$classification
  # HOG5's nine non-significant cross pairs now count against cross_max.
  want[res0$hog == "HOG5"] <- "underpowered"
  # HOG4's clique is one species short; that species was tested and
  # failed, but at power 0.1 the failure is no rejection, so it counts as a
  # gap and the clique is partial_present rather than unclassified.
  want[res0$hog == "HOG4"] <- "partial_present"
  expect_equal(res1$classification, want)
  # HOG6's outside species are absent, not tested: still lineage_specific.
  expect_equal(res1$classification[res1$hog == "HOG6"], "lineage_specific")
  expect_equal(res1$n_underpowered_cross[res1$hog == "HOG5"], 9L)
})


test_that("differentiated must survive underpowered cross pairs", {
  sp <- gcg_six
  run <- function(e) {
    cl <- gene_clique_graph(e, alpha_graph = Inf)
    classify_gene_cliques(cl, e, sp, lineage = gcg_lin)
  }
  # cross_max is 4: five underpowered pairs could hold a fifth
  # significant one, four cannot.
  five <- run(gcg_up_diff(c(rep(0.1, 5), rep(0.95, 4))))
  expect_equal(five$classification, "underpowered")
  expect_equal(five$n_underpowered_cross, 5L)
  expect_equal(five$n_sig_cross, 0L)

  four <- run(gcg_up_diff(c(rep(0.1, 4), rep(0.95, 5))))
  expect_equal(four$classification, "differentiated")
  expect_equal(four$n_underpowered_cross, 4L)

  # One significant cross pair plus four underpowered ones exceeds it.
  e <- gcg_up_diff(c(rep(0.1, 4), rep(0.95, 5)))
  cross <- which(e$power == 0.95)[1]
  e$q.value[cross] <- 0.01
  mixed <- run(e)
  expect_equal(mixed$n_sig_cross, 1L)
  expect_equal(mixed$classification, "underpowered")

  na_pw <- run(gcg_up_diff(NA_real_))
  expect_equal(na_pw$classification, "differentiated")
  expect_equal(na_pw$n_underpowered_cross, 0L)

  e_nl <- gcg_up_diff(0.1)
  cl_nl <- gene_clique_graph(e_nl, alpha_graph = Inf)
  no_lin <- classify_gene_cliques(cl_nl, e_nl, sp)
  expect_true(is.na(no_lin$n_underpowered_cross))
})


test_that("an underpowered missing species is a gap for partial_present", {
  # Five species form a significant clique; the sixth was compared against
  # every member and failed. Rejected (powered) it blocks partial_present;
  # unknown (underpowered) it does not.
  five <- function(pw) {
    h <- gcg_pairs(gcg_six[1:5], paste0(c("a", "b", "c", "d", "e"), 1),
                   "HOG1", 0.01)
    h$power <- 0.99
    miss <- data.frame(
      gene1 = paste0(c("a", "b", "c", "d", "e"), 1), gene2 = "f_g",
      species1 = gcg_six[1:5], species2 = gcg_six[6], hog = "HOG1",
      q.value = 0.95, effect_size = 1, power = pw,
      stringsAsFactors = FALSE
    )
    rbind(h, miss)
  }
  run <- function(pw) {
    e <- five(pw)
    cl <- gene_clique_graph(e, alpha_graph = 0.9)
    classify_gene_cliques(cl, e, gcg_six)
  }

  up <- run(0.1)
  expect_equal(up$missing_reason, "underpowered")
  expect_equal(up$classification, "partial_present")

  for (pw in list(0.95, NA_real_)) {
    blocked <- run(pw)
    expect_equal(blocked$missing_reason, "tested_ns")
    expect_false(blocked$classification == "partial_present")
  }
})


test_that("min_power is validated", {
  e <- gcg_up_lineage(0.1)
  cl <- gene_clique_graph(e, alpha_graph = 0.9)
  for (bad in list(-0.1, 1.5, NA_real_, c(0.5, 0.6), "0.8")) {
    expect_error(
      classify_gene_cliques(cl, e, gcg_six, min_power = bad),
      "min_power must be a single number in \\[0, 1\\]"
    )
  }
  empty <- classify_gene_cliques(cl[0, ], e, gcg_six)
  expect_true("n_underpowered_cross" %in% names(empty))
})
