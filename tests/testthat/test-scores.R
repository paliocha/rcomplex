# score (bits) and evalue on every edge table, and the clique score (WP13).

head_cols <- c(
  "gene1", "gene2", "hog", "score", "evalue", "q_value", "effect_size",
  "power"
)

test_that("score and evalue read the combined p on all three methods", {
  expect_scores <- function(e) {
    expect_gt(nrow(e), 0L)
    expect_identical(names(e)[seq_along(head_cols)], head_cols)
    expect_false(anyNA(e$p_value))
    expect_equal(e$evalue, e$n_tests * e$p_value)
    expect_equal(e$score, -log2(e$p_value))
  }
  f <- make_spec_nets()
  hy <- find_coexpressologs(f$networks, f$ortho, seed = 1L)
  rk <- find_coexpressologs(f$networks, f$ortho,
    method = "rank", null_networks = f$nulls
  )
  pm <- find_coexpressologs(f$networks, f$ortho,
    method = "permutation", seed = 1L
  )
  expect_scores(hy)
  expect_scores(rk)
  expect_scores(pm)

  # hypergeometric: the combined p is the larger directional p, and
  # n_tests counts the pairs of the direction with more tests
  cmp <- rcomplex:::compare_neighborhoods(
    f$networks$species1, f$networks$species2, f$ortho
  )
  expect_equal(
    hy$p_value, pmax(cmp$species1.p_value_con, cmp$species2.p_value_con)
  )
  expect_equal(unique(hy$n_tests), nrow(cmp))
  # "min" takes the smaller one
  hy_min <- find_coexpressologs(f$networks, f$ortho,
    pval_combine = "min", seed = 1L
  )
  expect_equal(
    hy_min$p_value, pmin(cmp$species1.p_value_con, cmp$species2.p_value_con)
  )
  # rank: the calibrated (empirical) p, not the raw rank p
  raw <- rcomplex:::compare_specificity(
    f$networks$species1, f$networks$species2, f$ortho
  )
  raw <- raw[match(rk$gene1, raw$gene1), ]
  expect_false(isTRUE(all.equal(
    rk$p_value, pmax(raw$species1.p_value, raw$species2.p_value)
  )))
  # permutation: the HOG p-value, joined onto its pairs
  hog <- withr::with_seed(1L, rcomplex:::permutation_hog_test(
    f$networks$species1, f$networks$species2, cmp
  ))
  expect_equal(pm$p_value, hog$p_value[match(pm$hog, hog$hog)])
})

test_that("an edge frame without p-value columns gets NA scores", {
  comp <- data.frame(
    gene1 = "A1", gene2 = "B1", hog = 1L,
    species1.effect_size = 2, species2.effect_size = 2,
    species1.q_value_con = 0.01, species2.q_value_con = 0.02
  )
  e <- rcomplex:::comparison_to_edges(comp, "SP_A", "SP_B")
  expect_identical(names(e)[seq_along(head_cols)], head_cols)
  expect_true(is.na(e$score) && is.na(e$evalue) && is.na(e$n_tests))
})

test_that("a null-network comparison has about one evalue below 1", {
  set.seed(7)
  n <- 300
  s <- 20
  load <- rep(1:10, each = 30)
  lat <- matrix(rnorm(10 * s), 10)
  mk <- function(prefix) {
    x <- lat[load, ] + matrix(rnorm(n * s), n)
    rownames(x) <- paste0(prefix, seq_len(n))
    x
  }
  xa <- mk("A")
  xb <- mk("B")
  na <- compute_network(xa, density = 0.05)
  nb <- compute_network(xb, density = 0.05)
  ortho <- data.frame(
    gene1 = rownames(xa), gene2 = rownames(xb),
    hog = paste0("HOG", seq_len(n))
  )
  real <- find_coexpressologs(list(A = na, B = nb), ortho, seed = 1L)
  null <- find_coexpressologs(
    list(A = na, B = null_network(xb, nb, seed = 2L)), ortho,
    seed = 1L
  )
  expect_gt(sum(real$evalue < 1), 50L)
  # E[#(evalue < 1)] <= 1 under the null; 3 leaves room for chance
  expect_lte(sum(null$evalue < 1), 3L)
})

test_that("a clique score is the sum of its member edges' scores", {
  fx <- make_clique_fixture_3sp()
  e <- fx$edges
  cons <- e[e$type == "conserved", ]
  sum_score <- function(hog, genes) {
    in_clique <- cons$hog == hog & cons$gene1 %in% genes &
      cons$gene2 %in% genes
    c(sum(cons$score[in_clique]), sum(in_clique))
  }

  cl <- find_cliques(e, fx$target_species)
  expect_gt(nrow(cl), 0L)
  for (i in seq_len(nrow(cl))) {
    got <- sum_score(cl$hog[i], unlist(cl[i, fx$target_species]))
    expect_equal(got[2], cl$n_edges[i])
    expect_equal(cl$score[i], got[1])
  }

  gg <- gene_clique_graph(e)
  expect_gt(nrow(gg), 0L)
  for (id in unique(gg$clique_id)) {
    m <- gg[gg$clique_id == id, ]
    got <- sum_score(m$hog[1], m$gene)
    expect_equal(got[2], m$n_edges[1])
    expect_equal(unique(m$score), got[1])
  }
})

test_that("the permutation path joins HOG results by id, not position", {
  f <- make_spec_nets()
  # integer ids in reverse order: positional indexing would pick the
  # wrong HOG's row of the p-sorted permutation table
  f$ortho$hog <- rev(seq_len(nrow(f$ortho)))
  pm <- find_coexpressologs(f$networks, f$ortho,
    method = "permutation", seed = 1L
  )
  cmp <- rcomplex:::compare_neighborhoods(
    f$networks$species1, f$networks$species2, f$ortho
  )
  hog <- withr::with_seed(1L, rcomplex:::permutation_hog_test(
    f$networks$species1, f$networks$species2, cmp
  ))
  at <- match(as.character(pm$hog), hog$hog)
  expect_false(anyNA(at))
  expect_equal(pm$p_value, hog$p_value[at])
  expect_equal(pm$q_value, hog$q_value[at])
})
