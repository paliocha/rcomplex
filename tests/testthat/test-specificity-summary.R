# summarize_specificity() on synthetic comparison frames, null_network()
# and .check_null_networks(). Nothing here runs compare_specificity().

spec_frame <- function(n = 60L, n2 = 200L, seed = 1L) {
  withr::with_seed(seed, {
    grid <- function() sample(seq_len(n2), n, replace = TRUE) / n2
    p1 <- grid()
    p2 <- grid()
    p1[c(3L, 17L)] <- NA
    p2[c(17L, 40L)] <- NA
    data.frame(
      Species1 = paste0("A_", seq_len(n)),
      Species2 = paste0("B_", seq_len(n)),
      hog = paste0("HOG", rep(seq_len(n / 2L), each = 2L)),
      Species1.neigh = 10L, Species1.mapped = 5L,
      Species1.auroc = stats::runif(n), Species1.p.val = p1,
      Species1.effect.size = stats::runif(n),
      Species1.jaccard = stats::runif(n),
      Species2.neigh = 8L, Species2.mapped = 4L,
      Species2.auroc = stats::runif(n), Species2.p.val = p2,
      Species2.effect.size = stats::runif(n),
      Species2.jaccard = stats::runif(n),
      stringsAsFactors = FALSE
    )
  })
}

spec_null <- function(n0 = 500L, n2 = 200L, seed = 2L) {
  withr::with_seed(seed, {
    grid <- function() sample(seq_len(n2), n0, replace = TRUE) / n2
    list(sp1 = grid(), sp2 = grid())
  })
}

edge_names <- c(
  "gene1", "gene2", "species1", "species2", "hog", "q.value",
  "effect_size", "jaccard", "power", "type"
)


test_that("empirical p is (1 + #null <= p) / (1 + n_null)", {
  cmp <- spec_frame()
  p0 <- spec_null()
  s <- summarize_specificity(cmp, null_p = p0)
  r <- s$results
  brute <- function(p, p0) {
    vapply(p, function(x) (1 + sum(p0 <= x)) / (1 + length(p0)), 1)
  }
  expect_equal(r$Species1.p.emp, brute(r$Species1.p.val, p0$sp1))
  expect_equal(r$Species2.p.emp, brute(r$Species2.p.val, p0$sp2))
  expect_equal(s$summary$n_null, c(sp1 = 500L, sp2 = 500L))

  # NA null draws are not draws
  p0_na <- list(sp1 = c(p0$sp1, NA, NA), sp2 = p0$sp2)
  expect_equal(summarize_specificity(cmp, null_p = p0_na)$results, r)
})


test_that("q-values are non-decreasing in p, with and without a null", {
  cmp <- spec_frame()
  mono <- function(q, p) all(diff(q[order(p)]) >= -1e-12)
  r <- summarize_specificity(cmp, null_p = spec_null())$results
  expect_true(mono(r$Species1.q.val.con, r$Species1.p.emp))
  expect_true(mono(r$Species2.q.val.con, r$Species2.p.emp))

  raw <- summarize_specificity(cmp)$results
  expect_false(any(grepl("p\\.emp$", names(raw))))
  expect_true(mono(raw$Species1.q.val.con, raw$Species1.p.val))
  expect_true(mono(raw$Species2.q.val.con, raw$Species2.p.val))
})


test_that("NA rows are dropped and counted", {
  cmp <- spec_frame()
  s <- summarize_specificity(cmp)
  expect_equal(nrow(s$results), 57L)
  expect_equal(s$summary$n_dropped, 3L)
  expect_equal(s$summary$n_null, c(sp1 = 0L, sp2 = 0L))
  expect_equal(s$summary$gene_pairs$total, 57L)
  expect_false(any(is.na(s$results$Species1.q.val.con)))
  expect_named(
    s$summary,
    c("gene_pairs", "genes", "orthogroups", "pi0", "n_null", "n_dropped")
  )
})


test_that("pi0_method = 'none' is Benjamini-Hochberg on the empirical p", {
  s <- summarize_specificity(spec_frame(),
    null_p = spec_null(), pi0_method = "none"
  )
  r <- s$results
  expect_equal(unname(s$summary$pi0), c(1, 1))
  expect_equal(r$Species1.q.val.con, p.adjust(r$Species1.p.emp, "BH"))
  expect_equal(r$Species2.q.val.con, p.adjust(r$Species2.p.emp, "BH"))
})


test_that("sp1/sp2 give an edge frame with NA power", {
  cmp <- spec_frame()
  expect_null(summarize_specificity(cmp)$edges)
  s <- summarize_specificity(cmp,
    null_p = spec_null(), sp1 = "SP_A", sp2 = "SP_B"
  )
  expect_named(s$edges, edge_names)
  expect_equal(nrow(s$edges), nrow(s$results))
  expect_true(all(is.na(s$edges$power)))
  expect_equal(
    s$edges,
    comparison_to_edges(s$results, "SP_A", "SP_B")
  )
  expect_error(summarize_specificity(cmp, sp1 = "SP_A"), "Both sp1 and sp2")
})


test_that("an all-NA comparison returns the empty structure", {
  cmp <- spec_frame()
  cmp$Species1.p.val <- NA_real_
  s <- summarize_specificity(cmp, null_p = spec_null(), sp1 = "A", sp2 = "B")
  expect_equal(nrow(s$results), 0L)
  expect_equal(s$summary$gene_pairs$total, 0L)
  expect_equal(s$summary$n_dropped, nrow(cmp))
  expect_equal(s$summary$n_null, c(sp1 = 500L, sp2 = 500L))
  expect_true(all(is.na(s$summary$pi0)))
  expect_named(s$edges, edge_names)
  expect_equal(nrow(s$edges), 0L)
})


test_that("summarize_specificity validates its input", {
  cmp <- spec_frame()
  expect_error(
    summarize_specificity(cmp[, setdiff(names(cmp), "Species2.p.val")]),
    "compare_specificity"
  )
  expect_error(summarize_specificity(cmp, null_p = 1:3), "null_p")
})


# ---- null_network() ------------------------------------------------------

null_fixture <- function() {
  x <- withr::with_seed(1L, matrix(stats::rnorm(480), 40L))
  rownames(x) <- paste0("g", seq_len(40L))
  list(x = x, net = compute_network(x, density = 0.1))
}


test_that("null_network keeps the universe and the parameters", {
  f <- null_fixture()
  nn <- null_network(f$x, f$net, seed = 1L)
  expect_identical(rownames(nn$network), rownames(f$net$network))
  expect_equal(nn$n_genes, f$net$n_genes)
  expect_identical(nn$params, f$net$params)
  expect_identical(nn$store_density, f$net$store_density)
  expect_s4_class(nn$network, "dgCMatrix")
  expect_false(identical(nn$network, f$net$network))
})


test_that("null_network is reproducible under seed and moves across seeds", {
  f <- null_fixture()
  a <- null_network(f$x, f$net, seed = 1L)$network
  expect_identical(null_network(f$x, f$net, seed = 1L)$network, a)
  expect_false(identical(null_network(f$x, f$net, seed = 2L)$network, a))
})


test_that(".check_null_networks validates coverage and universe", {
  f <- null_fixture()
  xb <- f$x
  rownames(xb) <- paste0("h", seq_len(40L))
  nets <- list(A = f$net, B = compute_network(xb, density = 0.1))
  null_a <- null_network(f$x, f$net, seed = 1L)
  null_b <- null_network(xb, nets$B, seed = 1L)
  sp <- c("A", "B")

  expect_null(.check_null_networks(NULL, nets, sp))
  expect_error(
    .check_null_networks(list(A = null_a), nets, sp),
    "missing: B.*null_network\\(\\)"
  )
  expect_error(
    .check_null_networks(list(A = null_a, B = null_a), nets, sp),
    "null network for B.*null_network\\(\\)"
  )
  expect_error(
    .check_null_networks(list(A = null_a, B = list()), nets, sp),
    "null_networks\\[\\['B'\\]\\].*null_network\\(\\)"
  )

  out <- .check_null_networks(
    list(A = null_a, B = list(null_b, null_b)), nets, sp
  )
  expect_named(out, sp)
  expect_identical(out$A, list(null_a))
  expect_identical(out$B, list(null_b, null_b))
})
