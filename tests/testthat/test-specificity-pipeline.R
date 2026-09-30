# method = "rank" through find_coexpressologs(), density_sweep() and
# the clique consumers. Fixtures: make_spec_nets(), sparse_net() from
# helper-reference.R.

spec_edges <- function(f, ...) {
  find_coexpressologs(f$networks, f$ortho,
    method = "rank",
    null_networks = f$nulls, ...
  )
}

test_that("rank edges have the hypergeometric shape", {
  f <- make_spec_nets()
  e <- spec_edges(f)
  a <- find_coexpressologs(f$networks, f$ortho)
  expect_identical(names(e), names(a))
  expect_gt(nrow(e), 0L)
  expect_true(all(is.na(e$power)))
  expect_false(anyNA(e$type))
})

test_that("conserved calls fall in the shared cliques only", {
  f <- make_spec_nets()
  e <- spec_edges(f)
  idx <- as.integer(sub("^A", "", e$gene1))
  expect_gt(sum(e$type == "conserved" & idx %in% f$shared), 0L)
  expect_equal(sum(e$type == "conserved" & idx %in% f$one_sided), 0L)
  # the one-sided clique was really tested in direction 1 -> 2 and failed
  cmp <- compare_specificity(
    f$networks$sp1, f$networks$sp2, f$ortho,
    directions = "1to2"
  )
  p1 <- cmp$Species1.p.val[f$one_sided]
  expect_false(anyNA(p1))
  expect_true(all(p1 > 0.5))
})

test_that("summarize_specificity reports the pair counts", {
  f <- make_spec_nets()
  cmp <- compare_specificity(f$networks$sp1, f$networks$sp2, f$ortho)
  s <- summarize_specificity(cmp, sp1 = "sp1", sp2 = "sp2")
  expect_gt(s$summary$gene_pairs$total, 0L)
  expect_gt(s$summary$gene_pairs$reciprocal, 0L)
})

test_that("specificity arguments are validated", {
  f <- make_spec_nets()
  expect_error(
    find_coexpressologs(f$networks, f$ortho, null_networks = f$nulls),
    "only used with method"
  )
  expect_error(
    find_coexpressologs(f$networks, f$ortho, method = "rank"),
    "null_network\\(\\)"
  )
  expect_error(spec_edges(f, alternative = "less"), "greater")
  f_part <- f
  f_part$nulls$sp2 <- NULL
  expect_error(spec_edges(f_part), "missing: sp2")
  f_bad <- f
  rownames(f_bad$nulls$sp1$network)[1] <- "X"
  expect_error(spec_edges(f_bad), "not on the genes")
})

test_that("density_sweep at multiplier 1 equals find_coexpressologs", {
  f <- make_spec_nets()
  sw <- suppressMessages(density_sweep(f$networks, f$ortho,
    multipliers = 1, method = "rank", null_networks = f$nulls
  ))
  expect_equal(sw$edges[[1]], spec_edges(f))
})

test_that("clique consumers run on specificity edges", {
  f <- make_spec_nets()
  e <- spec_edges(f)
  sp <- c("sp1", "sp2")
  cl <- find_cliques(e, sp)
  expect_gt(nrow(cl), 0L)
  expect_no_error(classify_cliques(e, sp, c(sp1 = "a", sp2 = "b")))
  gcl <- gene_clique_graph(e, min_size = 2L)
  expect_gt(nrow(gcl), 0L)
  expect_no_error(classify_gene_cliques(gcl, e, sp))
})

test_that("dense and sparse input give the same specificity edges", {
  f <- make_spec_nets()
  lo <- min(f$networks$sp1$network[f$networks$sp1$network > 0])
  fs <- f
  # same neighbourhoods: the dense threshold 1 sits below every nonzero
  fs$networks <- lapply(f$networks, function(n) {
    sparse_net(modifyList(n, list(threshold = lo)))
  })
  expect_equal(spec_edges(fs), spec_edges(f))
})

test_that("coexpressolog_null refuses the specificity path", {
  f <- make_spec_nets()
  expect_error(
    coexpressolog_null(f$networks, f$ortho, method = "rank"),
    "hypergeometric path only"
  )
})

test_that("specificity p-values are uniform under independence", {
  mk <- function(seed, prefix) {
    set.seed(seed)
    x <- matrix(rnorm(80 * 15), 80, 15)
    rownames(x) <- paste0(prefix, 1:80)
    compute_network(x, sparse = FALSE)
  }
  n1 <- mk(1, "A")
  n2 <- mk(2, "B")
  ortho <- data.frame(
    Species1 = paste0("A", 1:80), Species2 = paste0("B", 1:80),
    hog = paste0("H", 1:80)
  )
  p <- compare_specificity(n1, n2, ortho)$Species1.p.val
  expect_gte(sum(!is.na(p)), 40L)
  expect_lt(abs(mean(p, na.rm = TRUE) - 0.5), 0.1)
  lo <- mean(p <= 0.2, na.rm = TRUE)
  expect_gte(lo, 0.1)
  expect_lte(lo, 0.3)

  # positive control: the same expression under the partner's labels
  n2_same <- n1
  dimnames(n2_same$network) <- list(paste0("B", 1:80), paste0("B", 1:80))
  p_same <- compare_specificity(n1, n2_same, ortho)$Species1.p.val
  expect_lt(mean(p_same, na.rm = TRUE), 0.2)
})

test_that("\"analytical\" is still accepted as the hypergeometric arm", {
  fx <- make_spec_nets()
  old <- find_coexpressologs(fx$networks, fx$ortho,
                             method = "analytical", seed = 1L)
  new <- find_coexpressologs(fx$networks, fx$ortho,
                             method = "hypergeometric", seed = 1L)
  expect_identical(old, new)
  expect_gt(nrow(new), 0L)
})
