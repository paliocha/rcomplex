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
  expect_true(all(is.finite(e$power) & e$power >= 0 & e$power <= 1))
  # under the default reference rank every called edge has power >= 0.5
  expect_true(all(e$power[e$type == "conserved"] >= 0.5))
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

test_that("rank power falls with the reference rank and checks p0", {
  f <- make_spec_nets()
  nets <- f$networks
  cmp <- compare_specificity(nets$sp1, nets$sp2, f$ortho)
  null_p <- list(
    sp1 = compare_specificity(nets$sp1, f$nulls$sp2, f$ortho,
      directions = "1to2"
    )$Species1.p.val,
    sp2 = compare_specificity(f$nulls$sp1, nets$sp2, f$ortho,
      directions = "2to1"
    )$Species2.p.val
  )
  pw <- function(p0) {
    summarize_specificity(cmp, null_p,
      sp1 = "sp1", sp2 = "sp2", p0 = p0
    )$edges$power
  }
  lo <- pw(0.5)
  hi <- pw(0.01)
  expect_true(all(hi >= lo - 1e-12, na.rm = TRUE))
  expect_gt(mean(hi, na.rm = TRUE), mean(lo, na.rm = TRUE))
  expect_error(pw(0), "p0")
  expect_error(pw(c(0.1, 0.2)), "p0")
  # "min" never reports less power than "max" (either direction suffices)
  # at one fixed reference rank (the defaults differ between the modes)
  e_min <- summarize_specificity(cmp, null_p,
    sp1 = "sp1", sp2 = "sp2", pval_combine = "min", p0 = 0.01
  )$edges
  e_max <- summarize_specificity(cmp, null_p,
    sp1 = "sp1", sp2 = "sp2", p0 = 0.01
  )$edges
  expect_true(all(e_min$power >= e_max$power - 1e-12, na.rm = TRUE))
})

# A rank-test results frame with a known AUROC grid, for pinning
# .rank_power() by hand. Grid values G(f) = a - b * log(f) per row.
rank_frame <- function(q1, q2, p1, p2, a = 0.6, b = 0.01, t = 30L, n = 1000L) {
  gf <- rcomplex:::.rank_grid_frac
  grid <- function(a) {
    g <- matrix(rep(a - b * log(gf), each = length(q1)), length(q1))
    colnames(g) <- format(gf, scientific = FALSE, drop0trailing = TRUE)
    g
  }
  d <- data.frame(
    Species1.p.val = p1, Species1.q.val.con = q1,
    Species1.mapped = t, Species1.n.cand = n,
    Species2.p.val = p2, Species2.q.val.con = q2,
    Species2.mapped = t, Species2.n.cand = n
  )
  d$Species1.auroc.grid <- grid(a)
  d$Species2.auroc.grid <- grid(a)
  d
}

test_that(".rank_power() matches a hand computation", {
  rp <- rcomplex:::.rank_power
  hm <- function(a, t, m) {
    q1 <- a / (2 - a)
    q2 <- 2 * a^2 / (1 + a)
    sqrt((a * (1 - a) + (t - 1) * (q1 - a^2) + (m - 1) * (q2 - a^2)) /
           (t * m))
  }
  # rows 1-3 called both ways; raw p 0.001, 0.002, 0.004 -> cut 0.004
  d <- rank_frame(q1 = c(0.01, 0.02, 0.05, 0.5), q2 = c(0.01, 0.02, 0.05, 0.5),
                  p1 = c(0.001, 0.002, 0.004, 0.3),
                  p2 = c(0.001, 0.002, 0.004, 0.3))
  g <- function(f) 0.6 - 0.01 * log(f)
  a_ref <- g(0.002) # median raw p of called rows
  want <- stats::pnorm((a_ref - g(0.004)) / hm(a_ref, 30, 1000 - 1 - 30))
  expect_equal(rp(d, alpha = 0.1), rep(want, 4), tolerance = 1e-12)
  # a flat grid: reference and threshold need the same AUROC, power 0.5
  flat <- rank_frame(q1 = c(0.01, 0.5), q2 = c(0.01, 0.5),
                     p1 = c(0.001, 0.3), p2 = c(0.001, 0.3), b = 0)
  expect_equal(rp(flat, alpha = 0.1, p0 = 0.01), c(0.5, 0.5))
  # nothing called: nothing was detectable, power 0
  none <- rank_frame(q1 = c(0.5, 0.6), q2 = c(0.5, 0.6),
                     p1 = c(0.2, 0.3), p2 = c(0.2, 0.3))
  expect_identical(rp(none, alpha = 0.1), c(0, 0))
  # "min" with direction 2 never significant: direction 1's power
  one <- rank_frame(q1 = c(0.01, 0.02, 0.5), q2 = c(0.5, 0.6, 0.7),
                    p1 = c(0.001, 0.004, 0.3), p2 = c(0.2, 0.3, 0.4))
  a1 <- g(0.0025) # median of 0.001 and 0.004
  w1 <- stats::pnorm((a1 - g(0.004)) / hm(a1, 30, 1000 - 1 - 30))
  expect_equal(rp(one, alpha = 0.1, pval_combine = "min"), rep(w1, 3),
               tolerance = 1e-12)
  # under "max" the uncalled direction bounds it at 0
  expect_identical(rp(one, alpha = 0.1, pval_combine = "max"), c(0, 0, 0))
  # both directions call, but on disjoint rows: under "max" nothing is
  # called both ways, so each direction takes its reference from its own
  # calls -- here equal to its threshold, so power 0.5
  disjoint <- rank_frame(q1 = c(0.01, 0.5), q2 = c(0.5, 0.01),
                         p1 = c(0.001, 0.3), p2 = c(0.3, 0.001))
  expect_equal(rp(disjoint, alpha = 0.1, pval_combine = "max"), c(0.5, 0.5))
})

test_that("both routes to rank edges carry the same power", {
  f <- make_spec_nets()
  nets <- f$networks
  cmp <- compare_specificity(nets$sp1, nets$sp2, f$ortho)
  null_p <- list(
    sp1 = compare_specificity(nets$sp1, f$nulls$sp2, f$ortho,
      directions = "1to2"
    )$Species1.p.val,
    sp2 = compare_specificity(f$nulls$sp1, nets$sp2, f$ortho,
      directions = "2to1"
    )$Species2.p.val
  )
  sm <- summarize_specificity(cmp, null_p, sp1 = "sp1", sp2 = "sp2")
  two_step <- comparison_to_edges(sm$results, "sp1", "sp2")
  expect_false(anyNA(two_step$power))
  expect_identical(two_step$power, sm$edges$power)
})

test_that(".rank_power() branches: clamp, no room, weaker side, bad grid", {
  rp <- rcomplex:::.rank_power
  d <- rank_frame(q1 = c(0.01, 0.02, 0.5), q2 = c(0.01, 0.02, 0.5),
                  p1 = c(0.001, 0.004, 0.3), p2 = c(0.001, 0.004, 0.3))
  # a reference below the grid is read as its first knot, 1e-5
  expect_identical(rp(d, 0.1, p0 = 1e-7), rp(d, 0.1, p0 = 1e-5))
  # no candidates left besides the translated set: no power
  full <- d
  full$Species1.mapped <- full$Species1.n.cand - 1L
  expect_identical(rp(full, 0.1, p0 = 1e-3), c(0, 0, 0))
  # "max" takes the weaker direction when the directions differ
  steep <- rank_frame(q1 = c(0.01, 0.02, 0.5), q2 = c(0.01, 0.02, 0.5),
                      p1 = c(0.001, 0.004, 0.3), p2 = c(0.001, 0.004, 0.3),
                      b = 0.05)
  mixed <- d
  mixed$Species2.auroc.grid <- steep$Species2.auroc.grid
  w1 <- rp(d, 0.1, p0 = 1e-3)
  w2 <- rp(steep, 0.1, p0 = 1e-3)
  expect_false(isTRUE(all.equal(w1, w2)))
  expect_equal(rp(mixed, 0.1, p0 = 1e-3), pmin(w1, w2))
  # a grid that lost its fraction names gives no power (and says so)
  lost <- d
  colnames(lost$Species1.auroc.grid) <- NULL
  expect_warning(pw_lost <- rp(lost, 0.1), "fraction names")
  expect_true(all(is.na(pw_lost)))
  # so does a frame missing a supporting column
  short <- d
  short$Species2.n.cand <- NULL
  expect_warning(pw_short <- rp(short, 0.1), "Species2.n.cand")
  expect_true(all(is.na(pw_short)))
})

test_that("p0 reaches the rank power through find_coexpressologs()", {
  f <- make_spec_nets()
  lo <- spec_edges(f, p0 = 0.5)
  hi <- spec_edges(f, p0 = 1e-5)
  expect_false(isTRUE(all.equal(lo$power, hi$power)))
  expect_error(spec_edges(f, p0 = 2), "p0")
})

test_that("p0 reaches density_sweep() and is refused for other methods", {
  f <- make_spec_nets()
  sweep <- function(p0) {
    suppressMessages(density_sweep(f$networks, f$ortho,
      multipliers = 1, method = "rank", null_networks = f$nulls, p0 = p0
    ))$edges[[1]]$power
  }
  expect_false(isTRUE(all.equal(sweep(0.5), sweep(1e-5))))
  expect_error(sweep(2), "p0")
  expect_error(
    find_coexpressologs(f$networks, f$ortho, p0 = 0.1),
    "p0 is only used with method = \"rank\""
  )
})

test_that("a rank frame without its grid warns instead of going NA quietly", {
  f <- make_spec_nets()
  nets <- f$networks
  cmp <- compare_specificity(nets$sp1, nets$sp2, f$ortho)
  null_p <- list(
    sp1 = compare_specificity(nets$sp1, f$nulls$sp2, f$ortho,
      directions = "1to2"
    )$Species1.p.val,
    sp2 = compare_specificity(f$nulls$sp1, nets$sp2, f$ortho,
      directions = "2to1"
    )$Species2.p.val
  )
  res <- summarize_specificity(cmp, null_p)$results
  flat <- res[, !grepl("auroc\\.grid", names(res))]
  expect_warning(comparison_to_edges(flat, "sp1", "sp2"), "saveRDS")
  expect_no_warning(comparison_to_edges(res, "sp1", "sp2"))
})

test_that("comparison_to_edges() refuses p0 on a hypergeometric frame", {
  f <- make_spec_nets()
  hy <- summarize_comparison(
    compare_neighborhoods(f$networks$sp1, f$networks$sp2, f$ortho)
  )$results
  expect_error(comparison_to_edges(hy, "sp1", "sp2", p0 = 0.1), "rho0")
})

test_that("rho0 is refused on the rank path", {
  f <- make_spec_nets()
  expect_error(spec_edges(f, rho0 = 2), "rho0 is only used")
  nets <- f$networks
  cmp <- compare_specificity(nets$sp1, nets$sp2, f$ortho)
  res <- summarize_specificity(cmp)$results
  expect_error(comparison_to_edges(res, "sp1", "sp2", rho0 = 2), "p0 sets")
})
