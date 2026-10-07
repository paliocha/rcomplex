test_that("summarize_comparison returns correct structure", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:10),
    gene2 = paste0("B_", 1:10),
    hog = rep(1:5, each = 2),
    species1.neigh = rep(10, 10),
    species1.ortho.neigh = rep(5, 10),
    species1.neigh.overlap = c(3, 0, 2, 4, 1, 3, 2, 0, 1, 5),
    species1.p_value_con = c(
      0.001, 1, 0.01, 0.0001, 0.5, 0.005, 0.05, 1, 0.3, 0.0001
    ),
    species1.p_value_div = c(
      0.99, 0.01, 0.9, 0.999, 0.5, 0.99, 0.9, 0.01, 0.7, 0.999
    ),
    species1.effect_size = c(5, 1, 3, 8, 1, 4, 2, 1, 1.5, 10),
    species2.neigh = rep(8, 10),
    species2.ortho.neigh = rep(4, 10),
    species2.neigh.overlap = c(2, 0, 1, 3, 0, 2, 1, 0, 1, 4),
    species2.p_value_con = c(0.01, 1, 0.1, 0.001, 1, 0.01, 0.1, 1, 0.5, 0.0001),
    species2.p_value_div = c(
      0.9, 0.01, 0.8, 0.99, 0.01, 0.9, 0.8, 0.01, 0.5, 0.999
    ),
    species2.effect_size = c(4, 1, 2, 6, 1, 3, 1.5, 1, 1, 8),
    stringsAsFactors = FALSE
  )

  result <- rcomplex:::summarize_comparison(comparison, pi0_method = "storey")

  expect_type(result, "list")
  expect_named(result, c("results", "summary"))
  expect_s3_class(result$results, "data.frame")
  expect_type(result$summary, "list")
  expect_named(result$summary, c("gene_pairs", "genes", "orthogroups", "pi0"))
})

test_that("zero-overlap rows are filtered by default", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:5),
    gene2 = paste0("B_", 1:5),
    hog = 1:5,
    species1.neigh = rep(10, 5),
    species1.ortho.neigh = rep(5, 5),
    species1.neigh.overlap = c(3, 0, 2, 0, 1),
    species1.p_value_con = c(0.001, 1, 0.01, 1, 0.5),
    species1.p_value_div = c(0.99, 0.01, 0.9, 0.01, 0.5),
    species1.effect_size = c(5, 1, 3, 1, 2),
    species2.neigh = rep(8, 5),
    species2.ortho.neigh = rep(4, 5),
    species2.neigh.overlap = c(2, 0, 1, 0, 1),
    species2.p_value_con = c(0.01, 1, 0.1, 1, 0.5),
    species2.p_value_div = c(0.9, 0.01, 0.8, 0.01, 0.5),
    species2.effect_size = c(4, 1, 2, 1, 1),
    stringsAsFactors = FALSE
  )

  result <- rcomplex:::summarize_comparison(comparison, pi0_method = "storey")
  expect_equal(nrow(result$results), 3) # rows 2 and 4 filtered

  result_no_filter <- rcomplex:::summarize_comparison(comparison,
    pi0_method = "storey",
    filter_zero = FALSE
  )
  expect_equal(nrow(result_no_filter$results), 5)
})

test_that("q-values are computed", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:5),
    gene2 = paste0("B_", 1:5),
    hog = 1:5,
    species1.neigh = rep(10, 5),
    species1.ortho.neigh = rep(5, 5),
    species1.neigh.overlap = rep(2, 5),
    species1.p_value_con = c(0.001, 0.01, 0.02, 0.03, 0.04),
    species1.p_value_div = c(0.9, 0.8, 0.7, 0.6, 0.5),
    species1.effect_size = rep(3, 5),
    species2.neigh = rep(8, 5),
    species2.ortho.neigh = rep(4, 5),
    species2.neigh.overlap = rep(2, 5),
    species2.p_value_con = c(0.002, 0.02, 0.03, 0.04, 0.05),
    species2.p_value_div = c(0.9, 0.8, 0.7, 0.6, 0.5),
    species2.effect_size = rep(2, 5),
    stringsAsFactors = FALSE
  )

  result <- rcomplex:::summarize_comparison(comparison, pi0_method = "storey")

  # q-value columns should exist
  expect_true("species1.q_value_con" %in% names(result$results))
  expect_true("species2.q_value_con" %in% names(result$results))

  # q-values should be >= raw p-values
  expect_true(all(result$results$species1.q_value_con >=
                    result$results$species1.p_value_con))
  expect_true(all(result$results$species2.q_value_con >=
                    result$results$species2.p_value_con))

  # Raw p-values should be unchanged
  expect_equal(
    result$results$species1.p_value_con,
    c(0.001, 0.01, 0.02, 0.03, 0.04)
  )
  expect_equal(
    result$results$species2.p_value_con,
    c(0.002, 0.02, 0.03, 0.04, 0.05)
  )
})

test_that("summary counts are correct", {
  comparison <- data.frame(
    gene1 = c("A_1", "A_1", "A_2"),
    gene2 = c("B_1", "B_2", "B_2"),
    hog = c(1, 1, 2),
    species1.neigh = rep(10, 3),
    species1.ortho.neigh = rep(5, 3),
    species1.neigh.overlap = rep(5, 3),
    species1.p_value_con = c(0.001, 0.5, 0.001),
    species1.p_value_div = c(0.99, 0.5, 0.99),
    species1.effect_size = c(5, 1, 5),
    species2.neigh = rep(8, 3),
    species2.ortho.neigh = rep(4, 3),
    species2.neigh.overlap = rep(4, 3),
    species2.p_value_con = c(0.001, 0.5, 0.001),
    species2.p_value_div = c(0.99, 0.5, 0.99),
    species2.effect_size = c(4, 1, 4),
    stringsAsFactors = FALSE
  )

  result <- rcomplex:::summarize_comparison(
    comparison, pi0_method = "storey", alpha = 0.05
  )

  expect_equal(result$summary$gene_pairs$total, 3)
  expect_equal(result$summary$orthogroups$total, 2)
})

test_that("empty comparison handled gracefully", {
  comparison <- data.frame(
    gene1 = character(0),
    gene2 = character(0),
    hog = integer(0),
    species1.neigh = integer(0),
    species1.ortho.neigh = integer(0),
    species1.neigh.overlap = integer(0),
    species1.p_value_con = numeric(0),
    species1.p_value_div = numeric(0),
    species1.effect_size = numeric(0),
    species2.neigh = integer(0),
    species2.ortho.neigh = integer(0),
    species2.neigh.overlap = integer(0),
    species2.p_value_con = numeric(0),
    species2.p_value_div = numeric(0),
    species2.effect_size = numeric(0),
    stringsAsFactors = FALSE
  )

  result <- rcomplex:::summarize_comparison(comparison, pi0_method = "storey")
  expect_equal(nrow(result$results), 0)
  expect_equal(result$summary$gene_pairs$total, 0L)
})

test_that("alternative='less' uses divergence p-values", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:5),
    gene2 = paste0("B_", 1:5),
    hog = 1:5,
    species1.neigh = rep(10, 5),
    species1.ortho.neigh = rep(5, 5),
    species1.neigh.overlap = c(0, 0, 0, 3, 5),
    species1.p_value_con = c(1, 1, 1, 0.01, 0.001),
    species1.p_value_div = c(0.001, 0.01, 0.02, 0.9, 0.99),
    species1.effect_size = c(0, 0, 0, 3, 5),
    species2.neigh = rep(8, 5),
    species2.ortho.neigh = rep(4, 5),
    species2.neigh.overlap = c(0, 0, 0, 2, 4),
    species2.p_value_con = c(1, 1, 1, 0.01, 0.001),
    species2.p_value_div = c(0.001, 0.01, 0.02, 0.9, 0.99),
    species2.effect_size = c(0, 0, 0, 3, 5),
    stringsAsFactors = FALSE
  )

  # With alternative="less", should use .p_value_div for thresholding
  result <- rcomplex:::summarize_comparison(comparison,
    pi0_method = "storey",
    alternative = "less", alpha = 0.05
  )

  # filter_zero defaults to FALSE for "less"
  expect_equal(nrow(result$results), 5)

  # q-value columns for divergence should exist
  expect_true("species1.q_value_div" %in% names(result$results))
  expect_true("species2.q_value_div" %in% names(result$results))

  # Divergence q-values for first three rows should be significant
  expect_true(result$results$species1.q_value_div[1] < 0.05)
  expect_true(result$results$species1.q_value_div[2] < 0.05)
  # Rows 4 and 5 have high div p-values, should not be significant
  expect_true(result$results$species1.q_value_div[4] > 0.05)
})

test_that("alternative='less' disables zero-overlap filtering by default", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:3),
    gene2 = paste0("B_", 1:3),
    hog = 1:3,
    species1.neigh = rep(10, 3),
    species1.ortho.neigh = rep(5, 3),
    species1.neigh.overlap = c(0, 0, 2),
    species1.p_value_con = c(1, 1, 0.01),
    species1.p_value_div = c(0.001, 0.01, 0.9),
    species1.effect_size = c(0, 0, 3),
    species2.neigh = rep(8, 3),
    species2.ortho.neigh = rep(4, 3),
    species2.neigh.overlap = c(0, 0, 1),
    species2.p_value_con = c(1, 1, 0.01),
    species2.p_value_div = c(0.001, 0.01, 0.9),
    species2.effect_size = c(0, 0, 2),
    stringsAsFactors = FALSE
  )

  # Zero-overlap rows are kept for divergence (the strongest signal)
  result <- rcomplex:::summarize_comparison(
    comparison, pi0_method = "storey", alternative = "less"
  )
  expect_equal(nrow(result$results), 3)

  # But can be overridden
  result_filtered <- rcomplex:::summarize_comparison(comparison,
    pi0_method = "storey",
    alternative = "less",
    filter_zero = TRUE
  )
  expect_equal(nrow(result_filtered$results), 1)
})


test_that("summarize_comparison with sp1/sp2 returns $edges", {
  comparison <- data.frame(
    gene1 = paste0("A_", 1:10),
    gene2 = paste0("B_", 1:10),
    hog = rep(1:5, each = 2),
    species1.neigh.overlap = c(5, 3, 0, 4, 2, 1, 6, 0, 3, 4),
    species2.neigh.overlap = c(4, 2, 0, 3, 1, 2, 5, 0, 4, 3),
    species1.p_value_con = c(
      0.001, 0.05, 0.9, 0.01, 0.1, 0.2, 0.001, 0.8, 0.03, 0.01
    ),
    species2.p_value_con = c(
      0.002, 0.06, 0.8, 0.02, 0.15, 0.25, 0.002, 0.7, 0.04, 0.02
    ),
    species1.p_value_div = rep(0.99, 10),
    species2.p_value_div = rep(0.99, 10),
    species1.effect_size = c(3.0, 1.5, 1.0, 2.5, 1.2, 1.1, 3.5, 1.0, 2.0, 2.5),
    species2.effect_size = c(2.5, 1.3, 1.0, 2.0, 1.1, 1.2, 3.0, 1.0, 2.5, 2.0)
  )

  # Without sp1/sp2: no $edges
  result1 <- rcomplex:::summarize_comparison(comparison, pi0_method = "storey")
  expect_null(result1$edges)

  # With sp1/sp2: has $edges
  result2 <- rcomplex:::summarize_comparison(comparison,
    pi0_method = "storey",
    sp1 = "SP_A", sp2 = "SP_B"
  )
  expect_true(!is.null(result2$edges))
  expect_true(is.data.frame(result2$edges))
  expect_true(all(c(
    "gene1", "gene2", "species1", "species2",
    "hog", "q_value", "effect_size", "type"
  ) %in%
    names(result2$edges)))
  expect_true(all(result2$edges$species1 == "SP_A"))
  expect_true(all(result2$edges$species2 == "SP_B"))

  # $edges should match calling comparison_to_edges separately
  separate <- rcomplex:::comparison_to_edges(result2$results, "SP_A", "SP_B")
  expect_equal(result2$edges, separate)
})


test_that("summarize_comparison errors when only one of sp1/sp2 provided", {
  comparison <- data.frame(
    gene1 = "A_1", gene2 = "B_1", hog = 1,
    species1.neigh.overlap = 5, species2.neigh.overlap = 4,
    species1.p_value_con = 0.01, species2.p_value_con = 0.02,
    species1.p_value_div = 0.99, species2.p_value_div = 0.99,
    species1.effect_size = 3.0, species2.effect_size = 2.5
  )

  expect_error(
    rcomplex:::summarize_comparison(comparison, sp1 = "SP_A"),
    "Both sp1 and sp2"
  )
  expect_error(
    rcomplex:::summarize_comparison(comparison, sp2 = "SP_B"),
    "Both sp1 and sp2"
  )
})


test_that(
  "summarize_comparison with sp1/sp2 returns empty $edges on zero rows",
  {
    # All zero overlap -> filtered out with default filter_zero=TRUE
    comparison <- data.frame(
      gene1 = c("A_1", "A_2"), gene2 = c("B_1", "B_2"),
      hog = c(1, 2),
      species1.neigh.overlap = c(0, 0), species2.neigh.overlap = c(0, 0),
      species1.p_value_con = c(1, 1), species2.p_value_con = c(1, 1),
      species1.p_value_div = c(0.5, 0.5), species2.p_value_div = c(0.5, 0.5),
      species1.effect_size = c(1, 1), species2.effect_size = c(1, 1)
    )

    result <- rcomplex:::summarize_comparison(comparison,
      pi0_method = "storey",
      sp1 = "SP_A", sp2 = "SP_B"
    )
    expect_equal(nrow(result$results), 0)
    expect_true(!is.null(result$edges))
    expect_equal(nrow(result$edges), 0)
    expect_true(all(c(
      "gene1", "gene2", "species1", "species2",
      "hog", "q_value", "effect_size", "type"
    ) %in%
      names(result$edges)))
  }
)


# ---- pi0 from randomized p-values (D4) ----

test_that(
  "summarize_comparison default estimates pi0 from randomized p-values",
  {
    td <- make_graded_nets()
    cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)

    set.seed(5)
    s <- rcomplex:::summarize_comparison(cmp)
    expect_named(s$summary$pi0, c("sp1", "sp2"))
    expect_true(all(s$summary$pi0 > 0 & s$summary$pi0 <= 1))

    # q-values are the exact p-values' BH values scaled by the recorded pi0
    r <- s$results
    expect_equal(
      r$species1.q_value_con,
      s$summary$pi0[["sp1"]] * p.adjust(r$species1.p_value_con, "BH")
    )
    expect_equal(
      r$species2.q_value_con,
      s$summary$pi0[["sp2"]] * p.adjust(r$species2.p_value_con, "BH")
    )

    # reproducible under set.seed()
    set.seed(5)
    expect_identical(rcomplex:::summarize_comparison(cmp), s)

    # divergence direction uses the lower tail: (div - eq) + U * eq
    set.seed(6)
    d <- rcomplex:::summarize_comparison(cmp, alternative = "less")
    expect_named(d$summary$pi0, c("sp1", "sp2"))
    expect_equal(
      d$results$species1.q_value_div,
      d$summary$pi0[["sp1"]] *
        p.adjust(d$results$species1.p_value_div, "BH")
    )
  }
)


test_that("pi0_method = 'none' and 'storey' behave as documented", {
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)

  none <- rcomplex:::summarize_comparison(cmp, pi0_method = "none")
  expect_equal(unname(none$summary$pi0), c(1, 1))
  expect_equal(
    none$results$species1.q_value_con,
    p.adjust(none$results$species1.p_value_con, "BH")
  )
  expect_equal(
    none$results$species2.q_value_con,
    p.adjust(none$results$species2.p_value_con, "BH")
  )

  st <- rcomplex:::summarize_comparison(cmp, pi0_method = "storey")
  ref <- compute_qvalues(st$results$species1.p_value_con, pi0_method = "storey")
  expect_equal(st$results$species1.q_value_con, ref$qvalues)
  expect_equal(st$summary$pi0[["sp1"]], ref$pi0)
  # storey / none do not touch the RNG
  set.seed(8)
  u1 <- runif(1)
  set.seed(8)
  invisible(rcomplex:::summarize_comparison(cmp, pi0_method = "storey"))
  expect_identical(runif(1), u1)
})


test_that("randomized pi0 requires the p_value_gt / p_value_eq columns", {
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)
  old <- cmp[, !grepl("p_value_(gt|eq)$", names(cmp))]
  expect_error(rcomplex:::summarize_comparison(old), "p_value_gt")
  expect_silent(rcomplex:::summarize_comparison(old, pi0_method = "storey"))
  expect_silent(rcomplex:::summarize_comparison(old, pi0_method = "none"))
})


test_that("empty result records undefined pi0", {
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)
  s <- rcomplex:::summarize_comparison(cmp[0, ])
  expect_equal(nrow(s$results), 0L)
  expect_named(s$summary$pi0, c("sp1", "sp2"))
  expect_true(all(is.na(s$summary$pi0)))
})


# ---- seed: reproducibility as a property of the call --------------------
# The randomized-p draws behind pi0 used to come from the global RNG with
# no way to pin them from the call site, so q-values -- and every count
# thresholded on them downstream -- moved between runs unless the caller
# remembered set.seed(). `seed` makes that a property of the call.

test_that("summarize_comparison(seed = ) pins the randomized-p q-values", {
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)

  # Guard against a vacuous test: the draw has to actually move pi0 on
  # this fixture, or "different seeds differ" would prove nothing. It
  # does -- pi0 spans about 0.73 to 0.89 over seeds 1:6.
  pi0s <- vapply(1:6, function(s) {
    rcomplex:::summarize_comparison(cmp, seed = s)$summary$pi0[["sp1"]]
  }, numeric(1))
  expect_gt(diff(range(pi0s)), 0.05)

  # same seed, two calls, identical q-values
  a <- rcomplex:::summarize_comparison(cmp, seed = 99)
  expect_identical(rcomplex:::summarize_comparison(cmp, seed = 99), a)

  # a different seed moves pi0, and the q-values with it
  b <- rcomplex:::summarize_comparison(cmp, seed = 100)
  expect_false(identical(a$summary$pi0, b$summary$pi0))
  expect_false(identical(
    a$results$species1.q_value_con, b$results$species1.q_value_con
  ))

  # seed = NULL reproduces the old behaviour exactly: seeding the
  # session first and seeding the call give the same answer.
  set.seed(99)
  expect_identical(rcomplex:::summarize_comparison(cmp), a)
})


test_that("a seed does not change the pi0-free methods' answers", {
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)
  # storey and none draw nothing, so a seed can only pin the stream
  expect_identical(
    rcomplex:::summarize_comparison(cmp, pi0_method = "storey", seed = 3),
    rcomplex:::summarize_comparison(cmp, pi0_method = "storey")
  )
  expect_identical(
    rcomplex:::summarize_comparison(cmp, pi0_method = "none", seed = 3),
    rcomplex:::summarize_comparison(cmp, pi0_method = "none")
  )
})


test_that("a seeded call restores the caller's stream", {
  # Same contract as detect_modules() (test-module-determinism.R): a seed
  # buys a private stream, so a later draw that has no seed of its own
  # continues from the caller's own set.seed() rather than from wherever
  # B rounds of pi0est() happened to leave things -- or from a position
  # this function's seed decided.
  td <- make_graded_nets()
  cmp <- rcomplex:::compare_neighborhoods(td$net1, td$net2, td$ortho)

  restored <- function(f) {
    set.seed(7)
    before <- get(".Random.seed", envir = globalenv())
    f()
    identical(before, get(".Random.seed", envir = globalenv()))
  }
  expect_true(restored(function() rcomplex:::summarize_comparison(cmp, seed = 42)))
  expect_true(restored(function() {
    rcomplex:::summarize_comparison(cmp, pi0_method = "storey", seed = 42)
  }))

  # seed = NULL must still advance the stream, or two unseeded calls in
  # one session would silently share a pi0 draw.
  set.seed(7)
  before <- get(".Random.seed", envir = globalenv())
  invisible(rcomplex:::summarize_comparison(cmp))
  expect_false(identical(before, get(".Random.seed", envir = globalenv())))
})
