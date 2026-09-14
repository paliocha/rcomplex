# coexpressolog_null(): degree-preserving edge-swap null (P6).
#
# Uses the complex_py fixture data (sparse networks, ortholog-restricted
# gene universe) with n_perm = 19 and fixed seeds. Settings passed through
# `...` (pval_combine = "max", pi0_method = "none") are deterministic and
# reproduce the canonical fixture calls in the observed run.

null_fx <- function(...) testthat::test_path("fixtures", "complex_py", ...)

load_null_fixture <- local({
  cache <- NULL
  function() {
    skip_if_not(all(file.exists(null_fx(c(
      "ortho_pairs.tsv", "sp1_expr.tsv",
      "sp2_expr.tsv"
    )))))
    if (!is.null(cache)) {
      return(cache)
    }
    read_expr <- function(file) {
      d <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
      x <- as.matrix(d[, -1])
      rownames(x) <- d$Genes
      x
    }
    ortho <- read.delim(null_fx("ortho_pairs.tsv"), stringsAsFactors = FALSE)
    x1 <- read_expr(null_fx("sp1_expr.tsv"))
    x2 <- read_expr(null_fx("sp2_expr.tsv"))
    x1 <- x1[rownames(x1) %in% ortho$Species1, ]
    x2 <- x2[rownames(x2) %in% ortho$Species2, ]
    networks <- list(
      sp1 = compute_network(x1, density = 0.03, sparse = TRUE),
      sp2 = compute_network(x2, density = 0.03, sparse = TRUE)
    )
    cache <<- list(networks = networks, ortho = ortho)
    cache
  }
})

test_that("coexpressolog_null requires sparse networks", {
  d <- make_cmp_nets()
  nets_dense <- list(A = d$net1, B = d$net2)
  expect_error(
    coexpressolog_null(nets_dense, d$ortho, n_perm = 2L),
    "as_sparse_network"
  )
  # mixed dense/sparse is rejected too
  nets_mixed <- list(A = sparse_net(d$net1), B = d$net2)
  expect_error(
    coexpressolog_null(nets_mixed, d$ortho, n_perm = 2L),
    "as_sparse_network"
  )
})

test_that("the seed is validated up front, and only where it must be", {
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)

  # A seed that will not survive as.integer() has to be caught up front:
  # NA would otherwise reach set.seed(NA), and a length > 1 seed would
  # error on a condition of length > 1, not on the seed. Wrong length,
  # wrong type and wrong magnitude report apart, so no message names a
  # limit the value did not cross. A seed at the integer limit is legal:
  # permutation b derives its seed modulo 2^31 - 1, so nothing overflows.
  run_seed <- function(s) {
    suppressWarnings(coexpressolog_null(nets, d$ortho,
      n_perm = 5L, seed = s,
      pi0_method = "none", pval_combine = "max"
    ))
  }
  expect_s3_class(run_seed(.Machine$integer.max), "data.frame")
  expect_error(run_seed(c(1L, 2L)), "single value; got integer of length 2")
  expect_error(run_seed("seven"), "set.seed\\(\\) accepts; got character")
  expect_error(run_seed(NA), "set.seed\\(\\) accepts")
  expect_error(run_seed(NA_integer_), "set.seed\\(\\) accepts")
  # Coercion decides acceptance, not type: "7" and TRUE are legal seeds
  # everywhere else in the package, so they have to be legal here too.
  expect_equal(run_seed("7"), run_seed(7L))
  expect_equal(run_seed(TRUE), run_seed(1L))
  # Both signs overflow, and neither message may claim the value was
  # merely too large.
  expect_error(run_seed(3e9), "within \\+/- .Machine\\$integer.max")
  expect_error(run_seed(3e9), "got 3e\\+09")
  expect_error(run_seed(-3e9), "within \\+/- .Machine\\$integer.max")
  expect_error(run_seed(-3e9), "got -3e\\+09")
})

test_that("observed conserved calls exceed the rewired null", {
  d <- load_null_fixture()
  # n_perm = 19 is the fixture's standard permutation count; the floor
  # warning it triggers is covered by its own test, not this one
  res <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 19L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  ))

  expect_s3_class(res, "data.frame")
  # one row per species pair plus the total
  expect_identical(res$statistic, c("sp1~sp2", "total"))
  expect_identical(
    names(res),
    c(
      "statistic", "observed", "null_mean", "null_sd", "null_max",
      "fold", "p_emp", "n_ge", "null_se", "p_emp_lo", "p_emp_hi"
    )
  )

  total <- res[res$statistic == "total", ]
  expect_gt(total$observed, total$null_max)
  expect_equal(total$p_emp, 1 / 20)
  # two-species fixture: the pair row equals the total row
  expect_equal(res$observed[1], total$observed)
  expect_equal(res$p_emp[1], 1 / 20)

  null_mat <- attr(res, "null")
  expect_identical(dim(null_mat), c(19L, 2L))
  expect_identical(colnames(null_mat), c("sp1~sp2", "total"))
  expect_equal(res$null_mean, vapply(
    1:2, function(j) mean(null_mat[, j]),
    numeric(1)
  ))
  expect_equal(res$null_max, vapply(
    1:2, function(j) max(null_mat[, j]),
    numeric(1)
  ))
})

test_that("shuffled orthologs give a non-significant null", {
  d <- load_null_fixture()
  ortho_shuf <- d$ortho
  set.seed(99)
  ortho_shuf$Species2 <- sample(ortho_shuf$Species2)
  # n_perm = 19 is the fixture's standard permutation count; the floor
  # warning it triggers is covered by its own test, not this one
  res <- suppressWarnings(coexpressolog_null(
    d$networks, ortho_shuf,
    n_perm = 19L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  ))
  expect_gt(res$p_emp[res$statistic == "total"], 0.05)
})

test_that("n_cores = 2 reproduces the serial result", {
  skip_if(.Platform$OS.type != "unix")
  d <- load_null_fixture()
  # n_perm = 5 is below the 19 the p < 0.05 floor needs; that warning is
  # not what this test is about
  res1 <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 5L, seed = 7L, n_cores = 1L,
    pval_combine = "max", pi0_method = "none"
  ))
  res2 <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 5L, seed = 7L, n_cores = 2L,
    pval_combine = "max", pi0_method = "none"
  ))
  expect_equal(res1, res2)
})


test_that("the effective seed is recorded, and it replays the run", {
  d <- load_null_fixture()
  # n_perm = 19 is the fixture's standard permutation count; the floor
  # warning it triggers is covered by its own test, not this one
  seeded <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 19L, seed = 4L,
    pval_combine = "max", pi0_method = "none"
  ))
  expect_identical(attr(seeded, "seed"), 4L)

  # The reproducibility contract for a default call: the run is not fixed
  # across calls, but the seed it used comes back on the result and
  # reproduces it exactly.
  msgs <- capture_messages(
    unseeded <- suppressWarnings(coexpressolog_null(
      d$networks, d$ortho,
      n_perm = 19L,
      pval_combine = "max", pi0_method = "none"
    ))
  )
  drawn <- attr(unseeded, "seed")
  expect_true(is.numeric(drawn) && length(drawn) == 1L && !is.na(drawn))
  # the message has to name the seed that was actually used: a run that
  # errors before returning leaves nothing else to replay it with
  expect_match(paste(msgs, collapse = ""), paste0("\\b", drawn, "\\b"))

  replay <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 19L, seed = drawn,
    pval_combine = "max", pi0_method = "none"
  ))
  expect_equal(replay, unseeded)
})


test_that("the Monte Carlo columns are what they claim", {
  d <- load_null_fixture()
  # n_perm = 19 is the fixture's standard permutation count; the floor
  # warning it triggers is covered by its own test, not this one
  res <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 19L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  ))
  null_mat <- attr(res, "null")
  n_ge <- vapply(
    seq_along(res$statistic),
    function(j) sum(null_mat[, j] >= res$observed[j]), integer(1)
  )
  expect_identical(res$n_ge, n_ge)
  expect_equal(res$p_emp, (res$n_ge + 1) / 20)
  expect_equal(res$null_se, res$null_sd / sqrt(19))
  # no null draw reaches the observed value here, so the lower endpoint is
  # 0 by definition rather than by qbeta(), whose first shape would be 0
  expect_identical(res$n_ge, c(0L, 0L))
  expect_equal(res$p_emp_lo, c(0, 0))
  expect_equal(res$p_emp_hi, rep(stats::qbeta(0.975, 1, 19), 2))
})



test_that("a p_emp whose interval still covers 0.05 is flagged", {
  d <- load_null_fixture()
  # n_perm = 20 puts p_emp = 1/21 = 0.0476 just under 0.05 while the exact
  # interval on the exceedance probability reaches 0.168, so the call
  # rests on the seed rather than on the data
  expect_warning(
    coexpressolog_null(
      d$networks, d$ortho,
      n_perm = 20L, seed = 1L,
      pval_combine = "max", pi0_method = "none"
    ),
    "cannot be separated from non-significance"
  )
  # at n_perm = 19 the same run gives p_emp = 1/20 = 0.05 exactly, which
  # is not below 0.05 -- the smallest attainable p_emp can never satisfy
  # p < 0.05, so this falls into the "unreachable" warning, not "flagged"
  expect_warning(
    coexpressolog_null(
      d$networks, d$ortho,
      n_perm = 19L, seed = 1L,
      pval_combine = "max", pi0_method = "none"
    ),
    "unreachable for any signal \\(use n_perm >= 20\\)"
  )
})


test_that("neighbouring seeds do not share rewirings", {
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)
  # The conserved count is 0 in nearly every rewiring of these networks,
  # which is far too coarse to tell one rewiring from another. Summing the
  # Jaccard column gives a continuous statistic instead, so two runs
  # agreeing on it means they really did rewire the same way.
  total_jaccard <- function(edges) {
    if (is.null(edges) || nrow(edges) == 0L) {
      return(c(total = 0))
    }
    c(total = sum(edges$jaccard))
  }
  run <- function(s) {
    attr(suppressWarnings(coexpressolog_null(
      nets, d$ortho,
      statistic = total_jaccard, n_perm = 6L, seed = s,
      pval_combine = "max", pi0_method = "none"
    )), "null")
  }
  n1 <- run(1L)
  n2 <- run(2L)
  # set.seed(seed + b) made the null at seed 2 the null at seed 1 shifted
  # by one permutation, so "try another seed" reused all but one rewiring:
  # measured identical() TRUE on this fixture under that scheme, and equal
  # in 0 of 5 shifted rows under the per-task derivation
  expect_false(identical(n1[2:6, ], n2[1:5, ]))
  expect_false(any(n1[2:6, ] == n2[1:5, ]))
})


# Tiny matched-networks fixture for the missing-pair and RNG-state tests:
# a disjoint perfect matching gives every gene exactly one neighbour, so a
# rewired permutation in which no ortholog pair shares a mapped neighbour
# has zero overlap > 0 rows and the species pair is absent from the null
# edges entirely. With seed = 1 the first of 6 permutations is such a run.
make_match_nets <- function() {
  mk <- function(prefix) {
    n <- 8L
    g <- paste0(prefix, seq_len(n))
    m <- matrix(0, n, n, dimnames = list(g, g))
    for (e in list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))) {
      m[e[1], e[2]] <- m[e[2], e[1]] <- 10
    }
    list(network = m, threshold = 5)
  }
  list(
    networks = list(A = sparse_net(mk("A")), B = sparse_net(mk("B"))), # nolint
    ortho = data.frame(
      Species1 = paste0("A", 1:8),
      Species2 = paste0("B", 1:8),
      hog = paste0("H", 1:8),
      stringsAsFactors = FALSE
    )
  )
}


test_that("the interval endpoints hold when every draw exceeds", {
  d <- make_match_nets()
  # matched networks give 0 conserved calls observed and in every
  # permutation, so n_ge == n_perm and the upper endpoint is 1
  res <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 6L, seed = 1L,
    pval_combine = "max", pi0_method = "none"
  ))
  expect_identical(res$n_ge, c(6L, 6L))
  expect_equal(res$p_emp, rep(1, 2))
  expect_equal(res$p_emp_hi, rep(1, 2))
  expect_equal(res$p_emp_lo, rep(stats::qbeta(0.025, 6, 1), 2))
})


test_that("n_perm below 19 warns that p < 0.05 is unreachable", {
  d <- make_match_nets()
  expect_warning(
    coexpressolog_null(
      d$networks, d$ortho,
      n_perm = 6L, seed = 1L,
      pval_combine = "max", pi0_method = "none"
    ),
    "smallest attainable p-value of 0.143"
  )
})


test_that(
  "null runs missing a species pair record 0 for the built-in statistic",
  {
    d <- make_match_nets()
    res <- suppressWarnings(coexpressolog_null(
      d$networks, d$ortho,
      n_perm = 6L, seed = 1L,
      pval_combine = "max", pi0_method = "none"
    ))
    expect_identical(res$statistic, c("A~B", "total"))
    null_mat <- attr(res, "null")
    expect_true(all(is.finite(null_mat)))
    # the pair-less rewiring is a null observation of 0, not an abort
    expect_equal(unname(null_mat[1L, "A~B"]), 0)
  }
)


test_that("a user statistic missing a name still errors", {
  d <- make_match_nets()
  per_pair <- function(edges) {
    if (is.null(edges) || nrow(edges) == 0L) {
      return(c(total = 0))
    }
    pair <- paste(edges$species1, edges$species2, sep = "~")
    counts <- vapply(
      split(edges$type == "conserved", pair), sum,
      numeric(1)
    )
    c(counts, total = sum(counts))
  }
  expect_error(
    coexpressolog_null(d$networks, d$ortho,
      statistic = per_pair,
      n_perm = 6L, seed = 1L,
      pval_combine = "max", pi0_method = "none"
    ),
    "statistic is missing"
  )
})


test_that("the serial path restores the caller's RNG state", {
  d <- make_match_nets()
  set.seed(11)
  before <- runif(5)
  set.seed(11)
  invisible(suppressWarnings(coexpressolog_null(
    d$networks, d$ortho,
    n_perm = 2L, seed = 3L,
    pval_combine = "max", pi0_method = "none"
  )))
  after <- runif(5)
  # the per-permutation set.seed() inside the serial loop must not leak:
  # the caller's stream continues exactly where the observed run left it
  expect_equal(after, before)
})


test_that("an NA null statistic reports NA instead of aborting the run", {
  # A user statistic may be undefined on some permutation -- mean() over a
  # rewiring that called nothing conserved is NaN. The Monte Carlo
  # diagnostics must carry that through as NA; before they were NA-safe,
  # `if (any(weak))` saw NA and the whole run died with
  # "missing value where TRUE/FALSE needed".
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)
  stat_na <- function(edges) {
    keep <- edges$type == "conserved"
    c(total = mean(edges$jaccard[keep]))
  }
  res <- expect_silent(
    suppressWarnings(
      coexpressolog_null(nets, d$ortho,
        n_perm = 19L, statistic = stat_na, seed = 11L
      )
    )
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("n_ge", "null_se", "p_emp_lo", "p_emp_hi") %in%
                    names(res)))
  # This fixture must actually exercise the NA path -- otherwise the
  # assertions below pass vacuously and the regression they guard against
  # (NA reaching `if (any(weak))`) goes unexercised.
  na_rows <- is.na(res$n_ge)
  expect_true(any(na_rows))
  # An NA row must stay NA on every derived column rather than being
  # scored as significant.
  expect_true(all(is.na(res$p_emp[na_rows])))
  expect_true(all(is.na(res$p_emp_lo[na_rows])))
  expect_true(all(is.na(res$p_emp_hi[na_rows])))
})
