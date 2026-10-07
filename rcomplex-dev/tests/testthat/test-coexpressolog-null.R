# coexpressolog_null(): degree-preserving edge-swap null (P6).
#
# Uses the complex_py fixture data (sparse networks, ortholog-restricted
# gene universe) with 100 permutations and fixed seeds.

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
    names(ortho)[1:2] <- c("gene1", "gene2")  # fixture predates gene1/gene2
    x1 <- read_expr(null_fx("sp1_expr.tsv"))
    x2 <- read_expr(null_fx("sp2_expr.tsv"))
    x1 <- x1[rownames(x1) %in% ortho$gene1, ]
    x2 <- x2[rownames(x2) %in% ortho$gene2, ]
    networks <- list(
      species1 = compute_network(x1, density = 0.03, sparse = TRUE),
      species2 = compute_network(x2, density = 0.03, sparse = TRUE)
    )
    cache <<- list(networks = networks, ortho = ortho)
    cache
  }
})

test_that("coexpressolog_null requires sparse networks", {
  d <- make_cmp_nets()
  nets_dense <- list(A = d$net1, B = d$net2)
  expect_error(
    coexpressolog_null(nets_dense, d$ortho),
    "as_sparse_network"
  )
  # mixed dense/sparse is rejected too
  nets_mixed <- list(A = sparse_net(d$net1), B = d$net2)
  expect_error(
    coexpressolog_null(nets_mixed, d$ortho),
    "as_sparse_network"
  )
})

test_that("swap_factor must be a finite positive number", {
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)
  # Each of these used to rewire nothing and return the observed graph as
  # its own null, silently.
  for (sf in list(NA_real_, NaN, -5, 0, Inf, c(1, 2), "10")) {
    expect_error(
      coexpressolog_null(nets, d$ortho, swap_factor = sf,
                         seed = 1L),
      "swap_factor must be a single finite number > 0"
    )
  }
})

test_that("every network is validated", {
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)
  wide <- nets$B
  wide$network <- cbind(wide$network, wide$network[, 1:5])
  # C never enters the observed run, so find_coexpressologs() never checks
  # it, but the null still rewires it: a wide matrix would index past the
  # kernel's nrow x nrow bit matrix.
  expect_error(
    coexpressolog_null(c(nets, list(C = wide)), d$ortho, seed = 1L),
    "square"
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
      seed = s, pval_combine = "max"
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
  res <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 1L,
    pval_combine = "max"
  ))

  expect_s3_class(res, "data.frame")
  # one row per species pair plus the total
  expect_identical(res$statistic, c("species1~species2", "total"))
  expect_identical(
    names(res),
    c(
      "statistic", "observed", "null_mean", "null_sd", "null_max",
      "fold", "p_emp", "n_ge", "null_se", "p_emp_lo", "p_emp_hi"
    )
  )

  total <- res[res$statistic == "total", ]
  expect_gt(total$observed, total$null_max)
  expect_equal(total$p_emp, 1 / 101)
  # two-species fixture: the pair row equals the total row
  expect_equal(res$observed[1], total$observed)
  expect_equal(res$p_emp[1], 1 / 101)

  null_mat <- attr(res, "null")
  expect_identical(dim(null_mat), c(100L, 2L))
  expect_identical(colnames(null_mat), c("species1~species2", "total"))
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
  ortho_shuf$gene2 <- sample(ortho_shuf$gene2)
  res <- suppressWarnings(coexpressolog_null(
    d$networks, ortho_shuf, seed = 1L,
    pval_combine = "max"
  ))
  expect_gt(res$p_emp[res$statistic == "total"], 0.05)
})

test_that("n_cores = 2 reproduces the serial result", {
  skip_if(.Platform$OS.type != "unix")
  d <- load_null_fixture()
  res1 <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 7L, n_cores = 1L,
    pval_combine = "max"
  ))
  res2 <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 7L, n_cores = 2L,
    pval_combine = "max"
  ))
  expect_equal(res1, res2)
})


test_that("the effective seed is recorded, and it replays the run", {
  d <- load_null_fixture()
  seeded <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 4L,
    pval_combine = "max"
  ))
  expect_identical(attr(seeded, "seed"), 4L)

  # The reproducibility contract for a default call: the run is not fixed
  # across calls, but the seed it used comes back on the result and
  # reproduces it exactly.
  msgs <- capture_messages(
    unseeded <- suppressWarnings(coexpressolog_null(
      d$networks, d$ortho,
      pval_combine = "max"
    ))
  )
  drawn <- attr(unseeded, "seed")
  expect_true(is.numeric(drawn) && length(drawn) == 1L && !is.na(drawn))
  # the message has to name the seed that was actually used: a run that
  # errors before returning leaves nothing else to replay it with
  expect_match(paste(msgs, collapse = ""), paste0("\\b", drawn, "\\b"))

  replay <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = drawn,
    pval_combine = "max"
  ))
  expect_equal(replay, unseeded)
})


test_that("the Monte Carlo columns are what they claim", {
  d <- load_null_fixture()
  res <- suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 1L,
    pval_combine = "max"
  ))
  null_mat <- attr(res, "null")
  n_ge <- vapply(
    seq_along(res$statistic),
    function(j) sum(null_mat[, j] >= res$observed[j]), integer(1)
  )
  expect_identical(res$n_ge, n_ge)
  expect_equal(res$p_emp, (res$n_ge + 1) / 101)
  expect_equal(res$null_se, res$null_sd / sqrt(100))
  # no null draw reaches the observed value here, so the lower endpoint is
  # 0 by definition rather than by qbeta(), whose first shape would be 0
  expect_identical(res$n_ge, c(0L, 0L))
  expect_equal(res$p_emp_lo, c(0, 0))
  expect_equal(res$p_emp_hi, rep(stats::qbeta(0.975, 1, 100), 2))
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
      gene1 = paste0("A", 1:8),
      gene2 = paste0("B", 1:8),
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
    d$networks, d$ortho, seed = 1L,
    pval_combine = "max"
  ))
  expect_identical(res$n_ge, c(100L, 100L))
  expect_equal(res$p_emp, rep(1, 2))
  expect_equal(res$p_emp_hi, rep(1, 2))
  expect_equal(res$p_emp_lo, rep(stats::qbeta(0.025, 100, 1), 2))
})


test_that(
  "null runs missing a species pair record 0 for the built-in statistic",
  {
    d <- make_match_nets()
    res <- suppressWarnings(coexpressolog_null(
      d$networks, d$ortho, seed = 1L,
      pval_combine = "max"
    ))
    expect_identical(res$statistic, c("A~B", "total"))
    null_mat <- attr(res, "null")
    expect_true(all(is.finite(null_mat)))
    # the pair-less rewiring is a null observation of 0, not an abort
    expect_equal(unname(null_mat[1L, "A~B"]), 0)
  }
)


test_that("the serial path restores the caller's RNG state", {
  d <- make_match_nets()
  set.seed(11)
  before <- runif(5)
  set.seed(11)
  invisible(suppressWarnings(coexpressolog_null(
    d$networks, d$ortho, seed = 3L,
    pval_combine = "max"
  )))
  after <- runif(5)
  # the per-permutation set.seed() inside the serial loop must not leak:
  # the caller's stream continues exactly where the observed run left it
  expect_equal(after, before)
})


# .rewire_degseq(): the C++ degree-preserving swap kernel behind
# coexpressolog_null(). It replaced igraph::rewire(keeping_degseq()), so
# these tests pin the properties that made igraph's chain a valid null:
# degrees kept, simple output, and uniform sampling over the realizations
# of a degree sequence.

# Symmetric binary dgCMatrix with both triangles stored, from 1-based
# endpoint vectors.
sym_adj <- function(from, to, n, names = NULL) {
  Matrix::sparseMatrix(
    i = c(from, to), j = c(to, from), x = 1, dims = c(n, n),
    dimnames = if (is.null(names)) NULL else list(names, names)
  )
}

random_adj <- function(n, p) {
  ut <- which(upper.tri(diag(n)), arr.ind = TRUE)
  keep <- stats::runif(nrow(ut)) < p
  sym_adj(ut[keep, 1], ut[keep, 2], n)
}

edge_keys <- function(a) {
  s <- Matrix::summary(a)
  s <- s[s$i < s$j, ]
  sort(paste(s$i, s$j))
}


test_that(".rewire_degseq keeps every degree and returns a simple graph", {
  set.seed(1)
  n <- 200L
  a <- random_adj(n, 0.05)
  dimnames(a) <- list(paste0("g", seq_len(n)), paste0("g", seq_len(n)))
  r <- .rewire_degseq(a, swap_factor = 10)

  expect_s4_class(r, "dgCMatrix")
  expect_identical(dimnames(r), dimnames(a))
  expect_equal(Matrix::rowSums(r), Matrix::rowSums(a))
  expect_true(Matrix::isSymmetric(r))
  expect_true(all(Matrix::diag(r) == 0))
  # sparseMatrix() sums duplicate entries, so a multi-edge would read x = 2
  expect_true(all(r@x == 1))
  expect_false(identical(edge_keys(r), edge_keys(a)))
})


test_that(".rewire_degseq returns graphs without a legal swap unchanged", {
  set.seed(2)
  a <- random_adj(50L, 0.1)
  expect_identical(edge_keys(.rewire_degseq(a, swap_factor = 0)),
                   edge_keys(a))
  # one edge: nothing to pair it with
  one <- sym_adj(1, 2, 4L)
  expect_identical(edge_keys(.rewire_degseq(one, 10)), edge_keys(one))
  # star: every pair of edges shares the hub, so each trial is a no-op or
  # would create a loop
  star <- sym_adj(rep(1, 5), 2:6, 6L)
  expect_identical(edge_keys(.rewire_degseq(star, 10)), edge_keys(star))
  # complete graph: every swap would create a multi-edge
  k5 <- which(upper.tri(diag(5)), arr.ind = TRUE)
  full <- sym_adj(k5[, 1], k5[, 2], 5L)
  expect_identical(edge_keys(.rewire_degseq(full, 10)), edge_keys(full))
  # no edges at all
  empty <- Matrix::drop0(sym_adj(1, 2, 3L) * 0)
  expect_identical(Matrix::nnzero(.rewire_degseq(empty, 10)), 0L)
})


test_that(".rewire_degseq draws from R's RNG stream", {
  set.seed(3)
  a <- random_adj(100L, 0.08)
  # coexpressolog_null() relies on set.seed(.task_seed()) fixing each
  # rewiring, which only holds if the kernel draws from R's stream
  set.seed(10)
  r1 <- .rewire_degseq(a, 5)
  after1 <- stats::runif(1)
  set.seed(10)
  r2 <- .rewire_degseq(a, 5)
  after2 <- stats::runif(1)
  set.seed(11)
  r3 <- .rewire_degseq(a, 5)
  expect_identical(r1, r2)
  expect_identical(after1, after2)
  expect_false(identical(edge_keys(r1), edge_keys(r3)))
})


test_that(".rewire_degseq mixes away from the original edges", {
  set.seed(4)
  a <- random_adj(400L, 0.05)
  orig <- edge_keys(a)
  retained <- function(sf) {
    mean(edge_keys(.rewire_degseq(a, sf)) %in% orig)
  }
  expect_identical(retained(0), 1)
  # a uniformly random graph with these degrees keeps about one edge in
  # twenty; 0.2 leaves room for degree heterogeneity, not for a chain
  # that stalls near its start
  expect_lt(retained(10), 0.2)
})


test_that(".rewire_degseq samples degree-sequence realizations uniformly", {
  # Every simple graph on 6 labelled nodes with the start graph's degrees,
  # enumerated over all 2^15 edge subsets of K6. Counting rejected trials
  # is what makes the swap chain's stationary distribution uniform; a
  # kernel that retried on rejection, or never flipped the second edge,
  # would fail this.
  ut <- which(upper.tri(diag(6)), arr.ind = TRUE)
  # the 6-cycle 1-2-3-4-5-6-1, as rows of `ut`
  start <- c(1, 3, 6, 10, 15, 11)
  deg <- tabulate(c(ut[start, 1], ut[start, 2]), 6)

  subsets <- 0:(2^15 - 1)
  bits <- vapply(0:14, function(k) bitwAnd(subsets, 2^k) > 0, logical(2^15))
  inc <- matrix(0, 15, 6)
  inc[cbind(1:15, ut[, 1])] <- 1
  inc[cbind(1:15, ut[, 2])] <- 1
  node_deg <- bits %*% inc
  valid <- subsets[apply(node_deg, 1, function(x) all(x == deg))]
  # 2-regular on 6 labelled nodes: 60 hexagons plus 10 pairs of triangles
  expect_identical(length(valid), 70L)

  a <- sym_adj(ut[start, 1], ut[start, 2], 6L)
  pair_id <- (ut[, 2] - 1) * 6 + ut[, 1]
  set.seed(5)
  n_draw <- 3000L
  keys <- vapply(seq_len(n_draw), function(i) {
    s <- Matrix::summary(.rewire_degseq(a, 50))
    s <- s[s$i < s$j, ]
    sum(2^(match((s$j - 1) * 6 + s$i, pair_id) - 1))
  }, numeric(1))

  expect_true(all(keys %in% valid))
  counts <- table(factor(keys, levels = valid))
  expect_true(all(counts > 0))
  expect_gt(stats::chisq.test(as.vector(counts))$p.value, 1e-3)
})


test_that("the rewiring kernel rejects swap factors it cannot represent", {
  a <- sym_adj(c(1, 3), c(2, 4), 4L)
  for (sf in c(NA_real_, NaN, -1, Inf, 1e300)) {
    expect_error(rewire_degseq_cpp(a@p, a@i, a@x, sf),
                 "swap_factor must be finite and non-negative")
  }
})


test_that("the rewiring kernel validates the dgCMatrix slots it reads", {
  a <- sym_adj(c(1, 3, 5), c(2, 4, 6), 6L)
  bad <- a
  bad@i[1L] <- 100L # row index outside 6 x 6; slot assignment skips validity
  expect_error(.rewire_degseq(bad, 10), "row indices must be strictly")
  # ... and through coexpressolog_null(), for a network that .net_check()
  # alone lets through
  d <- make_cmp_nets()
  nets <- lapply(list(A = d$net1, B = d$net2), sparse_net)
  broken <- nets$B
  broken$network@i[which(broken$network@i > 0L)[1L]] <- 100000L
  expect_error(
    coexpressolog_null(c(nets, list(C = broken)), d$ortho, seed = 1L),
    "row indices must be strictly"
  )
})


test_that("a small swap_factor still makes at least one trial", {
  # 4 disjoint edges: any swap of two of them is legal, so one trial always
  # changes the graph. floor(0.1 * 4) = 0 trials used to return it intact.
  a <- sym_adj(c(1, 3, 5, 7), c(2, 4, 6, 8), 8L)
  set.seed(1)
  expect_false(identical(edge_keys(.rewire_degseq(a, 0.1)), edge_keys(a)))
})
