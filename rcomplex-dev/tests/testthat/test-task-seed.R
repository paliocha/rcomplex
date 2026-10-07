# Tests for .task_seed() (R/modules.R): the shared per-task RNG seed used
# by coexpressolog_null() and the module-detection sweeps. The property
# under test is the one the roborev history flagged: a direct
# root + k1 * stream + k2 * index combination is affine in root and index,
# so two roots exactly 40503 apart alias -- the whole permutation vector
# at one root reproduces the other's, shifted by one index. Hashing each
# argument before combining removes that.

test_that(".task_seed() output is always a legal set.seed() input", {
  roots <- c(
    -.Machine$integer.max, -1L, 0L, 1L,
    .Machine$integer.max, 1000L, 1000L + 40503L
  )
  for (root in roots) {
    for (stream in c(1L, 2L, 100L, 137L)) {
      for (index in c(1L, 2L, 40503L, 100000L)) {
        seed <- .task_seed(root, stream, index)
        expect_type(seed, "integer")
        expect_false(is.na(seed))
        expect_true(seed >= 0L && seed <= .Machine$integer.max)
        expect_silent(set.seed(seed))
      }
    }
  }
})


test_that(".task_seed() does not alias two roots 40503 apart", {
  # The historical bug: task_seed(root, s, idx) == task_seed(root + 40503,
  # s, idx - 1) for every idx, so the whole permutation stream at one root
  # was the other's, shifted by one index.
  root1 <- 1000L
  root2 <- root1 + 40503L
  s1 <- vapply(2:20, function(idx) .task_seed(root1, 1L, idx), integer(1))
  s2 <- vapply(1:19, function(idx) .task_seed(root2, 1L, idx), integer(1))
  expect_false(any(s1 == s2))
})


test_that(".task_seed() does not alias at other 40503 multiples", {
  root1 <- -500000L
  for (k in c(1L, 2L, 5L, -3L)) {
    root2 <- root1 + 40503L * k
    s1 <- vapply(1:15, function(idx) .task_seed(root1, 1L, idx), integer(1))
    s2 <- vapply(1:15, function(idx) {
      .task_seed(root2, 1L, idx - k)
    }, integer(1))
    expect_false(all(s1 == s2))
  }
})


test_that(".task_seed() is deterministic and varies with each argument", {
  base <- .task_seed(42L, 1L, 1L)
  expect_identical(.task_seed(42L, 1L, 1L), base)
  expect_false(.task_seed(43L, 1L, 1L) == base)
  expect_false(.task_seed(42L, 2L, 1L) == base)
  expect_false(.task_seed(42L, 1L, 2L) == base)
})


test_that(".task_seed() spreads outputs without gross collisions", {
  seeds <- vapply(seq_len(2000), function(index) {
    .task_seed(7L, 1L, index)
  }, integer(1))
  expect_identical(length(unique(seeds)), length(seeds))
})


test_that(".hash32() never returns NA or a negative value", {
  xs <- c(
    -.Machine$integer.max, -1L, 0L, 1L, .Machine$integer.max,
    2147483647L, 40503L * 1:5
  )
  hashed <- vapply(xs, .hash32, integer(1))
  expect_false(any(is.na(hashed)))
  expect_true(all(hashed >= 0L))
})


test_that(".hash32() does not collapse a root and its negation", {
  # The historical bug: `x %% p` discards sign before anything else runs,
  # so -.Machine$integer.max (== -p) and 0 both reduced to residue 0 and
  # hashed identically, even though they are distinct, individually legal
  # seeds -- set.seed(-.Machine$integer.max) and set.seed(0) are not the
  # same caller seed.
  expect_false(.hash32(-.Machine$integer.max) == .hash32(0L))
  expect_false(.hash32(-.Machine$integer.max) == .hash32(.Machine$integer.max))
})

test_that(".task_seed() does not alias root -.Machine$integer.max with 0", {
  root1 <- -.Machine$integer.max
  root2 <- 0L
  s1 <- vapply(1:10, function(idx) .task_seed(root1, 1L, idx), integer(1))
  s2 <- vapply(1:10, function(idx) .task_seed(root2, 1L, idx), integer(1))
  expect_false(any(s1 == s2))
})
