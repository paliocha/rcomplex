# The package-wide RNG contract, enforced once for every seeded entry
# point rather than re-asserted by hand per function.
#
# The contract (see .seed_scope() in R/rng.R): the caller's stream advances
# by exactly what the function drew from it, and by nothing else. With a
# seed the function draws from a private stream and restores the caller's;
# with seed = NULL it draws from the caller's stream and leaves it
# advanced.
#
# Before 0.3.0 three contracts coexisted: pinning the exit state at
# set.seed(seed), restoring the ambient state, and a bare set.seed(seed)
# with no restore at all. This file is what keeps a fourth from appearing.

ambient <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

# Built once at file scope: every case below reuses it, and the loops run
# each case several times.
rng_fx <- local({
  td <- make_graded_nets()
  nets <- list(SP_A = td$net1, SP_B = td$net2)
  list(
    td = td,
    nets = nets,
    cmp = compare_neighborhoods(td$net1, td$net2, td$ortho),
    sparse_nets = lapply(nets, sparse_net),
    mf = rng_module_fixture(),
    cf = make_clique_fixture(),
    mx = rng_matrix_classification(),
    tf = rng_tag_fixture()
  )
})


test_that("the contract table covers every seeded entry point", {
  # Fails the moment a new exported function grows a seed argument without
  # joining the table, which is what makes the checks below a contract
  # rather than a sample of one.
  covered <- vapply(
    rng_contract_cases(rng_fx), function(x) x$name, character(1)
  )
  expect_setequal(covered, rng_seeded_exports())
})


test_that("every seeded entry point routes its seed through .seed_scope", {
  ns <- asNamespace("rcomplex")
  for (nm in rng_seeded_exports()) {
    body_txt <- paste(deparse(body(get(nm, envir = ns))), collapse = " ")
    expect_true(grepl(".seed_scope(", body_txt, fixed = TRUE),
      label = paste0(nm, "() calls .seed_scope()")
    )
  }
})


test_that("a seeded call restores the caller's ambient stream", {
  skip_on_cran()
  for (case in rng_contract_cases(rng_fx)) {
    set.seed(7L)
    before <- ambient()
    invisible(case$call(42L))
    expect_identical(ambient(), before,
      label = paste0(case$name, "() left the ambient stream where it was")
    )
  }
})


test_that("a seeded call with no ambient stream leaves none behind", {
  skip_on_cran()
  for (case in rng_contract_cases(rng_fx)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(case$call(42L))
    # Leaving one behind would turn the session's next unseeded draw from
    # clock-and-PID entropy into a continuation of this call's seed.
    expect_false(
      exists(".Random.seed", envir = globalenv(), inherits = FALSE),
      label = paste0(case$name, "() created a .Random.seed")
    )
  }
  set.seed(1L)
})


test_that("a seed decides the result, the ambient stream does not", {
  skip_on_cran()
  for (case in rng_contract_cases(rng_fx)) {
    set.seed(1L)
    a <- case$call(42L)
    set.seed(2L)
    b <- case$call(42L)
    expect_equal(a, b,
      label = paste0(case$name, "() at seed 42 from two ambient states")
    )
  }
})


test_that("seed = NULL draws from the ambient stream and advances it", {
  skip_on_cran()
  for (case in rng_contract_cases(rng_fx)) {
    # Advancing is what keeps consecutive unseeded calls from silently
    # sharing one set of draws.
    set.seed(7L)
    before <- ambient()
    invisible(case$call(NULL))
    expect_false(identical(ambient(), before),
      label = paste0(case$name, "(seed = NULL) advanced the stream")
    )

    # And the caller's own set.seed() is what makes it reproducible.
    set.seed(3L)
    a <- case$call(NULL)
    set.seed(3L)
    b <- case$call(NULL)
    expect_equal(a, b,
      label = paste0(case$name, "(seed = NULL) under a caller's set.seed()")
    )
  }
})


test_that(".seed_scope restores through an error", {
  # The restore rides on on.exit(), so a stop() inside the seeded body must
  # not strand the caller on the private stream.
  f <- function(seed = NULL) {
    rcomplex:::.seed_scope(seed)
    stats::runif(3)
    stop("boom")
  }
  set.seed(7L)
  before <- ambient()
  expect_error(f(42L), "boom")
  expect_identical(ambient(), before)
})


test_that(".seed_scope is a no-op at seed = NULL", {
  set.seed(7L)
  before <- ambient()
  expect_false(rcomplex:::.seed_scope(NULL))
  expect_identical(ambient(), before)
})
