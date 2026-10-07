# Fixtures for the package-wide RNG contract test (test-rng-contract.R).
#
# Every case is a self-contained thunk of the form function(seed), so the
# contract assertions can be written once and looped over the table.

#' Every exported entry point that takes a `seed`
#'
#' Generics forward through `...`, so the `seed` formal lives on the
#' `.default` method; both are searched. A new seeded export appears here
#' automatically, and the contract test fails until the table covers it.
rng_seeded_exports <- function() {
  ns <- asNamespace("rcomplex")
  out <- character(0)
  for (nm in sort(getNamespaceExports("rcomplex"))) {
    f <- get(nm, envir = ns)
    if (!is.function(f)) next
    cand <- nm
    default <- paste0(nm, ".default")
    if (exists(default, envir = ns, inherits = FALSE)) {
      cand <- c(cand, default)
    }
    for (cn in cand) {
      g <- get(cn, envir = ns)
      if (is.function(g) && "seed" %in% names(formals(g))) {
        out <- c(out, cn)
      }
    }
  }
  unique(out)
}


#' The contract table: one entry per seeded entry point
#'
#' `call` runs the function at the given seed (possibly `NULL`) and returns
#' something comparable with `identical()`. Everything is sized for speed;
#' the assertions are about the RNG stream, not about the statistics.
#'
#' @param fx The fixture bundle built at the top of test-rng-contract.R.
rng_contract_cases <- function(fx) {
  td <- fx$td
  nets <- fx$nets
  cmp <- fx$cmp
  sparse_nets <- fx$sparse_nets
  list(
    list(
      name = "summarize_comparison",
      internal = TRUE,
      call = function(seed) {
        rcomplex:::summarize_comparison(cmp, seed = seed)$results
      }
    ),
    list(
      name = "find_coexpressologs.default",
      call = function(seed) find_coexpressologs(nets, td$ortho, seed = seed)
    ),
    list(
      name = "density_sweep.default",
      call = function(seed) {
        suppressMessages(density_sweep(nets, td$ortho,
          multipliers = 1, method = "hypergeometric", seed = seed
        ))$edges
      }
    ),
    list(
      name = "permutation_hog_test",
      internal = TRUE,
      call = function(seed) {
        rcomplex:::permutation_hog_test(td$net1, td$net2, cmp,
          min_exceedances = 3L, max_permutations = 40L, seed = seed
        )
      }
    ),
    list(
      name = "coexpressolog_null",
      variant = "pi0 none",
      # n_perm = 2 is far below the 19 permutations p < 0.05 needs, and an
      # unseeded call announces the seed it drew; neither is what the
      # contract is about, and the assertions are unaffected by both.
      call = function(seed) {
        suppressMessages(suppressWarnings(
          coexpressolog_null(sparse_nets, td$ortho,
            n_perm = 2L, swap_factor = 1L, seed = seed,
            pi0_method = "none", pval_combine = "max"
          )
        ))
      }
    ),
    list(
      name = "coexpressolog_null",
      variant = "pi0 randomized",
      # The scope covers the observed run, not just the permutation loop.
      # Under pi0_method = "none" the observed run draws nothing, so that
      # is the only case that would notice the scope sliding back below it.
      call = function(seed) {
        suppressMessages(suppressWarnings(
          coexpressolog_null(sparse_nets, td$ortho,
            n_perm = 2L, swap_factor = 1L, seed = seed,
            pi0_method = "randomized", pval_combine = "max"
          )
        ))
      }
    ),
    list(
      name = "null_network",
      call = function(seed) {
        null_network(fx$null_x, fx$null_net, seed = seed)$network
      }
    ),
    list(
      name = "null_network",
      variant = "block",
      call = function(seed) {
        null_network(fx$null_x, fx$null_net,
          seed = seed,
          block = rep_len(1:3, ncol(fx$null_x))
        )$network
      }
    )
  )
}
