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
  list(
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
