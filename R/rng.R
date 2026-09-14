#' The package RNG contract (internal)
#'
#' Every exported function in rcomplex that consumes randomness takes a
#' `seed` argument defaulting to `NULL`, and every one of them follows the
#' single rule implemented here: **the caller's stream advances by exactly
#' what the function drew from it, and by nothing else.**
#'
#' * `seed = NULL` -- the function draws from the ambient stream and leaves
#'   it advanced, exactly as [sample()] does. Consecutive unseeded calls
#'   therefore differ, and a `set.seed()` in the caller's own script makes
#'   the whole pipeline reproducible.
#' * `seed = <value>` -- the function draws from a private stream started at
#'   `seed`, and the ambient stream is restored on exit, byte for byte,
#'   including the case where `.Random.seed` did not exist beforehand. A
#'   seeded call is invisible to anything the caller draws afterwards.
#'
#' The rule replaces three contracts that coexisted before 0.3.0: pinning the
#' exit state at `set.seed(seed)` (which silently handed the caller's next
#' unseeded draw a stream determined by *this* function's seed), restoring
#' the ambient stream, and a bare `set.seed(seed)` that displaced the caller
#' by however much the function consumed.
#'
#' Note what the rule is *not*: it does not pin the exit state. Pinning made
#' a downstream unseeded draw -- `summarize_comparison()`'s randomized-p pi0,
#' say -- depend on the upstream function's seed rather than on the caller's,
#' so two pipelines differing only in their top-level `set.seed()` produced
#' identical q-values. Restoring makes that draw depend on the caller's own
#' seed, which is what a caller means by seeding a script.
#'
#' Parallel paths need no special handling: an `mclapply()` fork never
#' propagates `.Random.seed` back to the parent, so a forked worker cannot
#' displace the ambient stream at all. Reproducibility across core counts
#' comes from seeding each task from its own identity (see `.task_seed()`),
#' not from the ambient stream.
#'
#' @param seed A single seed, or `NULL` for the ambient stream.
#' @param envir The frame whose exit restores the stream. Defaults to the
#'   caller, which is what every call site wants.
#' @return `TRUE` invisibly when a private stream was started, `FALSE` when
#'   `seed` was `NULL`.
#' @noRd
.seed_scope <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) {
    return(invisible(FALSE))
  }
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  set.seed(seed)
  # on.exit() has to be registered from the frame it belongs to, so it is
  # called there rather than here. The saved state travels as a literal in
  # the registered call, which keeps the handler independent of this frame.
  do.call(
    base::on.exit,
    list(substitute(.seed_restore(OLD), list(OLD = old)), add = TRUE),
    envir = envir
  )
  invisible(TRUE)
}


#' Put the ambient RNG stream back (internal)
#'
#' `old = NULL` means there was no `.Random.seed` when the scope opened, so
#' the one the seeded call created is removed again. Leaving it behind would
#' turn a fresh session's first unseeded draw from clock-and-PID entropy into
#' a continuation of whatever seed the call used.
#'
#' @noRd
.seed_restore <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old, envir = globalenv())  # nolint
  }
  invisible(NULL)
}

#' Whether `mclapply()` forking is safe to use here (internal)
#'
#' `fork()`-based parallelism and gcov code-coverage instrumentation do not
#' mix: a forked child can inherit a coverage-counter file lock the parent
#' (or another thread) held at fork time and never release it, hanging the
#' process indefinitely rather than erroring. This is a documented `covr`
#' limitation (<https://github.com/r-lib/covr/issues/322>), not a bug in the
#' forked code itself, so the fix is to fall back to the serial path -- not
#' to fork with fewer cores -- whenever a coverage run is detected. `covr`
#' sets `R_COVR = "true"` in the child process it runs tests in, which is
#' the same signal `covr::in_covr()` uses; reading the env var directly
#' avoids taking a hard dependency on `covr` for a package that only
#' Suggests it.
#'
#' @param n_cores Requested core count.
#' @return `TRUE` if `mclapply()` forking is both requested and safe.
#' @noRd
.can_fork <- function(n_cores) {
  .Platform$OS.type == "unix" && n_cores > 1L &&
    !identical(Sys.getenv("R_COVR"), "true")
}
