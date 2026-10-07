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


#' Check the results of an `mclapply()` call
#'
#' `mclapply()` returns a `try-error` for a task that failed and `NULL` for
#' one whose worker died without reporting (e.g. a segfault), with only a
#' warning. Unchecked, either surfaces later as an unrelated error or as a
#' silently shortened result. Every fork site routes its results through
#' here.
#'
#' @param res List returned by `mclapply()` (or `lapply()`), one element
#'   per task.
#' @param labels Label of each task for the message (a resolution, a
#'   permutation index).
#' @param what What the tasks are, e.g. `"permutation"`.
#' @return `res`, unchanged, when every task returned a result. A failed
#'   task raises an error naming it, with the worker's condition as its
#'   parent and the worker's classes kept.
#' @noRd
.check_fork_results <- function(res, labels, what) {
  errs <- which(vapply(res, inherits, logical(1), "try-error"))
  if (length(errs)) {
    e <- res[[errs[1L]]]
    head <- paste0(what, " ", labels[errs[1L]], " failed")
    cnd <- attr(e, "condition")
    if (is.null(cnd)) {
      rlang::abort(paste0(head, ": ", trimws(as.character(e))), call = NULL)
    }
    # The worker's condition becomes the parent of the new error, so rlang
    # renders its message once whatever its shape (bullets, several lines,
    # a cnd_header() method) and records the backtrace; its own classes are
    # kept on the new error, so class-based handlers match on either path.
    rlang::abort(head,
      class = setdiff(
        class(cnd), c("rlang_error", "error", "condition", "simpleError")
      ),
      parent = cnd, call = NULL
    )
  }
  failed <- vapply(res, is.null, logical(1))
  if (any(failed)) {
    stop(
      "no result for ", what, " ", paste(labels[failed], collapse = ", "),
      " (on a forked run a worker may have crashed; rerun with n_cores = 1)",
      call. = FALSE
    )
  }
  res
}


#' Avalanche-mix a value into [0, 2^31 - 2] (internal)
#'
#' A Thomas Wang-style integer hash, reimplemented with `%%` so every
#' `bitwXor`/`bitwShiftR` argument stays inside 2^31 - 1 (never exceeding the
#' 32-bit signed range). The two internal multiplications do exceed 2^53 in
#' double precision and so round before the `%%`; this is a source of
#' rounding, not overflow (`x %% p` after the multiply always lands back in
#' `[0, p)`), and the avalanche property below does not depend on those
#' products being exact. Two rounds of xor-shift plus modular
#' multiplication are what make this a permutation with no useful linear
#' structure left in it: unlike a bare `root + k * index` term, there is no
#' fixed offset `d` with `.hash32(x + d) - .hash32(x)` constant across `x`,
#' which is the property `.task_seed()` needs from it.
#'
#' `x %% p` on a signed input discards the sign before anything else runs:
#' `-p %% p == 0 == 0 %% p`, so `.hash32(-.Machine$integer.max)` and
#' `.hash32(0)` used to be identical even though `.Machine$integer.max` is
#' the documented other end of the accepted seed domain. The sign is
#' folded in explicitly here (as a high bit XORed into the reduced
#' magnitude) before the avalanche runs, so a negative and a non-negative
#' input with the same magnitude no longer collide at this step. The
#' domain (`+/- .Machine$integer.max`, `2p + 1` values) is still larger
#' than the `p`-periodic range this folds into, so some collisions
#' further apart remain possible in principle; this closes the specific,
#' cheap-to-reach `x` vs `-x` case, not the pigeonhole bound.
#'
#' @noRd
.hash32 <- function(x) {
  p <- 2147483647
  sign_bit <- if (x < 0) 1073741824L else 0L # 2^30, well inside int32 range
  x <- as.integer(abs(as.numeric(x)) %% p)
  x <- bitwXor(x, sign_bit)
  x <- bitwXor(x, bitwShiftR(x, 15L))
  x <- as.integer((as.numeric(x) * 2246822519) %% p)
  x <- bitwXor(x, bitwShiftR(x, 13L))
  x <- as.integer((as.numeric(x) * 3266489917) %% p)
  bitwXor(x, bitwShiftR(x, 16L))
}

#' Deterministic per-task RNG seed (internal)
#'
#' Maps (root, stream, index) to a legal R seed, so every parallel task seeds
#' itself from its own identity rather than inheriting one. That is what makes
#' the result independent of how mclapply chunks the work: with the default
#' RNGkind a forked child runs parallel:::mc.set.stream(), whose non-L'Ecuyer
#' branch deletes .Random.seed, and the child then re-seeds from clock and PID
#' at its first draw. igraph::cluster_leiden() consumes the R stream, so
#' set.seed() in the parent reaches no worker at n_cores > 1.
#'
#' `root`, `stream` and `index` are each run through `.hash32()` before being
#' combined, rather than combined directly with fixed multipliers. A direct
#' `root + 2654435761 * stream + 40503 * index` form is affine in `root` and
#' `index`, so any two roots exactly 40503 apart (mod 2^31 - 1) alias: the
#' whole permutation vector at root `r + 40503` reproduces the one at root
#' `r`, shifted by one index -- the same failure `.task_seed()` was written
#' to remove from `seed_root + b`, just at a rarer, silent distance. Hashing
#' each argument first breaks that: two hashed roots differing by a fixed
#' amount no longer differ by that same amount at every index, so no root
#' pair can alias more than a single, coincidental index.
#'
#' The arithmetic is done in doubles and folded modulo 2^31 - 1 so it cannot
#' overflow integer range: a naive root + offset form gives NA for a seed near
#' the limit, and set.seed(NA) errors. Products stay exact in double
#' (2654435761 * 1e6 is well under 2^53).
#'
#' `.hash32(root)` used to fold `root %% p` directly before hashing, which
#' is periodic in `root` with period `p = 2^31 - 1` and, worse, collapses
#' sign: `-p %% p == 0 == p %% p == 0 %% p`, so root `-.Machine$integer.max`
#' (a legal seed callers such as `coexpressolog_null()` accept) hashed
#' identically to root `0`. `.hash32()` now folds the sign in separately
#' (see its own docs) so that specific collision is closed; the domain is
#' still larger than the range this folds into (`2p + 1` legal roots
#' against `p` residues), so some pair of roots further apart than `p`
#' can in principle still alias -- the risk this leaves is a coincidental
#' one at cryptographically-unlikely distances, not the cheap, guaranteed
#' `x` vs `-x` case that existed before.
#'
#' Streams: 1 = initial sweep, 2 = K=1 permutations, 100 + iter = consensus
#' sweep at that iteration. The consensus stream must vary with iter because
#' the same resolutions are re-run on a different graph each round.
#'
#' @noRd
.task_seed <- function(root, stream, index) {
  p <- 2147483647
  r <- .hash32(root)
  s <- .hash32((as.numeric(stream) * 2654435761) %% p)
  i <- .hash32((as.numeric(index) * 40503) %% p)
  as.integer((as.numeric(r) + as.numeric(s) + as.numeric(i)) %% p)
}
