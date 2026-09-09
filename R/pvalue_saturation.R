# Resolution diagnostics for permutation p-values.
#
# A permutation p-value is a rational with denominator n_perm + 1, so it can
# only take n_perm + 1 values no matter how strong the signal is. Once a
# batch of tests piles up on the smallest of them the p-value has stopped
# ordering those tests, and anything downstream that reads order rather than
# a threshold -- ranking, weighting, top-k selection, a colour scale on
# -log10(q) -- is reading tie-break noise. This file measures that
# resolution so it is seen rather than assumed.


#' How much resolution is left in a set of p-values
#'
#' Reports how many distinct values a set of p-values (or q-values) actually
#' takes, how many are tied at the smallest one and at 1, and -- when
#' `n_perm` is supplied -- whether that smallest value is the permutation
#' floor `1 / (n_perm + 1)` rather than a statement about the evidence.
#'
#' @details
#' A permutation p-value cannot go below `1 / (n_perm + 1)`, because the
#' observed labelling is itself one of the draws. Tests whose true p-values
#' lie anywhere below that floor all come back holding exactly it. They are
#' then indistinguishable from each other, and so are the q-values derived
#' from them: within a tied block, Benjamini-Hochberg maps one input value to
#' one output value.
#'
#' The distinction the diagnostic draws:
#'
#' \itemize{
#'   \item \strong{permutation-limited} -- the minimum equals the floor. The
#'     tie is an artefact of how many permutations were run, and more
#'     permutations would separate the tied tests.
#'   \item \strong{evidence-limited} -- the minimum sits above the floor.
#'     More permutations would not help; the data are what bind.
#'   \item \strong{below the floor} -- the minimum is smaller than the floor,
#'     so `p` cannot be raw permutation p-values from `n_perm` draws. The
#'     usual cause is passing calibrated or corrected values with the
#'     `n_perm` of the test they came from. Neither statement above holds,
#'     and the function says so rather than picking one.
#' }
#'
#' Either way the consequence for a tied block is the same, and it is the
#' point of this function: **among tied values the p-value carries no
#' ordering at all**, so ranking, weighting or thresholding on `p.value`,
#' `q.value` or `-log10(q)` is arbitrary inside that block. Rank on a
#' continuous effect size instead. For [module_preservation()] that is
#' `Zsummary_std`, which is standardised to unit null variance and keeps
#' separating modules long after the p-value has saturated: on the
#' eight-species Pooideae run, the 35 module-directions tied at the q-value
#' floor spanned `Zsummary_std` from 6.4 to 66.7. Reserve p and q for the
#' significance call.
#'
#' Report the tie count wherever a p-value is reported, so a reader can see
#' the resolution the test actually had.
#'
#' @param p Numeric vector of p-values or q-values in `[0, 1]`. `NA`s are
#'   counted and dropped.
#' @param n_perm Number of permutations behind `p`, if it came from a
#'   permutation test. Supplying it turns on the floor comparison; leave it
#'   `NULL` (default) for corrected or calibrated values, where the observed
#'   minimum is used as the effective floor instead.
#'
#' @return An object of class `"pvalue_resolution"`: a list with
#'   \describe{
#'     \item{n, n_missing}{Non-missing and missing value counts.}
#'     \item{n_distinct}{Number of distinct non-missing values, counted with
#'       the same relative tolerance as `n_at_min` so the two agree.}
#'     \item{min}{The observed minimum.}
#'     \item{n_at_min}{How many values are tied at that minimum.}
#'     \item{n_at_one}{How many values are at or above 1, within
#'       tolerance. Values above 1 are invalid input and are counted
#'       here rather than silently dropped.}
#'     \item{n_perm, floor}{The supplied `n_perm` and the implied floor
#'       `1 / (n_perm + 1)`; both `NULL`/`NA` when `n_perm` was not given.}
#'     \item{floor_status}{`"at"`, `"above"` or `"below"` -- where the
#'       minimum sits relative to the floor. `NA` when `n_perm` was not
#'       given.}
#'     \item{permutation_limited}{`TRUE` when the minimum equals the floor,
#'       `FALSE` when it is above it, `NA` when `n_perm` was not given or
#'       the minimum lies below the floor, where the question is ill-posed
#'       (see `floor_status`).}
#'     \item{n_off_grid}{How many values do not lie on the
#'       `k / (n_perm + 1)` grid, i.e. are not raw permutation p-values.
#'       `NA` when `n_perm` was not given.}
#'     \item{suggested_n_perm}{Permutations needed before the tied block
#'       could take distinct values -- a lower bound, since the true
#'       p-values may lie far further below the floor. `NA` when the minimum
#'       is untied, when the tie is not at the floor, and when `n_perm` was
#'       not given and the observed minimum is too large to be a floor.}
#'   }
#'   with a `print` method.
#'
#' @examples
#' # The eight-species Pooideae run: 511 module-directions from
#' # preservation_paired(), n_perm = 2000. These are the counts it reported.
#' q <- c(
#'   rep(0.00071, 35), rep(1, 20),
#'   rep(seq(0.002, 0.99, length.out = 170), length.out = 456)
#' )
#' pvalue_resolution(q)
#' # 511 values, 172 distinct, 35 tied at 0.00071, 20 at 1: among those 35
#' # the q-value orders nothing, and Zsummary_std spanned 6.4 to 66.7 there.
#'
#' # The same shape for the p-values behind them, where the floor is visible.
#' # The rest sit on the 1/2001 grid, as raw permutation p-values must.
#' p <- c(
#'   rep(1 / 2001, 35), rep(1, 20),
#'   rep(round(seq(4, 1981, length.out = 170)) / 2001, length.out = 456)
#' )
#' pvalue_resolution(p, n_perm = 2000)
#'
#' @seealso [module_preservation()] for `Zsummary_std`, the continuous
#'   effect size to rank on; [tag_permutation()], whose label space puts a
#'   floor under its p-value that no amount of sampling can lower.
#' @export
pvalue_resolution <- function(p, n_perm = NULL) {
  if (!is.numeric(p) || length(p) == 0L) {
    stop("p must be a non-empty numeric vector")
  }

  n_missing <- sum(is.na(p))
  v <- p[!is.na(p)]
  if (length(v) == 0L) {
    stop("p contains no non-missing values")
  }
  if (any(v < 0 | v > 1)) {
    stop("p must lie in [0, 1]")
  }

  if (!is.null(n_perm)) {
    ok <- is.numeric(n_perm) && length(n_perm) == 1L && !is.na(n_perm) &&
      n_perm >= 1 && n_perm == round(n_perm)
    if (!ok) {
      stop("n_perm must be a single positive whole number")
    }
    n_perm <- as.numeric(n_perm) # 1e9 permutations overflow integer
  }

  # 1/(n_perm + 1) is not representable in binary, so a p-value formed in
  # C++ and the same quantity recomputed here can differ in the last bits.
  # Every comparison against the floor is therefore relative, never `==`.
  tol <- sqrt(.Machine$double.eps)

  # One grouping feeds both counts, so a printout cannot claim "6 distinct
  # (100%)" and "4 tied at the minimum" about the same six values.
  u <- sort(unique(v))
  grp <- .tol_groups(u, tol)
  min_obs <- u[1L]
  n_at_min <- sum(v <= max(u[grp == 1L]))
  # >= 1, not a two-sided window: a value above 1 is invalid input and
  # must be reported rather than dropped out of the count. Matches
  # preservation_matrix_test()'s $saturation exactly.
  n_at_one <- sum(v >= 1 - tol)

  floor_p <- NA_real_
  floor_status <- NA_character_
  at_floor <- NA
  n_off_grid <- NA_integer_
  if (!is.null(n_perm)) {
    floor_p <- 1 / (n_perm + 1)
    floor_status <- if (abs(min_obs - floor_p) <= tol * floor_p) {
      "at"
    } else if (min_obs > floor_p) {
      "above"
    } else {
      "below"
    }
    # Below the floor, "permutation-limited" and "evidence-limited" are both
    # false and the flag must not read as either; floor_status carries it.
    at_floor <- switch(floor_status,
      at = TRUE,
      above = FALSE,
      below = NA
    )
    grid <- v * (n_perm + 1)
    n_off_grid <- sum(abs(grid - round(grid)) > tol * pmax(grid, 1))
    # Values under the floor cannot come from n_perm draws. The usual cause
    # is passing calibrated or corrected values together with the n_perm of
    # the test they came from -- the floor then describes the input to the
    # correction, not the vector in hand.
    if (identical(floor_status, "below")) {
      warning(
        "minimum (", signif(min_obs, 3), ") is below the ",
        "permutation floor 1/(n_perm + 1) = ", signif(floor_p, 3),
        "; p does not look like raw permutation p-values"
      )
    }
  }

  # A block of n_at_min values sharing one grid point needs at least that
  # many grid points below it before they can separate. Purely a lower
  # bound: their true p-values may sit orders of magnitude further down.
  suggested <- NA_real_
  if (n_at_min > 1L) {
    eff_floor <- if (is.null(n_perm)) min_obs else floor_p
    usable <- if (is.null(n_perm)) {
      # Without n_perm the observed minimum only stands in for a floor if it
      # is small enough to be one: a minimum of 0.2 would mean a run of four
      # permutations. Above that the tie is in the data and any permutation
      # count would be advice invented out of nothing.
      eff_floor > 0 && 1 / eff_floor - 1 >= .min_credible_n_perm
    } else {
      isTRUE(at_floor)
    }
    if (usable) {
      suggested <- .round_up_nice(n_at_min / eff_floor - 1)
    }
  }

  structure(
    list(
      n                   = length(v),
      n_missing           = n_missing,
      n_distinct          = grp[length(grp)],
      min                 = min_obs,
      n_at_min            = n_at_min,
      n_at_one            = n_at_one,
      n_perm              = n_perm,
      floor               = floor_p,
      floor_status        = floor_status,
      permutation_limited = at_floor,
      n_off_grid          = n_off_grid,
      suggested_n_perm    = suggested
    ),
    class = "pvalue_resolution"
  )
}


# A minimum this large cannot credibly be a permutation floor: 1/(99 + 1)
# is already a coarser run than anyone reports. Used only when n_perm was
# not supplied and the observed minimum has to stand in for the floor.
.min_credible_n_perm <- 99

# The relative tolerance every tie count in this package compares with.
# It lives here, once, because the counts drifted apart when the constant
# was written out at each site: two diagnostics of one vector must not
# disagree about which values are distinct.
.tie_tol <- function() sqrt(.Machine$double.eps)



#' Count values tied at the minimum, to floating-point tolerance
#'
#' Two doubles a few last bits apart are one value wobbled by arithmetic,
#' not two resolvable quantities. Exact `==` against the minimum
#' under-reports a tie, which matters because a tie count is reported
#' precisely to say how far a quantity can rank -- and because the same
#' vector is summarised by more than one diagnostic in this package, which
#' must not disagree about it.
#'
#' @param v Numeric vector; `NA` are dropped.
#' @param tol Relative tolerance.
#' @return List with `n_distinct`, `min` and `n_at_min`.
#' @noRd
.tol_min_ties <- function(v, tol = .tie_tol()) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) {
    return(list(n_distinct = 0L, min = NA_real_, n_at_min = 0L))
  }
  u <- sort(unique(v))
  g <- .tol_groups(u, tol)
  mn <- u[1L]
  list(
    n_distinct = length(unique(g)),
    min = mn,
    n_at_min = sum(v <= mn + tol * abs(mn))
  )
}


# Group sorted unique values that lie within a relative tolerance of the
# group's smallest member. Two doubles a few last bits apart are one value
# wobbled by arithmetic, not two resolvable p-values, and the distinct count
# and the tie count must not disagree about which it is.
.tol_groups <- function(u, tol) {
  n <- length(u)
  if (n <= 1L) {
    return(rep(1L, n))
  }
  if (all(diff(u) > tol * abs(u[-n]))) {
    return(seq_len(n)) # nothing near-equal: skip the loop
  }
  g <- integer(n)
  g[1L] <- 1L
  ref <- u[1L]
  k <- 1L
  for (i in seq.int(2L, n)) {
    if (u[i] > ref + tol * abs(ref)) {
      k <- k + 1L
      ref <- u[i]
    }
    g[i] <- k
  }
  g
}


# Round up to the next 1, 2 or 5 times a power of ten, so a suggested n_perm
# reads as a number someone would actually type rather than 70034.
.round_up_nice <- function(x) {
  if (!is.finite(x) || x <= 0) {
    return(NA_real_)
  }
  cand <- c(1, 2, 5, 10) * 10^floor(log10(x))
  cand[which(cand >= x)[1L]]
}


# Wrap a paragraph into the indented block the print method uses. Written
# out rather than cat(fill = TRUE) because fill ignores the prefix.
.cat_para <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(paste(strwrap(txt, width = 74, prefix = "  "), collapse = "\n"), "\n",
    sep = ""
  )
}


#' @export
print.pvalue_resolution <- function(x, ...) {
  fmt <- function(v) format(signif(v, 3))

  cat("p-value resolution\n")
  cat("  values          ", x$n,
    if (x$n_missing > 0L) paste0(" (", x$n_missing, " missing)") else "",
    "\n",
    sep = ""
  )
  cat("  distinct        ", x$n_distinct, " (",
    round(100 * x$n_distinct / x$n), "% of the values)\n",
    sep = ""
  )
  cat("  minimum         ", fmt(x$min),
    if (x$n_at_min == 1L) {
      " (unique)"
    } else {
      paste0(" (", x$n_at_min, " values tied there)")
    },
    "\n",
    sep = ""
  )
  cat("  at 1            ", x$n_at_one, "\n", sep = "")

  if (!is.null(x$n_perm)) {
    cat("  n_perm          ", format(x$n_perm, scientific = FALSE),
      ", floor 1/(n_perm + 1) = ", fmt(x$floor), "\n",
      sep = ""
    )
    # The object is the deliverable: a result reprinted long after the
    # construction-time warning has scrolled away must still say which of
    # the three positions the minimum holds, not guess between two.
    cat("  status          ",
      switch(x$floor_status,
        at = "permutation-limited: the minimum IS the floor",
        above = "evidence-limited: the minimum is above the floor",
        below = "below the floor: not raw permutation p-values"
      ),
      "\n",
      sep = ""
    )
    if (!is.na(x$n_off_grid) && x$n_off_grid > 0L) {
      cat("  off grid        ", x$n_off_grid,
        " values are not multiples of the floor\n",
        sep = ""
      )
    }
  }

  cat("\n")
  if (x$n_at_min > 1L) {
    .cat_para(
      x$n_at_min, " values are tied at ", fmt(x$min),
      " and carry NO ordering among themselves. Ranking, weighting ",
      "or top-k selection over that block reads tie-break noise, ",
      "not evidence."
    )
    .cat_para(
      "Rank on a continuous effect size instead -- Zsummary_std ",
      "from module_preservation() -- and keep the p-value for the ",
      "significance call alone."
    )
    if (!is.na(x$suggested_n_perm)) {
      .cat_para(
        "Separating the tie needs n_perm >= ",
        format(x$suggested_n_perm, scientific = FALSE),
        ". That is a lower bound: the true p-values may lie far ",
        "further below the floor."
      )
    } else if (identical(x$floor_status, "above")) {
      .cat_para(
        "More permutations will not help: the minimum is above the ",
        "floor, so the tie is in the data, not in the sampling."
      )
    } else if (identical(x$floor_status, "below")) {
      .cat_para(
        "The minimum sits BELOW the floor ", fmt(x$floor), ", so these ",
        "are not raw permutation p-values from ",
        format(x$n_perm, scientific = FALSE), " draws -- corrected or ",
        "calibrated values, most likely. Whether more permutations would ",
        "separate the tie cannot be read off this vector."
      )
    } else if (is.null(x$n_perm)) {
      .cat_para(
        "Whether more permutations would separate the tie cannot be ",
        "told without n_perm, and the minimum (", fmt(x$min), ") is too ",
        "large to stand in for a floor. Supply n_perm to find out."
      )
    }
  } else {
    .cat_para(
      "The minimum is untied, so the values still order the tests. ",
      "Check again after any change to n_perm or to the number of ",
      "tests."
    )
  }
  if (x$n_at_one > 0L) {
    .cat_para(
      x$n_at_one, " values sit at 1. That is a real result ",
      "-- no draw was less extreme than the observed statistic -- ",
      "not a truncation, but those tests are unordered too."
    )
  }

  invisible(x)
}
