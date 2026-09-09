#' All unordered species pairs as a `pairs` data frame
#'
#' Builds the `pairs` argument of [preservation_paired()] for every
#' `choose(n, 2)` contrast, so an all-pairs preservation matrix does not have
#' to be hand-written. Running all pairs rather than a designated few is what
#' makes [preservation_matrix_test()] worth doing: every extra contrast adds
#' its modules to the class means the statistic is built from, and any
#' species the designated few left out widens the relabelling null as well.
#' Eight species split 4/4 give `choose(8, 4) = 70` free labellings, against
#' the 16 that a within-genus null over four genera of two can reach.
#'
#' @param species Character vector of species identifiers. A *named* vector
#'   --- a `group` map, say --- is read as its names, since a bare species
#'   vector is not normally named.
#' @param sep Separator used to build `pair_name` (default `"."`), matching
#'   the `paste(sp1, sp2, sep = ".")` default of [preservation_paired()], so
#'   the two agree when `pair_name` is left implicit.
#'
#' @return A data frame with columns `sp1`, `sp2` and `pair_name`, one row per
#'   unordered pair. Both directions of each contrast are run by
#'   [preservation_paired()], so each pair appears once.
#'
#' @examples
#' all_species_pairs(c("BDIS", "BSYL", "HVUL"))
#'
#' @export
all_species_pairs <- function(species, sep = ".") {
  if (!is.null(names(species))) {
    species <- names(species)
  }
  species <- as.character(species)
  if (length(species) < 2L) {
    stop("species must name at least two species")
  }
  if (anyNA(species) || any(!nzchar(species))) {
    stop("species must not contain NA or empty identifiers")
  }
  if (anyDuplicated(species) > 0L) {
    stop(
      "species must be unique; repeated: ",
      paste(unique(species[duplicated(species)]), collapse = ", ")
    )
  }
  idx <- utils::combn(length(species), 2L)
  data.frame(
    sp1 = species[idx[1L, ]],
    sp2 = species[idx[2L, ]],
    pair_name = paste(species[idx[1L, ]], species[idx[2L, ]], sep = sep),
    stringsAsFactors = FALSE
  )
}


#' Distinct permutations of a label multiset (internal)
#'
#' Enumerates each arrangement once rather than the `n!` orderings of
#' indistinguishable labels: fixing each *distinct* value in turn as the head
#' and recursing on the remainder visits the multinomial space exactly.
#' Enumerating `n!` and deduplicating would be 576 times the work at eight
#' species split 4/4.
#'
#' @noRd
.pmt_perms <- function(x) {
  n <- length(x)
  if (n <= 1L) {
    return(list(x))
  }
  out <- list()
  for (v in sort(unique(x))) {
    i <- match(v, x)
    for (r in .pmt_perms(x[-i])) {
      out[[length(out) + 1L]] <- c(v, r)
    }
  }
  out
}


#' Size of a relabelling space (internal)
#'
#' The multinomial coefficient of the label counts, or its product over blocks
#' for the restricted null. Computed through `lgamma` and returned as a double
#' because the free space overflows integer range well before it stops being
#' worth counting.
#'
#' @noRd
.pmt_space <- function(labels, blocks) {
  one <- function(z) {
    cnt <- tabulate(z)
    exp(lgamma(length(z) + 1) - sum(lgamma(cnt + 1)))
  }
  if (is.null(blocks)) {
    return(round(one(labels)))
  }
  round(prod(vapply(split(labels, blocks), one, numeric(1))))
}


#' Enumerate a relabelling space as a matrix (internal)
#'
#' Rows are labellings, columns species, in the column order of `labels`. The
#' restricted null is the product over blocks, so it is built by walking a
#' grid of per-block permutation indices rather than by filtering the free
#' space, which would be astronomically larger.
#'
#' @noRd
.pmt_enumerate <- function(labels, blocks) {
  if (is.null(blocks)) {
    return(do.call(rbind, .pmt_perms(labels)))
  }
  idx_by_block <- split(seq_along(labels), blocks)
  per <- lapply(idx_by_block, function(idx) {
    do.call(
      rbind,
      .pmt_perms(labels[idx])
    )
  })
  grid <- expand.grid(lapply(per, function(m) seq_len(nrow(m))),
    KEEP.OUT.ATTRS = FALSE
  )
  out <- matrix(NA_integer_, nrow(grid), length(labels))
  for (b in seq_along(per)) {
    out[, idx_by_block[[b]]] <- per[[b]][grid[[b]], , drop = FALSE]
  }
  out
}


#' One sampled labelling (internal)
#' @noRd
.pmt_draw <- function(labels, blocks) {
  if (is.null(blocks)) {
    return(labels[sample.int(length(labels))])
  }
  out <- labels
  for (idx in split(seq_along(labels), blocks)) {
    out[idx] <- labels[idx][sample.int(length(idx))]
  }
  out
}


#' Class code of each module-direction (internal)
#'
#' `0` for a concordant row, whatever level the two species share, and
#' otherwise an id unique to the unordered pair of distinct levels. Pooling
#' the concordant classes is what makes the two-level case reduce exactly to
#' `mean(concordant) - mean(discordant)`; splitting them would give three
#' classes for a binary trait and no such reduction.
#'
#' The statistic and the reported `$class_means` both partition through this
#' one definition. Two copies of it would let the reported classes go on
#' looking right while the statistic pooled differently.
#'
#' @noRd
.pmt_class_code <- function(a, b, n_lev) {
  lo <- pmin(a, b)
  hi <- pmax(a, b)
  ifelse(lo == hi, 0L, lo + hi * n_lev)
}


#' Readable name of a class code (internal)
#'
#' Inverts `.pmt_class_code()`: `lo = code %% n_lev`, `hi = code %/% n_lev`,
#' which is exact because `lo` is at most `n_lev - 1` whenever `lo < hi`.
#'
#' @noRd
.pmt_class_name <- function(code, lev) {
  n_lev <- length(lev)
  out <- rep("concordant", length(code))
  # lev[0] is character(0), which would silently shorten a vectorised
  # subscript, so the concordant rows never index lev at all.
  disc <- code != 0L
  out[disc] <- paste(lev[code[disc] %% n_lev], lev[code[disc] %/% n_lev],
    sep = " vs "
  )
  out
}


#' The relabelling statistic as a closure over the fixed rows (internal)
#'
#' `ri` / `ti` index the reference and test species of each module-direction,
#' so a labelling enters only through the trait values it gives those two
#' species. Everything a labelling cannot change --- which rows exist, their
#' effect sizes, the block structure --- is captured once.
#'
#' Rows are partitioned by [.pmt_class_code()], the same helper the reported
#' `$class_means` uses.
#'
#' @noRd
.pmt_make_stat <- function(z, ri, ti, n_lev, form) {
  if (identical(form, "difference")) {
    return(function(lab) {
      conc <- lab[ri] == lab[ti]
      # A labelling with no rows on one side has no difference to report.
      # Dropping such draws silently would shrink the null without saying so.
      if (!any(conc) || all(conc)) {
        return(NA_real_)
      }
      mean(z[conc]) - mean(z[!conc])
    })
  }
  function(lab) {
    code <- .pmt_class_code(lab[ri], lab[ti], n_lev)
    parts <- split(z, code)
    if (length(parts) < 2L) {
      return(NA_real_)
    }
    m <- vapply(parts, mean, numeric(1))
    w <- lengths(parts)
    mbar <- sum(w * m) / sum(w)
    sqrt(sum(w * (m - mbar)^2) / sum(w))
  }
}


#' Run one relabelling null and summarise its resolution (internal)
#'
#' Enumerates while the space is small enough to be walked exactly, otherwise
#' samples. Enumeration is decided on the size of the space and not on
#' `n_perm`, so raising `n_perm` for precision cannot demote an exact null to
#' a sampled one.
#'
#' @noRd
.pmt_null <- function(stat_fun, obs, labels, blocks, n_perm, enum_max,
                      what, n_perm_supplied) {
  n_space <- .pmt_space(labels, blocks)
  exact <- n_space <= enum_max
  # A sampled null counts draws, not distinct labellings: the same labelling
  # can be drawn twice, so calling them labellings would overstate how much
  # of the space a reported count covers.
  unit <- if (exact) "labellings" else "draws"
  if (exact) {
    mat <- .pmt_enumerate(labels, blocks)
    null <- apply(mat, 1L, stat_fun)
    if (n_perm_supplied) {
      message(
        what, " null enumerated exactly over ", nrow(mat),
        " labellings; n_perm = ", n_perm, " ignored"
      )
    }
  } else {
    null <- vapply(
      seq_len(n_perm),
      function(i) stat_fun(.pmt_draw(labels, blocks)),
      numeric(1)
    )
  }
  n_all <- length(null)
  keep <- is.finite(null)
  if (!all(keep)) {
    warning(
      sum(!keep), " of ", n_all, " ", unit, " in the ", what,
      " null leave one side of the statistic empty and were dropped; ",
      "the null is over the ", sum(keep), " scorable ones"
    )
    null <- null[keep]
  }
  if (length(null) == 0L) {
    stop("no labelling in the ", what, " null scores a statistic")
  }
  # The statistic is a difference of means over sets that change with the
  # labelling, so two labellings that are mathematically tied need not be
  # bitwise equal. An exact >= would drop those from the tail and undercount
  # the tie at the maximum, which is the resolution diagnostic itself.
  tol <- 1e-9 * max(1, abs(obs), max(abs(null)))
  n_ge <- sum(null >= obs - tol)
  n_max <- sum(null >= max(null) - tol)
  n_draw <- length(null)
  # An enumerated null already contains the observed labelling, so no +1;
  # a sampled one does not, hence the usual (r + 1) / (m + 1).
  p_value <- if (exact) n_ge / n_draw else (n_ge + 1L) / (n_draw + 1L)
  p_min <- if (exact) 1 / n_draw else 1 / (n_draw + 1L)
  p_attainable <- if (exact) {
    n_max / n_draw
  } else {
    (n_max + 1L) / (n_draw + 1L)
  }
  floor_p <- max(p_min, p_attainable)
  if (floor_p > 0.05) {
    warning(
      "the ", what, " null over ", n_draw, " ", unit, " has a ",
      "smallest attainable p-value of ", signif(floor_p, 3),
      ", so p < 0.05 is unreachable for any signal",
      if (p_attainable > p_min) {
        paste0(
          " (", n_max, " ", unit, " share the maximum, so ",
          "separating them matters more than adding species)"
        )
      } else {
        " (add species, or relax the label counts)"
      }
    )
  }
  list(
    null_distribution = null,
    p_value = p_value,
    p_min = p_min,
    p_attainable = p_attainable,
    n_tied_max = n_max,
    n_distinct = length(unique(signif(null, 12))),
    rank = n_ge,
    exact = exact,
    n_labellings = n_space,
    n_scored = n_draw
  )
}


#' Reject a species-keyed map that names a species twice (internal)
#' @noRd
.pmt_check_unique <- function(nm, what) {
  if (anyDuplicated(nm) > 0L) {
    stop(
      what, " names a species more than once, so its trait would be read ",
      "from whichever entry comes first: ",
      paste(unique(nm[duplicated(nm)]), collapse = ", ")
    )
  }
  invisible(NULL)
}


#' Trait relabelling test on the all-pairs preservation matrix
#'
#' Asks whether modules are less well preserved between species that differ
#' in a trait than between species that share it, using every contrast in an
#' all-pairs [preservation_paired()] run. The null relabels species, holding
#' the label counts and the whole preservation matrix fixed.
#'
#' @section Why the effect size and not the q-value:
#' The statistic ranks contrasts by `Zsummary_std`, never by `p.value` or
#' `q.value`. A permutation p-value saturates: on the eight-species Pooideae
#' set with `n_perm = 2000` and 511 module-directions, `q.value` takes 172
#' distinct values, 35 modules sit tied at the floor of 0.00071 and 20 at
#' exactly 1. Among those 35 floored modules `Zsummary_std` spans 6.4 to
#' 66.7 --- a tenfold range of effect that the p-value cannot see at all.
#' Anything that ranks, weights or orders modules must therefore read the
#' standardised effect, which is continuous; p and q are for the significance
#' *call* only. `$saturation` reports how far the supplied q-values are
#' collapsed, so a caller who reaches for them anyway can see what resolution
#' is left.
#'
#' @section The statistic:
#' Each row of `classification` is one module in one direction of one
#' contrast, and its two species are trait-*concordant* when they share a
#' trait value and *discordant* otherwise. With two trait levels the
#' statistic is
#' \code{mean(Zsummary_std | concordant) - mean(Zsummary_std | discordant)},
#' so a **positive** value means discordant pairs are the less preserved
#' ones, which is the biological hypothesis, and the test is one-sided upward.
#'
#' With more than two levels a single difference no longer describes the
#' design, because discordance is no longer one thing. The classes are then
#' the concordant rows (pooled across levels) plus one class per unordered
#' pair of distinct levels, and the statistic is the row-weighted standard
#' deviation of the class means --- their between-class dispersion, on the
#' `Zsummary_std` scale. Larger is still more structure, so the test stays
#' one-sided upward, but it is no longer signed and does not say which
#' direction the difference runs; read `$class_means` for that. Because the
#' concordant classes are pooled, a two-level trait has exactly two classes
#' and the dispersion form is never used for it: the binary case is the
#' difference above, exactly.
#'
#' Every module-direction counts once, so a contrast with many modules
#' contributes more than one with few. That is deliberate --- the unit is the
#' module, and a species pair with more modules carries more evidence --- but
#' it means a single hugely partitioned species can dominate; check
#' `$rows_per_pair`.
#'
#' @section The two nulls:
#' `p_free` permutes trait labels over all species, preserving the label
#' counts (`choose(8, 4) = 70` labellings for a balanced binary trait on
#' eight species). `p_blocked` permutes only *within* each `block`, which
#' holds the phylogeny fixed and is the conservative, phylogenetically
#' controlled null (16 labellings for four genera of two). Both are reported
#' because their agreement is the diagnostic: a free p-value far smaller than
#' the blocked one means the larger label space bought its resolution by
#' breaking phylogenetic control rather than by measuring the trait.
#'
#' Neither floor is `1 / n_labellings`, because relabellings that leave every
#' pair's concordance alone reproduce the statistic exactly. Swapping the two
#' trait values everywhere is one such, so a binary design always carries at
#' least a twofold tie at the maximum: 2/70 = 0.029 free on eight species
#' split 4/4, and 2/16 = 0.125 within four blocks of two --- which is why the
#' restricted null of a four-genus paired design cannot reach alpha on its
#' own. With more levels the dispersion cannot tell the level *names* apart
#' either, so every renaming of the levels that preserves the label counts
#' ties as well. That is `g!` renamings only when the design is balanced;
#' in general it is `prod(factorial(table(table(labels))))`, since a renaming
#' has to map each level onto one of the same size. Balanced three levels on
#' six species tie 6 of 90 labellings, a floor of 0.067; an unbalanced 2/2/1
#' split on five species ties only 2 of 30 --- the same 0.067, but from two
#' ties rather than six. `p_attainable` reports whichever floor actually
#' binds and `n_tied_max` how many labellings share it; the function warns
#' when that floor exceeds 0.05, since no signal can then reach significance.
#'
#' Both floors describe an *enumerated* null. When a space is too large to
#' enumerate, `p_min`, `p_attainable` and `n_tied_max` are counted over the
#' `n_perm` draws instead: `p_min` is then the sampling floor
#' `1 / (n_perm + 1)`, which can sit far below the labelling-space floor the
#' design really imposes, and `n_tied_max` counts drawn labellings with
#' repeats rather than distinct ones. `exact` says which regime a null is in,
#' and only an exact one bounds the design.
#'
#' @section Multiplicity scope of the q-values:
#' [preservation_paired()] corrects within each contrast: `q.value` is
#' Benjamini-Hochberg over that contrast's own modules, a dozen or two tests.
#' An all-pairs matrix reframes those contrasts as one analysis --- every
#' module-direction tests the same hypothesis --- so the per-contrast
#' correction under-corrects, by roughly the number of contrasts. Whether a
#' global correction is available instead is a question about the
#' permutation floor, not about the data: with `n_tests` module-directions
#' and a smallest attainable p-value of `1 / (n_perm + 1)`, BH cannot put any
#' module below `n_tests / (n_perm + 1)`. For 511 tests that is 0.26 at
#' `n_perm = 2000`, which nothing can pass, and 0.026 at `n_perm = 20000`,
#' which is usable. Pass `n_perm_pres` and `$saturation$q_floor_global`
#' reports that number, with a warning when it exceeds 0.05. The relabelling
#' test itself is untouched --- it never reads a q-value --- but any
#' per-module significance call taken off the same matrix is.
#'
#' @section Why within-block pairs are excluded:
#' In a paired design such as the Pooideae one --- each genus contributing an
#' annual and a perennial --- every within-genus pair is trait-discordant
#' while every trait-concordant pair is between-genus. Trait status and
#' phylogenetic distance are then perfectly confounded within that subset,
#' and the confound runs *against* the hypothesis: close relatives are better
#' preserved, so the within-genus discordant pairs are pulled upward and mask
#' a real effect. `exclude_within_block = TRUE` (the default) drops those
#' rows. The exclusion is by block membership, which no relabelling changes,
#' so the same rows are excluded under every labelling and the null stays
#' valid.
#'
#' @param classification The `classification` data frame from
#'   [preservation_paired()], ideally over an all-pairs `pairs` table (see
#'   [all_species_pairs()]). Must carry `reference`, `test` and the effect
#'   column named by `statistic`; `q.value` is used only for the
#'   `$saturation` diagnostic. Rows whose effect is `NA` --- `"untested"`
#'   modules, and any module whose null correlation was unavailable --- are
#'   dropped, since nothing was measured for them.
#' @param group Named vector mapping species to trait values. Any number of
#'   levels is allowed. Species not appearing in `classification` are ignored,
#'   so the permuted label multiset is the one the matrix actually holds.
#' @param block Optional named vector mapping species to a phylogenetic group
#'   (a genus, say). Supplying it enables the restricted within-block null and
#'   the within-block exclusion. With `block = NULL` only the free null is
#'   run and `exclude_within_block` has nothing to act on.
#' @param statistic Which effect column to average. `"zsummary"` (default)
#'   uses `Zsummary_std`, which has unit variance under the permutation null
#'   and is therefore comparable across modules and contrasts; `"zsummary_raw"`
#'   uses `Zsummary`, whose null spread depends on the correlation between the
#'   two preservation statistics and so differs from module to module. There
#'   is deliberately no q-value option.
#' @param exclude_within_block Drop rows whose two species share a block
#'   (default `TRUE`). Ignored when `block` is `NULL`.
#' @param n_perm Number of labellings to draw when a space is too large to
#'   enumerate. `NULL` (default) means 10000 if it comes to that. Ignored, with
#'   a message when it was supplied explicitly, for a space that is enumerated.
#' @param enum_max Enumerate a null while its space holds at most this many
#'   labellings (default 50000). Lower it to cap runtime on wide designs; the
#'   free space grows as the multinomial coefficient, so 20 species split
#'   10/10 is already 184756. A sampled null reports its floor over draws,
#'   not over the labelling space --- see the last paragraph of *The two
#'   nulls*.
#' @param n_perm_pres The `n_perm` used in the [module_preservation()] runs
#'   behind `classification`, if known. Nothing in the test reads it; it lets
#'   `$saturation` report `q_floor_global`, the smallest q-value a *global*
#'   Benjamini-Hochberg correction over the whole matrix could reach,
#'   `n_tests / (n_perm_pres + 1)`. `NULL` (default) leaves that field `NA`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{The statistic under the true labelling.}
#'     \item{statistic}{Echo of the input.}
#'     \item{form}{`"difference"` for two trait levels, `"dispersion"`
#'       otherwise.}
#'     \item{class_means}{Data frame with `class`, `n` and `mean_z`: the
#'       observed class means the statistic is built from.}
#'     \item{rows_per_pair}{Data frame with `sp1`, `sp2`, `concordant`,
#'       `same_block` and `n`, one row per species pair contributing.}
#'     \item{free}{List for the free null: `null_distribution`, `p_value`,
#'       `p_min`, `p_attainable`, `n_tied_max`, `n_distinct`, `rank`,
#'       `exact`, `n_labellings`, `n_scored`. `p_attainable` is the tie-aware
#'       floor --- labellings sharing a statistic cannot be separated, so
#'       nothing scores below the tail of the maximum --- and `n_tied_max`
#'       says how many labellings that tail holds. When `exact` is `FALSE`
#'       both are counted over the `n_scored` draws rather than over the
#'       `n_labellings` of the space.}
#'     \item{blocked}{The same list for the within-block null, or `NULL` when
#'       `block` was not supplied.}
#'     \item{p_free,p_blocked}{The two p-values, lifted out for convenience.
#'       `p_blocked` is `NA` without a `block`.}
#'     \item{saturation}{List describing the resolution of the supplied
#'       q-values: `n_tests`, `n_distinct`, `q_floor`, `n_at_floor`,
#'       `n_at_one`, plus `p_min_pres` and `q_floor_global` (both `NA`
#'       without `n_perm_pres`). `n_tests` counts every row of the supplied
#'       `classification` carrying a finite q-value --- the multiplicity
#'       population is the whole matrix, so this is deliberately *not*
#'       `n_rows` below, which is the subset this test averages over after
#'       dropping unmeasured rows and within-block pairs. Reported next to
#'       the p-values so the reader can see how much of the input's
#'       significance scale is collapsed.}
#'     \item{n_rows,n_pairs,n_excluded}{Rows used, distinct species pairs
#'       behind them, and rows dropped for sharing a block.}
#'     \item{species,group,block}{The design actually tested.}
#'   }
#'
#' @examples
#' \dontrun{
#' pairs <- all_species_pairs(names(group))
#' res <- preservation_paired(mods, nets, ortho, pairs, group = group)
#' pmt <- preservation_matrix_test(res$classification, group,
#'   block = genus
#' )
#' c(pmt$observed, pmt$p_free, pmt$p_blocked)
#' }
#'
#' @seealso [preservation_paired()], [all_species_pairs()],
#'   [tag_permutation()]
#' @export
preservation_matrix_test <- function(classification, group, block = NULL,
                                     statistic = c(
                                       "zsummary",
                                       "zsummary_raw"
                                     ),
                                     exclude_within_block = TRUE,
                                     n_perm = NULL,
                                     enum_max = 50000L,
                                     n_perm_pres = NULL) {
  # Captured before n_perm is defaulted: missing() reports FALSE once an
  # argument has been written to, so this cannot be asked for later.
  n_perm_supplied <- !is.null(n_perm)
  statistic <- match.arg(statistic)
  z_col <- if (statistic == "zsummary") "Zsummary_std" else "Zsummary"

  if (!is.data.frame(classification)) {
    stop("classification must be a data frame from preservation_paired()")
  }
  need <- c("reference", "test", z_col)
  miss <- setdiff(need, names(classification))
  if (length(miss) > 0L) {
    stop(
      "classification missing columns: ", paste(miss, collapse = ", "),
      " (expected the classification table of preservation_paired())"
    )
  }
  if (is.null(names(group)) || length(group) == 0L) {
    stop("group must be a named vector mapping species to trait values")
  }
  if (!is.null(block) && is.null(names(block))) {
    stop("block must be a named vector mapping species to a group")
  }
  # group[species] takes the first match, so a repeated name would hand a
  # species whichever trait happened to be listed first -- silently, and with
  # the wrong labelling then permuted as if it were the design.
  .pmt_check_unique(names(group), "group")
  if (!is.null(block)) {
    .pmt_check_unique(names(block), "block")
  }
  if (is.null(n_perm)) {
    n_perm <- 10000L
  }
  bad_perm <- !is.numeric(n_perm) || length(n_perm) != 1L || is.na(n_perm)
  if (bad_perm || n_perm < 1) {
    stop("n_perm must be a single positive number or NULL")
  }
  n_perm <- as.integer(min(n_perm, .Machine$integer.max))
  bad_enum <- !is.numeric(enum_max) || length(enum_max) != 1L
  if (bad_enum || is.na(enum_max) || enum_max < 1) {
    stop("enum_max must be a single positive number")
  }
  if (!is.null(n_perm_pres)) {
    bad_pres <- !is.numeric(n_perm_pres) || length(n_perm_pres) != 1L
    if (bad_pres || is.na(n_perm_pres) || n_perm_pres < 1) {
      stop("n_perm_pres must be a single positive number or NULL")
    }
  }

  cls <- classification
  ref <- as.character(cls$reference)
  tst <- as.character(cls$test)
  z <- as.numeric(cls[[z_col]])

  # Saturation of the supplied q-values, reported rather than used: it is the
  # measurement that forces the statistic onto Zsummary_std in the first
  # place, so the reader should see it next to the p-values. Counted over the
  # rows as supplied, before the NA and within-block drops below: the
  # multiplicity population is every test the matrix ran, not the subset this
  # statistic averages.
  qv <- if ("q.value" %in% names(cls)) as.numeric(cls$q.value) else NA_real_
  qv_ok <- qv[is.finite(qv)]
  n_tests <- length(qv_ok)
  # Benjamini-Hochberg over the whole matrix cannot reach below
  # n_tests * p_min, so at a given n_perm the global correction the all-pairs
  # framing calls for may not be available at all.
  p_min_pres <- if (is.null(n_perm_pres)) NA_real_ else 1 / (n_perm_pres + 1)
  saturation <- list(
    n_tests = n_tests,
    n_distinct = length(unique(qv_ok)),
    q_floor = if (n_tests > 0L) min(qv_ok) else NA_real_,
    n_at_floor = if (n_tests > 0L) sum(qv_ok == min(qv_ok)) else 0L,
    n_at_one = sum(qv_ok >= 1),
    p_min_pres = p_min_pres,
    # No q-values means no multiplicity to correct, which is not the same as
    # a floor of zero.
    q_floor_global = if (n_tests == 0L) {
      NA_real_
    } else {
      min(1, n_tests * p_min_pres)
    }
  )
  if (isTRUE(saturation$q_floor_global > 0.05)) {
    warning(
      "a global Benjamini-Hochberg correction over the ", n_tests,
      " tests of this matrix cannot reach below q = ",
      signif(saturation$q_floor_global, 3), " at n_perm = ", n_perm_pres,
      "; the supplied q-values are corrected per contrast, which is the ",
      "wrong scope for an all-pairs matrix. Raise n_perm in ",
      "module_preservation()"
    )
  }

  # A module whose statistic is NA was never measured, and a self-contrast is
  # not a comparison; either would enter the mean as if it carried evidence.
  keep <- is.finite(z) & !is.na(ref) & !is.na(tst) & ref != tst
  keep <- keep & ref %in% names(group) & tst %in% names(group)
  if (!any(keep)) {
    stop(
      "no usable rows: every row has a missing ", z_col,
      ", compares a species with itself, or names a species outside group"
    )
  }
  no_group <- !(ref %in% names(group)) | !(tst %in% names(group))
  n_dropped_group <- sum(no_group & is.finite(z))
  if (n_dropped_group > 0L) {
    warning(
      n_dropped_group, " row(s) name a species with no group entry ",
      "and were dropped"
    )
  }
  ref <- ref[keep]
  tst <- tst[keep]
  z <- z[keep]

  n_excluded <- 0L
  if (!is.null(block)) {
    miss_bl <- setdiff(unique(c(ref, tst)), names(block))
    if (length(miss_bl) > 0L) {
      stop("block missing entries for: ", paste(miss_bl, collapse = ", "))
    }
    if (isTRUE(exclude_within_block)) {
      same <- as.character(block[ref]) == as.character(block[tst])
      n_excluded <- sum(same)
      ref <- ref[!same]
      tst <- tst[!same]
      z <- z[!same]
      if (length(z) == 0L) {
        stop(
          "excluding within-block rows leaves nothing to test; ",
          "set exclude_within_block = FALSE"
        )
      }
    }
  }

  species <- sort(unique(c(ref, tst)))
  grp <- as.character(group[species])
  if (anyNA(grp)) {
    stop(
      "group has NA trait values for: ",
      paste(species[is.na(grp)], collapse = ", ")
    )
  }
  levels_g <- sort(unique(grp))
  n_lev <- length(levels_g)
  if (n_lev < 2L) {
    stop(
      "group takes a single value over the tested species, so no pair ",
      "is trait-discordant and there is nothing to relabel"
    )
  }
  labels <- match(grp, levels_g)
  ri <- match(ref, species)
  ti <- match(tst, species)

  # Two levels give exactly two classes -- concordant (pooled over levels)
  # and the single discordant pair -- so the general dispersion would be the
  # unsigned half of the difference. Take the difference itself there: it is
  # signed, and the sign is the hypothesis.
  form <- if (n_lev == 2L) "difference" else "dispersion"
  stat_fun <- .pmt_make_stat(z, ri, ti, n_lev, form)
  obs <- stat_fun(labels)
  if (!is.finite(obs)) {
    stop(
      "the observed labelling leaves one side of the statistic empty; ",
      "there is no trait contrast to test"
    )
  }

  blocks_v <- if (is.null(block)) {
    NULL
  } else {
    as.integer(factor(as.character(block[species])))
  }
  free <- .pmt_null(
    stat_fun, obs, labels, NULL, n_perm, enum_max,
    "free", n_perm_supplied
  )
  blocked <- if (is.null(block)) {
    NULL
  } else {
    .pmt_null(
      stat_fun, obs, labels, blocks_v, n_perm, enum_max,
      "within-block", n_perm_supplied
    )
  }

  # Observed class means: the pieces the statistic is assembled from, and the
  # only place the direction of a dispersion statistic can be read.
  cls_name <- .pmt_class_name(
    .pmt_class_code(labels[ri], labels[ti], n_lev), levels_g
  )
  parts <- split(z, cls_name)
  class_means <- data.frame(
    class = names(parts),
    n = as.integer(lengths(parts)),
    mean_z = vapply(parts, mean, numeric(1)),
    stringsAsFactors = FALSE
  )
  class_means <- class_means[order(
    class_means$class != "concordant",
    class_means$class
  ), , drop = FALSE]
  rownames(class_means) <- NULL

  key <- paste(pmin(ref, tst), pmax(ref, tst), sep = "\x01")
  first <- !duplicated(key)
  rows_per_pair <- data.frame(
    sp1 = pmin(ref, tst)[first],
    sp2 = pmax(ref, tst)[first],
    concordant = (grp[ri] == grp[ti])[first],
    same_block = if (is.null(block)) {
      rep(NA, sum(first))
    } else {
      (as.character(block[ref]) == as.character(block[tst]))[first]
    },
    n = as.integer(table(key)[key[first]]),
    stringsAsFactors = FALSE
  )
  rownames(rows_per_pair) <- NULL

  list(
    observed = obs,
    statistic = statistic,
    form = form,
    class_means = class_means,
    rows_per_pair = rows_per_pair,
    free = free,
    blocked = blocked,
    p_free = free$p_value,
    p_blocked = if (is.null(blocked)) NA_real_ else blocked$p_value,
    saturation = saturation,
    n_rows = length(z),
    n_pairs = sum(first),
    n_excluded = n_excluded,
    species = species,
    group = stats::setNames(grp, species),
    block = if (is.null(block)) {
      NULL
    } else {
      stats::setNames(as.character(block[species]), species)
    }
  )
}
