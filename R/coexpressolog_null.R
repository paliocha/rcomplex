# Degree-preserving edge-swap null for co-expressolog statistics (P6).


#' Default statistic: conserved-call counts per species pair
#'
#' Counts rows with `type == "conserved"` per species pair (named
#' `paste(species1, species2, sep = "~")`) plus their sum as `"total"`.
#'
#' @param edges Edge data frame from `find_coexpressologs()`.
#' @return Named numeric vector.
#' @noRd
.coexpressolog_conserved_counts <- function(edges) {  # nolint
  if (is.null(edges) || nrow(edges) == 0L) {
    return(c(total = 0))
  }
  pair <- paste(edges$species1, edges$species2, sep = "~")
  conserved <- edges$type == "conserved"
  counts <- vapply(split(conserved, pair), sum, numeric(1))
  c(counts, total = sum(conserved))
}


#' Rewire a binary network by degree-preserving edge swaps
#'
#' @param a Symmetric binary `dgCMatrix` with both triangles stored.
#' @param swap_factor Swap trials as a multiple of the edge count.
#' @return `dgCMatrix` with the same dimensions, dimnames and degrees.
#' @noRd
.rewire_degseq <- function(a, swap_factor) {
  # the kernel validates the slots and reads each edge once from the upper
  # triangle, so nothing here indexes into @p / @i unchecked
  r <- rewire_degseq_cpp(a@p, a@i, a@x, swap_factor)
  Matrix::sparseMatrix(
    i = c(r$from, r$to), j = c(r$to, r$from), x = 1,
    dims = dim(a), dimnames = dimnames(a), index1 = FALSE
  )
}


#' Degree-preserving edge-swap null for co-expressolog statistics
#'
#' Tests whether an observed co-expressolog statistic exceeds what
#' network topology alone produces. Each species network is binarised at
#' its analysis threshold and rewired by degree-preserving edge swaps,
#' which keep every gene's degree but destroy the correspondence between
#' network neighbourhoods and the ortholog mapping.
#' \code{\link{find_coexpressologs}} then runs on the rewired networks
#' with exactly the same arguments (\code{...}) as the
#' observed run, and the statistic is compared against the resulting
#' null distribution.
#'
#' @details
#' Rewiring operates on the network thresholded at \code{net$threshold}
#' (the analysis density), not on the stored superset: sparse-store
#' entries with values at or above the threshold become edges, the rest
#' are dropped. The rewired networks are unweighted (every stored value
#' is 1, with \code{threshold = 1} and \code{store_threshold = 1}), so
#' only membership-based consumers are valid downstream; edge weights
#' carry no information after rewiring. A network with \code{m} edges gets
#' \code{ceiling(swap_factor * m)} swap trials, each the trial of igraph's
#' \code{keeping_degseq()} rewiring: two distinct edges are drawn uniformly,
#' a swap that would create a loop or a multi-edge is rejected, and the
#' rejected trial still counts. Counting it is what makes the chain sample
#' the realizations of a degree sequence uniformly. Networks stay simple (no
#' loops, no multi-edges). Adjacency is held as an \code{n x n} bit matrix
#' during rewiring, \code{n^2 / 8} bytes per network per worker (32 MB at
#' 16 000 genes).
#'
#' Permutation \code{b} runs in a worker that seeds itself from the base
#' seed and its own index (the package-wide per-task seed derivation), so
#' results are reproducible, independent of \code{n_cores}, and two
#' neighbouring base seeds share no rewiring. On Unix the permutations run
#' under \code{parallel::mclapply()}; on Windows they run serially. The seed
#' covers the observed run as well as the null runs, so \code{...}
#' settings that draw --- the default \code{pi0_method = "randomized"},
#' for one --- are pinned too. \code{method = "permutation"} in
#' \code{...} is allowed but slow (a full permutation test per rewired
#' network).
#'
#' @param networks Named list of sparse network objects
#'   (\code{compute_network(sparse = TRUE)} outputs), keyed by species
#'   abbreviation. Dense networks are rejected; convert them with
#'   \code{\link{as_sparse_network}}.
#' @param orthologs Data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog}, as for
#'   \code{\link{find_coexpressologs}}.
#' @param statistic Function mapping the \code{find_coexpressologs()}
#'   edge data frame to a named numeric vector, or \code{NULL} (default)
#'   for the count of rows with \code{type == "conserved"} per species
#'   pair (named \code{paste(species1, species2, sep = "~")}) plus
#'   \code{"total"}. For the built-in statistic, a null run whose
#'   rewiring leaves a species pair with no overlap > 0 rows records 0
#'   conserved calls for that pair; a user-supplied statistic missing a
#'   name from the observed run errors.
#' @param n_perm Number of rewired permutations (default 100).
#' @param swap_factor Swap trials per permutation, rejected ones
#'   included, as a multiple of the edge count of each thresholded network
#'   (default 10). Must be a single finite number > 0. The trial count is
#'   rounded up, so any positive factor makes at least one trial on a
#'   network with two or more edges; a network with fewer than two edges has
#'   no swap to make and is returned unchanged.
#' @param n_cores Number of parallel workers for the permutation loop
#'   (default 1).
#' @param seed Base seed for the run. \code{NULL} (default) draws one
#'   from the ambient RNG stream and leaves it advanced by that one draw,
#'   so consecutive unseeded calls give different nulls and a
#'   \code{set.seed()} in the caller's script reproduces the whole run.
#'   With a seed the call draws from a private stream and restores the
#'   caller's on exit --- the package-wide contract, see
#'   \code{\link{detect_modules}}.
#'
#'   What a default call guarantees is replayability, not a fixed null.
#'   The seed actually used --- drawn or supplied --- is recorded as
#'   \code{attr(result, "seed")}, and a drawn one is also announced in a
#'   message, so any run can be reproduced exactly by passing that value
#'   back as \code{seed}. Two default calls return different nulls; the
#'   size of that difference is what \code{null_se} and the
#'   \code{p_emp_lo} / \code{p_emp_hi} interval report. The default was
#'   \code{1L} before 0.3.0, which pinned the null of every default call
#'   while leaving the observed statistic free to drift, and hid the
#'   Monte Carlo error entirely.
#'
#'   Two seeds exactly \code{2^31 - 1} apart (both legal, since a seed may
#'   lie anywhere in \code{+/- .Machine$integer.max}) give every rewiring
#'   permutation the same seed, so the null is identical for the two runs
#'   even though the observed statistic is not: see the aliasing note under
#'   \code{.task_seed()}.
#' @param ... Passed unchanged to \code{\link{find_coexpressologs}} for
#'   both the observed and every null run (\code{species_pairs},
#'   \code{method}, \code{alternative}, \code{alpha},
#'   \code{pval_combine}, \code{pi0_method}, ...).
#'
#' @return Data frame with one row per statistic and columns
#'   \code{statistic}, \code{observed}, \code{null_mean},
#'   \code{null_sd}, \code{null_max}, \code{fold} (observed /
#'   null_mean; \code{NA} when the null mean is 0), \code{p_emp}
#'   (\code{(1 + n_ge) / (n_perm + 1)}), \code{n_ge} (null draws at or
#'   above the observed value), \code{null_se}
#'   (\code{null_sd / sqrt(n_perm)}, the Monte Carlo error of
#'   \code{null_mean}; the relative error of \code{fold} is
#'   \code{null_se / null_mean}, which is large whenever
#'   \code{null_mean} is near 0) and \code{p_emp_lo} / \code{p_emp_hi},
#'   an exact Clopper-Pearson 95\% interval for the exceedance
#'   probability that \code{p_emp} estimates --- an interval on that
#'   probability, not on \code{p_emp} itself. The \code{n_perm x k}
#'   matrix of null statistics is attached as \code{attr(, "null")}, and
#'   the base seed the run used as \code{attr(, "seed")}.
#'
#' @examples
#' \dontrun{
#' res <- coexpressolog_null(networks, orthologs,
#'   n_perm = 100L,
#'   n_cores = 4L
#' )
#' res[res$statistic == "total", ]
#'
#' # replay that exact run, whether or not it was seeded
#' same <- coexpressolog_null(networks, orthologs,
#'   n_perm = 100L,
#'   n_cores = 4L, seed = attr(res, "seed")
#' )
#' }
#'
#' @export
coexpressolog_null <- function(networks, orthologs, statistic = NULL,
                               n_perm = 100L, swap_factor = 10L,
                               n_cores = 1L, seed = NULL, ...) {
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  is_sparse <- vapply(networks, .net_is_sparse, logical(1))
  if (!all(is_sparse)) {
    stop(
      "coexpressolog_null() requires sparse networks; convert ",
      paste(names(networks)[!is_sparse], collapse = ", "),
      " with as_sparse_network()"
    )
  }
  # The rewiring kernel indexes an nrow x nrow bit matrix with column
  # indices, so every network must pass the square/dimnames validation that
  # find_coexpressologs() would apply -- including networks that the
  # species_pairs in `...` leave out of the observed run, which are rewired
  # all the same.
  for (net in networks) .net_check(net, net$threshold)
  n_perm <- as.integer(n_perm)
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a single integer >= 1")
  }
  # An NA, NaN or non-positive swap_factor used to reach the kernel as a
  # trial count below 1, which rewires nothing: the "null" was the observed
  # graph, and the run returned without a word.
  if (!is.numeric(swap_factor) || length(swap_factor) != 1L ||
        !is.finite(swap_factor) || swap_factor <= 0) {
    stop("swap_factor must be a single finite number > 0")
  }
  builtin_stat <- is.null(statistic)
  if (builtin_stat) {
    statistic <- .coexpressolog_conserved_counts
  }
  if (!is.function(statistic)) {
    stop("statistic must be NULL or a function(edges) -> named numeric")
  }

  # The base seed is drawn from the ambient stream when none was given, so
  # the caller sees exactly one draw and consecutive unseeded calls differ.
  # The scope covers the observed run too: with the default
  # pi0_method = "randomized" that run draws, and leaving it outside the
  # scope would make the observed statistic irreproducible under a seed.
  # A drawn root is announced here rather than at the return, so a run that
  # errors halfway through still names the seed that would reproduce it.
  # The root is recorded on the result either way; that is what makes an
  # unseeded run replayable at all.
  if (is.null(seed)) {
    seed_root <- sample.int(.Machine$integer.max, 1L)
    message(
      "coexpressolog_null(): drawn seed ", seed_root,
      "; replay this run with seed = ", seed_root
    )
  } else {
    # Anything that will not survive as.integer() -- a double outside the
    # integer range, a string that is not a number, NA, a vector -- has to
    # be caught before the bound check, which would otherwise let NA
    # through to set.seed(NA) and error on a length > 1 condition instead
    # of on the seed.
    #
    # Coercion decides acceptance, and only the message is chosen by type.
    # A rule like !is.numeric(seed) would reject seed = "7" and
    # seed = TRUE, which set.seed() and therefore .seed_scope() accept, and
    # would leave this the one seeded entry point in the package with its
    # own idea of a legal seed. Wrong length, wrong type and wrong
    # magnitude report apart, so no message names a limit the value did not
    # cross, and the magnitude message names both ends because -3e9 fails
    # at the low one.
    if (length(seed) != 1L) {
      stop(
        "seed must be NULL or a single value; got ", class(seed)[1L],
        " of length ", length(seed)
      )
    }
    seed_root <- suppressWarnings(as.integer(seed))
    if (is.na(seed_root)) {
      if (is.numeric(seed) && !is.na(seed)) {
        stop(
          "seed must lie within +/- .Machine$integer.max (",
          .Machine$integer.max, "); got ", format(seed)
        )
      }
      stop(
        "seed must be NULL or a value set.seed() accepts; got ",
        class(seed)[1L], " ", format(seed)
      )
    }
  }
  .seed_scope(seed_root)

  observed <- statistic(find_coexpressologs(networks, orthologs, ...))
  if (!is.numeric(observed) || is.null(names(observed))) {
    stop("statistic must return a named numeric vector")
  }

  rewire_net <- function(net) {
    a <- net$network
    a@x <- as.numeric(a@x >= net$threshold)
    a <- Matrix::drop0(a)
    modifyList(net, list(
      network = .rewire_degseq(a, swap_factor), threshold = 1,
      store_threshold = 1
    ))
  }

  one_perm <- function(b) {
    # Not set.seed(seed_root + b): that made the null at root r the null at
    # root r + 1 shifted by one permutation, so "try another seed" reused
    # all but one rewiring. .task_seed() also folds modulo 2^31 - 1, so a
    # root near the integer limit cannot overflow to set.seed(NA).
    set.seed(.task_seed(seed_root, 1L, b))
    nets_perm <- lapply(networks, rewire_net)
    statistic(find_coexpressologs(nets_perm, orthologs, ...))
  }

  # one_perm() calls set.seed() in the caller's session on the serial path
  # (no fork), which is why the scope opened above covers the whole
  # function rather than the observed run alone.
  use_mc <- .can_fork(n_cores)
  if (use_mc) {
    null_list <- parallel::mclapply(seq_len(n_perm), one_perm,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
  } else {
    null_list <- lapply(seq_len(n_perm), one_perm)
  }

  errs <- which(vapply(null_list, inherits, logical(1), "try-error"))
  if (length(errs)) {
    e <- null_list[[errs[1L]]]
    stop(
      "permutation ", errs[1L], " failed: ",
      conditionMessage(attr(e, "condition"))
    )
  }
  failed <- vapply(null_list, is.null, logical(1))
  if (any(failed)) {
    stop(
      "parallel workers returned NULL for permutations: ",
      paste(which(failed), collapse = ", ")
    )
  }

  nm <- names(observed)
  null_mat <- matrix(NA_real_,
    nrow = n_perm, ncol = length(nm),
    dimnames = list(NULL, nm)
  )
  for (b in seq_len(n_perm)) {
    s <- null_list[[b]]
    if (!is.numeric(s) || is.null(names(s))) {
      stop(
        "permutation ", b,
        " statistic did not return a named numeric vector"
      )
    }
    miss <- setdiff(nm, names(s))
    if (length(miss) > 0L) {
      if (builtin_stat) {
        # a rewired run can leave a species pair with no overlap > 0
        # rows at all; for the built-in conserved count that IS a null
        # observation of 0 conserved calls, not an error
        s[miss] <- 0
      } else {
        stop(
          "permutation ", b, " statistic is missing: ",
          paste(miss, collapse = ", ")
        )
      }
    }
    null_mat[b, ] <- s[nm]
  }

  null_mean <- vapply(
    seq_along(nm), function(j) mean(null_mat[, j]),
    numeric(1)
  )
  null_sd <- vapply(
    seq_along(nm), function(j) stats::sd(null_mat[, j]),
    numeric(1)
  )
  null_max <- vapply(
    seq_along(nm), function(j) max(null_mat[, j]),
    numeric(1)
  )
  n_ge <- vapply(seq_along(nm), function(j) {
    sum(null_mat[, j] >= observed[[j]])
  }, integer(1))
  p_emp <- (1 + n_ge) / (n_perm + 1)
  fold <- ifelse(null_mean == 0, NA_real_,
    unname(observed) / null_mean
  )
  # Monte Carlo error of the null mean, and an exact Clopper-Pearson 95%
  # interval on the exceedance probability that p_emp estimates. Both come
  # free from the null matrix, and together they say how much of the result
  # is the data and how much is this particular seed. The endpoints are
  # written out rather than left to qbeta(), whose shape arguments would be
  # 0 there.
  #
  # A user statistic can return NA or NaN on a permutation -- mean() over
  # a rewiring that called nothing conserved is NaN -- and then n_ge is
  # NA. Such a row reports NA rather than aborting the run, so the
  # subscripts below are NA-safe and the warning ignores those rows.
  null_se <- null_sd / sqrt(n_perm)
  ok <- !is.na(n_ge)
  p_emp_lo <- ifelse(ok, 0, NA_real_)
  pos <- ok & n_ge > 0L
  p_emp_lo[pos] <- stats::qbeta(0.025, n_ge[pos], n_perm - n_ge[pos] + 1L)
  p_emp_hi <- ifelse(ok, 1, NA_real_)
  lt <- ok & n_ge < n_perm
  p_emp_hi[lt] <- stats::qbeta(0.975, n_ge[lt] + 1L, n_perm - n_ge[lt])

  # Two mutually exclusive ways this n_perm cannot support the call it is
  # being asked to make: at or below 19 permutations the smallest
  # attainable p_emp, 1 / (n_perm + 1), is 0.05 or larger and so can never
  # be < 0.05, and above that a p_emp under 0.05 whose interval still
  # covers 0.05 is a call the next seed may not repeat.
  if (1 / (n_perm + 1) >= 0.05) {
    warning(
      "the edge-swap null over ", n_perm,
      " permutations has a smallest attainable p-value of ",
      signif(1 / (n_perm + 1), 3),
      ", so p < 0.05 is unreachable for any signal (use n_perm >= 20)"
    )
  } else {
    weak <- !is.na(p_emp) & !is.na(p_emp_hi) &
      p_emp < 0.05 & p_emp_hi >= 0.05
    if (any(weak)) {
      warning(
        "p_emp < 0.05 at n_perm = ", n_perm, " for ",
        paste(nm[weak], collapse = ", "),
        ", but the 95% interval on the exceedance probability reaches ",
        signif(max(p_emp_hi[weak]), 3),
        ": those calls cannot be separated from non-significance at this ",
        "n_perm"
      )
    }
  }

  out <- data.frame(
    statistic = nm,
    observed = unname(observed),
    null_mean = null_mean,
    null_sd = null_sd,
    null_max = null_max,
    fold = fold,
    p_emp = p_emp,
    n_ge = n_ge,
    null_se = null_se,
    p_emp_lo = p_emp_lo,
    p_emp_hi = p_emp_hi,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "null") <- null_mat
  attr(out, "seed") <- seed_root
  out
}
