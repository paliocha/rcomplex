# Degree-preserving edge-swap null for co-expressolog statistics (P6).


#' Default statistic: conserved-call counts per species pair
#'
#' Counts rows with `type == "conserved"` per species pair (named
#' `paste(species1, species2, sep = "~")`) plus their sum as `"total"`.
#'
#' @param edges Edge data frame from `find_coexpressologs()`.
#' @return Named numeric vector.
#' @noRd
.coexpressolog_conserved_counts <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(c(total = 0))
  }
  pair <- paste(edges$species1, edges$species2, sep = "~")
  conserved <- edges$type == "conserved"
  counts <- vapply(split(conserved, pair), sum, numeric(1))
  c(counts, total = sum(conserved))
}


#' Degree-preserving edge-swap null for co-expressolog statistics
#'
#' Tests whether an observed co-expressolog statistic exceeds what
#' network topology alone produces. Each species network is binarised at
#' its analysis threshold and rewired by degree-preserving edge swaps
#' (\code{igraph::keeping_degseq()}), which keeps every gene's degree but
#' destroys the correspondence between network neighbourhoods and the
#' ortholog mapping. \code{\link{find_coexpressologs}} then runs on the
#' rewired networks with exactly the same arguments (\code{...}) as the
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
#' carry no information after rewiring. Each network is rewired with
#' \code{swap_factor * igraph::ecount(g)} swaps and stays simple (no
#' loops, no multi-edges).
#'
#' Permutation \code{b} runs in a worker that calls
#' \code{set.seed(seed + b)} first, so results are reproducible and
#' independent of \code{n_cores}. On Unix the permutations run under
#' \code{parallel::mclapply()}; on Windows they run serially. The serial
#' path saves and restores the caller's RNG state around the loop, so
#' the ambient stream continues where the observed run left it. The
#' observed run uses the
#' ambient RNG: call \code{set.seed()} beforehand when \code{...}
#' includes settings that draw from it (e.g. the default
#' \code{pi0_method = "randomized"}). \code{method = "permutation"} in
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
#' @param swap_factor Edge swaps per permutation, as a multiple of the
#'   edge count of each thresholded network (default 10).
#' @param n_cores Number of parallel workers for the permutation loop
#'   (default 1).
#' @param seed Integer base seed; worker \code{b} uses
#'   \code{seed + b} (default 1).
#' @param ... Passed unchanged to \code{\link{find_coexpressologs}} for
#'   both the observed and every null run (\code{species_pairs},
#'   \code{method}, \code{alternative}, \code{alpha},
#'   \code{pval_combine}, \code{pi0_method}, ...).
#'
#' @return Data frame with one row per statistic and columns
#'   \code{statistic}, \code{observed}, \code{null_mean},
#'   \code{null_sd}, \code{null_max}, \code{fold} (observed /
#'   null_mean; \code{NA} when the null mean is 0) and \code{p_emp}
#'   (\code{(1 + sum(null >= observed)) / (n_perm + 1)}). The
#'   \code{n_perm x k} matrix of null statistics is attached as
#'   \code{attr(, "null")}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' res <- coexpressolog_null(networks, orthologs, n_perm = 100L,
#'                           n_cores = 4L)
#' res[res$statistic == "total", ]
#' }
#'
#' @export
coexpressolog_null <- function(networks, orthologs, statistic = NULL,
                               n_perm = 100L, swap_factor = 10L,
                               n_cores = 1L, seed = 1L, ...) {
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  is_sparse <- vapply(networks, .net_is_sparse, logical(1))
  if (!all(is_sparse)) {
    stop("coexpressolog_null() requires sparse networks; convert ",
         paste(names(networks)[!is_sparse], collapse = ", "),
         " with as_sparse_network()")
  }
  n_perm <- as.integer(n_perm)
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a single integer >= 1")
  }
  builtin_stat <- is.null(statistic)
  if (builtin_stat) {
    statistic <- .coexpressolog_conserved_counts
  }
  if (!is.function(statistic)) {
    stop("statistic must be NULL or a function(edges) -> named numeric")
  }

  observed <- statistic(find_coexpressologs(networks, orthologs, ...))
  if (!is.numeric(observed) || is.null(names(observed))) {
    stop("statistic must return a named numeric vector")
  }

  rewire_net <- function(net) {
    a <- net$network
    a@x <- as.numeric(a@x >= net$threshold)
    a <- Matrix::drop0(a)
    g <- igraph::graph_from_adjacency_matrix(a, mode = "undirected",
                                             diag = FALSE)
    g <- igraph::rewire(
      g,
      igraph::keeping_degseq(loops = FALSE,
                             niter = swap_factor * igraph::ecount(g))
    )
    a_perm <- igraph::as_adjacency_matrix(g, sparse = TRUE)
    dimnames(a_perm) <- dimnames(a)
    a_perm <- methods::as(
      methods::as(methods::as(a_perm, "dMatrix"), "generalMatrix"),
      "CsparseMatrix"
    )
    modifyList(net, list(network = a_perm, threshold = 1,
                         store_threshold = 1))
  }

  one_perm <- function(b) {
    set.seed(seed + b)
    nets_perm <- lapply(networks, rewire_net)
    statistic(find_coexpressologs(nets_perm, orthologs, ...))
  }

  use_mc <- .Platform$OS.type == "unix" && n_cores > 1L
  if (use_mc) {
    null_list <- parallel::mclapply(seq_len(n_perm), one_perm,
                                    mc.cores = n_cores,
                                    mc.preschedule = FALSE)
  } else {
    # one_perm() calls set.seed() in the caller's session here (no fork):
    # save/restore .Random.seed so the ambient stream continues where
    # the observed run left it, exactly as under mclapply()
    has_seed <- exists(".Random.seed", envir = globalenv(),
                       inherits = FALSE)
    old_seed <- if (has_seed) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(),
                        inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    null_list <- lapply(seq_len(n_perm), one_perm)
  }

  errs <- which(vapply(null_list, inherits, logical(1), "try-error"))
  if (length(errs)) {
    e <- null_list[[errs[1L]]]
    stop("permutation ", errs[1L], " failed: ",
         conditionMessage(attr(e, "condition")))
  }
  failed <- vapply(null_list, is.null, logical(1))
  if (any(failed)) {
    stop("parallel workers returned NULL for permutations: ",
         paste(which(failed), collapse = ", "))
  }

  nm <- names(observed)
  null_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(nm),
                     dimnames = list(NULL, nm))
  for (b in seq_len(n_perm)) {
    s <- null_list[[b]]
    if (!is.numeric(s) || is.null(names(s))) {
      stop("permutation ", b,
           " statistic did not return a named numeric vector")
    }
    miss <- setdiff(nm, names(s))
    if (length(miss) > 0L) {
      if (builtin_stat) {
        # a rewired run can leave a species pair with no overlap > 0
        # rows at all; for the built-in conserved count that IS a null
        # observation of 0 conserved calls, not an error
        s[miss] <- 0
      } else {
        stop("permutation ", b, " statistic is missing: ",
             paste(miss, collapse = ", "))
      }
    }
    null_mat[b, ] <- s[nm]
  }

  null_mean <- vapply(seq_along(nm), function(j) mean(null_mat[, j]),
                      numeric(1))
  null_sd <- vapply(seq_along(nm), function(j) stats::sd(null_mat[, j]),
                    numeric(1))
  null_max <- vapply(seq_along(nm), function(j) max(null_mat[, j]),
                     numeric(1))
  p_emp <- vapply(seq_along(nm), function(j) {
    (1 + sum(null_mat[, j] >= observed[[j]])) / (n_perm + 1)
  }, numeric(1))
  fold <- ifelse(null_mean == 0, NA_real_,
                 unname(observed) / null_mean)

  out <- data.frame(
    statistic = nm,
    observed = unname(observed),
    null_mean = null_mean,
    null_sd = null_sd,
    null_max = null_max,
    fold = fold,
    p_emp = p_emp,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "null") <- null_mat
  out
}
