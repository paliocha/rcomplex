# Sparse (dgCMatrix) network support.
#
# A network object is `list(network = , threshold = , ...)` where `network`
# is either a dense numeric matrix or a `dgCMatrix` holding both triangles
# of the thresholded co-expression matrix (column j = neighbours of gene j).
# Every C++ consumer of a network goes through `.net_cpp_args()`, so storage
# class validation and the store guard live in one place.


#' Is a network object stored as a sparse dgCMatrix?
#' @noRd
.net_is_sparse <- function(net) inherits(net$network, "dgCMatrix")


#' Validate a network's storage class and run the store guard
#'
#' Accepts a dense numeric matrix or a `dgCMatrix`. Any other `Matrix` class
#' (e.g. the `dsCMatrix` returned by `Matrix::Matrix(m, sparse = TRUE)` for
#' symmetric input) is rejected. A sparse network must be square with
#' identical, non-NULL row and column names (column j = neighbours of gene
#' j; the R wrappers map gene names to indices through the row names) and
#' must not store diagonal entries (the sparse effective-density formula
#' counts every stored entry).
#' Analysing at a threshold below what the store holds would silently miss
#' discarded edges, so that is an error: against `store_threshold` when the
#' network carries one, otherwise against the smallest stored value
#' (hand-built sparse networks). Exact equality passes.
#'
#' @param net Network object.
#' @param thr Analysis threshold.
#' @return `net$network`, invisibly.
#' @noRd
.net_check <- function(net, thr) {
  m <- net$network
  if (.net_is_sparse(net)) {
    if (nrow(m) != ncol(m)) {
      stop("network must be a square dgCMatrix (got ", nrow(m), " x ",
           ncol(m), "); see as_sparse_network()")
    }
    dn <- dimnames(m)
    if (is.null(dn[[1L]]) || is.null(dn[[2L]]) ||
        !identical(dn[[1L]], dn[[2L]])) {
      stop("network must have identical, non-NULL row and column names; ",
           "see as_sparse_network()")
    }
    dp <- diff(m@p)
    if (length(m@x) > 0L && length(dp) == ncol(m) && all(dp >= 0L) &&
        m@p[1L] == 0L && m@p[length(m@p)] == length(m@i) &&
        length(m@i) == length(m@x) &&
        any(m@i == rep.int(seq_len(ncol(m)) - 1L, dp))) {
      # malformed slots skip this check and fall through to the C++
      # validator in neighbor_lists_sparse() (canonical error messages)
      stop("network must not store diagonal entries; ",
           "see as_sparse_network()")
    }
    if (!is.null(net$store_threshold)) {
      if (thr < net$store_threshold) {
        stop("threshold ", thr, " is below the stored threshold ",
             net$store_threshold,
             if (!is.null(net$store_density)) {
               paste0(" (store_density = ", net$store_density, ")")
             },
             "; recompute with a larger store_density")
      }
    } else if (length(m@x) > 0L && thr < min(m@x)) {
      stop("threshold ", thr, " is below the smallest stored value ",
           min(m@x), "; the sparse network cannot represent this density ",
           "(use threshold >= that value or rebuild with a larger ",
           "store_density)")
    }
  } else if (methods::is(m, "Matrix")) {
    stop("network must be a dgCMatrix; see as_sparse_network()")
  }
  invisible(m)
}


#' Build the C++ argument list for a network object
#'
#' Runs `.net_check()` and returns the arguments for the matching C++ entry
#' point: `list(net = <matrix>, thr = thr)` for dense networks,
#' `list(p = , i = , x = , thr = thr)` (the `dgCMatrix` slots, passed
#' without conversion) for sparse ones.
#'
#' @param net Network object.
#' @param thr Analysis threshold to pass to C++.
#' @noRd
.net_cpp_args <- function(net, thr) {
  m <- .net_check(net, thr)
  if (.net_is_sparse(net)) {
    list(p = m@p, i = m@i, x = m@x, thr = thr)
  } else {
    list(net = m, thr = thr)
  }
}


#' Storage class of a network pair for a two-species C++ call
#'
#' Both networks must use the same storage class (the C++ entry points take
#' either two dense matrices or two sets of `dgCMatrix` slots).
#'
#' @return `TRUE` if both sparse, `FALSE` if both dense; error if mixed.
#' @noRd
.net_pair_sparse <- function(net1, net2) {
  sparse1 <- .net_is_sparse(net1)
  if (sparse1 != .net_is_sparse(net2)) {
    stop("net1 and net2 must be both dense or both sparse; ",
         "see as_sparse_network()")
  }
  sparse1
}


#' Convert a dense network matrix to a dgCMatrix
#'
#' Keeps off-diagonal entries `>= thr`, both triangles, with dimnames.
#' Intended for tests and small matrices: `row()`/`col()` allocate two
#' n x n integer matrices. Never use `Matrix::Matrix(m, sparse = TRUE)` for
#' this - it returns a `dsCMatrix` for symmetric input.
#'
#' @param m Dense numeric matrix (n x n).
#' @param thr Store threshold.
#' @return A `dgCMatrix`.
#' @noRd
dense_to_dgc <- function(m, thr) {
  ij <- which(m >= thr & row(m) != col(m), arr.ind = TRUE)
  Matrix::sparseMatrix(
    i = ij[, 1L], j = ij[, 2L], x = m[ij],
    dims = dim(m), dimnames = dimnames(m)
  )
}


#' Convert a dense network object to the sparse representation
#'
#' Thresholds the dense co-expression matrix at the `store_density`
#' quantile and repacks the surviving off-diagonal entries (both
#' triangles) as a `Matrix::dgCMatrix`, exactly as
#' `compute_network(sparse = TRUE)` does. All other fields of `net` are
#' carried over; `store_density`, `store_threshold` and
#' `params$store_density` are added. Values below `store_threshold` are
#' discarded, so downstream analyses at a threshold below it are refused.
#'
#' @param net Dense network object: output of
#'   `compute_network(sparse = FALSE)`, or a hand-built
#'   `list(network = <matrix>, threshold = )` with gene dimnames.
#' @param store_density Fraction of top edges to keep in the store
#'   (default 0.05). Must be in (0, 1) and, when the network records its
#'   density in `params$density`, at least that density (the store must
#'   contain every analysis edge).
#' @return The network object with `network` as a `dgCMatrix` plus
#'   `store_density` and `store_threshold` fields.
#'
#' @examples
#' \dontrun{
#' net <- compute_network(x, density = 0.03, sparse = FALSE)
#' net_sparse <- as_sparse_network(net, store_density = 0.05)
#' }
#'
#' @export
as_sparse_network <- function(net, store_density = 0.05) {
  if (!is.list(net) || is.null(net$network) || is.null(net$threshold)) {
    stop("net must be a network object from compute_network()")
  }
  if (.net_is_sparse(net)) {
    stop("net is already sparse")
  }
  m <- net$network
  if (!is.matrix(m) || !is.numeric(m)) {
    stop("network must be a dense numeric matrix; got ",
         paste(class(m), collapse = "/"))
  }
  if (nrow(m) != ncol(m) || is.null(rownames(m))) {
    stop("network must be a square matrix with gene names")
  }
  if (!is.numeric(store_density) || length(store_density) != 1L ||
      store_density <= 0 || store_density >= 1) {
    stop("store_density must be a single number in (0, 1) (exclusive)")
  }
  d <- net$params$density
  if (!is.null(d) && store_density < d) {
    stop("store_density (", store_density, ") must be >= the network's ",
         "density (", d, "); the store must contain every analysis edge")
  }

  store_thr <- density_threshold_cpp(m, store_density)
  slots <- extract_sparse_cpp(m, store_thr, 1L)
  spnet <- methods::new(
    "dgCMatrix", i = slots$i, p = slots$p, x = slots$x,
    Dim = dim(m), Dimnames = dimnames(m)
  )
  modifyList(net, list(
    network = spnet,
    params = c(net$params, list(store_density = store_density)),
    store_density = store_density,
    store_threshold = store_thr
  ))
}


#' One-line storage description of a network object (for print methods)
#' @noRd
.net_describe <- function(net) {
  if (.net_is_sparse(net)) {
    paste0("dgCMatrix, ", length(net$network@x), " stored entries",
           if (!is.null(net$store_density)) {
             paste0(", store_density = ", net$store_density)
           })
  } else {
    paste0("dense matrix, ", nrow(net$network), " x ", ncol(net$network))
  }
}
