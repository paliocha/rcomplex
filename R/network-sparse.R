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
#' symmetric input) is rejected. For sparse networks that carry a
#' `store_threshold`, analysing at a threshold below it would silently miss
#' discarded edges, so that is an error. Hand-built sparse networks without
#' `store_threshold` are not guarded.
#'
#' @param net Network object.
#' @param thr Analysis threshold.
#' @return `net$network`, invisibly.
#' @noRd
.net_check <- function(net, thr) {
  m <- net$network
  if (.net_is_sparse(net)) {
    if (!is.null(net$store_threshold) && thr < net$store_threshold) {
      stop("threshold ", thr, " is below the stored threshold ",
           net$store_threshold, " (store_density = ", net$store_density,
           "); recompute with a larger store_density")
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
