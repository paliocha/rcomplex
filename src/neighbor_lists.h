// neighbor_lists.h
// Per-column neighbour lists from a thresholded co-expression network.
//
// Both hot C++ paths (neighborhood comparison, HOG permutation) only need
// neighbour lists, so list construction is decoupled from the storage format:
//
//   neighbor_lists_dense  — arma::mat (column-major, colptr access)
//   neighbor_lists_sparse — dgCMatrix slots (p, i, x); n = p.size() - 1
//
// Both keep entries with value >= thr and row != col, and return one
// ascending-sorted vector per column. Column j of a symmetric network lists
// the neighbours of gene j. Sparse lists inherit their order from the
// dgCMatrix invariant (row indices sorted within each column), which
// neighbor_lists_sparse() enforces in a serial validation pre-pass.
//
// Header-only. No std::unordered_map (Homebrew clang ABI issue).

#ifndef RCOMPLEX_NEIGHBOR_LISTS_H
#define RCOMPLEX_NEIGHBOR_LISTS_H

#include <RcppArmadillo.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

inline std::vector<std::vector<int>> neighbor_lists_dense(
    const arma::mat& m, double thr, int n_cores
) {
    const int n = static_cast<int>(m.n_rows);
    std::vector<std::vector<int>> neighbors(n);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int i = 0; i < n; ++i) {
        const double* col_i = m.colptr(i);
        for (int j = 0; j < n; ++j) {
            if (i != j && col_i[j] >= thr) {
                neighbors[i].push_back(j);
            }
        }
    }
    return neighbors;
}

inline std::vector<std::vector<int>> neighbor_lists_sparse(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x,
    double thr,
    int n_cores
) {
    // Serial validation pre-pass. Runs before the parallel region so that
    // Rcpp::stop() is never called from inside it. Rejects slots the loop
    // below would otherwise read out of bounds or silently miscount: empty
    // p, p not starting at 0, p / i / x length mismatch, non-monotone p, and
    // row indices that are not strictly increasing within a column or not in
    // [0, n) (duplicates, out-of-range rows, non-square input).
    if (p.size() < 1) {
        Rcpp::stop("dgCMatrix slot p must have length >= 1 (got 0)");
    }
    const int n = static_cast<int>(p.size()) - 1;
    if (p[0] != 0) {
        Rcpp::stop("dgCMatrix slot p must start at 0 (got %d)", p[0]);
    }
    if (static_cast<R_xlen_t>(p[n]) != x.size() || i.size() != x.size()) {
        Rcpp::stop("dgCMatrix slots are inconsistent: p[n] = %d, "
                   "length(i) = %d, length(x) = %d",
                   p[n], static_cast<int>(i.size()),
                   static_cast<int>(x.size()));
    }
    // Verify ALL of p before any row-index scan: the scan for column c
    // reads i[p[c]..p[c+1]) and is safe only once every pointer is known
    // to be non-decreasing (with p[0] = 0 and p[n] = length(i) already
    // checked, that bounds each p[c+1] by length(i)). Checking p per
    // column while scanning would read i[] out of bounds inside the
    // validator itself for a later inflated pointer, e.g. p = {0, 100, 3}
    // with length(i) = 3.
    for (int c = 0; c < n; ++c) {
        if (p[c] > p[c + 1]) {
            Rcpp::stop("dgCMatrix slot p must be non-decreasing (column %d)",
                       c + 1);
        }
    }
    for (int c = 0; c < n; ++c) {
        int prev = -1;
        for (int k = p[c]; k < p[c + 1]; ++k) {
            const int r = i[k];
            if (r <= prev || r >= n) {
                Rcpp::stop("dgCMatrix row indices must be strictly increasing "
                           "within each column and in [0, n) (column %d)",
                           c + 1);
            }
            prev = r;
        }
    }

    std::vector<std::vector<int>> neighbors(n);

    // Raw pointers: no R API calls inside the parallel region
    const int* pp = p.begin();
    const int* ii = i.begin();
    const double* xx = x.begin();

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int c = 0; c < n; ++c) {
        for (int k = pp[c]; k < pp[c + 1]; ++k) {
            const int r = ii[k];
            if (r != c && xx[k] >= thr) {
                neighbors[c].push_back(r);
            }
        }
    }
    return neighbors;
}

#endif // RCOMPLEX_NEIGHBOR_LISTS_H
