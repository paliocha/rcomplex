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
// dgCMatrix invariant (row indices sorted within each column).
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
    const int n = static_cast<int>(p.size()) - 1;
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
