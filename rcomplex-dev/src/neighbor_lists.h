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

// Serial validation pre-pass for dgCMatrix slots. Runs before any parallel
// region so that Rcpp::stop() is never called from inside one. Rejects slots
// the readers would otherwise read out of bounds or silently miscount: empty
// p, p not starting at 0, p / i / x length mismatch, non-monotone p, and row
// indices that are not strictly increasing within a column or not in [0, n)
// (duplicates, out-of-range rows, non-square input). Returns n.
inline int validate_dgc_slots(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x
) {
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
    return n;
}

inline std::vector<std::vector<int>> neighbor_lists_sparse(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x,
    double thr,
    int n_cores
) {
    const int n = validate_dgc_slots(p, i, x);

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

// ---------------------------------------------------------------------------
// Weighted induced subgraphs (module preservation)
//
// Preservation statistics read edge weights, not just membership, and only
// ever touch the ortholog-mappable genes. Both builders below therefore return
// the subgraph induced on `keep` in local index space (0 .. keep.size() - 1),
// with weights divided by `scale`.
//
// `keep` must be strictly ascending. That is what makes each returned
// neighbour list ascending too: the dense loop and the dgCMatrix column layout
// both walk rows in ascending global order, and an ascending `keep` makes the
// global-to-local map monotone. Triangle enumeration relies on that ordering
// to intersect two lists in linear time.
// ---------------------------------------------------------------------------

struct WeightedNeighbors {
    std::vector<std::vector<int>> idx;
    std::vector<std::vector<double>> w;
    double scale = 1.0;  // divisor applied to every weight
};

// Build the global -> local index map for an induced subgraph. Validates that
// `keep` is strictly ascending and in range.
inline std::vector<int> induced_local_map(
    const Rcpp::IntegerVector& keep, int n
) {
    std::vector<int> local(n, -1);
    int prev = -1;
    for (R_xlen_t u = 0; u < keep.size(); ++u) {
        const int g = keep[u];
        if (g <= prev || g < 0 || g >= n) {
            Rcpp::stop("induced subgraph indices must be strictly increasing "
                       "and in [0, n) (position %d)", static_cast<int>(u) + 1);
        }
        local[g] = static_cast<int>(u);
        prev = g;
    }
    return local;
}

// Largest thresholded off-diagonal weight within the induced subgraph, used to
// rescale adjacencies onto [0, 1]. Returns 1.0 when the subgraph has no edges
// so callers never divide by zero. Note the Z-scores are invariant to this
// choice (a common factor cancels in (obs - mean) / sd); it only puts the
// reported observed statistics on a comparable footing across species.
inline double induced_max_weight(const WeightedNeighbors& g) {
    double mx = 0.0;
    for (const auto& row : g.w) {
        for (const double v : row) {
            if (v > mx) mx = v;
        }
    }
    return (mx > 0.0) ? mx : 1.0;
}

inline void rescale_weights(WeightedNeighbors& g, double scale) {
    g.scale = scale;
    for (auto& row : g.w) {
        for (double& v : row) v /= scale;
    }
}

inline WeightedNeighbors induced_weighted_lists_dense(
    const arma::mat& m, double thr, const Rcpp::IntegerVector& keep
) {
    const int n = static_cast<int>(m.n_rows);
    const std::vector<int> local = induced_local_map(keep, n);
    const int mm = static_cast<int>(keep.size());

    WeightedNeighbors out;
    out.idx.resize(mm);
    out.w.resize(mm);

    for (int u = 0; u < mm; ++u) {
        const int g = keep[u];
        const double* col_g = m.colptr(g);
        for (int r = 0; r < n; ++r) {
            if (r == g) continue;
            const int lr = local[r];
            if (lr < 0) continue;
            const double v = col_g[r];
            if (v >= thr) {
                out.idx[u].push_back(lr);
                out.w[u].push_back(v);
            }
        }
    }
    return out;
}

inline WeightedNeighbors induced_weighted_lists_sparse(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x,
    double thr,
    const Rcpp::IntegerVector& keep
) {
    const int n = validate_dgc_slots(p, i, x);
    const std::vector<int> local = induced_local_map(keep, n);
    const int mm = static_cast<int>(keep.size());

    WeightedNeighbors out;
    out.idx.resize(mm);
    out.w.resize(mm);

    const int* pp = p.begin();
    const int* ii = i.begin();
    const double* xx = x.begin();

    for (int u = 0; u < mm; ++u) {
        const int g = keep[u];
        for (int k = pp[g]; k < pp[g + 1]; ++k) {
            const int r = ii[k];
            if (r == g) continue;
            const int lr = local[r];
            if (lr < 0) continue;
            if (xx[k] >= thr) {
                out.idx[u].push_back(lr);
                out.w[u].push_back(xx[k]);
            }
        }
    }
    return out;
}

#endif // RCOMPLEX_NEIGHBOR_LISTS_H
