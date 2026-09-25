// network_block.cpp
// Sparse mutual-rank network built block by block of columns, never
// holding the n x n matrix; reproduces mutual_rank_inplace_cpp() +
// density_threshold_cpp() + extract_sparse_cpp() exactly.
//
// A pair (i, j) scores sqrt(R_ij * R_ji) and both ranks are at most n, so
// a pair at or above threshold T has min(R_ij, R_ji) >= T^2 / n (raw MR).
// Pass 1 records, per column, the correlation cut at a grid of rank
// fractions f. Pass 2 keeps only pairs in the top fraction f of BOTH
// columns (exact rank in its own column, loosened cut in the other), joins
// the two orientations, and takes the density thresholds among them. If
// the store threshold is too low to prove every qualifying pair was a
// candidate, f widens; at f = 1 every pair is a candidate.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <vector>

#include "density_k.h"
#include "rank_column.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct Trip {
    int lo, hi, col;
    double r;
};

// Same clamp / abs / NaN handling as mutual_rank_inplace_cpp(); returns
// false on NaN (column left partially transformed).
bool prep_column(double* col, R_xlen_t n, bool abs_cor) {
    for (R_xlen_t r = 0; r < n; ++r) {
        if (std::isnan(col[r])) return false;
        const double v = std::clamp(col[r], -1.0, 1.0);
        col[r] = abs_cor ? std::fabs(v) : v;
    }
    return true;
}

}  // namespace

//' Blockwise sparse mutual-rank network
//'
//' @param zt Standardised expression (samples x genes) such that
//'   `crossprod(zt)` is the correlation matrix.
//' @param log_transform,abs_cor As in [mutual_rank_inplace_cpp()].
//' @param density,store_density Analysis and store densities.
//' @param block_size Columns per correlation block.
//' @param n_cores Number of OpenMP threads.
//' @return List with dgCMatrix slots `i`, `p`, `x`, the `threshold` and
//'   `store_threshold`, the number of joined candidate pairs
//'   `n_candidates`, and the grid `fraction` that sufficed.
//'
//' @keywords internal
// [[Rcpp::export]]
List mr_block_network_cpp(const arma::mat& zt, bool log_transform,
                          bool abs_cor, double density,
                          double store_density, int block_size,
                          int n_cores) {
    const R_xlen_t n = zt.n_cols;
    if (n < 3) stop("Fewer than 3 genes");
    if (!(density > 0.0 && density <= store_density && store_density < 1.0)) {
        stop("need 0 < density <= store_density < 1");
    }
    if (block_size < 1) stop("block_size must be >= 1");

    std::vector<double> grid;
    for (double m : {2.0, 4.0, 8.0, 16.0}) {
        const double f = m * store_density;
        if (f < 1.0 && (grid.empty() || f > grid.back())) grid.push_back(f);
    }
    if (grid.empty() || grid.back() < 0.5) grid.push_back(0.5);
    grid.push_back(1.0);
    const int n_grid = static_cast<int>(grid.size());
    const double dn = static_cast<double>(n);
    const double log_n = std::log(dn);
    const double eps = 64 * DBL_EPSILON;

    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif
    std::vector<std::vector<double>> bufs(max_threads,
                                          std::vector<double>(n));
    std::vector<std::vector<R_xlen_t>> idxs(max_threads,
                                            std::vector<R_xlen_t>(n));
    std::vector<char> nan_seen(max_threads, 0);

    // Ascending position of the per-column cut for each grid fraction.
    std::vector<R_xlen_t> cut_pos(n_grid);
    for (int g = 0; g < n_grid; ++g) {
        if (log_transform) {
            const auto d = std::max<R_xlen_t>(
                0, static_cast<R_xlen_t>(std::floor(dn * grid[g])) - 1);
            cut_pos[g] = n - 1 - d;
        } else {
            cut_pos[g] = std::max<R_xlen_t>(
                0, static_cast<R_xlen_t>(std::ceil(dn * (1.0 - grid[g]))) - 1);
        }
    }

    // Pass 1: per-column cuts (n x n_grid).
    arma::mat cuts(n, n_grid);
    for (R_xlen_t c0 = 0; c0 < n; c0 += block_size) {
        const R_xlen_t c1 = std::min<R_xlen_t>(n, c0 + block_size);
        arma::mat C = zt.t() * zt.cols(c0, c1 - 1);
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
        for (R_xlen_t b = 0; b < c1 - c0; ++b) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            double* col = C.colptr(b);
            if (!prep_column(col, n, abs_cor)) {
                nan_seen[tid] = 1;
                continue;
            }
            std::vector<double>& buf = bufs[tid];
            std::copy(col, col + n, buf.begin());
            for (int g = 0; g < n_grid; ++g) {
                std::nth_element(buf.begin(), buf.begin() + cut_pos[g],
                                 buf.end());
                cuts(c0 + b, g) = buf[cut_pos[g]];
            }
        }
        if (std::any_of(nan_seen.begin(), nan_seen.end(),
                        [](char x) { return x != 0; })) {
            stop("sim contains NaN");
        }
    }

    const std::size_t tri = static_cast<std::size_t>(n) * (n - 1) / 2;
    const std::size_t k_s = density_k(store_density, tri);
    const std::size_t k_d = density_k(density, tri);

    std::vector<int> plo, phi;
    std::vector<double> pv;
    double t_s = 0.0, t_d = 0.0;
    int g = 0;
    for (; g < n_grid; ++g) {
        const double f = grid[g];
        const double lim = log_transform ? dn * f : dn * (1.0 - f);
        std::vector<std::vector<Trip>> trips(max_threads);

        // Pass 2: exact ranks, emit candidate orientations.
        for (R_xlen_t c0 = 0; c0 < n; c0 += block_size) {
            const R_xlen_t c1 = std::min<R_xlen_t>(n, c0 + block_size);
            arma::mat C = zt.t() * zt.cols(c0, c1 - 1);
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
            for (R_xlen_t b = 0; b < c1 - c0; ++b) {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                const int i = static_cast<int>(c0 + b);
                double* col = C.colptr(b);
                prep_column(col, n, abs_cor);
                std::vector<double>& rk = bufs[tid];
                std::copy(col, col + n, rk.begin());
                rank_column_inplace(rk.data(), n, !log_transform, idxs[tid]);
                for (R_xlen_t j = 0; j < n; ++j) {
                    if (j == i) continue;
                    const bool in_rank = log_transform ? rk[j] <= lim
                                                       : rk[j] >= lim;
                    if (in_rank && col[j] >= cuts(j, g) - eps) {
                        const int jj = static_cast<int>(j);
                        trips[tid].push_back({std::min(i, jj),
                                              std::max(i, jj), i, rk[j]});
                    }
                }
            }
        }

        std::vector<Trip> all;
        for (auto& t : trips) {
            all.insert(all.end(), t.begin(), t.end());
            std::vector<Trip>().swap(t);
        }
        std::sort(all.begin(), all.end(), [](const Trip& a, const Trip& b) {
            if (a.lo != b.lo) return a.lo < b.lo;
            if (a.hi != b.hi) return a.hi < b.hi;
            return a.col < b.col;
        });

        // Join: both orientations present -> (lo, hi) pair, lower first.
        plo.clear();
        phi.clear();
        pv.clear();
        for (std::size_t a = 0; a + 1 < all.size(); ++a) {
            if (all[a].lo != all[a + 1].lo || all[a].hi != all[a + 1].hi) {
                continue;
            }
            double v = std::sqrt(all[a].r * all[a + 1].r);
            if (log_transform) {
                v = std::clamp(1.0 - std::log(v) / log_n, 0.0, 1.0);
            }
            plo.push_back(all[a].lo);
            phi.push_back(all[a].hi);
            pv.push_back(v);
            ++a;
        }
        if (pv.size() < k_s) continue;

        std::vector<double> tmp(pv);
        std::nth_element(tmp.begin(), tmp.end() - k_s, tmp.end());
        t_s = *(tmp.end() - k_s);
        std::nth_element(tmp.begin(), tmp.end() - k_d, tmp.end());
        t_d = *(tmp.end() - k_d);

        if (g == n_grid - 1) break;
        // ponytail: 1e-12 slack guards the rounding of T^2 / n^(2(1-T));
        // exact rational bounds if a boundary case ever bites.
        const bool valid = log_transform
            ? std::pow(dn, 2.0 * (1.0 - t_s)) * (1.0 + 1e-12) <= lim
            : t_s * t_s / dn * (1.0 - 1e-12) >= lim;
        if (valid) break;
    }

    // Store: kept pairs, both orientations, CSC. Pairs are sorted by
    // (lo, hi), so each column receives its rows in ascending order.
    IntegerVector p(n + 1);
    std::size_t nnz = 0;
    for (std::size_t a = 0; a < pv.size(); ++a) {
        if (pv[a] >= t_s) {
            ++p[plo[a] + 1];
            ++p[phi[a] + 1];
            nnz += 2;
        }
    }
    if (nnz > static_cast<std::size_t>(INT_MAX)) {
        stop("stored network exceeds 2^31 - 1 entries");
    }
    for (R_xlen_t c = 0; c < n; ++c) p[c + 1] += p[c];
    IntegerVector ri(nnz);
    NumericVector rx(nnz);
    std::vector<int> fill(p.begin(), p.end() - 1);
    for (std::size_t a = 0; a < pv.size(); ++a) {
        if (pv[a] >= t_s) {
            int k = fill[plo[a]]++;
            ri[k] = phi[a];
            rx[k] = pv[a];
            k = fill[phi[a]]++;
            ri[k] = plo[a];
            rx[k] = pv[a];
        }
    }

    return List::create(
        Named("i") = ri, Named("p") = p, Named("x") = rx,
        Named("threshold") = t_d, Named("store_threshold") = t_s,
        Named("n_candidates") = static_cast<double>(pv.size()),
        Named("fraction") = grid[std::min(g, n_grid - 1)]);
}
