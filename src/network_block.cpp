// network_block.cpp
// Sparse mutual-rank network built block by block of columns, never
// holding the n x n matrix; reproduces mutual_rank_inplace_cpp() +
// density_threshold_cpp() + extract_sparse_cpp() exactly.
//
// One correlation pass ranks every column exactly and keeps, per gene, the
// partners whose rank lies in the top fraction f of that column. A pair is
// a candidate when each gene is on the other's list. A pair (i, j) scores
// sqrt(R_ij * R_ji) and both ranks are at most n, so a pair at or above
// threshold T has min(R_ij, R_ji) >= T^2 / n (raw MR). If the store
// threshold is too low to prove every qualifying pair was a candidate, f
// widens by 1.5x; at f = 1 every pair is a candidate.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
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

// Average ranks are multiples of 0.5, exact in a float while n <= 2^23.
struct Entry {
    int j;
    float r;
};
using List_t = std::vector<Entry>;

// Rank of i in list l (sorted by j), or -1 if absent.
double find_rank(const List_t& l, int i) {
    auto it = std::lower_bound(
        l.begin(), l.end(), i,
        [](const Entry& e, int x) { return e.j < x; });
    return (it != l.end() && it->j == i) ? it->r : -1.0;
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
//'   `n_candidates`, and the rank `fraction` that sufficed.
//'
//' @keywords internal
// [[Rcpp::export]]
List mr_block_network_cpp(const arma::mat& zt, bool log_transform,
                          bool abs_cor, double density,
                          double store_density, int block_size,
                          int n_cores) {
    const R_xlen_t n = zt.n_cols;
    if (n < 3) stop("Fewer than 3 genes");
    if (n > (R_xlen_t{1} << 23)) stop("more than 2^23 genes");
    if (!(density > 0.0 && density <= store_density && store_density < 1.0)) {
        stop("need 0 < density <= store_density < 1");
    }
    if (block_size < 1) stop("block_size must be >= 1");

    const double dn = static_cast<double>(n);
    const double log_n = std::log(dn);
    const int ni = static_cast<int>(n);

    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif
    std::vector<std::vector<R_xlen_t>> idxs(max_threads,
                                            std::vector<R_xlen_t>(n));
    std::vector<char> nan_seen(max_threads, 0);

    const std::size_t tri = static_cast<std::size_t>(n) * (n - 1) / 2;
    const std::size_t k_s = density_k(store_density, tri);
    const std::size_t k_d = density_k(density, tri);

    // MR of (i, j) from list ranks, as mutual_rank_inplace_cpp() computes
    // it; -1 when the pair is not on both lists.
    std::vector<List_t> lists(n);
    auto mr = [&](int i, const Entry& e) {
        const double rji = find_rank(lists[e.j], i);
        if (rji < 0) return -1.0;
        double v = std::sqrt(static_cast<double>(e.r) * rji);
        if (log_transform) {
            v = std::clamp(1.0 - std::log(v) / log_n, 0.0, 1.0);
        }
        return v;
    };

    // At store_density 0.05 on BDIS leaf data the list fraction needed was
    // ~0.11, so the 3x start usually passes first time.
    double f = std::min(1.0, 3.0 * store_density);
    double t_s = 0.0, t_d = 0.0;
    std::size_t n_cand = 0;
    for (;;) {
        const double lim = log_transform ? dn * f : dn * (1.0 - f);

        // Single pass: exact ranks per column, keep the top-f partners.
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
                const R_xlen_t i = c0 + b;
                double* col = C.colptr(b);
                bool bad = false;
                for (R_xlen_t r = 0; r < n; ++r) {
                    if (std::isnan(col[r])) {
                        bad = true;
                        break;
                    }
                    const double v = std::clamp(col[r], -1.0, 1.0);
                    col[r] = abs_cor ? std::fabs(v) : v;
                }
                if (bad) {
                    nan_seen[tid] = 1;
                    continue;
                }
                rank_column_inplace(col, n, !log_transform, idxs[tid]);
                auto keep = [&](R_xlen_t j) {
                    return j != i &&
                        (log_transform ? col[j] <= lim : col[j] >= lim);
                };
                std::size_t cnt = 0;
                for (R_xlen_t j = 0; j < n; ++j) cnt += keep(j);
                List_t& l = lists[i];
                l.reserve(cnt);
                for (R_xlen_t j = 0; j < n; ++j) {
                    if (keep(j)) {
                        l.push_back({static_cast<int>(j),
                                     static_cast<float>(col[j])});
                    }
                }
            }
            if (std::any_of(nan_seen.begin(), nan_seen.end(),
                            [](char x) { return x != 0; })) {
                stop("sim contains NaN");
            }
        }

        // Join: candidate MR values of the upper triangle, counted first.
        std::vector<std::size_t> off(n + 1, 0);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 64) num_threads(n_cores) if(n_cores > 1)
#endif
        for (int i = 0; i < ni; ++i) {
            std::size_t cnt = 0;
            for (const Entry& e : lists[i]) {
                if (e.j > i && mr(i, e) >= 0) ++cnt;
            }
            off[i + 1] = cnt;
        }
        for (R_xlen_t i = 0; i < n; ++i) off[i + 1] += off[i];
        n_cand = off[n];

        if (n_cand >= k_s) {
            std::vector<double> cand(n_cand);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 64) num_threads(n_cores) if(n_cores > 1)
#endif
            for (int i = 0; i < ni; ++i) {
                std::size_t k = off[i];
                for (const Entry& e : lists[i]) {
                    if (e.j <= i) continue;
                    const double v = mr(i, e);
                    if (v >= 0) cand[k++] = v;
                }
            }
            std::nth_element(cand.begin(), cand.end() - k_s, cand.end());
            t_s = *(cand.end() - k_s);
            std::nth_element(cand.begin(), cand.end() - k_d, cand.end());
            t_d = *(cand.end() - k_d);

            // ponytail: 1e-12 slack guards the rounding of T^2 / n^(2(1-T));
            // exact rational bounds if a boundary case ever bites.
            const bool valid = f >= 1.0 || (log_transform
                ? std::pow(dn, 2.0 * (1.0 - t_s)) * (1.0 + 1e-12) <= lim
                : t_s * t_s / dn * (1.0 - 1e-12) >= lim);
            if (valid) break;
        }
        for (List_t& l : lists) List_t().swap(l);
        f = std::min(1.0, f * 1.5);
    }

    // Store: one upper-triangle sweep keeps the pairs at or above T_s, the
    // lists are freed, and the kept pairs are laid out in both orientations
    // as CSC. Walking i upward, column j receives its rows below j in
    // ascending order, then its own kept pairs (rows above j, sorted).
    struct Kept {
        int j;
        double v;
    };
    std::vector<std::vector<Kept>> kept(n);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int i = 0; i < ni; ++i) {
        for (const Entry& e : lists[i]) {
            if (e.j <= i) continue;
            const double v = mr(i, e);
            if (v >= t_s) kept[i].push_back({e.j, v});
        }
    }
    std::vector<List_t>().swap(lists);

    IntegerVector p(n + 1);
    std::size_t nnz = 0;
    for (int i = 0; i < ni; ++i) {
        for (const Kept& k : kept[i]) {
            ++p[i + 1];
            ++p[k.j + 1];
        }
        nnz += 2 * kept[i].size();
    }
    if (nnz > static_cast<std::size_t>(INT_MAX)) {
        stop("stored network exceeds 2^31 - 1 entries");
    }
    for (R_xlen_t c = 0; c < n; ++c) p[c + 1] += p[c];
    IntegerVector ri(nnz);
    NumericVector rx(nnz);
    std::vector<int> fill(p.begin(), p.end() - 1);
    for (int i = 0; i < ni; ++i) {
        for (const Kept& k : kept[i]) {
            int a = fill[i]++;
            ri[a] = k.j;
            rx[a] = k.v;
            a = fill[k.j]++;
            ri[a] = i;
            rx[a] = k.v;
        }
        std::vector<Kept>().swap(kept[i]);
    }

    return List::create(
        Named("i") = ri, Named("p") = p, Named("x") = rx,
        Named("threshold") = t_d, Named("store_threshold") = t_s,
        Named("n_candidates") = static_cast<double>(n_cand),
        Named("fraction") = f);
}
