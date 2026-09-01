// hog_permutation.cpp
// Permutation-based HOG-level test for co-expression conservation
//
// Tests whether genes in a Hierarchical Ortholog Group (HOG) show greater
// cross-species co-expression conservation than expected by chance, using a
// gene-identity permutation null:
//
//   H0: Replacing this HOG's genes with randomly chosen genes from the
//       same networks yields equally large neighborhood overlap.
//
// Test statistic: Sum of fold-enrichments across all pair x direction
//   T = sum_{i,j} [x1_ij/E1_ij + x2_ij/E2_ij]
// where x is the observed overlap and E = m * k / (n - 1) is the
// hypergeometric expectation under the self-excluded urn: the anchor gene
// leaves the ortholog-reachable set (k) and the population (n - 1).
//
// Intersection modes:
// - Bit-vector with popcount when max(n1, n2) <= 100,000
// - Flag-vector (sparse set/count/clear) when max(n1, n2) > 100,000
//
// Adaptive stopping: Besag & Clifford (1991) sequential Monte Carlo p-values.
// Permutations stop once min_exceedances permutation statistics exceed T_obs,
// yielding p = (n_exceed + 1) / (n_perm + 1).
//
// References:
// - Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
//   Biometrika, 78(2), 301-304. doi:10.1093/biomet/78.2.301

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <numeric>
#include <span>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "neighbor_lists.h"
#include "sample_k_distinct.h"

using namespace Rcpp;


// ---- Bit-vector helpers ----

// Row view into a flat bit-vector array: row i of width nw words.
static std::span<const uint64_t> bv_row(const std::vector<uint64_t>& bv,
                                        int i, int nw) {
    return {&bv[static_cast<size_t>(i) * nw], static_cast<size_t>(nw)};
}

static std::span<uint64_t> bv_row_mut(std::vector<uint64_t>& bv,
                                      int i, int nw) {
    return {&bv[static_cast<size_t>(i) * nw], static_cast<size_t>(nw)};
}

static int bv_and_popcount(std::span<const uint64_t> a,
                           std::span<const uint64_t> b) {
    int count = 0;
    for (size_t w = 0; w < a.size(); ++w) {
        count += std::popcount(a[w] & b[w]);
    }
    return count;
}

static inline void bv_set(std::span<uint64_t> vec, int bit) {
    vec[bit >> 6] |= (1ULL << (bit & 63));
}


// ---- Fold-enrichment computation (bit-vector mode) ----

static double compute_T_bitvec(
    const std::vector<int>& sp1_genes,
    const std::vector<int>& sp2_genes,
    const std::vector<uint64_t>& neigh1_bv,
    const std::vector<uint64_t>& reach1_bv,
    const std::vector<uint64_t>& neigh2_bv,
    const std::vector<uint64_t>& reach2_bv,
    const std::vector<int>& neigh1_sz,
    const std::vector<int>& reach1_sz,
    const std::vector<int>& neigh2_sz,
    const std::vector<int>& reach2_sz,
    int n1, int n2, int n1w, int n2w
) {
    double T = 0.0;

    // Direction 1: anchor = net1. The anchor a leaves the reachable set
    // (k1 drops by one when a is reachable) and the population (n1 - 1);
    // the overlap is unaffected because a is never its own neighbour.
    for (int b : sp2_genes) {
        int k1_all = reach1_sz[b];
        if (k1_all == 0) continue;
        auto r1 = bv_row(reach1_bv, b, n1w);
        for (int a : sp1_genes) {
            int m1 = neigh1_sz[a];
            if (m1 == 0) continue;
            int k1 = k1_all - static_cast<int>((r1[a >> 6] >> (a & 63)) & 1ULL);
            if (k1 == 0) continue;
            int x1 = bv_and_popcount(bv_row(neigh1_bv, a, n1w), r1);
            double E1 = static_cast<double>(m1) * k1 / (n1 - 1);
            T += x1 / E1;
        }
    }

    // Direction 2: anchor = net2
    for (int a : sp1_genes) {
        int k2_all = reach2_sz[a];
        if (k2_all == 0) continue;
        auto r2 = bv_row(reach2_bv, a, n2w);
        for (int b : sp2_genes) {
            int m2 = neigh2_sz[b];
            if (m2 == 0) continue;
            int k2 = k2_all - static_cast<int>((r2[b >> 6] >> (b & 63)) & 1ULL);
            if (k2 == 0) continue;
            int x2 = bv_and_popcount(bv_row(neigh2_bv, b, n2w), r2);
            double E2 = static_cast<double>(m2) * k2 / (n2 - 1);
            T += x2 / E2;
        }
    }

    return T;
}


// ---- Fold-enrichment computation (flag-vector mode) ----

static double compute_T_flags(
    const std::vector<int>& sp1_genes,
    const std::vector<int>& sp2_genes,
    const std::vector<std::vector<int>>& neighbors1,
    const std::vector<std::vector<int>>& reachable1,
    const std::vector<std::vector<int>>& neighbors2,
    const std::vector<std::vector<int>>& reachable2,
    int n1, int n2,
    std::vector<char>& flags1,
    std::vector<char>& flags2
) {
    double T = 0.0;

    // Direction 1: anchor = net1 (self-excluded urn, see compute_T_bitvec)
    for (int b : sp2_genes) {
        const auto& r1 = reachable1[b];
        int k1_all = static_cast<int>(r1.size());
        if (k1_all == 0) continue;

        for (int x : r1) flags1[x] = 1;

        for (int a : sp1_genes) {
            int m1 = static_cast<int>(neighbors1[a].size());
            if (m1 == 0) continue;
            int k1 = k1_all - (flags1[a] ? 1 : 0);
            if (k1 == 0) continue;
            int x1 = 0;
            for (int x : neighbors1[a]) {
                if (flags1[x]) ++x1;
            }
            double E1 = static_cast<double>(m1) * k1 / (n1 - 1);
            T += x1 / E1;
        }

        for (int x : r1) flags1[x] = 0;
    }

    // Direction 2: anchor = net2
    for (int a : sp1_genes) {
        const auto& r2 = reachable2[a];
        int k2_all = static_cast<int>(r2.size());
        if (k2_all == 0) continue;

        for (int x : r2) flags2[x] = 1;

        for (int b : sp2_genes) {
            int m2 = static_cast<int>(neighbors2[b].size());
            if (m2 == 0) continue;
            int k2 = k2_all - (flags2[b] ? 1 : 0);
            if (k2 == 0) continue;
            int x2 = 0;
            for (int x : neighbors2[b]) {
                if (flags2[x]) ++x2;
            }
            double E2 = static_cast<double>(m2) * k2 / (n2 - 1);
            T += x2 / E2;
        }

        for (int x : r2) flags2[x] = 0;
    }

    return T;
}


// Shared body for both entry points: everything after neighbour-list
// construction. neighbors1[i] / neighbors2[j] are ascending 0-based lists
// (see neighbor_lists.h).
static Rcpp::DataFrame hog_permutation_test_core(
    const std::vector<std::vector<int>>& neighbors1,
    const std::vector<std::vector<int>>& neighbors2,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    const Rcpp::List& hog_sp1_list,
    const Rcpp::List& hog_sp2_list,
    bool test_greater,
    int min_exceedances,
    int max_permutations,
    int n_cores
) {
    const int n1 = static_cast<int>(neighbors1.size());
    const int n2 = static_cast<int>(neighbors2.size());
    const int n_ortho = ortho_sp1_idx.size();
    const int n_hogs = hog_sp1_list.size();

    // ---- Build ortholog mappings ----
    std::vector<std::vector<int>> sp2_to_sp1(n2);
    std::vector<std::vector<int>> sp1_to_sp2(n1);

    for (int i = 0; i < n_ortho; ++i) {
        int s1 = ortho_sp1_idx[i];
        int s2 = ortho_sp2_idx[i];
        if (s1 >= 0 && s1 < n1 && s2 >= 0 && s2 < n2) {
            sp2_to_sp1[s2].push_back(s1);
            sp1_to_sp2[s1].push_back(s2);
        }
    }

    // ---- Compute reachable sets ----
    // reachable1[b] = sp1 genes reachable from sp2 gene b's net2 neighbors
    // reachable2[a] = sp2 genes reachable from sp1 gene a's net1 neighbors
    std::vector<std::vector<int>> reachable1(n2);
    std::vector<std::vector<int>> reachable2(n1);

    {
        std::vector<char> seen(std::max(n1, n2), 0);
        for (int b = 0; b < n2; ++b) {
            for (int nb : neighbors2[b]) {
                for (int a : sp2_to_sp1[nb]) {
                    if (!seen[a]) {
                        seen[a] = 1;
                        reachable1[b].push_back(a);
                    }
                }
            }
            for (int a : reachable1[b]) seen[a] = 0;
        }
        for (int a = 0; a < n1; ++a) {
            for (int nb : neighbors1[a]) {
                for (int b : sp1_to_sp2[nb]) {
                    if (!seen[b]) {
                        seen[b] = 1;
                        reachable2[a].push_back(b);
                    }
                }
            }
            for (int b : reachable2[a]) seen[b] = 0;
        }
    }

    // ---- Precompute sizes ----
    std::vector<int> neigh1_sz(n1), neigh2_sz(n2);
    std::vector<int> reach1_sz(n2), reach2_sz(n1);
    for (int i = 0; i < n1; ++i) neigh1_sz[i] = static_cast<int>(neighbors1[i].size());
    for (int j = 0; j < n2; ++j) neigh2_sz[j] = static_cast<int>(neighbors2[j].size());
    for (int b = 0; b < n2; ++b) reach1_sz[b] = static_cast<int>(reachable1[b].size());
    for (int a = 0; a < n1; ++a) reach2_sz[a] = static_cast<int>(reachable2[a].size());

    // ---- Choose intersection mode and build bit-vectors ----
    bool use_bitvec = (std::max(n1, n2) <= 100000);
    if (!use_bitvec) {
        REprintf("Network size > 100K genes: using flag-vector mode "
                 "(slower than bit-vector; expected for large genomes)\n");
    }
    int n1w = (n1 + 63) / 64;
    int n2w = (n2 + 63) / 64;
    std::vector<uint64_t> neigh1_bv, reach1_bv, neigh2_bv, reach2_bv;

    if (use_bitvec) {
        neigh1_bv.assign(static_cast<size_t>(n1) * n1w, 0);
        reach1_bv.assign(static_cast<size_t>(n2) * n1w, 0);
        neigh2_bv.assign(static_cast<size_t>(n2) * n2w, 0);
        reach2_bv.assign(static_cast<size_t>(n1) * n2w, 0);

        for (int i = 0; i < n1; ++i)
            for (int j : neighbors1[i])
                bv_set(bv_row_mut(neigh1_bv, i, n1w), j);

        for (int b = 0; b < n2; ++b)
            for (int a : reachable1[b])
                bv_set(bv_row_mut(reach1_bv, b, n1w), a);

        for (int j = 0; j < n2; ++j)
            for (int k : neighbors2[j])
                bv_set(bv_row_mut(neigh2_bv, j, n2w), k);

        for (int a = 0; a < n1; ++a)
            for (int b : reachable2[a])
                bv_set(bv_row_mut(reach2_bv, a, n2w), b);
    }

    // ---- Prepare per-thread resources ----
    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif

    std::vector<std::mt19937> thread_rng(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        thread_rng[t].seed(
            static_cast<uint32_t>(R::runif(0.0, 4294967296.0)));
    }

    std::vector<std::vector<int>> thread_perm1(max_threads);
    std::vector<std::vector<int>> thread_perm2(max_threads);
    std::vector<std::vector<char>> thread_f1(max_threads);
    std::vector<std::vector<char>> thread_f2(max_threads);
    if (!use_bitvec) {
        for (int t = 0; t < max_threads; ++t) {
            thread_f1[t].assign(n1, 0);
            thread_f2[t].assign(n2, 0);
        }
    }

    // Extract HOG gene lists from R List objects
    std::vector<std::vector<int>> hog_sp1(n_hogs);
    std::vector<std::vector<int>> hog_sp2(n_hogs);
    for (int h = 0; h < n_hogs; ++h) {
        Rcpp::IntegerVector v1 = hog_sp1_list[h];
        Rcpp::IntegerVector v2 = hog_sp2_list[h];
        hog_sp1[h].assign(v1.begin(), v1.end());
        hog_sp2[h].assign(v2.begin(), v2.end());
    }

    // ---- HOG-level permutation tests ----
    // Sort HOGs by descending total gene count so large HOGs (most work)
    // are scheduled first, reducing tail latency with guided scheduling
    std::vector<int> hog_order(n_hogs);
    std::iota(hog_order.begin(), hog_order.end(), 0);
    std::sort(hog_order.begin(), hog_order.end(), [&](int a, int b) {
        return (hog_sp1[a].size() + hog_sp2[a].size()) >
               (hog_sp1[b].size() + hog_sp2[b].size());
    });

    Rcpp::NumericVector out_T_obs(n_hogs);
    Rcpp::IntegerVector out_n_perm(n_hogs);
    Rcpp::IntegerVector out_n_exceed(n_hogs);
    Rcpp::NumericVector out_p_value(n_hogs);

#ifdef _OPENMP
    #pragma omp parallel for schedule(guided) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int hi = 0; hi < n_hogs; ++hi) {
        int h = hog_order[hi];
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& rng = thread_rng[tid];
        auto& perm_sp1 = thread_perm1[tid];
        auto& perm_sp2 = thread_perm2[tid];

        int M = static_cast<int>(hog_sp1[h].size());
        int N_sp2 = static_cast<int>(hog_sp2[h].size());

        if (M == 0 || N_sp2 == 0 || M > n1 || N_sp2 > n2) {
            out_T_obs[h] = 0.0;
            out_n_perm[h] = 0;
            out_n_exceed[h] = 0;
            out_p_value[h] = 1.0;
            continue;
        }

        double T_obs;
        if (use_bitvec) {
            T_obs = compute_T_bitvec(
                hog_sp1[h], hog_sp2[h],
                neigh1_bv, reach1_bv, neigh2_bv, reach2_bv,
                neigh1_sz, reach1_sz, neigh2_sz, reach2_sz,
                n1, n2, n1w, n2w);
        } else {
            T_obs = compute_T_flags(
                hog_sp1[h], hog_sp2[h],
                neighbors1, reachable1, neighbors2, reachable2,
                n1, n2, thread_f1[tid], thread_f2[tid]);
        }

        if (test_greater && T_obs <= 0.0) {
            out_T_obs[h] = 0.0;
            out_n_perm[h] = 0;
            out_n_exceed[h] = 0;
            out_p_value[h] = 1.0;
            continue;
        }

        // Besag & Clifford adaptive stopping
        int n_exceed = 0;
        int n_perm = 0;

        while (n_perm < max_permutations && n_exceed < min_exceedances) {
            sample_k_distinct(n1, M, rng, perm_sp1);
            sample_k_distinct(n2, N_sp2, rng, perm_sp2);

            double T_perm;
            if (use_bitvec) {
                T_perm = compute_T_bitvec(
                    perm_sp1, perm_sp2,
                    neigh1_bv, reach1_bv, neigh2_bv, reach2_bv,
                    neigh1_sz, reach1_sz, neigh2_sz, reach2_sz,
                    n1, n2, n1w, n2w);
            } else {
                T_perm = compute_T_flags(
                    perm_sp1, perm_sp2,
                    neighbors1, reachable1, neighbors2, reachable2,
                    n1, n2, thread_f1[tid], thread_f2[tid]);
            }

            if (test_greater) {
                if (T_perm >= T_obs) ++n_exceed;
            } else {
                if (T_perm <= T_obs) ++n_exceed;
            }
            ++n_perm;
        }

        out_T_obs[h] = T_obs;
        out_n_perm[h] = n_perm;
        out_n_exceed[h] = n_exceed;
        out_p_value[h] = static_cast<double>(n_exceed + 1) / (n_perm + 1);
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("T_obs") = out_T_obs,
        Rcpp::Named("n_perm") = out_n_perm,
        Rcpp::Named("n_exceed") = out_n_exceed,
        Rcpp::Named("p_value") = out_p_value
    );
}


//' Permutation-based HOG-level conservation test
//'
//' Tests each HOG for co-expression conservation using a gene-identity
//' permutation null with adaptive stopping (Besag & Clifford, 1991).
//'
//' @param net1 Co-expression network matrix for species 1 (n1 x n1)
//' @param net2 Co-expression network matrix for species 2 (n2 x n2)
//' @param thr1 Co-expression threshold for species 1
//' @param thr2 Co-expression threshold for species 2
//' @param ortho_sp1_idx 0-based net1 indices for full ortholog table
//' @param ortho_sp2_idx 0-based net2 indices for full ortholog table
//' @param hog_sp1_list List of integer vectors: unique 0-based sp1 indices per HOG
//' @param hog_sp2_list List of integer vectors: unique 0-based sp2 indices per HOG
//' @param test_greater If TRUE, test conservation (T >= T_obs); if FALSE, divergence
//' @param min_exceedances Besag-Clifford stopping parameter (default 50)
//' @param max_permutations Maximum permutations per HOG (default 10000)
//' @param n_cores Number of OpenMP threads (default 1)
//' @return DataFrame with T_obs, n_perm, n_exceed, p_value per HOG
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame hog_permutation_test_cpp(
    const arma::mat& net1,
    const arma::mat& net2,
    double thr1,
    double thr2,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    const Rcpp::List& hog_sp1_list,
    const Rcpp::List& hog_sp2_list,
    bool test_greater,
    int min_exceedances,
    int max_permutations,
    int n_cores
) {
    return hog_permutation_test_core(
        neighbor_lists_dense(net1, thr1, n_cores),
        neighbor_lists_dense(net2, thr2, n_cores),
        ortho_sp1_idx, ortho_sp2_idx, hog_sp1_list, hog_sp2_list,
        test_greater, min_exceedances, max_permutations, n_cores);
}


//' Permutation-based HOG-level conservation test (sparse networks)
//'
//' Same as [hog_permutation_test_cpp()] but takes the slots of a
//' `dgCMatrix` (column-compressed, both triangles stored) for each network
//' instead of a dense matrix. Column j lists the neighbours of gene j.
//'
//' @param p1 `@p` slot of net1 (column pointers, length n1 + 1)
//' @param i1 `@i` slot of net1 (0-based row indices)
//' @param x1 `@x` slot of net1 (stored values)
//' @param thr1 Co-expression threshold for species 1
//' @param p2 `@p` slot of net2
//' @param i2 `@i` slot of net2
//' @param x2 `@x` slot of net2
//' @param thr2 Co-expression threshold for species 2
//' @param ortho_sp1_idx 0-based net1 indices for full ortholog table
//' @param ortho_sp2_idx 0-based net2 indices for full ortholog table
//' @param hog_sp1_list List of integer vectors: unique 0-based sp1 indices per HOG
//' @param hog_sp2_list List of integer vectors: unique 0-based sp2 indices per HOG
//' @param test_greater If TRUE, test conservation (T >= T_obs); if FALSE, divergence
//' @param min_exceedances Besag-Clifford stopping parameter (default 50)
//' @param max_permutations Maximum permutations per HOG (default 10000)
//' @param n_cores Number of OpenMP threads (default 1)
//' @return DataFrame with T_obs, n_perm, n_exceed, p_value per HOG
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame hog_permutation_test_sparse_cpp(
    const Rcpp::IntegerVector& p1,
    const Rcpp::IntegerVector& i1,
    const Rcpp::NumericVector& x1,
    double thr1,
    const Rcpp::IntegerVector& p2,
    const Rcpp::IntegerVector& i2,
    const Rcpp::NumericVector& x2,
    double thr2,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    const Rcpp::List& hog_sp1_list,
    const Rcpp::List& hog_sp2_list,
    bool test_greater,
    int min_exceedances,
    int max_permutations,
    int n_cores
) {
    return hog_permutation_test_core(
        neighbor_lists_sparse(p1, i1, x1, thr1, n_cores),
        neighbor_lists_sparse(p2, i2, x2, thr2, n_cores),
        ortho_sp1_idx, ortho_sp2_idx, hog_sp1_list, hog_sp2_list,
        test_greater, min_exceedances, max_permutations, n_cores);
}
