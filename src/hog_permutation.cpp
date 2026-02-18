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
// where x is the observed overlap and E is the hypergeometric expectation.
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
#include <vector>
#include <random>
#include <cstdint>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// ---- Portable popcount ----

#ifdef _MSC_VER
#include <intrin.h>
static inline int popcnt64(uint64_t x) {
    return static_cast<int>(__popcnt64(x));
}
#else
static inline int popcnt64(uint64_t x) {
    return __builtin_popcountll(x);
}
#endif


// ---- Bit-vector operations ----

static int bv_and_popcount(const uint64_t* a, const uint64_t* b, int n_words) {
    int count = 0;
    for (int w = 0; w < n_words; ++w) {
        count += popcnt64(a[w] & b[w]);
    }
    return count;
}

static inline void bv_set(uint64_t* vec, int bit) {
    vec[bit >> 6] |= (1ULL << (bit & 63));
}


// ---- Sampling ----

// Sample k distinct integers from [0, n). Efficient when k << n.
static void sample_k_distinct(int n, int k, std::mt19937& rng,
                               std::vector<int>& out) {
    out.resize(k);
    std::uniform_int_distribution<int> dist(0, n - 1);
    for (int i = 0; i < k; ++i) {
        int x;
        bool dup;
        do {
            x = dist(rng);
            dup = false;
            for (int j = 0; j < i; ++j) {
                if (out[j] == x) { dup = true; break; }
            }
        } while (dup);
        out[i] = x;
    }
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

    // Direction 1: anchor = net1. For each sp2 gene, set reachable once.
    for (int b : sp2_genes) {
        int k1 = reach1_sz[b];
        if (k1 == 0) continue;
        const uint64_t* r1 = &reach1_bv[static_cast<size_t>(b) * n1w];
        for (int a : sp1_genes) {
            int m1 = neigh1_sz[a];
            if (m1 == 0) continue;
            const uint64_t* nb1 = &neigh1_bv[static_cast<size_t>(a) * n1w];
            int x1 = bv_and_popcount(nb1, r1, n1w);
            double E1 = static_cast<double>(m1) * k1 / n1;
            T += x1 / E1;
        }
    }

    // Direction 2: anchor = net2. For each sp1 gene, set reachable once.
    for (int a : sp1_genes) {
        int k2 = reach2_sz[a];
        if (k2 == 0) continue;
        const uint64_t* r2 = &reach2_bv[static_cast<size_t>(a) * n2w];
        for (int b : sp2_genes) {
            int m2 = neigh2_sz[b];
            if (m2 == 0) continue;
            const uint64_t* nb2 = &neigh2_bv[static_cast<size_t>(b) * n2w];
            int x2 = bv_and_popcount(nb2, r2, n2w);
            double E2 = static_cast<double>(m2) * k2 / n2;
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

    // Direction 1: for each sp2 gene, set reachable flags once, test all sp1
    for (int b : sp2_genes) {
        const auto& r1 = reachable1[b];
        int k1 = static_cast<int>(r1.size());
        if (k1 == 0) continue;

        for (int x : r1) flags1[x] = 1;

        for (int a : sp1_genes) {
            int m1 = static_cast<int>(neighbors1[a].size());
            if (m1 == 0) continue;
            int x1 = 0;
            for (int x : neighbors1[a]) {
                if (flags1[x]) ++x1;
            }
            double E1 = static_cast<double>(m1) * k1 / n1;
            T += x1 / E1;
        }

        for (int x : r1) flags1[x] = 0;
    }

    // Direction 2: for each sp1 gene, set reachable flags once, test all sp2
    for (int a : sp1_genes) {
        const auto& r2 = reachable2[a];
        int k2 = static_cast<int>(r2.size());
        if (k2 == 0) continue;

        for (int x : r2) flags2[x] = 1;

        for (int b : sp2_genes) {
            int m2 = static_cast<int>(neighbors2[b].size());
            if (m2 == 0) continue;
            int x2 = 0;
            for (int x : neighbors2[b]) {
                if (flags2[x]) ++x2;
            }
            double E2 = static_cast<double>(m2) * k2 / n2;
            T += x2 / E2;
        }

        for (int x : r2) flags2[x] = 0;
    }

    return T;
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
    const int n1 = static_cast<int>(net1.n_rows);
    const int n2 = static_cast<int>(net2.n_rows);
    const int n_ortho = ortho_sp1_idx.size();
    const int n_hogs = hog_sp1_list.size();

    // ---- Phase 1: Build ortholog mappings ----
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

    // ---- Phase 2: Compute neighbor lists ----
    // Use colptr() for cache-friendly access. Since the network matrices are
    // symmetric, column i has the same values as row i. Column access is
    // sequential in Armadillo's column-major layout.
    std::vector<std::vector<int>> neighbors1(n1);
    std::vector<std::vector<int>> neighbors2(n2);

#ifdef _OPENMP
    if (n_cores > 1) omp_set_num_threads(n_cores);
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int i = 0; i < n1; ++i) {
        const double* col_i = net1.colptr(i);
        for (int j = 0; j < n1; ++j) {
            if (i != j && col_i[j] >= thr1) {
                neighbors1[i].push_back(j);
            }
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int i = 0; i < n2; ++i) {
        const double* col_i = net2.colptr(i);
        for (int j = 0; j < n2; ++j) {
            if (i != j && col_i[j] >= thr2) {
                neighbors2[i].push_back(j);
            }
        }
    }

    // ---- Phase 3: Compute reachable sets ----
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

    // ---- Phase 4: Precompute sizes ----
    std::vector<int> neigh1_sz(n1), neigh2_sz(n2);
    std::vector<int> reach1_sz(n2), reach2_sz(n1);
    for (int i = 0; i < n1; ++i) neigh1_sz[i] = static_cast<int>(neighbors1[i].size());
    for (int j = 0; j < n2; ++j) neigh2_sz[j] = static_cast<int>(neighbors2[j].size());
    for (int b = 0; b < n2; ++b) reach1_sz[b] = static_cast<int>(reachable1[b].size());
    for (int a = 0; a < n1; ++a) reach2_sz[a] = static_cast<int>(reachable2[a].size());

    // ---- Phase 5: Choose intersection mode and build bit-vectors ----
    bool use_bitvec = (std::max(n1, n2) <= 100000);
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
                bv_set(&neigh1_bv[static_cast<size_t>(i) * n1w], j);

        for (int b = 0; b < n2; ++b)
            for (int a : reachable1[b])
                bv_set(&reach1_bv[static_cast<size_t>(b) * n1w], a);

        for (int j = 0; j < n2; ++j)
            for (int k : neighbors2[j])
                bv_set(&neigh2_bv[static_cast<size_t>(j) * n2w], k);

        for (int a = 0; a < n1; ++a)
            for (int b : reachable2[a])
                bv_set(&reach2_bv[static_cast<size_t>(a) * n2w], b);
    }

    // ---- Phase 6: Prepare per-thread resources ----
    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif

    std::vector<uint32_t> seeds(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        seeds[t] = static_cast<uint32_t>(R::runif(0.0, 4294967296.0));
    }

    std::vector<std::mt19937> thread_rng(max_threads);
    for (int t = 0; t < max_threads; ++t) thread_rng[t].seed(seeds[t]);

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

    // ---- Phase 7: HOG-level permutation tests ----
    Rcpp::NumericVector out_T_obs(n_hogs);
    Rcpp::IntegerVector out_n_perm(n_hogs);
    Rcpp::IntegerVector out_n_exceed(n_hogs);
    Rcpp::NumericVector out_p_value(n_hogs);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int h = 0; h < n_hogs; ++h) {
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

        // Compute observed statistic
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

        // Besag & Clifford (1991) adaptive stopping permutation loop
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
