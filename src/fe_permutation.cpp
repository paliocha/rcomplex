// fe_permutation.cpp
// Permutation-based HOG test using precomputed fold-enrichment matrix
//
// Takes a precomputed combined FE matrix where combined[a, b] sums the
// fold-enrichments from both network directions:
//   combined[a, b] = |N1(a) ∩ R1(b)| / E1(a,b) + |N2(b) ∩ R2(a)| / E2(a,b)
//
// Each permutation reduces to a submatrix sum — O(M×N) per permutation
// instead of O(M×N×n_words) with bit-vector intersections. The FE matrix
// is computed upstream (typically on GPU via torch matmuls).

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sample_k_distinct.h"

using namespace Rcpp;


// Submatrix sum: sum of combined[rows, cols].
// Iterates columns in outer loop for sequential access in column-major layout.
static double submatrix_sum(const arma::mat& combined,
                            const std::vector<int>& rows,
                            const std::vector<int>& cols) {
    double s = 0.0;
    for (int b : cols) {
        const double* col = combined.colptr(b);
        for (int a : rows) {
            s += col[a];
        }
    }
    return s;
}


//' Permutation test using precomputed fold-enrichment matrix
//'
//' @param combined Precomputed combined FE matrix (n1 x n2), where
//'   combined\[a, b\] = FE1(b->a) + FE2(a->b).
//' @param hog_sp1_list List of 0-based sp1 gene indices per HOG
//' @param hog_sp2_list List of 0-based sp2 gene indices per HOG
//' @param test_greater If TRUE, test conservation (T >= T_obs)
//' @param min_exceedances Besag-Clifford stopping parameter
//' @param max_permutations Maximum permutations per HOG
//' @param n_cores Number of OpenMP threads
//' @return DataFrame with T_obs, n_perm, n_exceed, p_value per HOG
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame fe_hog_permutation_test_cpp(
    const arma::mat& combined,
    const Rcpp::List& hog_sp1_list,
    const Rcpp::List& hog_sp2_list,
    bool test_greater,
    int min_exceedances,
    int max_permutations,
    int n_cores
) {
    const int n1 = static_cast<int>(combined.n_rows);
    const int n2 = static_cast<int>(combined.n_cols);
    const int n_hogs = hog_sp1_list.size();

    // Extract HOG gene lists
    std::vector<std::vector<int>> hog_sp1(n_hogs), hog_sp2(n_hogs);
    for (int h = 0; h < n_hogs; ++h) {
        Rcpp::IntegerVector v1 = hog_sp1_list[h];
        Rcpp::IntegerVector v2 = hog_sp2_list[h];
        hog_sp1[h].assign(v1.begin(), v1.end());
        hog_sp2[h].assign(v2.begin(), v2.end());
    }

    // Per-thread RNG seeded from R's RNG (before parallel region)
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

    Rcpp::NumericVector out_T_obs(n_hogs);
    Rcpp::IntegerVector out_n_perm(n_hogs);
    Rcpp::IntegerVector out_n_exceed(n_hogs);
    Rcpp::NumericVector out_p_value(n_hogs);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int h = 0; h < n_hogs; ++h) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& rng = thread_rng[tid];
        auto& perm1 = thread_perm1[tid];
        auto& perm2 = thread_perm2[tid];

        const auto& sp1 = hog_sp1[h];
        const auto& sp2 = hog_sp2[h];
        int M = static_cast<int>(sp1.size());
        int N = static_cast<int>(sp2.size());

        if (M == 0 || N == 0 || M > n1 || N > n2) {
            out_T_obs[h] = 0.0;
            out_n_perm[h] = 0;
            out_n_exceed[h] = 0;
            out_p_value[h] = 1.0;
            continue;
        }

        double T_obs = submatrix_sum(combined, sp1, sp2);

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
            sample_k_distinct(n1, M, rng, perm1);
            sample_k_distinct(n2, N, rng, perm2);

            double T_perm = submatrix_sum(combined, perm1, perm2);

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
        Named("T_obs") = out_T_obs,
        Named("n_perm") = out_n_perm,
        Named("n_exceed") = out_n_exceed,
        Named("p_value") = out_p_value
    );
}
