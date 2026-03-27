// coclassification.cpp
// Build co-classification matrix and compute adaptive threshold for
// consensus module detection (Jeub et al., 2018).
//
// For each of K resolution partitions, counts how often each gene pair
// is assigned to the same module. Returns normalized co-classification
// frequencies and per-resolution expected co-classification under random
// assignment (for adaptive thresholding).

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

//' Build co-classification matrix and compute adaptive threshold
//'
//' For each of K resolution partitions, counts how often each gene pair
//' is assigned to the same module. Returns the co-classification frequency
//' matrix and per-resolution expected co-classification under random
//' assignment (for adaptive thresholding per Jeub et al. 2018).
//'
//' @param memberships List of K IntegerVectors (1-based module IDs from
//'   igraph, each length N)
//' @param n_genes Number of genes (N)
//' @return List with: coclassification (N x N matrix, entries in
//'   0 to 1 range), expected (length-K numeric vector of per-resolution
//'   expected co-classification), adaptive_threshold (scalar = mean of
//'   expected)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List build_coclassification_cpp(
    const Rcpp::List& memberships,
    int n_genes)
{
    int K = memberships.size();
    arma::mat C(n_genes, n_genes, arma::fill::zeros);
    Rcpp::NumericVector expected_vec(K);

    double inv_n = 1.0 / n_genes;

    for (int k = 0; k < K; ++k) {
        if (k % 10 == 0) Rcpp::checkUserInterrupt();

        Rcpp::IntegerVector mem = memberships[k];

        // Group genes by module (1-based IDs from igraph)
        int max_mod = 0;
        for (int i = 0; i < n_genes; ++i)
            if (mem[i] > max_mod) max_mod = mem[i];

        std::vector<std::vector<int>> modules(max_mod + 1);
        for (int i = 0; i < n_genes; ++i)
            modules[mem[i]].push_back(i);

        // Expected co-classification: sum_m (s_m/N)^2
        double exp_k = 0.0;
        for (auto& mod : modules) {
            if (mod.empty()) continue;
            double frac = static_cast<double>(mod.size()) * inv_n;
            exp_k += frac * frac;
        }
        expected_vec[k] = exp_k;

        // Update co-classification matrix (column-major access)
        for (auto& mod : modules) {
            int s = static_cast<int>(mod.size());
            if (s < 2) {
                // Singleton: only diagonal
                if (s == 1) C(mod[0], mod[0]) += 1.0;
                continue;
            }
            for (int bi = 0; bi < s; ++bi) {
                double* col_b = C.colptr(mod[bi]);
                for (int ai = 0; ai < s; ++ai) {
                    col_b[mod[ai]] += 1.0;
                }
            }
        }
    }

    // Normalize to [0, 1]
    C /= static_cast<double>(K);

    // Adaptive threshold = mean expected co-classification
    double adaptive_tau = 0.0;
    for (int k = 0; k < K; ++k) adaptive_tau += expected_vec[k];
    adaptive_tau /= K;

    return Rcpp::List::create(
        Named("coclassification") = C,
        Named("expected") = expected_vec,
        Named("adaptive_threshold") = adaptive_tau
    );
}
