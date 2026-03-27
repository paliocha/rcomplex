// coclassification.cpp
// Co-classification matrix with per-pair null subtraction for iterative
// consensus module detection (Jeub et al., 2018).

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

//' Build co-classification matrix with per-pair excess
//'
//' For each of K resolution partitions, counts how often each gene pair
//' is assigned to the same module. Also accumulates the per-pair expected
//' co-classification under random assignment (Jeub et al. 2018):
//' for same-module pairs at resolution k, E_k(i,j) = (s_m/N)^2 where
//' s_m is the size of the module containing both i and j.
//'
//' @param memberships List of K IntegerVectors (1-based module IDs from
//'   igraph, each length N)
//' @param n_genes Number of genes (N)
//' @return List with: coclassification (N x N raw co-classification,
//'   values in 0 to 1), excess_coclassification (N x N, max(C - E, 0)),
//'   expected (length-K numeric vector of per-resolution expected scalars)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List build_coclassification_cpp(
    const Rcpp::List& memberships,
    int n_genes)
{
    int K = memberships.size();
    arma::mat C(n_genes, n_genes, arma::fill::zeros);
    arma::mat E(n_genes, n_genes, arma::fill::zeros);
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

        // Update C (co-classification) and E (per-pair expected) matrices,
        // and accumulate the scalar expected value for resolution_scan.
        double exp_k = 0.0;
        for (auto& mod : modules) {
            int s = static_cast<int>(mod.size());
            if (s < 1) continue;
            double frac = static_cast<double>(s) * inv_n;
            double frac_sq = frac * frac;
            exp_k += frac_sq;
            // Singletons only affect the diagonal, which is zeroed in R
            if (s < 2) continue;
            for (int bi = 0; bi < s; ++bi) {
                double* col_c = C.colptr(mod[bi]);
                double* col_e = E.colptr(mod[bi]);
                for (int ai = 0; ai < s; ++ai) {
                    col_c[mod[ai]] += 1.0;
                    col_e[mod[ai]] += frac_sq;
                }
            }
        }
        expected_vec[k] = exp_k;
    }

    // Normalize to [0, 1]
    double inv_K = 1.0 / static_cast<double>(K);
    C *= inv_K;
    E *= inv_K;

    // Excess co-classification: reuse E to store max(C - E, 0)
    E = arma::clamp(C - E, 0.0, 1.0);

    return Rcpp::List::create(
        Named("coclassification") = C,
        Named("excess_coclassification") = E,
        Named("expected") = expected_vec
    );
}
