// coclassification.cpp
// Co-classification matrix with per-pair null subtraction for iterative
// consensus module detection (Jeub et al., 2018).

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

//' Build co-classification matrix with per-pair excess
//'
//' @param memberships List of K IntegerVectors (1-based module IDs)
//' @param n_genes Number of genes (N)
//' @param return_excess If TRUE, return excess (adaptive path).
//'   If FALSE, return raw co-classification (fixed threshold path).
//'   Only one N x N matrix is returned, saving ~4.6 GB for N = 24k.
//' @return List with one of coclassification or excess_coclassification,
//'   plus expected (length-K per-resolution expected scalars)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List build_coclassification_cpp(
    const Rcpp::List& memberships,
    int n_genes,
    bool return_excess = true)
{
    int K = memberships.size();
    arma::mat C(n_genes, n_genes, arma::fill::zeros);
    arma::mat E;
    if (return_excess) E.zeros(n_genes, n_genes);
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
                for (int ai = 0; ai < s; ++ai)
                    col_c[mod[ai]] += 1.0;
                if (return_excess) {
                    double* col_e = E.colptr(mod[bi]);
                    for (int ai = 0; ai < s; ++ai)
                        col_e[mod[ai]] += frac_sq;
                }
            }
        }
        expected_vec[k] = exp_k;
    }

    double inv_K = 1.0 / static_cast<double>(K);
    C *= inv_K;

    if (return_excess) {
        E *= inv_K;
        // Excess = max(C - E, 0) in-place, no temporary
        C -= E;
        E.reset();
        C.transform([](double val) { return val < 0.0 ? 0.0 : val; });
        return Rcpp::List::create(
            Named("excess_coclassification") = C,
            Named("expected") = expected_vec
        );
    } else {
        return Rcpp::List::create(
            Named("coclassification") = C,
            Named("expected") = expected_vec
        );
    }
}
