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


//' Sparse co-classification restricted to original network edges
//'
//' Computes co-classification and per-pair excess only for gene pairs
//' connected in the thresholded co-expression network, reducing memory
//' from O(N^2) to O(|E|).
//'
//' @param memberships List of K IntegerVectors (1-based module IDs)
//' @param n_genes Number of genes (N)
//' @param edges |E| x 2 IntegerMatrix of 0-based vertex indices
//' @return List with from, to (0-based), coclassification, excess,
//'   and expected (length-K per-resolution expected scalars)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List build_sparse_coclassification_cpp(
    const Rcpp::List& memberships,
    int n_genes,
    const Rcpp::IntegerMatrix& edges)
{
    int K = memberships.size();
    int n_edges = edges.nrow();

    std::vector<double> C_vec(n_edges, 0.0);
    std::vector<double> E_vec(n_edges, 0.0);
    Rcpp::NumericVector expected_vec(K);

    double inv_n = 1.0 / n_genes;

    for (int k = 0; k < K; ++k) {
        if (k % 10 == 0) Rcpp::checkUserInterrupt();

        Rcpp::IntegerVector mem = memberships[k];

        // Build mod_of and mod_size (1-based module IDs from igraph)
        int max_mod = 0;
        for (int i = 0; i < n_genes; ++i)
            if (mem[i] > max_mod) max_mod = mem[i];

        std::vector<int> mod_of(n_genes);
        std::vector<int> mod_size(max_mod + 1, 0);
        for (int i = 0; i < n_genes; ++i) {
            mod_of[i] = mem[i];
            mod_size[mem[i]]++;
        }

        // Per-resolution expected scalar: sum_m (s_m / N)^2
        double exp_k = 0.0;
        for (int m = 0; m <= max_mod; ++m) {
            if (mod_size[m] > 0) {
                double frac = static_cast<double>(mod_size[m]) * inv_n;
                exp_k += frac * frac;
            }
        }
        expected_vec[k] = exp_k;

        for (int idx = 0; idx < n_edges; ++idx) {
            int i = edges(idx, 0);
            int j = edges(idx, 1);
            if (mod_of[i] == mod_of[j]) {
                C_vec[idx] += 1.0;
                double frac = static_cast<double>(mod_size[mod_of[i]]) * inv_n;
                E_vec[idx] += frac * frac;
            }
        }
    }

    double inv_K = 1.0 / static_cast<double>(K);
    Rcpp::NumericVector coclassification(n_edges);
    Rcpp::NumericVector excess(n_edges);
    for (int idx = 0; idx < n_edges; ++idx) {
        coclassification[idx] = C_vec[idx] * inv_K;
        double e = E_vec[idx] * inv_K;
        double ex = coclassification[idx] - e;
        excess[idx] = ex > 0.0 ? ex : 0.0;
    }

    Rcpp::IntegerVector from_out(n_edges);
    Rcpp::IntegerVector to_out(n_edges);
    for (int idx = 0; idx < n_edges; ++idx) {
        from_out[idx] = edges(idx, 0);
        to_out[idx] = edges(idx, 1);
    }

    return Rcpp::List::create(
        Named("from") = from_out,
        Named("to") = to_out,
        Named("coclassification") = coclassification,
        Named("excess") = excess,
        Named("expected") = expected_vec
    );
}


//' Leading eigenvalue of sparse excess co-classification matrix
//'
//' Computes the spectral norm (largest eigenvalue) of the excess
//' co-classification matrix restricted to the original edge set.
//' Used for the K = 1 community structure test.
//'
//' @param memberships List of K IntegerVectors (1-based module IDs)
//' @param n_genes Number of genes (N)
//' @param edges |E| x 2 IntegerMatrix of 0-based vertex indices
//' @return Leading eigenvalue (non-negative scalar)
//' @keywords internal
// [[Rcpp::export]]
double sparse_excess_spectral_norm_cpp(
    const Rcpp::List& memberships,
    int n_genes,
    const Rcpp::IntegerMatrix& edges)
{
    int K = memberships.size();
    int n_edges = edges.nrow();

    std::vector<double> C_vec(n_edges, 0.0);
    std::vector<double> E_vec(n_edges, 0.0);

    double inv_n = 1.0 / n_genes;

    for (int k = 0; k < K; ++k) {
        if (k % 10 == 0) Rcpp::checkUserInterrupt();

        Rcpp::IntegerVector mem = memberships[k];

        int max_mod = 0;
        for (int i = 0; i < n_genes; ++i)
            if (mem[i] > max_mod) max_mod = mem[i];

        std::vector<int> mod_of(n_genes);
        std::vector<int> mod_size(max_mod + 1, 0);
        for (int i = 0; i < n_genes; ++i) {
            mod_of[i] = mem[i];
            mod_size[mem[i]]++;
        }

        for (int idx = 0; idx < n_edges; ++idx) {
            int i = edges(idx, 0);
            int j = edges(idx, 1);
            if (mod_of[i] == mod_of[j]) {
                C_vec[idx] += 1.0;
                double frac = static_cast<double>(mod_size[mod_of[i]]) * inv_n;
                E_vec[idx] += frac * frac;
            }
        }
    }

    double inv_K = 1.0 / static_cast<double>(K);

    std::vector<arma::uword> row_idx;
    std::vector<arma::uword> col_idx;
    std::vector<double> vals;

    for (int idx = 0; idx < n_edges; ++idx) {
        double c = C_vec[idx] * inv_K;
        double e = E_vec[idx] * inv_K;
        double ex = c - e;
        if (ex > 0.0) {
            arma::uword i = static_cast<arma::uword>(edges(idx, 0));
            arma::uword j = static_cast<arma::uword>(edges(idx, 1));
            row_idx.push_back(i);
            col_idx.push_back(j);
            vals.push_back(ex);
            row_idx.push_back(j);
            col_idx.push_back(i);
            vals.push_back(ex);
        }
    }

    if (vals.empty()) return 0.0;

    arma::uword nnz = vals.size();
    arma::umat locations(2, nnz);
    arma::vec values(nnz);
    for (arma::uword k = 0; k < nnz; ++k) {
        locations(0, k) = row_idx[k];
        locations(1, k) = col_idx[k];
        values(k) = vals[k];
    }

    arma::sp_mat X(locations, values, n_genes, n_genes);

    arma::vec eigval;
    arma::mat eigvec;
    bool ok = arma::eigs_sym(eigval, eigvec, X, 1);

    if (ok && eigval.n_elem > 0) {
        return eigval(0);
    }

    if (n_genes < 5000) {
        arma::mat X_dense(X);
        arma::vec eigval_dense = arma::eig_sym(X_dense);
        return eigval_dense(eigval_dense.n_elem - 1);
    }

    Rcpp::warning(
        "eigs_sym failed to converge and matrix too large for dense fallback");
    return 0.0;
}
