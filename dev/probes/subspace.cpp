// was rcomplex::top_eigs_sym_cpp (subspace_preservation kernel) until 0.4.0; see git history (main @ b245ac0)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Top-k eigenpairs (largest algebraic) of a sparse symmetric matrix given
// as CSC slots, eigenvalues in decreasing order. ARPACK (arma::eigs_sym)
// on large inputs; dense eig_sym when n is small, when k is too close to
// n for ARPACK, or when ARPACK fails to converge and n allows it.
// [[Rcpp::export]]
Rcpp::List top_eigs_sym_cpp(const arma::uvec& p, const arma::uvec& i,
                            const arma::vec& x, int n, int k) {
    if (k < 1 || k > n) Rcpp::stop("k must be in 1..n");
    arma::sp_mat X(i, p, x, n, n);
    arma::vec val;
    arma::mat vec;
    bool ok = false;
    if (n > 1000 && k < n - 1) {
        ok = arma::eigs_sym(val, vec, X, k, "la");
    }
    if (!ok) {
        if (n > 20000) {
            Rcpp::stop("eigs_sym failed to converge; matrix too large "
                       "for the dense fallback");
        }
        arma::vec dval;
        arma::mat dvec;
        if (!arma::eig_sym(dval, dvec, arma::mat(X))) {
            Rcpp::stop("eig_sym failed");
        }
        // ascending; keep the last k
        val = dval.tail(k);
        vec = dvec.tail_cols(k);
    }
    arma::uvec ord = arma::sort_index(val, "descend");
    return Rcpp::List::create(
        Rcpp::Named("values") = arma::vec(val.elem(ord)),
        Rcpp::Named("vectors") = arma::mat(vec.cols(ord)));
}
