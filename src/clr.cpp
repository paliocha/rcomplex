// clr.cpp
// CLR (Context Likelihood Ratio) transformation for correlation matrices
//
// Adapted from coexpr package's apply_clr_cpp (mutual_information.cpp)
// Applied directly to correlation matrices rather than MI matrices.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//' Apply CLR transformation to a correlation matrix
//'
//' Normalizes a correlation matrix using the Context Likelihood Ratio (CLR)
//' transformation (Faith et al. 2007). Computes row z-scores, clamps negatives
//' to zero, and combines row and column z-scores.
//'
//' @param cor_matrix Symmetric correlation matrix (n x n)
//' @param n_cores Number of OpenMP threads (default: 1)
//' @return CLR-transformed matrix (n x n), symmetric, non-negative
//'
//' @details
//' The CLR transformation:
//' \preformatted{
//' z_row(i,j) = (cor(i,j) - mean(row_i)) / sd(row_i)
//' z_col(i,j) = (cor(i,j) - mean(col_j)) / sd(col_j)
//' CLR(i,j) = sqrt(max(0, z_row)^2 + max(0, z_col)^2)
//' }
//'
//' This matches the original RComPlEx Rmd CLR implementation which uses
//' R's scale() function (column-wise z-scoring on the symmetric matrix,
//' equivalent to row-wise z-scoring due to symmetry).
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat apply_clr_to_cor_cpp(const arma::mat& cor_matrix, int n_cores = 1) {
    const auto n = cor_matrix.n_rows;

    if (n != cor_matrix.n_cols) {
        stop("cor_matrix must be a square matrix");
    }

    arma::vec row_means = arma::mean(cor_matrix, 1);
    arma::vec row_sds = arma::stddev(cor_matrix, 0, 1);
    row_sds.elem(arma::find(row_sds < 1e-10)).fill(1.0);

    // Vectorised row z-scoring and clamping (replaces per-element arithmetic)
    arma::mat Z = cor_matrix;
    Z.each_col() -= row_means;
    Z.each_col() /= row_sds;
    Z.transform([](double x) { return std::max(0.0, x); });

    // Combine row and column z-scores: CLR(i,j) = hypot(Z(i,j), Z(j,i))
    // Each (i,j)/(j,i) pair read+written by exactly one iteration — in-place safe
#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
    #pragma omp parallel for schedule(static, 1) if(n_cores > 1)
#endif
    for (arma::uword i = 0; i < n; ++i) {
        for (arma::uword j = i + 1; j < n; ++j) {
            double val = std::hypot(Z(i, j), Z(j, i));
            Z(i, j) = val;
            Z(j, i) = val;
        }
        Z(i, i) = 0.0;
    }

    return Z;
}
