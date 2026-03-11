// density_threshold.cpp
// Efficient density-based threshold computation using partial sort
//
// Uses std::ranges::nth_element for O(n) partial sort instead of full
// O(n log n) sort.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <ranges>
#include <vector>

using namespace Rcpp;

//' Compute density threshold from a symmetric matrix
//'
//' Extracts the upper triangle (excluding diagonal) and finds the value at
//' the given density percentile using O(n) partial sort.
//'
//' @param mat Symmetric matrix (n x n)
//' @param density Fraction of top edges to keep (e.g., 0.03 for 3%)
//' @return The threshold value such that `density` fraction of upper triangle
//'   values are >= the threshold
//'
//' @keywords internal
// [[Rcpp::export]]
double density_threshold_cpp(const arma::mat& mat, double density) {
    const auto n = mat.n_rows;

    if (n != mat.n_cols) {
        stop("mat must be a square matrix");
    }
    if (density <= 0.0 || density >= 1.0) {
        stop("density must be between 0 and 1 (exclusive)");
    }

    auto tri_size = n * (n - 1) / 2;
    std::vector<double> values(tri_size);

    size_t idx = 0;
    for (arma::uword j = 1; j < n; ++j) {
        for (arma::uword i = 0; i < j; ++i) {
            values[idx++] = mat(i, j);
        }
    }

    auto k = static_cast<decltype(tri_size)>(
        std::round(density * static_cast<double>(tri_size)));
    if (k == 0) k = 1;
    if (k >= tri_size) k = tri_size - 1;

    auto pos = tri_size - k;
    std::ranges::nth_element(values, values.begin() + static_cast<ptrdiff_t>(pos));

    return values[pos];
}
