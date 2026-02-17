// density_threshold.cpp
// Efficient density-based threshold computation using partial sort
//
// Uses std::nth_element for O(n) partial sort instead of full O(n log n) sort.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

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
    const uword n = mat.n_rows;

    if (n != mat.n_cols) {
        stop("mat must be a square matrix");
    }
    if (density <= 0.0 || density >= 1.0) {
        stop("density must be between 0 and 1 (exclusive)");
    }

    // Extract upper triangle values
    uword tri_size = n * (n - 1) / 2;
    std::vector<double> values(tri_size);

    uword idx = 0;
    for (uword j = 1; j < n; ++j) {
        for (uword i = 0; i < j; ++i) {
            values[idx++] = mat(i, j);
        }
    }

    // Find the threshold at the density percentile
    // We want the top `density` fraction, so the index from the top
    uword k = static_cast<uword>(std::round(density * tri_size));
    if (k == 0) k = 1;
    if (k >= tri_size) k = tri_size - 1;

    // nth_element partitions so that element at k is what would be there in sorted order
    // We want the k-th largest, so we use (tri_size - k) position in ascending order
    uword pos = tri_size - k;
    std::nth_element(values.begin(), values.begin() + pos, values.end());

    return values[pos];
}
