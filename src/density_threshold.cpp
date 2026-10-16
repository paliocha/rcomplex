// density_threshold.cpp
// Efficient density-based threshold computation using partial sort
//
// Uses std::ranges::nth_element for O(n) partial sort instead of full
// O(n log n) sort.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <vector>

#include "density_k.h"

// A 32-bit uword wraps the element count of an n x n matrix above 65,536
// genes. src/Makevars sets ARMA_64BIT_WORD; this fails the build without it.
static_assert(sizeof(arma::uword) == 8,
              "rcomplex needs ARMA_64BIT_WORD: n x n matrices above "
              "65,536 genes wrap a 32-bit index");

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
    const std::size_t n = mat.n_rows;

    if (n != mat.n_cols) {
        stop("mat must be a square matrix");
    }
    if (density <= 0.0 || density >= 1.0) {
        stop("density must be between 0 and 1 (exclusive)");
    }

    const std::size_t tri_size = n * (n - 1) / 2;
    std::vector<double> values(tri_size);

    std::size_t idx = 0;
    for (std::size_t j = 1; j < n; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            values[idx++] = mat(i, j);
        }
    }

    const std::size_t k = density_k(density, tri_size);

    const std::size_t pos = tri_size - k;
    std::ranges::nth_element(
        values, values.begin() + static_cast<std::ptrdiff_t>(pos));

    return values[pos];
}
