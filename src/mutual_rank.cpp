// mutual_rank.cpp
// Mutual rank transformation for co-expression network normalization
//
// Supports two modes:
// - log_transform=true: Obayashi & Kinoshita (2009) formula
//   S = 1 - log(sqrt(R_ij * R_ji)) / log(n), descending ranks, values in [0,1]
// - log_transform=false (default): Raw mutual rank (as in original RComPlEx Rmd)
//   MR = sqrt(R_ij * R_ji), ascending ranks, unbounded
//
// Adapted from coexpr package by Martin Paliocha

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Average-rank computation with projection-based sort.
// ascending=true  -> lowest value gets rank 1
// ascending=false -> highest value gets rank 1
static arma::vec compute_ranks_impl(const arma::vec& x, bool ascending) {
    const auto n = x.n_elem;
    arma::vec ranks(n);

    std::vector<arma::uword> indices(n);
    std::iota(indices.begin(), indices.end(), arma::uword{0});

    auto proj = [&x](arma::uword i) { return x(i); };
    if (ascending) {
        std::ranges::sort(indices, std::ranges::less{}, proj);
    } else {
        std::ranges::sort(indices, std::ranges::greater{}, proj);
    }

    arma::uword i = 0;
    while (i < n) {
        arma::uword j = i;
        while (j < n - 1 && x(indices[j]) == x(indices[j + 1])) {
            ++j;
        }
        double avg_rank = static_cast<double>(i + j + 2) / 2.0;
        for (arma::uword k = i; k <= j; ++k) {
            ranks(indices[k]) = avg_rank;
        }
        i = j + 1;
    }

    return ranks;
}

//' Cached mutual rank transformation
//'
//' Transforms a correlation matrix using mutual rank normalization.
//'
//' @param sim Symmetric correlation/similarity matrix (n x n)
//' @param log_transform If FALSE (default), uses raw mutual rank with ascending
//'   ranks (original Rmd formula). If TRUE, applies Obayashi & Kinoshita (2009)
//'   log-normalized formula with descending ranks (values in 0 to 1 range).
//' @param n_cores Number of OpenMP threads (default: 1)
//' @return Mutual rank normalized matrix (n x n)
//'
//' @details
//' When `log_transform = TRUE`:
//' \deqn{S_{ij} = 1 - \log(\sqrt{R_{ij} \cdot R_{ji}}) / \log(n)}
//' where ranks are descending (highest correlation = rank 1).
//'
//' When `log_transform = FALSE`:
//' \deqn{MR_{ij} = \sqrt{R_{ij} \cdot R_{ji}}}
//' where ranks are ascending (lowest correlation = rank 1), matching the
//' original RComPlEx R Markdown formula.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat mutual_rank_transform_cached_cpp(const arma::mat& sim,
                                           bool log_transform = false,
                                           int n_cores = 1) {
    const auto n = sim.n_rows;

    if (n != sim.n_cols) {
        stop("sim must be a square matrix");
    }
    if (n < 2) {
        stop("Matrix must have at least 2 rows/columns");
    }

    const double log_n = std::log(static_cast<double>(n));
    const bool ascending = !log_transform;

    arma::mat row_ranks(n, n);
    arma::mat result(n, n, arma::fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (arma::uword i = 0; i < n; ++i) {
        row_ranks.col(i) = compute_ranks_impl(sim.col(i), ascending);
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1) num_threads(n_cores) if(n_cores > 1)
#endif
    for (arma::uword i = 0; i < n; ++i) {
        for (arma::uword j = i + 1; j < n; ++j) {
            double mutual_rank = std::sqrt(row_ranks(j, i) * row_ranks(i, j));

            double val;
            if (log_transform) {
                val = std::clamp(1.0 - std::log(mutual_rank) / log_n, 0.0, 1.0);
            } else {
                val = mutual_rank;
            }

            result(i, j) = val;
            result(j, i) = val;
        }
        result(i, i) = log_transform ? 1.0 : 0.0;
    }

    return result;
}
