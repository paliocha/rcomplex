// sparse_extract.cpp
// Dense -> dgCMatrix slot extraction for compute_network(sparse = TRUE).
//
// Two passes over the dense matrix (column-major, colptr): count the kept
// entries per column, prefix-sum into p, then fill i and x. Row indices come
// out ascending within each column (the dgCMatrix invariant) because rows
// are scanned in order. OpenMP over columns; both passes write disjoint
// ranges and make no R API calls inside the parallel regions.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//' Extract dgCMatrix slots from a dense matrix at a store threshold
//'
//' Keeps the off-diagonal entries `>= thr` of both triangles and returns
//' the column-compressed slots of the corresponding `dgCMatrix`: 0-based
//' row indices `i` (ascending within each column), column pointers `p`
//' (length n + 1) and values `x`.
//'
//' @param m Dense numeric matrix (n x n)
//' @param thr Store threshold; entries below it are discarded
//' @param n_cores Number of threads for parallel computation (default 1)
//' @return List with integer vectors `i` and `p` and numeric vector `x`
//'
//' @keywords internal
// [[Rcpp::export]]
List extract_sparse_cpp(const arma::mat& m, double thr, int n_cores = 1) {
    if (m.n_rows != m.n_cols) {
        stop("m must be a square matrix");
    }
    const int n = static_cast<int>(m.n_rows);

    IntegerVector p(n + 1);  // zero-initialised by Rcpp
    int* pp = p.begin();

    // Pass 1: count kept entries per column
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int c = 0; c < n; ++c) {
        const double* col = m.colptr(c);
        int cnt = 0;
        for (int r = 0; r < n; ++r) {
            if (r != c && col[r] >= thr) ++cnt;
        }
        pp[c + 1] = cnt;
    }

    // Serial prefix sum
    for (int c = 0; c < n; ++c) {
        pp[c + 1] += pp[c];
    }
    const int nnz = pp[n];

    IntegerVector i(nnz);
    NumericVector x(nnz);
    int* ii = i.begin();
    double* xx = x.begin();

    // Pass 2: fill; rows scanned ascending -> sorted within each column
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int c = 0; c < n; ++c) {
        const double* col = m.colptr(c);
        int k = pp[c];
        for (int r = 0; r < n; ++r) {
            if (r != c && col[r] >= thr) {
                ii[k] = r;
                xx[k] = col[r];
                ++k;
            }
        }
    }

    return List::create(Named("i") = i, Named("p") = p, Named("x") = x);
}
