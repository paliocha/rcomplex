// module_auroc.cpp
// Set-level neighbour voting (EGAD, Ballouz et al. 2017) for
// module_auroc(). Per gene set T (genes set_gene[set_ptr[s] ..
// set_ptr[s + 1]), 0-based, fold labels 1..n_fold in set_fold): with the
// binary adjacency W (CSC slots p, i; symmetric, no diagonal) and
// a_f = W 1_{fold f}, a gene j scored with fold f held out gets
// (sum_g a_g[j] - a_f[j]) / deg(j), i.e. its edges to the training part of
// T over its degree. The AUROC of the fold-f genes against every gene
// outside T is computed by binary search against the sorted held-out
// scores (rank-sum, ties counted one half); the set's value is the mean
// over non-empty folds.
//
// Integer indices only, per-thread buffers, no R API inside the OpenMP
// region; sets are independent, so the result does not depend on n_cores.

// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector module_auroc_cpp(IntegerVector p, IntegerVector i,
                               IntegerVector set_ptr, IntegerVector set_gene,
                               IntegerVector set_fold, int n_fold,
                               int n_cores) {
    const int n = p.size() - 1;
    const int ns = set_ptr.size() - 1;
    if (n < 0 || ns < 0 || n_fold < 1) stop("invalid kernel arguments");
    if (set_gene.size() != set_fold.size() ||
        set_ptr[ns] != set_gene.size()) {
        stop("set_ptr, set_gene and set_fold are inconsistent");
    }
    for (R_xlen_t k = 0; k < set_gene.size(); ++k) {
        if (set_gene[k] < 0 || set_gene[k] >= n ||
            set_fold[k] < 1 || set_fold[k] > n_fold) {
            stop("set gene index or fold label out of range");
        }
    }
    const int* pp = p.begin();
    const int* ii = i.begin();
    const int* sp = set_ptr.begin();
    const int* sg = set_gene.begin();
    const int* sf = set_fold.begin();
    std::vector<double> out(ns, NA_REAL);
    double* op = out.data();

#ifdef _OPENMP
    #pragma omp parallel num_threads(n_cores) if(n_cores > 1)
#endif
    {
        std::vector<int> a(static_cast<size_t>(n_fold) * n, 0);
        std::vector<int> tot(n, 0);
        std::vector<char> in_t(n, 0);
        std::vector<double> test;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for (int s = 0; s < ns; ++s) {
            const int lo = sp[s], hi = sp[s + 1];
            for (int k = lo; k < hi; ++k) {
                const int g = sg[k];
                in_t[g] = 1;
                int* af = &a[static_cast<size_t>(sf[k] - 1) * n];
                for (int e = pp[g]; e < pp[g + 1]; ++e) {
                    ++af[ii[e]];
                    ++tot[ii[e]];
                }
            }
            double acc = 0.0;
            int nf = 0;
            for (int f = 0; f < n_fold; ++f) {
                const int* af = &a[static_cast<size_t>(f) * n];
                test.clear();
                for (int k = lo; k < hi; ++k) {
                    if (sf[k] - 1 != f) continue;
                    const int g = sg[k];
                    const int d = pp[g + 1] - pp[g];
                    test.push_back(d > 0 ? double(tot[g] - af[g]) / d : 0.0);
                }
                if (test.empty()) continue;
                std::sort(test.begin(), test.end());
                double sum = 0.0, nneg = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (in_t[j]) continue;
                    const int d = pp[j + 1] - pp[j];
                    const double v = d > 0 ? double(tot[j] - af[j]) / d : 0.0;
                    auto lb = std::lower_bound(test.begin(), test.end(), v);
                    auto ub = std::upper_bound(lb, test.end(), v);
                    sum += double(test.end() - ub) + 0.5 * double(ub - lb);
                    nneg += 1.0;
                }
                if (nneg == 0.0) continue;
                acc += sum / (double(test.size()) * nneg);
                ++nf;
            }
            if (nf) op[s] = acc / nf;
            // reset only what this set touched
            for (int k = lo; k < hi; ++k) {
                const int g = sg[k];
                in_t[g] = 0;
                int* af = &a[static_cast<size_t>(sf[k] - 1) * n];
                for (int e = pp[g]; e < pp[g + 1]; ++e) {
                    af[ii[e]] = 0;
                    tot[ii[e]] = 0;
                }
            }
        }
    }
    return NumericVector(out.begin(), out.end());
}
