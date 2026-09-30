// specificity.cpp
// Neighbourhood specificity: for each ortholog pair (i in species a, j* in
// species b) and direction a -> b, the anchor's neighbourhood N_a(i) is
// mapped through the orthologs to a species-b gene set T_i (own-HOG
// orthologs removed) and every species-b gene j is scored by the AUROC of
// T_i \ {j} against the ranks of column j of network b. The reported
// p-value is the rank of j* among all n_b candidates on the 1 / n_b grid.
//
// Sparse only. Ranks come from a per-network "rank store": the stored
// entries of column j are its top d_j values and take ranks
// (n - d_j) .. (n - 1); every unstored entry ties at the bottom with
// rank (n - d_j) / 2. rk_sym[e] holds the rank of column c within column
// i[e] (the transpose entry) so that one streaming pass over the stored
// rows of every k in T_i scatters rank_j(k) into S[j].
//
// Integer indices only (string mapping in R), no R API calls inside the
// OpenMP regions.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "neighbor_lists.h"
#include "rank_column.h"

using namespace Rcpp;

namespace {

struct RankStore {
    std::vector<double> rk_sym;  // rank of column c within column i[e]
    std::vector<double> mid;     // rank of every unstored entry, per column
};

RankStore build_rank_store(const IntegerVector& p, const IntegerVector& i,
                           const NumericVector& x, int n_cores) {
    const int n = validate_dgc_slots(p, i, x);
    const int* pp = p.begin();
    const int* ii = i.begin();
    std::vector<double> rk(x.begin(), x.end());
    double* rkp = rk.data();
    RankStore s;
    s.mid.resize(n);
    double* midp = s.mid.data();

#ifdef _OPENMP
    #pragma omp parallel num_threads(n_cores) if(n_cores > 1)
#endif
    {
        std::vector<R_xlen_t> buf;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 64)
#endif
        for (int c = 0; c < n; ++c) {
            const R_xlen_t d = pp[c + 1] - pp[c];
            buf.resize(d);
            rank_column_inplace(rkp + pp[c], d, true, buf);
            const double shift = n - 1.0 - static_cast<double>(d);
            for (int e = pp[c]; e < pp[c + 1]; ++e) rkp[e] += shift;
            midp[c] = (n - static_cast<double>(d)) / 2.0;
        }
    }

    // Transpose cursor pass: columns visited ascending and rows ascending
    // within a column, so cur[r] walks column r in step with row r's
    // appearances. Any mismatch means the pattern is not symmetric.
    s.rk_sym.resize(rk.size());
    std::vector<int> cur(pp, pp + n);
    for (int c = 0; c < n; ++c) {
        for (int e = pp[c]; e < pp[c + 1]; ++e) {
            const int r = ii[e];
            const int e2 = cur[r]++;
            if (e2 >= pp[r + 1] || ii[e2] != c) {
                stop("network pattern is not symmetric");
            }
            s.rk_sym[e] = rk[e2];
        }
    }
    return s;
}

// Same quantity as compute_direction()'s r.jaccard in
// neighborhood_comparison.cpp: other-network neighbours mapped back to the
// anchor network (anchor excluded), intersected with the anchor's own.
double jaccard(const std::vector<int>& anchor_neigh,
               const std::vector<int>& other_neigh,
               const std::vector<std::vector<int>>& mapping,
               int anchor, std::vector<char>& flags) {
    int k = 0, x = 0;
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) {
            if (idx == anchor || flags[idx]) continue;
            flags[idx] = 1;
            ++k;
        }
    }
    for (int nb : anchor_neigh) x += flags[nb];
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) flags[idx] = 0;
    }
    const int u = static_cast<int>(anchor_neigh.size()) + k - x;
    return u > 0 ? static_cast<double>(x) / u : 0.0;
}

// One direction: anchors in species a, candidates in species b.
void specificity_direction(
    const std::vector<std::vector<int>>& neigh_a,
    const std::vector<std::vector<int>>& neigh_b,
    const std::vector<std::vector<int>>& a_to_b,
    const std::vector<std::vector<int>>& b_to_a,
    const IntegerVector& p_b, const IntegerVector& i_b,
    const RankStore& store_b,
    const IntegerVector& pair_a, const IntegerVector& pair_b,
    int n_cores,
    int* out_neigh, int* out_mapped, double* out_auroc, double* out_p,
    double* out_jaccard
) {
    const int n_a = static_cast<int>(neigh_a.size());
    const int n_b = static_cast<int>(neigh_b.size());
    const int n_pairs = static_cast<int>(pair_a.size());

    // Tested pair rows grouped by anchor (counting sort)
    std::vector<int> aptr(n_a + 1, 0), rows(n_pairs);
    for (int q = 0; q < n_pairs; ++q) ++aptr[pair_a[q] + 1];
    for (int a = 0; a < n_a; ++a) aptr[a + 1] += aptr[a];
    {
        std::vector<int> fill(aptr.begin(), aptr.end() - 1);
        for (int q = 0; q < n_pairs; ++q) rows[fill[pair_a[q]]++] = q;
    }
    std::vector<int> anchors;
    for (int a = 0; a < n_a; ++a) {
        if (aptr[a + 1] > aptr[a]) anchors.push_back(a);
    }
    const int n_anchor = static_cast<int>(anchors.size());

    const int* pp = p_b.begin();
    const int* ii = i_b.begin();
    const double* rks = store_b.rk_sym.data();
    const double* mid = store_b.mid.data();
    const int* pb = pair_b.begin();

#ifdef _OPENMP
    #pragma omp parallel num_threads(n_cores) if(n_cores > 1)
#endif
    {
        std::vector<double> S(n_b), A(n_b);
        std::vector<int> cnt(n_b), T;
        std::vector<char> inT(n_b, 0), flags(n_a, 0);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif
        for (int u = 0; u < n_anchor; ++u) {
            const int i = anchors[u];
            // T_i: orthologs of the neighbours, minus own-HOG orthologs (2)
            for (int j : a_to_b[i]) inT[j] = 2;
            T.clear();
            for (int k : neigh_a[i]) {
                for (int j : a_to_b[k]) {
                    if (inT[j] == 0) {
                        inT[j] = 1;
                        T.push_back(j);
                    }
                }
            }
            const int t = static_cast<int>(T.size());
            if (t > 0) {
                std::fill(S.begin(), S.end(), 0.0);
                std::fill(cnt.begin(), cnt.end(), 0);
                for (int k : T) {
                    for (int e = pp[k]; e < pp[k + 1]; ++e) {
                        S[ii[e]] += rks[e];
                        ++cnt[ii[e]];
                    }
                }
                for (int j = 0; j < n_b; ++j) {
                    const int tj = t - (inT[j] == 1);
                    if (tj == 0 || tj >= n_b - 1) {
                        A[j] = NAN;
                        continue;
                    }
                    const double s = S[j] + (tj - cnt[j]) * mid[j];
                    A[j] = (s - tj * (tj + 1.0) / 2.0) /
                        (tj * (n_b - 1.0 - tj));
                }
            }
            for (int q = aptr[i]; q < aptr[i + 1]; ++q) {
                const int row = rows[q];
                const int js = pb[row];
                out_neigh[row] = static_cast<int>(neigh_a[i].size());
                out_mapped[row] = t;
                out_jaccard[row] = jaccard(neigh_a[i], neigh_b[js], b_to_a,
                                           i, flags);
                if (t == 0 || std::isnan(A[js])) {
                    out_auroc[row] = NA_REAL;
                    out_p[row] = NA_REAL;
                    continue;
                }
                int ge = 0;
                for (int j = 0; j < n_b; ++j) {
                    if (j != js && A[j] >= A[js]) ++ge;  // NaN never counts
                }
                out_auroc[row] = A[js];
                out_p[row] = (1.0 + ge) / n_b;
            }
            for (int j : T) inT[j] = 0;
            for (int j : a_to_b[i]) inT[j] = 0;
        }
    }
}

}  // namespace


//' Neighbourhood specificity of ortholog pairs (sparse networks)
//'
//' Takes the `dgCMatrix` slots of both networks (both triangles stored, no
//' diagonal, symmetric pattern) as [compare_neighborhoods_sparse_cpp()]
//' does. `pair_*` are the tested pairs, `ortho_*` the full ortholog
//' mapping; all indices 0-based.
//'
//' @inheritParams compare_neighborhoods_sparse_cpp
//' @param do_12,do_21 Compute direction 1 -> 2 / 2 -> 1.
//' @return DataFrame with, per requested direction, `Species1.neigh`,
//'   `Species1.mapped`, `Species1.auroc`, `Species1.p.val`,
//'   `Species1.jaccard` (and the `Species2.*` set), one row per pair.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame specificity_sparse_cpp(
    const Rcpp::IntegerVector& p1, const Rcpp::IntegerVector& i1,
    const Rcpp::NumericVector& x1, double thr1,
    const Rcpp::IntegerVector& p2, const Rcpp::IntegerVector& i2,
    const Rcpp::NumericVector& x2, double thr2,
    const Rcpp::IntegerVector& pair_sp1_idx,
    const Rcpp::IntegerVector& pair_sp2_idx,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    bool do_12, bool do_21, int n_cores
) {
    const auto neigh1 = neighbor_lists_sparse(p1, i1, x1, thr1, n_cores);
    const auto neigh2 = neighbor_lists_sparse(p2, i2, x2, thr2, n_cores);
    const int n1 = static_cast<int>(neigh1.size());
    const int n2 = static_cast<int>(neigh2.size());
    const int n_pairs = static_cast<int>(pair_sp1_idx.size());
    if (pair_sp2_idx.size() != n_pairs ||
        ortho_sp1_idx.size() != ortho_sp2_idx.size()) {
        stop("pair / ortholog index vectors must have equal lengths");
    }
    for (int q = 0; q < n_pairs; ++q) {
        if (pair_sp1_idx[q] < 0 || pair_sp1_idx[q] >= n1 ||
            pair_sp2_idx[q] < 0 || pair_sp2_idx[q] >= n2) {
            stop("pair index out of range (row %d)", q + 1);
        }
    }

    std::vector<std::vector<int>> sp2_to_sp1(n2), sp1_to_sp2(n1);
    for (R_xlen_t k = 0; k < ortho_sp1_idx.size(); ++k) {
        const int s1 = ortho_sp1_idx[k], s2 = ortho_sp2_idx[k];
        if (s1 >= 0 && s1 < n1 && s2 >= 0 && s2 < n2) {
            sp2_to_sp1[s2].push_back(s1);
            sp1_to_sp2[s1].push_back(s2);
        }
    }

    Rcpp::List out;
    if (do_12) {
        IntegerVector neigh(n_pairs), mapped(n_pairs);
        NumericVector auroc(n_pairs), pv(n_pairs), jac(n_pairs);
        const RankStore s2 = build_rank_store(p2, i2, x2, n_cores);
        specificity_direction(
            neigh1, neigh2, sp1_to_sp2, sp2_to_sp1, p2, i2, s2,
            pair_sp1_idx, pair_sp2_idx, n_cores,
            neigh.begin(), mapped.begin(), auroc.begin(), pv.begin(),
            jac.begin());
        out["Species1.neigh"] = neigh;
        out["Species1.mapped"] = mapped;
        out["Species1.auroc"] = auroc;
        out["Species1.p.val"] = pv;
        out["Species1.jaccard"] = jac;
    }
    if (do_21) {
        IntegerVector neigh(n_pairs), mapped(n_pairs);
        NumericVector auroc(n_pairs), pv(n_pairs), jac(n_pairs);
        const RankStore s1 = build_rank_store(p1, i1, x1, n_cores);
        specificity_direction(
            neigh2, neigh1, sp2_to_sp1, sp1_to_sp2, p1, i1, s1,
            pair_sp2_idx, pair_sp1_idx, n_cores,
            neigh.begin(), mapped.begin(), auroc.begin(), pv.begin(),
            jac.begin());
        out["Species2.neigh"] = neigh;
        out["Species2.mapped"] = mapped;
        out["Species2.auroc"] = auroc;
        out["Species2.p.val"] = pv;
        out["Species2.jaccard"] = jac;
    }
    return Rcpp::DataFrame(out);
}
