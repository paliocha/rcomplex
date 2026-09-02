// neighborhood_comparison.cpp
// Fast neighborhood comparison for cross-species co-expression analysis
//
// The main computational bottleneck: for each ortholog pair, tests overlap of
// co-expression neighborhoods in both directions using hypergeometric tests.
//
// Key optimizations:
// - All gene lookups use integer indices (string->int mapping done in R)
// - Precomputed neighbor lists (O(1) lookup vs O(n) per pair)
// - Bit-set or hash-set intersection for fast overlap counting
// - OpenMP parallelization over ortholog pairs
// - No R object allocation in the hot loop

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "neighbor_lists.h"

using namespace Rcpp;


// Per-direction hypergeometric test result
struct DirectionResult {
    int neigh = 0;
    int ortho_neigh = 0;
    int overlap = 0;
    double pval_con = 1.0;
    double pval_div = 1.0;
    double pval_gt = 0.0;   // P(X > x), ungated
    double pval_eq = 1.0;   // P(X = x); gt = 0, eq = 1 encodes p = 1 with a
                            // valid randomized-p decomposition (gt + U * eq
                            // is U at 1), matching compute_direction's
                            // degenerate m = k = x = 0 case
    double effect_size = 1.0;
    double jaccard = 0.0;
};

// Compute one direction of the neighborhood comparison.
//   anchor_neigh:  neighbors of the anchor gene in its own network
//   other_neigh:   neighbors of the paired gene in the other network
//   mapping:       other-network index -> vector of anchor-network indices
//   N:             total genes in the anchor network
//   anchor:        index of the anchor gene in its own network
//   anchor_flags:  scratch boolean vector (size N), must be pre-cleared
//
// Self-excluded urn: the anchor gene is never its own neighbour, so it
// leaves the ortholog-mapped set (k) and the hypergeometric population is
// the other N - 1 genes. The reported conservation p-value keeps the
// canonical x > 1 gate; pval_gt / pval_eq are ungated so that
// pval_gt + U * pval_eq is exactly uniform under H0 (randomized p-value).
static DirectionResult compute_direction(
    const std::vector<int>& anchor_neigh,
    const std::vector<int>& other_neigh,
    const std::vector<std::vector<int>>& mapping,
    int N,
    int anchor,
    std::vector<char>& anchor_flags
) {
    DirectionResult r{.neigh = static_cast<int>(anchor_neigh.size())};

    // Map other-network neighbors back to anchor-network via orthologs
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) {
            if (idx == anchor) continue;
            if (anchor_flags[idx] == 0) {
                anchor_flags[idx] = 1;
                ++r.ortho_neigh;
            }
        }
    }

    // Count intersection
    for (int nb : anchor_neigh) {
        if (anchor_flags[nb] != 0) ++r.overlap;
    }

    // Clear flags for reuse (the anchor flag was never set)
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) anchor_flags[idx] = 0;
    }

    const int m = r.neigh, k = r.ortho_neigh, x = r.overlap;
    const int Np = N - 1;   // population without the anchor gene
    r.pval_gt = R::phyper(x, m, Np - m, k, 0, 0);
    r.pval_eq = R::dhyper(x, m, Np - m, k, 0);
    if (x > 1) {
        r.pval_con = R::phyper(x - 1, m, Np - m, k, 0, 0);
    }
    if (k > 0 && m > 0) {
        r.effect_size = (static_cast<double>(x) / k) / (static_cast<double>(m) / Np);
        r.pval_div = R::phyper(x, m, Np - m, k, 1, 0);
    }

    int union_size = r.neigh + r.ortho_neigh - r.overlap;
    r.jaccard = (union_size > 0) ? static_cast<double>(r.overlap) / union_size : 0.0;

    return r;
}


// Shared body for both entry points: everything after neighbour-list
// construction. neighbors1[i] / neighbors2[j] are ascending 0-based lists
// (see neighbor_lists.h).
static Rcpp::DataFrame compare_neighborhoods_core(
    const std::vector<std::vector<int>>& neighbors1,
    const std::vector<std::vector<int>>& neighbors2,
    const Rcpp::IntegerVector& pair_sp1_idx,
    const Rcpp::IntegerVector& pair_sp2_idx,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    int n_cores
) {
    const int n_pairs = pair_sp1_idx.size();
    const int n1 = static_cast<int>(neighbors1.size());
    const int n2 = static_cast<int>(neighbors2.size());
    const int n_ortho = ortho_sp1_idx.size();

    // Build ortholog mapping: sp2 idx -> vector of sp1 indices, and vice versa
    // sp2_to_sp1[j] = list of sp1 indices that are orthologs of sp2 gene j
    std::vector<std::vector<int>> sp2_to_sp1(n2);
    std::vector<std::vector<int>> sp1_to_sp2(n1);

    for (int i = 0; i < n_ortho; ++i) {
        int s1 = ortho_sp1_idx[i];
        int s2 = ortho_sp2_idx[i];
        if (s1 >= 0 && s1 < n1 && s2 >= 0 && s2 < n2) {
            sp2_to_sp1[s2].push_back(s1);
            sp1_to_sp2[s1].push_back(s2);
        }
    }

    // Output vectors (9 per direction)
    Rcpp::IntegerVector sp1_neigh(n_pairs),  sp2_neigh(n_pairs);
    Rcpp::IntegerVector sp1_ortho_neigh(n_pairs), sp2_ortho_neigh(n_pairs);
    Rcpp::IntegerVector sp1_overlap(n_pairs), sp2_overlap(n_pairs);
    Rcpp::NumericVector sp1_pval(n_pairs),   sp2_pval(n_pairs);
    Rcpp::NumericVector sp1_pval_div(n_pairs), sp2_pval_div(n_pairs);
    Rcpp::NumericVector sp1_pval_gt(n_pairs), sp2_pval_gt(n_pairs);
    Rcpp::NumericVector sp1_pval_eq(n_pairs), sp2_pval_eq(n_pairs);
    Rcpp::NumericVector sp1_effect(n_pairs), sp2_effect(n_pairs);
    Rcpp::NumericVector sp1_jaccard(n_pairs), sp2_jaccard(n_pairs);

    // Pre-allocate per-thread scratch flags (cleared by compute_direction)
    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif
    std::vector<std::vector<char>> thread_flags1(max_threads, std::vector<char>(n1, 0));
    std::vector<std::vector<char>> thread_flags2(max_threads, std::vector<char>(n2, 0));

    // Parallel comparison loop
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (int p = 0; p < n_pairs; ++p) {
        int g1_idx = pair_sp1_idx[p];
        int g2_idx = pair_sp2_idx[p];

        if (g1_idx < 0 || g1_idx >= n1 || g2_idx < 0 || g2_idx >= n2) {
            // unreachable from R (the wrapper filters unmapped genes);
            // gt = 0, eq = 1 keeps the randomized-p decomposition valid
            sp1_neigh[p] = 0; sp1_ortho_neigh[p] = 0; sp1_overlap[p] = 0;
            sp1_pval[p] = 1.0; sp1_pval_div[p] = 1.0; sp1_effect[p] = 1.0;
            sp1_pval_gt[p] = 0.0; sp1_pval_eq[p] = 1.0;
            sp1_jaccard[p] = 0.0;
            sp2_neigh[p] = 0; sp2_ortho_neigh[p] = 0; sp2_overlap[p] = 0;
            sp2_pval[p] = 1.0; sp2_pval_div[p] = 1.0; sp2_effect[p] = 1.0;
            sp2_pval_gt[p] = 0.0; sp2_pval_eq[p] = 1.0;
            sp2_jaccard[p] = 0.0;
            continue;
        }

        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif

        // Direction 1: Species 1 -> Species 2
        DirectionResult d1 = compute_direction(
            neighbors1[g1_idx], neighbors2[g2_idx], sp2_to_sp1, n1, g1_idx,
            thread_flags1[tid]);
        sp1_neigh[p] = d1.neigh;   sp1_ortho_neigh[p] = d1.ortho_neigh;
        sp1_overlap[p] = d1.overlap; sp1_pval[p] = d1.pval_con;
        sp1_pval_div[p] = d1.pval_div; sp1_effect[p] = d1.effect_size;
        sp1_pval_gt[p] = d1.pval_gt; sp1_pval_eq[p] = d1.pval_eq;
        sp1_jaccard[p] = d1.jaccard;

        // Direction 2: Species 2 -> Species 1
        DirectionResult d2 = compute_direction(
            neighbors2[g2_idx], neighbors1[g1_idx], sp1_to_sp2, n2, g2_idx,
            thread_flags2[tid]);
        sp2_neigh[p] = d2.neigh;   sp2_ortho_neigh[p] = d2.ortho_neigh;
        sp2_overlap[p] = d2.overlap; sp2_pval[p] = d2.pval_con;
        sp2_pval_div[p] = d2.pval_div; sp2_effect[p] = d2.effect_size;
        sp2_pval_gt[p] = d2.pval_gt; sp2_pval_eq[p] = d2.pval_eq;
        sp2_jaccard[p] = d2.jaccard;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("Species1.neigh") = sp1_neigh,
        Rcpp::Named("Species1.ortho.neigh") = sp1_ortho_neigh,
        Rcpp::Named("Species1.neigh.overlap") = sp1_overlap,
        Rcpp::Named("Species1.p.val.con") = sp1_pval,
        Rcpp::Named("Species1.p.val.div") = sp1_pval_div,
        Rcpp::Named("Species1.p.val.gt") = sp1_pval_gt,
        Rcpp::Named("Species1.p.val.eq") = sp1_pval_eq,
        Rcpp::Named("Species1.effect.size") = sp1_effect,
        Rcpp::Named("Species1.jaccard") = sp1_jaccard,
        Rcpp::Named("Species2.neigh") = sp2_neigh,
        Rcpp::Named("Species2.ortho.neigh") = sp2_ortho_neigh,
        Rcpp::Named("Species2.neigh.overlap") = sp2_overlap,
        Rcpp::Named("Species2.p.val.con") = sp2_pval,
        Rcpp::Named("Species2.p.val.div") = sp2_pval_div,
        Rcpp::Named("Species2.p.val.gt") = sp2_pval_gt,
        Rcpp::Named("Species2.p.val.eq") = sp2_pval_eq,
        Rcpp::Named("Species2.effect.size") = sp2_effect,
        Rcpp::Named("Species2.jaccard") = sp2_jaccard
    );
}


//' Compare co-expression neighborhoods across species (integer-indexed)
//'
//' For each ortholog pair, tests the overlap of co-expression neighborhoods
//' in both directions (sp1->sp2 and sp2->sp1) using hypergeometric tests.
//' All gene identifiers are 0-based integer indices (string mapping done in R).
//'
//' @param net1 Co-expression network for species 1 (n1 x n1 matrix)
//' @param net2 Co-expression network for species 2 (n2 x n2 matrix)
//' @param thr1 Co-expression threshold for species 1
//' @param thr2 Co-expression threshold for species 2
//' @param pair_sp1_idx 0-based index into net1 for each ortholog pair
//' @param pair_sp2_idx 0-based index into net2 for each ortholog pair
//' @param ortho_sp1_idx 0-based net1 indices for full ortholog table
//' @param ortho_sp2_idx 0-based net2 indices for full ortholog table
//' @param n_cores Number of OpenMP threads (default: 1)
//' @return DataFrame with comparison results for each ortholog pair. The
//'   hypergeometric urn excludes the anchor gene (population n - 1, anchor
//'   dropped from the ortholog-mapped set); `*.p.val.gt` / `*.p.val.eq`
//'   are the ungated upper tail P(X > x) and point mass P(X = x).
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame compare_neighborhoods_cpp(
    const arma::mat& net1,
    const arma::mat& net2,
    double thr1,
    double thr2,
    const Rcpp::IntegerVector& pair_sp1_idx,
    const Rcpp::IntegerVector& pair_sp2_idx,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    int n_cores = 1
) {
    return compare_neighborhoods_core(
        neighbor_lists_dense(net1, thr1, n_cores),
        neighbor_lists_dense(net2, thr2, n_cores),
        pair_sp1_idx, pair_sp2_idx, ortho_sp1_idx, ortho_sp2_idx, n_cores);
}


//' Compare co-expression neighborhoods across species (sparse networks)
//'
//' Same as [compare_neighborhoods_cpp()] but takes the slots of a
//' `dgCMatrix` (column-compressed, both triangles stored) for each network
//' instead of a dense matrix. Column j lists the neighbours of gene j.
//'
//' @param p1 `@p` slot of net1 (column pointers, length n1 + 1)
//' @param i1 `@i` slot of net1 (0-based row indices)
//' @param x1 `@x` slot of net1 (stored values)
//' @param thr1 Co-expression threshold for species 1
//' @param p2 `@p` slot of net2
//' @param i2 `@i` slot of net2
//' @param x2 `@x` slot of net2
//' @param thr2 Co-expression threshold for species 2
//' @param pair_sp1_idx 0-based index into net1 for each ortholog pair
//' @param pair_sp2_idx 0-based index into net2 for each ortholog pair
//' @param ortho_sp1_idx 0-based net1 indices for full ortholog table
//' @param ortho_sp2_idx 0-based net2 indices for full ortholog table
//' @param n_cores Number of OpenMP threads (default: 1)
//' @return DataFrame with comparison results for each ortholog pair
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame compare_neighborhoods_sparse_cpp(
    const Rcpp::IntegerVector& p1,
    const Rcpp::IntegerVector& i1,
    const Rcpp::NumericVector& x1,
    double thr1,
    const Rcpp::IntegerVector& p2,
    const Rcpp::IntegerVector& i2,
    const Rcpp::NumericVector& x2,
    double thr2,
    const Rcpp::IntegerVector& pair_sp1_idx,
    const Rcpp::IntegerVector& pair_sp2_idx,
    const Rcpp::IntegerVector& ortho_sp1_idx,
    const Rcpp::IntegerVector& ortho_sp2_idx,
    int n_cores = 1
) {
    return compare_neighborhoods_core(
        neighbor_lists_sparse(p1, i1, x1, thr1, n_cores),
        neighbor_lists_sparse(p2, i2, x2, thr2, n_cores),
        pair_sp1_idx, pair_sp2_idx, ortho_sp1_idx, ortho_sp2_idx, n_cores);
}
