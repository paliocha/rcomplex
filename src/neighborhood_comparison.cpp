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

using namespace Rcpp;


// Per-direction hypergeometric test result
struct DirectionResult {
    int neigh;
    int ortho_neigh;
    int overlap;
    double pval_con;
    double pval_div;
    double effect_size;
};

// Compute one direction of the neighborhood comparison.
//   anchor_neigh:  neighbors of the anchor gene in its own network
//   other_neigh:   neighbors of the paired gene in the other network
//   mapping:       other-network index -> vector of anchor-network indices
//   N:             total genes in the anchor network
//   anchor_flags:  scratch boolean vector (size N), must be pre-cleared
static DirectionResult compute_direction(
    const std::vector<int>& anchor_neigh,
    const std::vector<int>& other_neigh,
    const std::vector<std::vector<int>>& mapping,
    int N,
    std::vector<char>& anchor_flags
) {
    DirectionResult r;
    r.neigh = static_cast<int>(anchor_neigh.size());
    r.pval_con = 1.0;
    r.pval_div = 1.0;
    r.effect_size = 1.0;

    // Map other-network neighbors back to anchor-network via orthologs
    r.ortho_neigh = 0;
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) {
            if (anchor_flags[idx] == 0) {
                anchor_flags[idx] = 1;
                ++r.ortho_neigh;
            }
        }
    }

    // Count intersection
    r.overlap = 0;
    for (int nb : anchor_neigh) {
        if (anchor_flags[nb] != 0) ++r.overlap;
    }

    // Clear flags for reuse
    for (int nb : other_neigh) {
        for (int idx : mapping[nb]) anchor_flags[idx] = 0;
    }

    int m = r.neigh, k = r.ortho_neigh, x = r.overlap;
    if (x > 1) {
        r.pval_con = R::phyper(x - 1, m, N - m, k, 0, 0);
    }
    if (k > 0 && m > 0) {
        r.effect_size = (static_cast<double>(x) / k) / (static_cast<double>(m) / N);
        r.pval_div = R::phyper(x, m, N - m, k, 1, 0);
    }

    return r;
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
//' @return DataFrame with comparison results for each ortholog pair
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
    const int n_pairs = pair_sp1_idx.size();
    const int n1 = static_cast<int>(net1.n_rows);
    const int n2 = static_cast<int>(net2.n_rows);
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

#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

    // Precompute neighbor lists for each gene in each network
    std::vector<std::vector<int>> neighbors1(n1);
    std::vector<std::vector<int>> neighbors2(n2);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            if (i != j && net1(i, j) >= thr1) {
                neighbors1[i].push_back(j);
            }
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (i != j && net2(i, j) >= thr2) {
                neighbors2[i].push_back(j);
            }
        }
    }

    // Output vectors (6 per direction)
    Rcpp::IntegerVector sp1_neigh(n_pairs),  sp2_neigh(n_pairs);
    Rcpp::IntegerVector sp1_ortho_neigh(n_pairs), sp2_ortho_neigh(n_pairs);
    Rcpp::IntegerVector sp1_overlap(n_pairs), sp2_overlap(n_pairs);
    Rcpp::NumericVector sp1_pval(n_pairs),   sp2_pval(n_pairs);
    Rcpp::NumericVector sp1_pval_div(n_pairs), sp2_pval_div(n_pairs);
    Rcpp::NumericVector sp1_effect(n_pairs), sp2_effect(n_pairs);

    // Parallel comparison loop
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int p = 0; p < n_pairs; ++p) {
        int g1_idx = pair_sp1_idx[p];
        int g2_idx = pair_sp2_idx[p];

        if (g1_idx < 0 || g1_idx >= n1 || g2_idx < 0 || g2_idx >= n2) {
            sp1_neigh[p] = 0; sp1_ortho_neigh[p] = 0; sp1_overlap[p] = 0;
            sp1_pval[p] = 1.0; sp1_pval_div[p] = 1.0; sp1_effect[p] = 1.0;
            sp2_neigh[p] = 0; sp2_ortho_neigh[p] = 0; sp2_overlap[p] = 0;
            sp2_pval[p] = 1.0; sp2_pval_div[p] = 1.0; sp2_effect[p] = 1.0;
            continue;
        }

        // Scratch flag vectors (one per network size, reused across directions)
        std::vector<char> flags1(n1, 0);
        std::vector<char> flags2(n2, 0);

        // Direction 1: Species 1 -> Species 2
        DirectionResult d1 = compute_direction(
            neighbors1[g1_idx], neighbors2[g2_idx], sp2_to_sp1, n1, flags1);
        sp1_neigh[p] = d1.neigh;   sp1_ortho_neigh[p] = d1.ortho_neigh;
        sp1_overlap[p] = d1.overlap; sp1_pval[p] = d1.pval_con;
        sp1_pval_div[p] = d1.pval_div; sp1_effect[p] = d1.effect_size;

        // Direction 2: Species 2 -> Species 1
        DirectionResult d2 = compute_direction(
            neighbors2[g2_idx], neighbors1[g1_idx], sp1_to_sp2, n2, flags2);
        sp2_neigh[p] = d2.neigh;   sp2_ortho_neigh[p] = d2.ortho_neigh;
        sp2_overlap[p] = d2.overlap; sp2_pval[p] = d2.pval_con;
        sp2_pval_div[p] = d2.pval_div; sp2_effect[p] = d2.effect_size;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("Species1.neigh") = sp1_neigh,
        Rcpp::Named("Species1.ortho.neigh") = sp1_ortho_neigh,
        Rcpp::Named("Species1.neigh.overlap") = sp1_overlap,
        Rcpp::Named("Species1.p.val.con") = sp1_pval,
        Rcpp::Named("Species1.p.val.div") = sp1_pval_div,
        Rcpp::Named("Species1.effect.size") = sp1_effect,
        Rcpp::Named("Species2.neigh") = sp2_neigh,
        Rcpp::Named("Species2.ortho.neigh") = sp2_ortho_neigh,
        Rcpp::Named("Species2.neigh.overlap") = sp2_overlap,
        Rcpp::Named("Species2.p.val.con") = sp2_pval,
        Rcpp::Named("Species2.p.val.div") = sp2_pval_div,
        Rcpp::Named("Species2.effect.size") = sp2_effect
    );
}
