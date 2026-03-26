// find_cliques.cpp
// Fast EVOTREE-style clique finding via two-level decomposition
//
// Level 1: Species-level maximal cliques (Bron-Kerbosch on <= 16 species)
// Level 2: Best gene assignment per species clique (backtracking with pruning)
//
// This avoids the combinatorial explosion of igraph::max_cliques() on
// gene-level graphs with many paralogs and low min_species thresholds.
//
// The core pipeline is in find_cliques_common.h (shared with
// find_cliques_stability.cpp).  This file provides the Rcpp-exported wrapper.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "find_cliques_common.h"

using namespace Rcpp;

//' Find co-expression cliques using two-level decomposition
//'
//' @param edge_hog 0-based HOG index per edge
//' @param edge_g1  0-based gene index for gene 1
//' @param edge_g2  0-based gene index for gene 2
//' @param edge_sp1 0-based species index for gene 1
//' @param edge_sp2 0-based species index for gene 2
//' @param edge_qval q-value per edge
//' @param edge_eff Effect size per edge
//' @param n_target_species Number of target species
//' @param min_species Minimum species for a valid clique
//' @param n_hogs Total number of distinct HOGs
//' @param n_genes Total number of distinct genes
//' @param max_genes_per_sp Maximum genes per species per HOG (default 10)
//' @return List with: hog_idx (0-based), genes (matrix, 0-based or NA),
//'   n_species, mean_q, max_q, mean_effect_size, n_edges
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List find_cliques_cpp(
    IntegerVector edge_hog,
    IntegerVector edge_g1,
    IntegerVector edge_g2,
    IntegerVector edge_sp1,
    IntegerVector edge_sp2,
    NumericVector edge_qval,
    NumericVector edge_eff,
    int n_target_species,
    int min_species,
    int n_hogs,
    int n_genes,
    int max_genes_per_sp = 10)
{
    if (n_target_species > 64)
        Rcpp::stop("n_target_species must be <= 64");

    int ne = edge_hog.size();

    // Extract raw pointers for fast access
    const int* hog_ptr = edge_hog.begin();
    const int* g1_ptr  = edge_g1.begin();
    const int* g2_ptr  = edge_g2.begin();
    const int* sp1_ptr = edge_sp1.begin();
    const int* sp2_ptr = edge_sp2.begin();
    const double* qval_ptr = edge_qval.begin();
    const double* eff_ptr = edge_eff.begin();

    // Group edges by HOG
    std::vector<std::vector<int>> hog_edges(n_hogs);
    for (int i = 0; i < ne; i++) {
        hog_edges[hog_ptr[i]].push_back(i);
    }

    // Full species mask: all species active
    uint64_t full_mask = (n_target_species == 64) ? ~0ULL : (1ULL << n_target_species) - 1;

    // Caller-provided flag vector for gene deduplication (reused across HOGs)
    std::vector<bool> seen_genes(n_genes, false);

    // Result accumulators
    std::vector<int> res_hog;
    std::vector<std::vector<int>> res_genes;
    std::vector<int> res_n_sp;
    std::vector<double> res_mean_q, res_max_q, res_mean_eff;
    std::vector<int> res_n_edges;

    for (int h = 0; h < n_hogs; h++) {
        if (hog_edges[h].empty()) continue;

        // Check for user interrupt every 500 HOGs
        if (h % 500 == 0) Rcpp::checkUserInterrupt();

        // Run the full per-HOG clique pipeline
        std::vector<CliqueResult> cliques = find_cliques_for_hog(
            h, hog_edges[h],
            g1_ptr, g2_ptr, sp1_ptr, sp2_ptr, qval_ptr, eff_ptr,
            full_mask, n_target_species, min_species, max_genes_per_sp,
            seen_genes);

        // Collect results
        for (auto& cr : cliques) {
            res_hog.push_back(cr.hog_idx);
            res_genes.push_back(std::move(cr.genes));
            res_n_sp.push_back(cr.n_species);
            res_mean_q.push_back(cr.mean_q);
            res_max_q.push_back(cr.max_q);
            res_mean_eff.push_back(cr.mean_effect);
            res_n_edges.push_back(cr.n_edges);
        }
    }

    // --- Convert to R structures ---
    int nr = static_cast<int>(res_hog.size());
    IntegerVector out_hog(nr);
    IntegerMatrix out_genes(nr, n_target_species);
    IntegerVector out_n_sp(nr);
    NumericVector out_mean_q(nr), out_max_q(nr), out_mean_eff(nr);
    IntegerVector out_n_edges(nr);

    std::fill(out_genes.begin(), out_genes.end(), NA_INTEGER);

    for (int i = 0; i < nr; i++) {
        out_hog[i] = res_hog[i];
        for (int j = 0; j < n_target_species; j++) {
            if (res_genes[i][j] >= 0) {
                out_genes(i, j) = res_genes[i][j];
            }
        }
        out_n_sp[i] = res_n_sp[i];
        out_mean_q[i] = res_mean_q[i];
        out_max_q[i] = res_max_q[i];
        out_mean_eff[i] = res_mean_eff[i];
        out_n_edges[i] = res_n_edges[i];
    }

    return List::create(
        Named("hog_idx") = out_hog,
        Named("genes") = out_genes,
        Named("n_species") = out_n_sp,
        Named("mean_q") = out_mean_q,
        Named("max_q") = out_max_q,
        Named("mean_effect_size") = out_mean_eff,
        Named("n_edges") = out_n_edges
    );
}
