// module_jaccard_permutation.cpp
// Permutation-based Jaccard test for cross-species module comparison
//
// Tests whether pairs of co-expression modules share more ortholog-mapped
// genes than expected by chance, using a null that permutes the ortholog
// mapping (which sp1 genes map to which sp2 genes).
//
// Batched permutation: one Fisher-Yates shuffle per iteration, shared across
// all active pairs.  Each pair has its own Besag-Clifford exceedance counter,
// so pairs with clearly null Jaccard accumulate exceedances quickly and stop
// early while significant pairs keep running.
//
// p-value = (n_exceed + 1) / (n_perm + 1)
//
// Design constraints (matching hog_permutation.cpp):
// - Integer indices only (no strings — Homebrew clang ABI issue)
// - OpenMP parallelization over module pairs within each iteration
// - Single std::mt19937 seeded from R::runif for sequential shuffle

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <random>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// Fisher-Yates shuffle of a vector
static void fisher_yates_shuffle(std::vector<int>& vec, std::mt19937& rng) {
    for (auto i = std::ssize(vec) - 1; i > 0; --i) {
        std::uniform_int_distribution<int> dist(0, static_cast<int>(i));
        std::swap(vec[i], vec[dist(rng)]);
    }
}


//' Permutation-based Jaccard test for module comparison (batched)
//'
//' Tests each module pair for significant Jaccard overlap using ortholog-mapping
//' permutation with Besag-Clifford adaptive stopping.  Uses batched permutation:
//' one shuffle per iteration shared across all active pairs.
//'
//' @param ortho_sp1_gene 0-based index of each ortholog's sp1 unique gene
//' @param ortho_sp2_gene 0-based sp2 gene index for each ortholog row
//' @param n_sp1_unique Number of unique sp1 genes in ortholog table
//' @param n_sp2_universe Total number of sp2 genes in index space
//' @param mod1_sp1_genes List of integer vectors: sp1 unique gene indices per sp1 module
//' @param mod_sp2_sets List of integer vectors: sp2 indices per sp2 module
//' @param mod_i_idx 0-based index into mod1_sp1_genes for each pair
//' @param mod_j_idx 0-based index into mod_sp2_sets for each pair
//' @param obs_jaccard Observed Jaccard for each pair
//' @param min_exceedances Besag-Clifford stopping parameter
//' @param max_permutations Maximum permutations
//' @param n_cores Number of OpenMP threads
//' @return DataFrame with n_perm, n_exceed, p_value per pair
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame module_jaccard_permutation_cpp(
    const Rcpp::IntegerVector& ortho_sp1_gene,
    const Rcpp::IntegerVector& ortho_sp2_gene,
    int n_sp1_unique,
    int n_sp2_universe,
    const Rcpp::List& mod1_sp1_genes,
    const Rcpp::List& mod_sp2_sets,
    const Rcpp::IntegerVector& mod_i_idx,
    const Rcpp::IntegerVector& mod_j_idx,
    const Rcpp::NumericVector& obs_jaccard,
    int min_exceedances,
    int max_permutations,
    int n_cores
) {
    const int n_ortho = ortho_sp1_gene.size();
    const int n_pairs = mod_i_idx.size();
    const int n_mod1 = mod1_sp1_genes.size();
    const int n_mod2 = mod_sp2_sets.size();

    // ---- Extract module sets from R lists ----
    std::vector<std::vector<int>> cpp_mod1_sp1(n_mod1);
    for (int i = 0; i < n_mod1; ++i) {
        Rcpp::IntegerVector v = mod1_sp1_genes[i];
        cpp_mod1_sp1[i].assign(v.begin(), v.end());
    }

    std::vector<std::vector<int>> cpp_mod2(n_mod2);
    for (int j = 0; j < n_mod2; ++j) {
        Rcpp::IntegerVector v = mod_sp2_sets[j];
        cpp_mod2[j].assign(v.begin(), v.end());
    }

    // ---- Per sp1 unique gene: which ortho rows map from it ----
    std::vector<std::vector<int>> sp1_gene_to_ortho_rows(n_sp1_unique);
    for (int r = 0; r < n_ortho; ++r) {
        int g = ortho_sp1_gene[r];
        if (g >= 0 && g < n_sp1_unique) {
            sp1_gene_to_ortho_rows[g].push_back(r);
        }
    }

    // ---- Base sp2 mapping for shuffling ----
    std::vector<int> perm_sp2(ortho_sp2_gene.begin(), ortho_sp2_gene.end());

    // ---- Module j membership flags for fast intersection ----
    std::vector<std::vector<char>> mod2_flags(n_mod2);
    for (int j = 0; j < n_mod2; ++j) {
        mod2_flags[j].assign(n_sp2_universe, 0);
        for (int g : cpp_mod2[j]) {
            if (g >= 0 && g < n_sp2_universe) {
                mod2_flags[j][g] = 1;
            }
        }
    }

    // ---- Pre-compute unique sp1 module indices (deduplicated) ----
    std::vector<int> unique_mod_i;
    {
        std::vector<char> mod_i_seen(n_mod1, 0);
        for (int p = 0; p < n_pairs; ++p) {
            int mi = mod_i_idx[p];
            if (!mod_i_seen[mi]) {
                mod_i_seen[mi] = 1;
                unique_mod_i.push_back(mi);
            }
        }
    }
    auto n_unique_mod_i = static_cast<int>(unique_mod_i.size());

    std::vector<int> mod_i_to_unique(n_mod1, -1);
    for (int u = 0; u < n_unique_mod_i; ++u) {
        mod_i_to_unique[unique_mod_i[u]] = u;
    }

    // ---- Single RNG for sequential shuffle ----
    std::mt19937 rng(
        static_cast<uint32_t>(R::runif(0.0, 4294967296.0)));

    // ---- Per-pair output and state ----
    std::vector<int> out_n_perm(n_pairs, 0);
    std::vector<int> out_n_exceed(n_pairs, 0);
    std::vector<char> pair_active(n_pairs, 1);

    int n_active = 0;
    for (int p = 0; p < n_pairs; ++p) {
        if (obs_jaccard[p] <= 0.0) {
            pair_active[p] = 0;
        } else {
            ++n_active;
        }
    }

#ifdef _OPENMP
#endif

    // Seen buffer for sequential remap (main thread only)
    std::vector<char> seen(n_sp2_universe, 0);
    // Per-iteration remapped sp2 sets, computed sequentially after shuffle
    std::vector<std::vector<int>> mod_i_remapped(n_unique_mod_i);

    // ---- Batched permutation loop ----
    for (int iter = 0; iter < max_permutations && n_active > 0; ++iter) {
        if (iter % 128 == 0) Rcpp::checkUserInterrupt();

        fisher_yates_shuffle(perm_sp2, rng);

        // Remap each unique sp1 module through shuffled orthologs
        for (int u = 0; u < n_unique_mod_i; ++u) {
            auto& remapped = mod_i_remapped[u];
            remapped.clear();

            for (int sp1_g : cpp_mod1_sp1[unique_mod_i[u]]) {
                for (int row : sp1_gene_to_ortho_rows[sp1_g]) {
                    int new_sp2 = perm_sp2[row];
                    if (new_sp2 >= 0 && new_sp2 < n_sp2_universe &&
                        !seen[new_sp2]) {
                        seen[new_sp2] = 1;
                        remapped.push_back(new_sp2);
                    }
                }
            }
            for (int g : remapped) seen[g] = 0;
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(n_cores) if(n_cores > 1)
#endif
        for (int p = 0; p < n_pairs; ++p) {
            if (!pair_active[p]) continue;

            int u = mod_i_to_unique[mod_i_idx[p]];
            int mj = mod_j_idx[p];

            const auto& remapped = mod_i_remapped[u];
            const auto& m2flags = mod2_flags[mj];

            int perm_overlap = 0;
            for (int g : remapped) {
                if (m2flags[g]) ++perm_overlap;
            }
            int perm_union = static_cast<int>(remapped.size()) +
                             static_cast<int>(cpp_mod2[mj].size()) - perm_overlap;
            double perm_jaccard = (perm_union > 0) ?
                static_cast<double>(perm_overlap) / perm_union : 0.0;

            if (perm_jaccard >= obs_jaccard[p]) ++out_n_exceed[p];
            ++out_n_perm[p];
        }

        n_active = 0;
        for (int p = 0; p < n_pairs; ++p) {
            if (pair_active[p]) {
                if (out_n_exceed[p] >= min_exceedances) {
                    pair_active[p] = 0;
                } else {
                    ++n_active;
                }
            }
        }
    }

    // ---- Compute p-values ----
    Rcpp::IntegerVector r_n_perm(n_pairs);
    Rcpp::IntegerVector r_n_exceed(n_pairs);
    Rcpp::NumericVector r_p_value(n_pairs);

    for (int p = 0; p < n_pairs; ++p) {
        r_n_perm[p] = out_n_perm[p];
        r_n_exceed[p] = out_n_exceed[p];
        r_p_value[p] = (out_n_perm[p] > 0) ?
            static_cast<double>(out_n_exceed[p] + 1) / (out_n_perm[p] + 1) :
            1.0;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("n_perm") = r_n_perm,
        Rcpp::Named("n_exceed") = r_n_exceed,
        Rcpp::Named("p_value") = r_p_value
    );
}
