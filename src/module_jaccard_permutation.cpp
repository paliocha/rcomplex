// module_jaccard_permutation.cpp
// Permutation-based Jaccard test for cross-species module comparison
//
// Tests whether pairs of co-expression modules share more ortholog-mapped
// genes than expected by chance, using a null that permutes the ortholog
// mapping (which sp1 genes map to which sp2 genes).
//
// Adaptive stopping: Besag & Clifford (1991) sequential Monte Carlo p-values.
// p-value = (n_exceed + 1) / (n_perm + 1)
//
// Design constraints (matching hog_permutation.cpp):
// - Integer indices only (no strings — Homebrew clang ABI issue)
// - OpenMP parallelization over module pairs
// - Per-thread std::mt19937 seeded from R::runif before parallel region

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>
#include <random>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// Fisher-Yates shuffle of a vector
static void fisher_yates_shuffle(std::vector<int>& vec, std::mt19937& rng) {
    int n = static_cast<int>(vec.size());
    for (int i = n - 1; i > 0; --i) {
        std::uniform_int_distribution<int> dist(0, i);
        std::swap(vec[i], vec[dist(rng)]);
    }
}


//' Permutation-based Jaccard test for module comparison
//'
//' Tests each module pair for significant Jaccard overlap using ortholog-mapping
//' permutation with Besag-Clifford adaptive stopping.
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
    std::vector<int> base_sp2(n_ortho);
    for (int r = 0; r < n_ortho; ++r) {
        base_sp2[r] = ortho_sp2_gene[r];
    }

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

    // ---- Per-thread resources ----
    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif

    std::vector<uint32_t> seeds(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        seeds[t] = static_cast<uint32_t>(R::runif(0.0, 4294967296.0));
    }

    std::vector<std::mt19937> thread_rng(max_threads);
    for (int t = 0; t < max_threads; ++t) thread_rng[t].seed(seeds[t]);

    std::vector<std::vector<int>> thread_sp2(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        thread_sp2[t] = base_sp2;
    }

    std::vector<std::vector<char>> thread_seen(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        thread_seen[t].assign(n_sp2_universe, 0);
    }

    std::vector<std::vector<int>> thread_mapped(max_threads);

    // ---- Output vectors ----
    Rcpp::IntegerVector out_n_perm(n_pairs);
    Rcpp::IntegerVector out_n_exceed(n_pairs);
    Rcpp::NumericVector out_p_value(n_pairs);

    // ---- Permutation loop ----
#ifdef _OPENMP
    if (n_cores > 1) omp_set_num_threads(n_cores);
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (int p = 0; p < n_pairs; ++p) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& rng = thread_rng[tid];
        auto& perm_sp2 = thread_sp2[tid];
        auto& seen = thread_seen[tid];
        auto& perm_mapped = thread_mapped[tid];

        int mi = mod_i_idx[p];
        int mj = mod_j_idx[p];
        double obs_j = obs_jaccard[p];

        if (obs_j <= 0.0) {
            out_n_perm[p] = 0;
            out_n_exceed[p] = 0;
            out_p_value[p] = 1.0;
            continue;
        }

        const auto& sp1_genes = cpp_mod1_sp1[mi];
        int mod2_size = static_cast<int>(cpp_mod2[mj].size());
        const auto& m2flags = mod2_flags[mj];

        int n_exceed = 0;
        int n_perm = 0;

        while (n_perm < max_permutations && n_exceed < min_exceedances) {
            fisher_yates_shuffle(perm_sp2, rng);

            // Remap sp1 module genes through shuffled orthologs
            perm_mapped.clear();
            for (int sp1_g : sp1_genes) {
                for (int row : sp1_gene_to_ortho_rows[sp1_g]) {
                    int new_sp2 = perm_sp2[row];
                    if (new_sp2 >= 0 && new_sp2 < n_sp2_universe &&
                        !seen[new_sp2]) {
                        seen[new_sp2] = 1;
                        perm_mapped.push_back(new_sp2);
                    }
                }
            }

            // Compute Jaccard with module j
            int perm_overlap = 0;
            for (int g : perm_mapped) {
                if (m2flags[g]) ++perm_overlap;
            }
            int perm_union = static_cast<int>(perm_mapped.size()) +
                             mod2_size - perm_overlap;
            double perm_jaccard = (perm_union > 0) ?
                static_cast<double>(perm_overlap) / perm_union : 0.0;

            for (int g : perm_mapped) seen[g] = 0;

            if (perm_jaccard >= obs_j) ++n_exceed;
            ++n_perm;
        }

        out_n_perm[p] = n_perm;
        out_n_exceed[p] = n_exceed;
        out_p_value[p] = static_cast<double>(n_exceed + 1) / (n_perm + 1);
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("n_perm") = out_n_perm,
        Rcpp::Named("n_exceed") = out_n_exceed,
        Rcpp::Named("p_value") = out_p_value
    );
}
