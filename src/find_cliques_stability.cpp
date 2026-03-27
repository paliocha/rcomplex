// find_cliques_stability.cpp
// Leave-k-out jackknife stability analysis for trait-exclusive cliques
//
// Subsets are drawn from all_species (the full analysis universe).
// Reduced cliques are found among remaining target species only.
// A clique is testable if at least 2 of its species remain active.
//
// OpenMP: parallelizes over species subsets within each k-level.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include "find_cliques_common.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


//' Leave-k-out jackknife stability for trait-exclusive cliques
//'
//' @param edge_hog,edge_g1,edge_g2,edge_sp1,edge_sp2 0-based edge vectors
//' @param edge_qval,edge_eff per-edge q-value and effect size
//' @param n_all_species total species in the analysis universe
//' @param n_hogs,n_genes total HOGs and genes
//' @param species_trait trait value per species (length n_all_species)
//' @param is_target 1/0 per species (length n_all_species): target membership
//' @param full_cliques output of find_cliques_cpp on the full dataset
//' @param max_k,max_genes_per_sp,jaccard_threshold,n_cores tuning parameters
//' @return List with stability, clique_disruption, stability_class, novel_cliques
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List find_cliques_stability_cpp(
    IntegerVector edge_hog, IntegerVector edge_g1, IntegerVector edge_g2,
    IntegerVector edge_sp1, IntegerVector edge_sp2,
    NumericVector edge_qval, NumericVector edge_eff,
    int n_all_species, int n_hogs, int n_genes,
    IntegerVector species_trait,
    IntegerVector is_target,
    Rcpp::List full_cliques,
    int max_k = 3, int max_genes_per_sp = 10,
    double jaccard_threshold = 0.8, int n_cores = 1)
{
    // ---- Validation ----
    if (n_all_species > 64)
        Rcpp::stop("n_all_species must be <= 64 (bitmask limit)");
    if (max_k < 1)
        Rcpp::stop("max_k must be >= 1");
    if (max_k >= n_all_species)
        Rcpp::stop("max_k must be < n_all_species");

    // ---- Build target_mask ----
    uint64_t target_mask = 0;
    int n_target = 0;
    for (int s = 0; s < n_all_species; ++s) {
        if (is_target[s]) {
            target_mask |= (1ULL << s);
            ++n_target;
        }
    }
    if (n_target < 2)
        Rcpp::stop("Need at least 2 target species");

    int ne = edge_hog.size();

    // ---- Raw pointers ----
    const int* hog_ptr  = edge_hog.begin();
    const int* g1_ptr   = edge_g1.begin();
    const int* g2_ptr   = edge_g2.begin();
    const int* sp1_ptr  = edge_sp1.begin();
    const int* sp2_ptr  = edge_sp2.begin();
    const double* qval_ptr = edge_qval.begin();
    const double* eff_ptr  = edge_eff.begin();
    const int* trait_ptr   = species_trait.begin();

    // ---- Group edges by HOG ----
    std::vector<std::vector<int>> hog_edges(n_hogs);
    for (int i = 0; i < ne; ++i)
        hog_edges[hog_ptr[i]].push_back(i);

    // ---- Unpack full_cliques ----
    IntegerVector fc_hog_idx = full_cliques["hog_idx"];
    IntegerMatrix fc_genes   = full_cliques["genes"];
    int n_full = fc_hog_idx.size();

    uint64_t full_mask = (n_all_species == 64) ? ~0ULL
                                               : (1ULL << n_all_species) - 1;

    std::vector<int> full_hog(n_full);
    std::vector<std::vector<int>> full_genes(n_full);
    std::vector<int> original_trait(n_full);

    for (int c = 0; c < n_full; ++c) {
        full_hog[c] = fc_hog_idx[c];
        full_genes[c].resize(n_all_species);
        for (int s = 0; s < n_all_species; ++s) {
            int g = fc_genes(c, s);
            full_genes[c][s] = (g == NA_INTEGER) ? -1 : g;
        }
        original_trait[c] = annotate_trait(
            full_genes[c], trait_ptr, target_mask, n_all_species);
    }

    // ---- HOG -> clique index ----
    std::vector<std::vector<int>> hog_to_cliques(n_hogs);
    for (int c = 0; c < n_full; ++c)
        hog_to_cliques[full_hog[c]].push_back(c);

    // ---- Identify exclusive cliques ----
    std::vector<int> exclusive_indices;
    std::vector<char> hog_has_excl(n_hogs, 0);

    for (int c = 0; c < n_full; ++c) {
        if (original_trait[c] >= 0) {
            exclusive_indices.push_back(c);
            hog_has_excl[full_hog[c]] = 1;
        }
    }

    std::vector<int> hogs_with_targets;
    for (int h = 0; h < n_hogs; ++h)
        if (hog_has_excl[h]) hogs_with_targets.push_back(h);

    int n_excl = static_cast<int>(exclusive_indices.size());

    // Early return if no exclusive cliques
    if (n_excl == 0) {
        return Rcpp::List::create(
            Named("stability") = Rcpp::DataFrame::create(
                Named("clique_idx")      = IntegerVector(0),
                Named("trait_value")     = IntegerVector(0),
                Named("k")              = IntegerVector(0),
                Named("n_subsets")      = IntegerVector(0),
                Named("n_stable")       = IntegerVector(0),
                Named("stability_score") = NumericVector(0),
                Named("sole_rep")       = Rcpp::LogicalVector(0)
            ),
            Named("clique_disruption") = Rcpp::DataFrame::create(
                Named("species_idx")          = IntegerVector(0),
                Named("trait_value")          = IntegerVector(0),
                Named("n_cliques_disrupted")  = IntegerVector(0)
            ),
            Named("stability_class") = IntegerVector(0),
            Named("novel_cliques")   = 0
        );
    }

    // ---- Reverse map: full index -> exclusive local index ----
    std::vector<int> fc_to_excl(n_full, -1);
    for (int e = 0; e < n_excl; ++e)
        fc_to_excl[exclusive_indices[e]] = e;

    // ---- sole_rep: count only TARGET species per trait ----
    int max_trait_val = 0;
    for (int s = 0; s < n_all_species; ++s)
        if (trait_ptr[s] > max_trait_val) max_trait_val = trait_ptr[s];

    std::vector<int> trait_species_count(max_trait_val + 1, 0);
    for (int s = 0; s < n_all_species; ++s) {
        if (((target_mask >> s) & 1) && trait_ptr[s] >= 0)
            trait_species_count[trait_ptr[s]]++;
    }

    std::vector<bool> sole_rep(n_excl, false);
    for (int e = 0; e < n_excl; ++e) {
        int t = original_trait[exclusive_indices[e]];
        if (t >= 0 && t <= max_trait_val && trait_species_count[t] == 1)
            sole_rep[e] = true;
    }

    // ---- Removal subsets from ALL species ----
    std::vector<std::vector<uint64_t>> all_subsets(max_k + 1);
    for (int k = 1; k <= max_k; ++k)
        all_subsets[k] = generate_subsets(n_all_species, k);

    // ---- Result accumulators ----
    std::vector<int> n_tested(n_excl * max_k, 0);
    std::vector<int> n_stable(n_excl * max_k, 0);
    std::vector<int> disruption_count(n_all_species, 0);
    int total_novel = 0;

    int max_threads = 1;
#ifdef _OPENMP
    if (n_cores > 1) max_threads = n_cores;
#endif

    // ---- Main loop over k-levels ----
    for (int k = 1; k <= max_k; ++k) {
        int n_subsets_k = static_cast<int>(all_subsets[k].size());

        std::vector<std::vector<bool>> thread_seen(
            max_threads, std::vector<bool>(n_genes, false));
        std::vector<std::vector<int>> thread_n_tested(
            max_threads, std::vector<int>(n_excl, 0));
        std::vector<std::vector<int>> thread_n_stable(
            max_threads, std::vector<int>(n_excl, 0));
        std::vector<int> thread_novel(max_threads, 0);

        std::vector<std::vector<int>> thread_disruption;
        if (k == 1) {
            thread_disruption.assign(
                max_threads, std::vector<int>(n_all_species, 0));
        }

#ifdef _OPENMP
        if (n_cores > 1) omp_set_num_threads(n_cores);
        #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
        for (int s = 0; s < n_subsets_k; ++s) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif

            uint64_t removal_mask = all_subsets[k][s];
            uint64_t active_mask = full_mask & ~removal_mask;
            uint64_t active_target = active_mask & target_mask;
            int n_active_target = __builtin_popcountll(active_target);

            // Need at least 2 target species for a meaningful clique
            if (n_active_target < 2) continue;

            for (int h : hogs_with_targets) {
                if (hog_edges[h].empty()) continue;

                // Early check: any exclusive clique has >= 2 active species?
                bool any_testable = false;
                for (int fc_idx : hog_to_cliques[h]) {
                    if (original_trait[fc_idx] < 0) continue;
                    int ca = 0;
                    for (int sp = 0; sp < n_all_species; ++sp) {
                        if (full_genes[fc_idx][sp] >= 0 &&
                            ((active_target >> sp) & 1))
                            ++ca;
                    }
                    if (ca >= 2) { any_testable = true; break; }
                }
                if (!any_testable) continue;

                // Find cliques among active target species (min 2)
                auto reduced = find_cliques_for_hog(
                    h, hog_edges[h], g1_ptr, g2_ptr, sp1_ptr, sp2_ptr,
                    qval_ptr, eff_ptr, active_target, n_all_species,
                    2, max_genes_per_sp, thread_seen[tid]);

                std::vector<char> rc_matched(reduced.size(), 0);

                // Match each full exclusive clique to a reduced clique
                for (int fc_idx : hog_to_cliques[h]) {
                    if (original_trait[fc_idx] < 0) continue;

                    // Testable: at least 2 of this clique's species remain
                    int clique_active = 0;
                    for (int sp = 0; sp < n_all_species; ++sp) {
                        if (full_genes[fc_idx][sp] >= 0 &&
                            ((active_target >> sp) & 1))
                            ++clique_active;
                    }
                    if (clique_active < 2) continue;

                    int excl_local = fc_to_excl[fc_idx];

                    bool found_stable = false;
                    for (size_t ri = 0; ri < reduced.size(); ++ri) {
                        double jacc = jaccard_remaining(
                            full_genes[fc_idx], reduced[ri].genes,
                            removal_mask, n_all_species);
                        if (jacc >= jaccard_threshold) {
                            int reduced_trait = annotate_trait(
                                reduced[ri].genes, trait_ptr,
                                active_target, n_all_species);
                            if (reduced_trait == original_trait[fc_idx]) {
                                found_stable = true;
                                rc_matched[ri] = 1;
                                break;
                            }
                        }
                    }

                    thread_n_tested[tid][excl_local]++;
                    if (found_stable) thread_n_stable[tid][excl_local]++;

                    if (k == 1 && !found_stable)
                        thread_disruption[tid][__builtin_ctzll(removal_mask)]++;
                }

                for (size_t r = 0; r < reduced.size(); ++r) {
                    if (!rc_matched[r]) thread_novel[tid]++;
                }
            }
        }

        Rcpp::checkUserInterrupt();

        // Merge thread results
        for (int t = 0; t < max_threads; ++t) {
            for (int e = 0; e < n_excl; ++e) {
                int flat = e * max_k + (k - 1);
                n_tested[flat] += thread_n_tested[t][e];
                n_stable[flat] += thread_n_stable[t][e];
                thread_n_tested[t][e] = 0;
                thread_n_stable[t][e] = 0;
            }
            total_novel += thread_novel[t];
            thread_novel[t] = 0;
            if (k == 1) {
                for (int sp = 0; sp < n_all_species; ++sp)
                    disruption_count[sp] += thread_disruption[t][sp];
            }
        }
    }

    // ---- Build output: stability DataFrame ----
    std::vector<int> out_clique_idx, out_trait_value, out_k;
    std::vector<int> out_n_subsets, out_n_stab;
    std::vector<double> out_stability_score;
    std::vector<int> out_sole_rep;

    for (int e = 0; e < n_excl; ++e) {
        int fc_idx = exclusive_indices[e];
        for (int k = 1; k <= max_k; ++k) {
            int flat = e * max_k + (k - 1);
            if (n_tested[flat] == 0) continue;
            out_clique_idx.push_back(fc_idx);
            out_trait_value.push_back(original_trait[fc_idx]);
            out_k.push_back(k);
            out_n_subsets.push_back(n_tested[flat]);
            out_n_stab.push_back(n_stable[flat]);
            out_stability_score.push_back(
                static_cast<double>(n_stable[flat]) / n_tested[flat]);
            out_sole_rep.push_back(sole_rep[e] ? 1 : 0);
        }
    }

    int n_stab_rows = static_cast<int>(out_clique_idx.size());
    Rcpp::IntegerVector r_clique_idx(out_clique_idx.begin(), out_clique_idx.end());
    Rcpp::IntegerVector r_trait_value(out_trait_value.begin(), out_trait_value.end());
    Rcpp::IntegerVector r_k(out_k.begin(), out_k.end());
    Rcpp::IntegerVector r_n_subsets(out_n_subsets.begin(), out_n_subsets.end());
    Rcpp::IntegerVector r_n_stable(out_n_stab.begin(), out_n_stab.end());
    Rcpp::NumericVector r_stability_score(out_stability_score.begin(),
                                          out_stability_score.end());
    Rcpp::LogicalVector r_sole_rep(n_stab_rows);
    for (int i = 0; i < n_stab_rows; ++i)
        r_sole_rep[i] = (out_sole_rep[i] == 1);

    Rcpp::DataFrame stability_df = Rcpp::DataFrame::create(
        Named("clique_idx")      = r_clique_idx,
        Named("trait_value")     = r_trait_value,
        Named("k")              = r_k,
        Named("n_subsets")      = r_n_subsets,
        Named("n_stable")       = r_n_stable,
        Named("stability_score") = r_stability_score,
        Named("sole_rep")       = r_sole_rep
    );

    // ---- clique_disruption: one row per ALL species ----
    Rcpp::IntegerVector r_sp_idx(n_all_species);
    Rcpp::IntegerVector r_sp_trait(n_all_species);
    Rcpp::IntegerVector r_disrupted(n_all_species);
    for (int sp = 0; sp < n_all_species; ++sp) {
        r_sp_idx[sp]    = sp;
        r_sp_trait[sp]  = trait_ptr[sp];
        r_disrupted[sp] = disruption_count[sp];
    }

    Rcpp::DataFrame disruption_df = Rcpp::DataFrame::create(
        Named("species_idx")         = r_sp_idx,
        Named("trait_value")         = r_sp_trait,
        Named("n_cliques_disrupted") = r_disrupted
    );

    // ---- stability_class: highest k where all tested subsets are stable ----
    Rcpp::IntegerVector stability_class_vec(n_excl);
    for (int e = 0; e < n_excl; ++e) {
        int highest_stable_k = 0;
        for (int k = 1; k <= max_k; ++k) {
            int flat = e * max_k + (k - 1);
            if (n_tested[flat] == 0) break;
            if (n_stable[flat] == n_tested[flat])
                highest_stable_k = k;
            else
                break;
        }
        stability_class_vec[e] = highest_stable_k;
    }

    return Rcpp::List::create(
        Named("stability")        = stability_df,
        Named("clique_disruption") = disruption_df,
        Named("stability_class")  = stability_class_vec,
        Named("novel_cliques")    = total_novel
    );
}
