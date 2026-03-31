// find_cliques_common.h
// Shared helpers for two-level clique decomposition:
//   Level 1: Species-level maximal cliques (Bron-Kerbosch on <= 64 species)
//   Level 2: Best gene assignment per species clique (backtracking with pruning)
//
// All functions are inline so this header can be included from multiple TUs
// (find_cliques.cpp and find_cliques_stability.cpp).
//
// Design choices:
// - No std::unordered_map / std::unordered_set (ABI safety on Homebrew clang)
// - Sorted vector + binary search for edge lookups (maps are tiny per HOG)
// - Flag-vector (caller-provided std::vector<bool>) instead of hash sets

#ifndef RCOMPLEX_FIND_CLIQUES_COMMON_H
#define RCOMPLEX_FIND_CLIQUES_COMMON_H

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Canonical pair key for an unordered gene pair
// ---------------------------------------------------------------------------
// Packs two 32-bit gene indices into a single 64-bit key with the smaller
// index in the upper 32 bits.  Ensures edge_key(a, b) == edge_key(b, a).
inline int64_t edge_key(int g1, int g2) {
    if (g1 > g2) std::swap(g1, g2);
    return (static_cast<int64_t>(g1) << 32) | static_cast<int64_t>(static_cast<unsigned int>(g2));
}

// ---------------------------------------------------------------------------
// Sorted-vector lookup (replaces std::unordered_map for ABI safety)
// ---------------------------------------------------------------------------
// The sorted vectors are tiny (~tens of entries per HOG), so binary search
// overhead is negligible compared to hash-map construction.
inline double lookup_sorted(const std::vector<std::pair<int64_t, double>>& sorted_map,
                            int64_t key, bool& found) {
    auto it = std::lower_bound(sorted_map.begin(), sorted_map.end(), key,
        [](const std::pair<int64_t, double>& p, int64_t k) { return p.first < k; });
    found = (it != sorted_map.end() && it->first == key);
    return found ? it->second : 0.0;
}

// ---------------------------------------------------------------------------
// AssignResult: output of gene-assignment backtracking for one species clique
// ---------------------------------------------------------------------------
struct AssignResult {
    std::vector<int> genes;   // one gene per species in the clique
    double sum_q;
    double max_q;
    double sum_effect;
    int n_edges;
    int n_missing;            // number of missing edges in this assignment
    bool found;
};

// ---------------------------------------------------------------------------
// Bron-Kerbosch with pivoting for species-level maximal cliques
// ---------------------------------------------------------------------------
// Operates on a small graph (<= 64 nodes). Finds all maximal cliques
// with at least min_size vertices.
inline void bron_kerbosch(
    std::vector<int>& R,
    std::vector<int> P,
    std::vector<int> X,
    const std::vector<std::vector<bool>>& adj,
    int min_size,
    std::vector<std::vector<int>>& out_cliques)
{
    if (P.empty() && X.empty()) {
        if (static_cast<int>(R.size()) >= min_size)
            out_cliques.push_back(R);
        return;
    }

    // Pivot: vertex in P U X maximizing |N(u) intersect P|
    int pivot = -1, best_count = -1;
    for (int u : P) {
        int c = 0;
        for (int v : P) if (adj[u][v]) c++;
        if (c > best_count) { best_count = c; pivot = u; }
    }
    for (int u : X) {
        int c = 0;
        for (int v : P) if (adj[u][v]) c++;
        if (c > best_count) { best_count = c; pivot = u; }
    }

    // Candidates = P \ N(pivot)
    std::vector<int> cands;
    for (int v : P) {
        if (pivot < 0 || !adj[pivot][v])
            cands.push_back(v);
    }

    for (int v : cands) {
        R.push_back(v);
        std::vector<int> newP, newX;
        for (int u : P) if (adj[v][u]) newP.push_back(u);
        for (int u : X) if (adj[v][u]) newX.push_back(u);
        bron_kerbosch(R, newP, newX, adj, min_size, out_cliques);
        R.pop_back();
        P.erase(std::find(P.begin(), P.end(), v));
        X.push_back(v);
    }
}

// ---------------------------------------------------------------------------
// Gene assignment backtracking
// ---------------------------------------------------------------------------
// For a given species-level clique, finds the gene assignment (one gene per
// species) with the lowest mean FDR across all present edges.
// When max_missing == 0 (default): prunes if ANY edge is missing.
// When max_missing > 0: tolerates up to max_missing missing edges,
// preferring assignments with fewer missing edges, then lower mean q-value.
//
// Uses sorted vectors + binary search instead of std::unordered_map.
inline void assign_bt(
    int depth,
    int k,
    const std::vector<int>& sp_order,
    const std::vector<std::vector<int>>& sp_genes,
    const std::vector<std::pair<int64_t, double>>& qval_sorted,
    const std::vector<std::pair<int64_t, double>>& eff_sorted,
    std::vector<int>& current,
    double cur_sum_q,
    double cur_max_q,
    double cur_sum_eff,
    int cur_missing,
    int max_missing,
    AssignResult& best,
    int& iterations,
    int max_iterations)
{
    if (++iterations > max_iterations) return;

    if (depth == k) {
        int total_possible = k * (k - 1) / 2;
        int present_edges = total_possible - cur_missing;
        if (present_edges == 0) return;  // no edges at all — skip
        double mean_q = cur_sum_q / present_edges;
        // Prefer fewer missing edges first, then lower mean q-value
        bool is_better = false;
        if (!best.found) {
            is_better = true;
        } else if (cur_missing < best.n_missing) {
            is_better = true;
        } else if (cur_missing == best.n_missing) {
            double best_mean = best.sum_q / best.n_edges;
            is_better = (mean_q < best_mean);
        }
        if (is_better) {
            best.genes = current;
            best.sum_q = cur_sum_q;
            best.max_q = cur_max_q;
            best.sum_effect = cur_sum_eff;
            best.n_edges = present_edges;
            best.n_missing = cur_missing;
            best.found = true;
        }
        return;
    }

    int sp = sp_order[depth];
    for (int gene : sp_genes[sp]) {
        double new_sum_q = cur_sum_q;
        double new_max_q = cur_max_q;
        double new_sum_eff = cur_sum_eff;
        int new_missing = cur_missing;
        bool over_budget = false;

        for (int j = 0; j < depth; j++) {
            int64_t key = edge_key(current[j], gene);
            bool qval_found = false;
            double qval = lookup_sorted(qval_sorted, key, qval_found);
            if (!qval_found) {
                new_missing++;
                if (new_missing > max_missing) { over_budget = true; break; }
            } else {
                new_sum_q += qval;
                if (qval > new_max_q) new_max_q = qval;
                bool eff_found = false;
                double eff_val = lookup_sorted(eff_sorted, key, eff_found);
                if (eff_found) new_sum_eff += eff_val;
            }
        }
        if (over_budget) continue;

        current[depth] = gene;
        assign_bt(depth + 1, k, sp_order, sp_genes, qval_sorted, eff_sorted,
                  current, new_sum_q, new_max_q, new_sum_eff,
                  new_missing, max_missing,
                  best, iterations, max_iterations);

        if (iterations > max_iterations) return;
    }
}

// ---------------------------------------------------------------------------
// CliqueResult: output of per-HOG clique detection
// ---------------------------------------------------------------------------
struct CliqueResult {
    int hog_idx;
    std::vector<int> genes;   // length n_target_species, -1 = absent
    int n_species;
    double mean_q, max_q, mean_effect;
    int n_edges;
    int n_missing;            // number of missing edges (0 when max_missing_edges=0)
};

// ---------------------------------------------------------------------------
// find_cliques_for_hog: full per-HOG clique pipeline
// ---------------------------------------------------------------------------
// Encapsulates the entire two-level decomposition for one HOG:
//   1. Collect genes per species from edges, filtering by species_mask
//   2. Build species adjacency matrix (bool)
//   3. Cap genes per species (keep most-connected)
//   4. Run Bron-Kerbosch on species graph -> species-level cliques
//   5. For each species clique: sort by gene count, run assign_bt
//   6. Return vector of CliqueResult
//
// The seen_genes flag vector is caller-provided (size n_genes, all false on
// entry). Used entries are cleared before return, following the flag-vector
// pattern from hog_permutation.cpp.
//
// species_mask (uint64_t bitmask) controls which species are active. Only edges where
// BOTH endpoint species have their bit set in the mask are considered.
inline std::vector<CliqueResult> find_cliques_for_hog(
    int hog_idx,
    const std::vector<int>& hog_edge_indices,
    const int* edge_g1, const int* edge_g2,
    const int* edge_sp1, const int* edge_sp2,
    const double* edge_qval, const double* edge_eff,
    uint64_t species_mask,
    int n_target_species,
    int min_species,
    int max_genes_per_sp,
    int max_missing_edges,
    std::vector<bool>& seen_genes,   // caller-provided, size n_genes
    int max_iterations = 2000000)
{
    std::vector<CliqueResult> results;

    // --- 1. Collect genes per species and build species adjacency ---
    std::vector<std::vector<int>> sp_genes(n_target_species);

    // Track which genes we've added (to restore seen_genes on exit)
    std::vector<int> genes_touched;

    // Species adjacency
    std::vector<std::vector<bool>> sp_adj(
        n_target_species, std::vector<bool>(n_target_species, false));

    // Collect raw (key, q-value, effect) triples from filtered edges.
    // FDR and effect are paired by index so we can deduplicate correctly.
    std::vector<std::pair<int64_t, double>> qval_pairs;
    std::vector<std::pair<int64_t, double>> eff_pairs;

    for (int ei : hog_edge_indices) {
        int g1 = edge_g1[ei], g2 = edge_g2[ei];
        int s1 = edge_sp1[ei], s2 = edge_sp2[ei];

        // Filter by species_mask: both endpoint species must be active
        if (!((species_mask >> s1) & 1) || !((species_mask >> s2) & 1)) continue;
        // Skip same-species edges
        if (s1 == s2) continue;

        sp_adj[s1][s2] = sp_adj[s2][s1] = true;

        if (!seen_genes[g1]) {
            sp_genes[s1].push_back(g1);
            seen_genes[g1] = true;
            genes_touched.push_back(g1);
        }
        if (!seen_genes[g2]) {
            sp_genes[s2].push_back(g2);
            seen_genes[g2] = true;
            genes_touched.push_back(g2);
        }

        int64_t key = edge_key(g1, g2);
        qval_pairs.push_back({key, edge_qval[ei]});
        eff_pairs.push_back({key, edge_eff[ei]});
    }

    // Clear seen_genes for all touched genes (restore to false)
    for (int g : genes_touched) {
        seen_genes[g] = false;
    }

    // If no edges passed the filter, return empty
    if (qval_pairs.empty()) return results;

    // --- Deduplicate edge maps: sort by key, keep min-FDR entry ---
    // We need FDR and effect paired, so sort by permuted index.
    int n_edges_raw = static_cast<int>(qval_pairs.size());
    std::vector<int> order(n_edges_raw);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return qval_pairs[a].first < qval_pairs[b].first;
    });

    // Build deduplicated sorted vectors: keep min FDR and its associated effect
    std::vector<std::pair<int64_t, double>> qval_sorted;
    std::vector<std::pair<int64_t, double>> eff_sorted;
    qval_sorted.reserve(n_edges_raw);
    eff_sorted.reserve(n_edges_raw);

    for (int i = 0; i < n_edges_raw; ++i) {
        int idx = order[i];
        int64_t key = qval_pairs[idx].first;
        double qv = qval_pairs[idx].second;
        double eff = eff_pairs[idx].second;

        if (qval_sorted.empty() || qval_sorted.back().first != key) {
            qval_sorted.push_back({key, qv});
            eff_sorted.push_back({key, eff});
        } else {
            // Duplicate key: keep the entry with smaller q-value
            if (qv < qval_sorted.back().second) {
                qval_sorted.back().second = qv;
                eff_sorted.back().second = eff;
            }
        }
    }

    // --- 2. Count present species and check minimum ---
    std::vector<int> present_sp;
    for (int s = 0; s < n_target_species; s++) {
        if (!sp_genes[s].empty())
            present_sp.push_back(s);
    }
    int n_present = static_cast<int>(present_sp.size());
    if (n_present < min_species) return results;

    // --- 3. Cap genes per species (keep most-connected) ---
    for (int s = 0; s < n_target_species; s++) {
        if (static_cast<int>(sp_genes[s].size()) <= max_genes_per_sp) continue;

        auto& genes = sp_genes[s];

        // Count cross-species edges per gene; keep most-connected
        // Mark genes of this species in seen_genes for fast lookup
        for (int g : genes) {
            seen_genes[g] = true;
        }

        // Linear scan of qval_sorted to count edges per gene.
        // genes.size() <= max_genes_per_sp (small), so inner scan is bounded.
        std::vector<int> counts(genes.size(), 0);

        for (const auto& kv : qval_sorted) {
            int ga = static_cast<int>(kv.first >> 32);
            int gb = static_cast<int>(kv.first & 0xFFFFFFFF);
            if (seen_genes[ga]) {
                for (int i = 0; i < static_cast<int>(genes.size()); ++i) {
                    if (genes[i] == ga) { counts[i]++; break; }
                }
            }
            if (seen_genes[gb]) {
                for (int i = 0; i < static_cast<int>(genes.size()); ++i) {
                    if (genes[i] == gb) { counts[i]++; break; }
                }
            }
        }

        // Clear seen_genes
        for (int g : genes) {
            seen_genes[g] = false;
        }

        // Sort genes by descending edge count
        std::vector<int> idx(genes.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return counts[a] > counts[b];
        });

        std::vector<int> top_genes;
        top_genes.reserve(max_genes_per_sp);
        for (int i = 0; i < max_genes_per_sp; ++i) {
            top_genes.push_back(genes[idx[i]]);
        }
        genes = std::move(top_genes);
    }

    // --- 4. Species-level clique/subset enumeration ---
    std::vector<std::vector<int>> sp_cliques;

    // Helper: run BK on the present-species adjacency graph
    auto run_bk = [&]() {
        std::vector<std::vector<bool>> local_adj(
            n_present, std::vector<bool>(n_present, false));
        for (int i = 0; i < n_present; i++) {
            for (int j = i + 1; j < n_present; j++) {
                if (sp_adj[present_sp[i]][present_sp[j]]) {
                    local_adj[i][j] = local_adj[j][i] = true;
                }
            }
        }
        std::vector<int> R;
        std::vector<int> P(n_present);
        std::iota(P.begin(), P.end(), 0);
        std::vector<int> X;
        bron_kerbosch(R, P, X, local_adj, min_species, sp_cliques);
    };

    if (max_missing_edges == 0) {
        run_bk();
    } else {
        // Tolerant: enumerate all species subsets of size >= min_species
        // where number of missing species-pair edges <= max_missing_edges.
        // Practical limit: ~25 species (C(25,12) ~ 5M subsets).
        // For larger n_present, the total enumeration becomes infeasible.
        static constexpr int MAX_TOLERANT_SPECIES = 25;
        if (n_present > MAX_TOLERANT_SPECIES) {
            // REprintf (not Rcpp::warning) because this function is called
            // inside omp parallel for in find_cliques_stability.cpp.
            // REprintf is fd-write-based and safe from threads.
            REprintf("Note: max_missing_edges tolerance ignored for HOG "
                     "with %d species (limit %d); using exact BK\n",
                     n_present, MAX_TOLERANT_SPECIES);
            run_bk();
        } else {
            for (int sz = n_present; sz >= min_species; sz--) {
                if (sz == 0) continue;
                // Gosper's hack: enumerate all C(n_present, sz) subsets.
                // Safe shift: n_present <= MAX_TOLERANT_SPECIES < 64.
                uint64_t mask = (1ULL << sz) - 1;
                uint64_t limit = 1ULL << n_present;
                while (mask < limit) {
                    int missing = 0;
                    bool feasible = true;
                    std::vector<int> subset;
                    subset.reserve(sz);
                    for (int i = 0; i < n_present; i++) {
                        if ((mask >> i) & 1) subset.push_back(i);
                    }
                    for (int i = 0; i < sz && feasible; i++) {
                        for (int j = i + 1; j < sz; j++) {
                            if (!sp_adj[present_sp[subset[i]]][present_sp[subset[j]]]) {
                                missing++;
                                if (missing > max_missing_edges) {
                                    feasible = false;
                                    break;
                                }
                            }
                        }
                    }
                    if (feasible) {
                        sp_cliques.push_back(std::move(subset));
                    }
                    // Next subset (Gosper's hack)
                    uint64_t c = mask & (~mask + 1);
                    uint64_t r = mask + c;
                    mask = (((r ^ mask) >> 2) / c) | r;
                }
            }
        }
    }

    if (sp_cliques.empty()) return results;

    // --- 5. Gene assignment for each species clique/subset ---
    for (auto& sp_cl : sp_cliques) {
        int k = static_cast<int>(sp_cl.size());

        // Map local species indices back to global
        std::vector<int> global_sp(k);
        for (int i = 0; i < k; i++) global_sp[i] = present_sp[sp_cl[i]];

        // Sort species by gene count (fewest first = better pruning)
        std::sort(global_sp.begin(), global_sp.end(),
                  [&](int a, int b) {
                      return sp_genes[a].size() < sp_genes[b].size();
                  });

        // Backtrack to find best gene assignment
        std::vector<int> current(k, -1);
        AssignResult best;
        best.found = false;
        best.sum_q = 0.0;
        best.max_q = 0.0;
        best.sum_effect = 0.0;
        best.n_edges = 0;
        best.n_missing = 0;
        int iterations = 0;

        assign_bt(0, k, global_sp, sp_genes, qval_sorted, eff_sorted,
                  current, 0.0, 0.0, 0.0, 0, max_missing_edges,
                  best, iterations, max_iterations);

        if (best.found) {
            CliqueResult cr;
            cr.hog_idx = hog_idx;
            cr.genes.assign(n_target_species, -1);
            for (int i = 0; i < k; i++) {
                cr.genes[global_sp[i]] = best.genes[i];
            }
            cr.n_species = k;
            cr.n_edges = best.n_edges;
            cr.n_missing = best.n_missing;
            cr.mean_q = (best.n_edges > 0) ? best.sum_q / best.n_edges : 1.0;
            cr.max_q = best.max_q;
            cr.mean_effect = (best.n_edges > 0) ? best.sum_effect / best.n_edges : 0.0;
            results.push_back(std::move(cr));
        }
    }

    return results;
}

// ---------------------------------------------------------------------------
// Trait annotation: check if all active species in a clique share a trait value
// ---------------------------------------------------------------------------
// Returns the common trait value if all active species with genes share it,
// or -1 if the trait is not exclusive (mixed values among present species).
inline int annotate_trait(
    const std::vector<int>& genes,
    const int* species_trait,
    uint64_t active_mask,
    int n_target_species)
{
    int first_trait = -1;
    for (int s = 0; s < n_target_species; ++s) {
        if (!((active_mask >> s) & 1)) continue;
        if (genes[s] < 0) continue;
        int t = species_trait[s];
        if (first_trait < 0) first_trait = t;
        else if (t != first_trait) return -1;
    }
    return first_trait;
}

// ---------------------------------------------------------------------------
// Jaccard similarity on genes from non-removed species
// ---------------------------------------------------------------------------
// Compares two gene assignments, ignoring species whose bit is set in
// removal_mask.  Used by the stability pipeline to measure how much a
// clique changes when species are removed.
inline double jaccard_remaining(
    const std::vector<int>& full_genes,
    const std::vector<int>& reduced_genes,
    uint64_t removal_mask,
    int n_target_species)
{
    int intersect = 0, union_size = 0;
    for (int s = 0; s < n_target_species; ++s) {
        if ((removal_mask >> s) & 1) continue;   // skip removed species
        int fg = full_genes[s], rg = reduced_genes[s];
        bool f_has = (fg >= 0), r_has = (rg >= 0);
        if (f_has || r_has) {
            ++union_size;
            if (f_has && r_has && fg == rg) ++intersect;
        }
    }
    return (union_size > 0) ? static_cast<double>(intersect) / union_size : 0.0;
}

// ---------------------------------------------------------------------------
// Generate all C(n,k) bitmasks with exactly k bits set (Gosper's hack)
// ---------------------------------------------------------------------------
// For n <= 64, k <= n.  Returns a vector of uint64_t bitmasks where each
// has exactly k of the lowest n bits set.
inline std::vector<uint64_t> generate_subsets(int n, int k) {
    std::vector<uint64_t> result;
    if (k == 0) {
        result.push_back(0);
        return result;
    }
    if (k > n) return result;

    // Guard: 1ULL << 64 is undefined behavior
    if (k == n) {
        // Only one subset: all n bits set
        result.push_back((n == 64) ? ~0ULL : (1ULL << n) - 1);
        return result;
    }
    // Now k < n, so k <= 63 and (1ULL << k) is safe
    uint64_t mask = (1ULL << k) - 1;
    uint64_t limit = (n < 64) ? (1ULL << n) : 0ULL;  // 0 = unused for n==64

    while (true) {
        // For n < 64: standard termination on limit
        if (n < 64 && mask >= limit) break;
        result.push_back(mask);

        // Gosper's hack: next k-subset in lexicographic order
        uint64_t c = mask & (~mask + 1);      // lowest set bit
        uint64_t r = mask + c;                 // carry into next higher bit
        uint64_t diff = ((r ^ mask) >> 2) / c;
        uint64_t next = diff | r;

        // For n == 64: detect overflow wrap-around (next <= mask)
        if (next <= mask) break;
        mask = next;
    }

    return result;
}

#endif // RCOMPLEX_FIND_CLIQUES_COMMON_H
