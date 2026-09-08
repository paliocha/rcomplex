// module_preservation.cpp
// Adjacency-based module preservation statistics (Langfelder et al. 2011)
//
// Reference modules are projected onto a test species through an ortholog map,
// and each projected module is scored on six statistics computed from the test
// network alone plus the reference network's per-gene vectors:
//
//   density      meanAdj, meanClusterCoeff, meanMAR
//   connectivity cor.kIM, cor.clusterCoeff, cor.MAR
//
// The null shuffles gene identities while holding edges constant (the WGCNA
// convention, and the same one Mahler et al. 2017 use), handing each module a
// contiguous block of the shuffled mappable genes. Reference vectors stay
// fixed, so the connectivity correlations collapse toward zero under the null.
//
// Cost: the permutation loop never iterates per module. One shuffle labels
// every mappable gene at once, so a single pass over the induced edge list
// serves all modules -- O(E) per permutation for the whole partition. Triangles
// are then enumerated on the within-module lists, whose length is roughly
// deg / n_modules, rather than on the full neighbour lists.
//
// Design constraints (matching hog_permutation.cpp):
// - Integer indices only (no strings -- Homebrew clang ABI issue)
// - No std::unordered_map anywhere; sorted vectors + linear intersection
// - One std::mt19937 per permutation, all seeded from R::runif before the
//   parallel region, so results do not depend on n_cores

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "neighbor_lists.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

constexpr int kNStats = 6;

// Scratch buffers reused across permutations by one thread.
struct Scratch {
    std::vector<int> label;                   // module id per gene, -1 = none
    std::vector<std::vector<int>> w_idx;      // within-module neighbours
    std::vector<std::vector<double>> w_val;   // ... and their weights
    std::vector<double> k_im;                 // intramodular connectivity
    std::vector<double> sq;                   // sum of squared adjacencies
    std::vector<double> tri;                  // weighted triangle sum
    std::vector<int> perm;                    // shuffled gene order

    explicit Scratch(int m)
        : label(m, -1), w_idx(m), w_val(m),
          k_im(m, 0.0), sq(m, 0.0), tri(m, 0.0), perm(m) {}

    void reset() {
        const auto m = label.size();
        for (std::size_t u = 0; u < m; ++u) {
            label[u] = -1;
            w_idx[u].clear();
            w_val[u].clear();
            k_im[u] = 0.0;
            sq[u] = 0.0;
            tri[u] = 0.0;
        }
    }
};

// Build the within-module adjacency lists and accumulate kIM / sum-of-squares.
// One pass over the induced edge list, all modules at once.
void build_within_module(const WeightedNeighbors& g, Scratch& s) {
    const auto m = static_cast<int>(g.idx.size());
    for (int u = 0; u < m; ++u) {
        const int lu = s.label[u];
        if (lu < 0) continue;
        const std::vector<int>& nb = g.idx[u];
        const std::vector<double>& wt = g.w[u];
        const auto deg = static_cast<int>(nb.size());
        for (int e = 0; e < deg; ++e) {
            const int v = nb[e];
            if (s.label[v] != lu) continue;
            const double a = wt[e];
            s.w_idx[u].push_back(v);
            s.w_val[u].push_back(a);
            s.k_im[u] += a;
            s.sq[u] += a * a;
        }
    }
}

// Weighted triangle sums. For each within-module edge (u, v) with u < v, every
// common neighbour t contributes a_ut * a_vt * a_uv to tri[t] -- and only to
// tri[t], so each triangle is counted exactly once per vertex overall. Both
// lists are ascending (induced_weighted_lists_* guarantees it), so the
// intersection is linear.
void accumulate_triangles(Scratch& s) {
    const auto m = static_cast<int>(s.label.size());
    for (int u = 0; u < m; ++u) {
        if (s.label[u] < 0) continue;
        const std::vector<int>& nu = s.w_idx[u];
        const std::vector<double>& wu = s.w_val[u];
        const auto du = static_cast<int>(nu.size());

        for (int e = 0; e < du; ++e) {
            const int v = nu[e];
            if (v <= u) continue;
            const double a_uv = wu[e];
            const std::vector<int>& nv = s.w_idx[v];
            const std::vector<double>& wv = s.w_val[v];
            const auto dv = static_cast<int>(nv.size());

            int a = 0;
            int b = 0;
            while (a < du && b < dv) {
                if (nu[a] < nv[b]) {
                    ++a;
                } else if (nu[a] > nv[b]) {
                    ++b;
                } else {
                    s.tri[nu[a]] += wu[a] * wv[b] * a_uv;
                    ++a;
                    ++b;
                }
            }
        }
    }
}

inline double clustering_coeff(double k, double sq, double tri) {
    const double denom = k * k - sq;
    return (denom > 0.0) ? (2.0 * tri) / denom : 0.0;
}

inline double max_adjacency_ratio(double k, double sq) {
    return (k > 0.0) ? sq / k : 0.0;
}

// Pearson correlation; NA when it would be degenerate. Two points always
// correlate at +/-1, so fewer than three genes yields NA rather than a
// meaningless +/-1.
double safe_cor(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n = static_cast<int>(x.size());
    if (n < 3) return NA_REAL;

    double mx = 0.0;
    double my = 0.0;
    for (int i = 0; i < n; ++i) {
        mx += x[i];
        my += y[i];
    }
    mx /= n;
    my /= n;

    double sxy = 0.0;
    double sxx = 0.0;
    double syy = 0.0;
    for (int i = 0; i < n; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if (sxx <= 0.0 || syy <= 0.0) return NA_REAL;
    return sxy / std::sqrt(sxx * syy);
}

// Score every module under the labelling currently in `s`. `blocks[k]` lists
// the test genes of module k in the order that pairs them with the reference
// vectors. Writes kNStats values per module into `out`.
void score_modules(const WeightedNeighbors& g, Scratch& s,
                   const std::vector<std::vector<int>>& blocks,
                   const std::vector<std::vector<double>>& ref_k,
                   const std::vector<std::vector<double>>& ref_cc,
                   const std::vector<std::vector<double>>& ref_mar,
                   std::vector<double>& out) {
    build_within_module(g, s);
    accumulate_triangles(s);

    const auto n_mod = static_cast<int>(blocks.size());
    std::vector<double> tk;
    std::vector<double> tc;
    std::vector<double> tm;

    for (int k = 0; k < n_mod; ++k) {
        const std::vector<int>& members = blocks[k];
        const auto sz = static_cast<int>(members.size());

        tk.clear();
        tc.clear();
        tm.clear();
        tk.reserve(sz);
        tc.reserve(sz);
        tm.reserve(sz);

        double sum_k = 0.0;
        double sum_cc = 0.0;
        double sum_mar = 0.0;

        for (int j = 0; j < sz; ++j) {
            const int u = members[j];
            const double kv = s.k_im[u];
            const double cc = clustering_coeff(kv, s.sq[u], s.tri[u]);
            const double mar = max_adjacency_ratio(kv, s.sq[u]);
            tk.push_back(kv);
            tc.push_back(cc);
            tm.push_back(mar);
            sum_k += kv;
            sum_cc += cc;
            sum_mar += mar;
        }

        const int base = k * kNStats;
        out[base + 0] = (sz > 1)
            ? sum_k / (static_cast<double>(sz) * (sz - 1))
            : 0.0;
        out[base + 1] = (sz > 0) ? sum_cc / sz : 0.0;
        out[base + 2] = (sz > 0) ? sum_mar / sz : 0.0;
        out[base + 3] = safe_cor(ref_k[k], tk);
        out[base + 4] = safe_cor(ref_cc[k], tc);
        out[base + 5] = safe_cor(ref_mar[k], tm);
    }
}

std::vector<std::vector<int>> as_int_lists(const List& x) {
    const auto n = static_cast<int>(x.size());
    std::vector<std::vector<int>> out(n);
    for (int k = 0; k < n; ++k) {
        IntegerVector v = x[k];
        out[k].assign(v.begin(), v.end());
    }
    return out;
}

std::vector<std::vector<double>> as_num_lists(const List& x) {
    const auto n = static_cast<int>(x.size());
    std::vector<std::vector<double>> out(n);
    for (int k = 0; k < n; ++k) {
        NumericVector v = x[k];
        out[k].assign(v.begin(), v.end());
    }
    return out;
}

void fisher_yates(std::vector<int>& v, std::mt19937& rng) {
    for (auto i = static_cast<int>(v.size()) - 1; i > 0; --i) {
        std::uniform_int_distribution<int> dist(0, i);
        std::swap(v[i], v[dist(rng)]);
    }
}

// Shared driver: observed pass plus n_perm null passes.
// Replace every weight with 1. On a hard-thresholded network the surviving
// weights are near-constant (a max/min ratio of ~1.04 for MR at density 0.03),
// so binarising costs almost nothing and makes avg.weight exactly the module
// edge density -- the interpretation NetRep documents for unweighted networks.
void binarise(WeightedNeighbors& g) {
    for (auto& row : g.w) {
        for (double& v : row) v = 1.0;
    }
    g.scale = 1.0;
}

List run_preservation(WeightedNeighbors& g,
                      const List& module_members,
                      const List& ref_kIM,
                      const List& ref_CC,
                      const List& ref_MAR,
                      int n_perm, int n_cores, bool binary) {
    const auto m = static_cast<int>(g.idx.size());
    const auto n_mod = static_cast<int>(module_members.size());

    if (n_mod == 0) {
        stop("module_members must contain at least one module");
    }
    if (ref_kIM.size() != n_mod || ref_CC.size() != n_mod ||
        ref_MAR.size() != n_mod) {
        stop("ref_kIM, ref_CC and ref_MAR must each have one entry per "
             "module (got %d, %d, %d for %d modules)",
             static_cast<int>(ref_kIM.size()),
             static_cast<int>(ref_CC.size()),
             static_cast<int>(ref_MAR.size()), n_mod);
    }

    double scale = 1.0;
    if (binary) {
        binarise(g);
    } else {
        scale = induced_max_weight(g);
        rescale_weights(g, scale);
    }

    const std::vector<std::vector<int>> members = as_int_lists(module_members);
    const std::vector<std::vector<double>> ref_k = as_num_lists(ref_kIM);
    const std::vector<std::vector<double>> ref_c = as_num_lists(ref_CC);
    const std::vector<std::vector<double>> ref_m = as_num_lists(ref_MAR);

    int total = 0;
    for (int k = 0; k < n_mod; ++k) {
        const auto sz = static_cast<int>(members[k].size());
        if (static_cast<int>(ref_k[k].size()) != sz ||
            static_cast<int>(ref_c[k].size()) != sz ||
            static_cast<int>(ref_m[k].size()) != sz) {
            stop("module %d: reference vectors must match the module size "
                 "(%d genes)", k + 1, sz);
        }
        for (const int u : members[k]) {
            if (u < 0 || u >= m) {
                stop("module %d: gene index out of range", k + 1);
            }
        }
        total += sz;
    }
    if (total > m) {
        stop("modules hold %d genes but only %d are mappable", total, m);
    }

    // ---- Observed ----
    std::vector<double> observed(static_cast<std::size_t>(n_mod) * kNStats,
                                NA_REAL);
    {
        Scratch s(m);
        for (int k = 0; k < n_mod; ++k) {
            for (const int u : members[k]) s.label[u] = k;
        }
        score_modules(g, s, members, ref_k, ref_c, ref_m, observed);
    }

    // ---- Null ----
    const auto cells = static_cast<std::size_t>(n_mod) * kNStats;
    std::vector<double> sum(cells, 0.0);
    std::vector<double> sumsq(cells, 0.0);
    std::vector<int> count(cells, 0);
    std::vector<int> exceed(cells, 0);

    if (n_perm > 0) {
        int n_threads = (n_cores > 1) ? n_cores : 1;
#ifdef _OPENMP
        n_threads = std::min(n_threads, omp_get_max_threads());
#else
        n_threads = 1;
#endif
        // One seed per permutation, drawn from R's RNG before the parallel
        // region (no R API calls are legal inside it). Seeding per iteration
        // rather than per thread is what makes the result independent of
        // n_cores: iteration i draws the same permutation however the
        // iterations are distributed across threads.
        std::vector<uint32_t> seeds(n_perm);
        for (int t = 0; t < n_perm; ++t) {
            seeds[t] = static_cast<uint32_t>(R::runif(0.0, 4294967296.0));
        }

        std::vector<std::vector<double>> t_sum(
            n_threads, std::vector<double>(cells, 0.0));
        std::vector<std::vector<double>> t_sumsq(
            n_threads, std::vector<double>(cells, 0.0));
        std::vector<std::vector<int>> t_count(
            n_threads, std::vector<int>(cells, 0));
        std::vector<std::vector<int>> t_exceed(
            n_threads, std::vector<int>(cells, 0));

#ifdef _OPENMP
        #pragma omp parallel num_threads(n_threads) if(n_threads > 1)
#endif
        {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            Scratch s(m);
            std::vector<double> stats(cells, NA_REAL);
            std::vector<std::vector<int>> blocks(n_mod);

#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (int iter = 0; iter < n_perm; ++iter) {
                std::mt19937 rng(seeds[iter]);
                s.reset();
                // Restart from the identity so the shuffle is a pure function
                // of this iteration's seed, not of whatever this thread did
                // on its previous iteration.
                for (int i2 = 0; i2 < m; ++i2) s.perm[i2] = i2;
                fisher_yates(s.perm, rng);

                int offset = 0;
                for (int k = 0; k < n_mod; ++k) {
                    const auto sz = static_cast<int>(members[k].size());
                    blocks[k].assign(s.perm.begin() + offset,
                                     s.perm.begin() + offset + sz);
                    for (const int u : blocks[k]) s.label[u] = k;
                    offset += sz;
                }

                score_modules(g, s, blocks, ref_k, ref_c, ref_m, stats);

                for (std::size_t c = 0; c < cells; ++c) {
                    const double v = stats[c];
                    if (!ISNAN(v)) {
                        t_sum[tid][c] += v;
                        t_sumsq[tid][c] += v * v;
                        t_count[tid][c] += 1;
                        // One-sided "greater": preservation means the observed
                        // statistic exceeds what random gene sets achieve.
                        if (!ISNAN(observed[c]) && v >= observed[c]) {
                            t_exceed[tid][c] += 1;
                        }
                    }
                }
            }
        }

        for (int t = 0; t < n_threads; ++t) {
            for (std::size_t c = 0; c < cells; ++c) {
                sum[c] += t_sum[t][c];
                sumsq[c] += t_sumsq[t][c];
                count[c] += t_count[t][c];
                exceed[c] += t_exceed[t][c];
            }
        }
    }

    NumericMatrix obs_out(n_mod, kNStats);
    NumericMatrix mean_out(n_mod, kNStats);
    NumericMatrix sd_out(n_mod, kNStats);
    IntegerMatrix n_out(n_mod, kNStats);
    IntegerMatrix exceed_out(n_mod, kNStats);
    NumericMatrix p_out(n_mod, kNStats);

    for (int k = 0; k < n_mod; ++k) {
        for (int j = 0; j < kNStats; ++j) {
            const std::size_t c = static_cast<std::size_t>(k) * kNStats + j;
            obs_out(k, j) = observed[c];
            // An NA observed statistic can never be exceeded, so reporting the
            // counts would let any consumer recompute (0 + 1) / (n + 1) -- the
            // most significant value attainable -- for a statistic that could
            // not be computed at all. All three outputs agree on NA.
            // Gate the counts on the observed value alone. A finite observed
            // statistic whose permutations were all NA is a different state --
            // the null was uncomputable, not the statistic -- and a usable
            // count of 0 is the only thing that records it.
            const bool have_obs = !ISNAN(observed[c]);
            n_out(k, j) = have_obs ? count[c] : NA_INTEGER;
            exceed_out(k, j) = have_obs ? exceed[c] : NA_INTEGER;
            p_out(k, j) = (have_obs && count[c] > 0)
                ? static_cast<double>(exceed[c] + 1) / (count[c] + 1)
                : NA_REAL;
            if (count[c] > 1) {
                const double mu = sum[c] / count[c];
                // Unbiased variance, clamped: catastrophic cancellation on a
                // constant statistic can leave a tiny negative value.
                double var = (sumsq[c] - count[c] * mu * mu) / (count[c] - 1);
                if (var < 0.0) var = 0.0;
                mean_out(k, j) = mu;
                sd_out(k, j) = std::sqrt(var);
            } else if (count[c] == 1) {
                mean_out(k, j) = sum[c];
                sd_out(k, j) = NA_REAL;
            } else {
                mean_out(k, j) = NA_REAL;
                sd_out(k, j) = NA_REAL;
            }
        }
    }

    return List::create(
        Named("observed") = obs_out,
        Named("perm_mean") = mean_out,
        Named("perm_sd") = sd_out,
        Named("n_perm_used") = n_out,
        Named("n_exceed") = exceed_out,
        Named("p_value") = p_out,
        Named("scale") = scale
    );
}

// Per-gene statistics for one labelling, used for the reference network's
// fixed vectors (and reusable for intramodular connectivity elsewhere).
List run_gene_stats(WeightedNeighbors& g, const List& module_members,
                    bool binary) {
    const auto m = static_cast<int>(g.idx.size());
    const auto n_mod = static_cast<int>(module_members.size());

    double scale = 1.0;
    if (binary) {
        binarise(g);
    } else {
        scale = induced_max_weight(g);
        rescale_weights(g, scale);
    }

    const std::vector<std::vector<int>> members = as_int_lists(module_members);

    Scratch s(m);
    for (int k = 0; k < n_mod; ++k) {
        for (const int u : members[k]) {
            if (u < 0 || u >= m) {
                stop("module %d: gene index out of range", k + 1);
            }
            s.label[u] = k;
        }
    }

    build_within_module(g, s);
    accumulate_triangles(s);

    NumericVector k_im(m);
    NumericVector cc(m);
    NumericVector mar(m);
    for (int u = 0; u < m; ++u) {
        k_im[u] = s.k_im[u];
        cc[u] = clustering_coeff(s.k_im[u], s.sq[u], s.tri[u]);
        mar[u] = max_adjacency_ratio(s.k_im[u], s.sq[u]);
    }

    return List::create(
        Named("kIM") = k_im,
        Named("CC") = cc,
        Named("MAR") = mar,
        Named("scale") = scale
    );
}

}  // namespace


//' Module preservation permutation engine (dense)
//'
//' @param net Dense network matrix.
//' @param thr Analysis threshold.
//' @param keep Ascending 0-based indices of the ortholog-mappable genes.
//' @param module_members List of integer vectors: local gene indices per
//'   module, ordered to match the reference vectors.
//' @param ref_kIM,ref_CC,ref_MAR Lists of numeric vectors: the reference
//'   network's per-gene statistics, one vector per module.
//' @param n_perm Number of permutations.
//' @param n_cores Number of OpenMP threads.
//' @param binary Treat every surviving edge as weight 1.
//' @return List with observed, perm_mean, perm_sd, n_perm_used, n_exceed,
//'   p_value, scale.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List module_preservation_dense_cpp(
    const arma::mat& net, double thr,
    const Rcpp::IntegerVector& keep,
    const Rcpp::List& module_members,
    const Rcpp::List& ref_kIM,
    const Rcpp::List& ref_CC,
    const Rcpp::List& ref_MAR,
    int n_perm, int n_cores, bool binary
) {
    WeightedNeighbors g = induced_weighted_lists_dense(net, thr, keep);
    return run_preservation(g, module_members, ref_kIM, ref_CC, ref_MAR,
                            n_perm, n_cores, binary);
}


//' Module preservation permutation engine (sparse)
//'
//' @param p,i,x dgCMatrix slots.
//' @param thr Analysis threshold.
//' @param keep Ascending 0-based indices of the ortholog-mappable genes.
//' @param module_members List of integer vectors: local gene indices per
//'   module, ordered to match the reference vectors.
//' @param ref_kIM,ref_CC,ref_MAR Lists of numeric vectors: the reference
//'   network's per-gene statistics, one vector per module.
//' @param n_perm Number of permutations.
//' @param n_cores Number of OpenMP threads.
//' @param binary Treat every surviving edge as weight 1.
//' @return List with observed, perm_mean, perm_sd, n_perm_used, n_exceed,
//'   p_value, scale.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List module_preservation_sparse_cpp(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x,
    double thr,
    const Rcpp::IntegerVector& keep,
    const Rcpp::List& module_members,
    const Rcpp::List& ref_kIM,
    const Rcpp::List& ref_CC,
    const Rcpp::List& ref_MAR,
    int n_perm, int n_cores, bool binary
) {
    WeightedNeighbors g = induced_weighted_lists_sparse(p, i, x, thr, keep);
    return run_preservation(g, module_members, ref_kIM, ref_CC, ref_MAR,
                            n_perm, n_cores, binary);
}


//' Per-gene intramodular statistics (dense)
//'
//' @param net Dense network matrix.
//' @param thr Analysis threshold.
//' @param keep Ascending 0-based indices of the genes to induce on.
//' @param module_members List of integer vectors: local gene indices per
//'   module.
//' @param binary Treat every surviving edge as weight 1.
//' @return List with kIM, CC, MAR (one value per induced gene) and scale.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List module_gene_stats_dense_cpp(
    const arma::mat& net, double thr,
    const Rcpp::IntegerVector& keep,
    const Rcpp::List& module_members, bool binary
) {
    WeightedNeighbors g = induced_weighted_lists_dense(net, thr, keep);
    return run_gene_stats(g, module_members, binary);
}


//' Per-gene intramodular statistics (sparse)
//'
//' @param p,i,x dgCMatrix slots.
//' @param thr Analysis threshold.
//' @param keep Ascending 0-based indices of the genes to induce on.
//' @param module_members List of integer vectors: local gene indices per
//'   module.
//' @param binary Treat every surviving edge as weight 1.
//' @return List with kIM, CC, MAR (one value per induced gene) and scale.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List module_gene_stats_sparse_cpp(
    const Rcpp::IntegerVector& p,
    const Rcpp::IntegerVector& i,
    const Rcpp::NumericVector& x,
    double thr,
    const Rcpp::IntegerVector& keep,
    const Rcpp::List& module_members, bool binary
) {
    WeightedNeighbors g = induced_weighted_lists_sparse(p, i, x, thr, keep);
    return run_gene_stats(g, module_members, binary);
}
