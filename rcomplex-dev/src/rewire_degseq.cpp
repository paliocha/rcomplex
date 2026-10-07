// rewire_degseq.cpp
// Degree-preserving edge swaps on a simple undirected graph: the null model
// of coexpressolog_null().
//
// Each trial is igraph_rewire()'s with loops = FALSE: two distinct edges are
// drawn uniformly, the second one's endpoints are flipped with probability
// 1/2, a swap that would be a no-op or create a loop or a multi-edge is
// rejected, and the rejected trial still counts. Counting rejections is what
// makes the chain's stationary distribution uniform over the degree
// sequence's realizations; retrying until a swap succeeds would bias it.
//
// Adjacency is an n x n bit matrix (n^2 / 8 bytes: 32 MB at n = 16 000), so
// the multi-edge check is one bit test where igraph edits a graph. Draws
// come from R's RNG, so set.seed() in the caller fixes the rewiring.

#include "neighbor_lists.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

// Takes the slots of a symmetric binary dgCMatrix (both triangles stored)
// and returns each rewired undirected edge once, 0-based.
// [[Rcpp::export]]
Rcpp::List rewire_degseq_cpp(const Rcpp::IntegerVector& p,
                             const Rcpp::IntegerVector& i,
                             const Rcpp::NumericVector& x,
                             double swap_factor) {
    // Every endpoint below comes from slots this has checked, so the bit
    // lookups cannot leave the n x n matrix.
    const int n = validate_dgc_slots(p, i, x);

    std::vector<int> u, v;
    for (int c = 0; c < n; ++c) {
        for (int k = p[c]; k < p[c + 1]; ++k) {
            if (i[k] < c) {
                u.push_back(i[k]);
                v.push_back(c);
            }
        }
    }
    const std::size_t m = u.size();

    // Rounded up, so any swap_factor > 0 makes at least one trial: truncating
    // a small factor on a small graph to zero trials would hand back the
    // observed graph as its own null. 2^64 is the first count a uint64_t
    // cannot hold.
    const double trials = std::ceil(swap_factor * static_cast<double>(m));
    if (!std::isfinite(swap_factor) || swap_factor < 0 ||
        trials >= 18446744073709551616.0)
        Rcpp::stop("swap_factor must be finite and non-negative, with "
                   "swap_factor * edge count below 2^64");

    if (m >= 2 && trials >= 1) {
        const std::uint64_t nn = static_cast<std::uint64_t>(n);
        std::vector<std::uint64_t> bits((nn * nn + 63) / 64, 0);
        auto bit = [nn](int a, int b) {
            return static_cast<std::uint64_t>(a) * nn +
                   static_cast<std::uint64_t>(b);
        };
        auto has = [&](int a, int b) {
            const std::uint64_t k = bit(a, b);
            return (bits[k >> 6] >> (k & 63)) & 1ULL;
        };
        auto set = [&](int a, int b, bool on) {
            for (std::uint64_t k : {bit(a, b), bit(b, a)}) {
                if (on) bits[k >> 6] |= 1ULL << (k & 63);
                else    bits[k >> 6] &= ~(1ULL << (k & 63));
            }
        };
        for (std::size_t e = 0; e < m; ++e) set(u[e], v[e], true);

        const double dm = static_cast<double>(m);
        const std::uint64_t n_trials = static_cast<std::uint64_t>(trials);
        for (std::uint64_t t = 0; t < n_trials; ++t) {
            if ((t & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();

            const std::size_t e1 = static_cast<std::size_t>(unif_rand() * dm);
            std::size_t e2;
            do {
                e2 = static_cast<std::size_t>(unif_rand() * dm);
            } while (e2 == e1);

            const int a = u[e1], b = v[e1];
            int c = u[e2], d = v[e2];
            if (unif_rand() < 0.5) std::swap(c, d);

            // the input is simple, so a != b and c != d
            if (a == c || b == d) continue;         // no-op swap
            if (a == d || b == c) continue;         // would create a loop
            if (has(a, d) || has(c, b)) continue;   // would create a multi-edge

            set(a, b, false);
            set(c, d, false);
            set(a, d, true);
            set(c, b, true);
            v[e1] = d;
            u[e2] = c;
            v[e2] = b;
        }
    }

    return Rcpp::List::create(Rcpp::Named("from") = u,
                              Rcpp::Named("to") = v);
}
