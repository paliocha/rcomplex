// rewire_degseq.cpp
// Degree-preserving edge swaps on a simple undirected graph: the null model
// of coexpressolog_null().
//
// Each trial is igraph_rewire()'s with loops = FALSE: draw two distinct
// edges uniformly, flip the second one's endpoints with probability 1/2,
// reject a swap that would be a no-op or create a loop or a multi-edge, and
// count the rejected trial anyway. Counting rejections is what makes the
// chain's stationary distribution uniform over the degree sequence's
// realizations; retrying until a swap succeeds would bias it.
//
// Adjacency is an n x n bit matrix (n^2 / 8 bytes: 32 MB at n = 16 000), so
// the multi-edge check is one bit test where igraph edits a graph. Draws
// come from R's RNG, so set.seed() in the caller fixes the rewiring.

#include <Rcpp.h>
#include <cstdint>
#include <vector>

// [[Rcpp::export]]
Rcpp::List rewire_degseq_cpp(Rcpp::IntegerVector from,
                             Rcpp::IntegerVector to,
                             int n, double niter) {
    const std::size_t m = from.size();
    std::vector<int> u(from.begin(), from.end());
    std::vector<int> v(to.begin(), to.end());

    if (m >= 2 && niter >= 1) {
        const std::uint64_t nn = static_cast<std::uint64_t>(n);
        std::vector<std::uint64_t> bits((nn * nn + 63) / 64, 0);
        auto bit = [nn](int x, int y) {
            return static_cast<std::uint64_t>(x) * nn +
                   static_cast<std::uint64_t>(y);
        };
        auto has = [&](int x, int y) {
            const std::uint64_t k = bit(x, y);
            return (bits[k >> 6] >> (k & 63)) & 1ULL;
        };
        auto set = [&](int x, int y, bool on) {
            for (std::uint64_t k : {bit(x, y), bit(y, x)}) {
                if (on) bits[k >> 6] |= 1ULL << (k & 63);
                else    bits[k >> 6] &= ~(1ULL << (k & 63));
            }
        };
        for (std::size_t e = 0; e < m; ++e) set(u[e], v[e], true);

        const double dm = static_cast<double>(m);
        const std::uint64_t n_trials = static_cast<std::uint64_t>(niter);
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
