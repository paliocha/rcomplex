// sample_k_distinct.h
// Sample k distinct integers from [0, n) using rejection sampling.
// Efficient when k << n (typical: HOGs have 2-10 genes, networks 5K-30K).

#ifndef RCOMPLEX_SAMPLE_K_DISTINCT_H
#define RCOMPLEX_SAMPLE_K_DISTINCT_H

#include <random>
#include <vector>

inline void sample_k_distinct(int n, int k, std::mt19937& rng,
                               std::vector<int>& out) {
    out.resize(k);
    std::uniform_int_distribution<int> dist(0, n - 1);
    for (int i = 0; i < k; ++i) {
        int x;
        bool dup;
        do {
            x = dist(rng);
            dup = false;
            for (int j = 0; j < i; ++j) {
                if (out[j] == x) { dup = true; break; }
            }
        } while (dup);
        out[i] = x;
    }
}

#endif // RCOMPLEX_SAMPLE_K_DISTINCT_H
