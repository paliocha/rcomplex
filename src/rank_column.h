// rank_column.h
// In-place average ranking of one column, shared by the dense and the
// blockwise mutual-rank kernels.

#ifndef RCOMPLEX_RANK_COLUMN_H
#define RCOMPLEX_RANK_COLUMN_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <ranges>
#include <vector>

// In-place average ranks of one column; same tie handling as
// compute_ranks_impl. `indices` is a caller-owned buffer of size n (reused
// across columns). Overwriting col[] while walking the sorted order is safe:
// every position is written exactly once, only after every comparison that
// reads it has been made (tie groups are contiguous in `indices`).
inline void rank_column_inplace(double* col, const R_xlen_t n,
                                const bool ascending,
                                std::vector<R_xlen_t>& indices) {
    std::iota(indices.begin(), indices.end(), R_xlen_t{0});

    auto proj = [col](R_xlen_t i) { return col[i]; };
    if (ascending) {
        std::ranges::sort(indices, std::ranges::less{}, proj);
    } else {
        std::ranges::sort(indices, std::ranges::greater{}, proj);
    }

    R_xlen_t i = 0;
    while (i < n) {
        R_xlen_t j = i;
        while (j < n - 1 && col[indices[j]] == col[indices[j + 1]]) {
            ++j;
        }
        const double avg_rank = static_cast<double>(i + j + 2) / 2.0;
        for (R_xlen_t k = i; k <= j; ++k) {
            col[indices[k]] = avg_rank;
        }
        i = j + 1;
    }
}

#endif  // RCOMPLEX_RANK_COLUMN_H
