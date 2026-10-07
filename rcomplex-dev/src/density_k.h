// density_k.h
// Number of strict-upper-triangle entries kept at a given density, shared
// by density_threshold_cpp() and the blockwise mutual-rank kernel.

#ifndef RCOMPLEX_DENSITY_K_H
#define RCOMPLEX_DENSITY_K_H

#include <cmath>
#include <cstddef>

inline std::size_t density_k(double density, std::size_t tri_size) {
    auto k = static_cast<std::size_t>(
        std::round(density * static_cast<double>(tri_size)));
    if (k == 0) k = 1;
    if (k >= tri_size) k = tri_size - 1;
    return k;
}

#endif  // RCOMPLEX_DENSITY_K_H
