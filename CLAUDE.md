# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

rcomplex is an R package for comparative co-expression network analysis across species. It maps orthologous genes (via PLAZA ortholog groups / HOGs), builds co-expression networks independently per species, then tests whether network neighborhoods are significantly preserved.

Two levels of testing:
- **Pair-level**: Hypergeometric tests with FDR correction and effect sizes (`compare_neighborhoods()` + `summarize_comparison()`)
- **HOG-level**: Gene-identity permutation test with adaptive stopping (`permutation_hog_test()`, recommended for multi-copy gene families)

Based on [Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106).

## Build & Test

```bash
# Regenerate Rcpp exports after modifying C++ [[Rcpp::export]] signatures
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate NAMESPACE and man/ from roxygen
Rscript -e 'devtools::document()'

# Install the package locally
R CMD INSTALL .

# Run tests (231 tests)
Rscript -e 'devtools::test()'

# Full R CMD check (expect 1 WARNING from R_ext/Boolean.h — unfixable R header issue with Homebrew clang)
R CMD check .
```

## Package Architecture

### R layer
| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()` — parse PLAZA ortholog group files |
| `R/network.R` | `compute_network()` — correlation + MR/CLR normalization + density threshold |
| `R/comparison.R` | `compare_neighborhoods()` — pair-level hypergeometric tests |
| `R/summary.R` | `summarize_comparison()` (pair-level FDR), `permutation_hog_test()` (HOG-level permutation) |
| `R/cliques.R` | `find_conserved_cliques()` — conserved clique detection via igraph |
| `R/rcomplex-package.R` | Package-level roxygen, namespace imports |

### C++ layer (src/, Armadillo + OpenMP)
| File | Purpose |
|------|---------|
| `src/mutual_rank.cpp` | MR normalization with cache-friendly column access |
| `src/clr.cpp` | CLR normalization |
| `src/density_threshold.cpp` | Quantile-based density thresholding |
| `src/neighborhood_comparison.cpp` | Pairwise neighborhood overlap (hypergeometric) |
| `src/hog_permutation.cpp` | Permutation engine: bit-vector/flag-vector intersections, Besag & Clifford stopping |

### Tests
| File | Covers |
|------|--------|
| `tests/testthat/test-network.R` | Network construction, MR/CLR, density threshold |
| `tests/testthat/test-comparison.R` | Neighborhood comparison, effect sizes |
| `tests/testthat/test-summary.R` | Pair-level FDR correction |
| `tests/testthat/test-permutation.R` | Permutation HOG test (correctness, adaptive stopping, p-value formula) |
| `tests/testthat/test-cliques.R` | Conserved clique detection |
| `tests/testthat/helper-reference.R` | Pure-R reference implementations for cross-checking C++ |

## Key Design Decisions

### Column-major memory access
Armadillo stores matrices column-major. Correlation and network matrices are symmetric, so `sim.col(i) == sim.row(i)` in value. All hot loops use column access (`sim.col(i)`, `net.colptr(i)`) for sequential reads instead of strided row access. This gives ~3x speedup on MR normalization and neighbor list extraction.

### Integer indices in C++, string mapping in R
Homebrew clang 21 has ABI issues with `std::unordered_map<std::string, ...>` (causes `__hash_memory` symbol not found at dlopen). All C++ functions use integer indices; the R wrappers handle string-to-int mapping via `match()`.

### HOG-level testing uses permutation, not Fisher's method
Fisher's method for combining pair-level p-values within HOGs is anti-conservative because pair-level tests are correlated (shared neighborhoods and ortholog mappings). `permutation_hog_test()` avoids this by permuting gene identities — the correlation structure is present in both observed and null distributions, so no independence assumption is needed.

### Build system
- `Makevars`: Uses `$(SHLIB_OPENMP_CXXFLAGS)` for portable OpenMP (not `-fopenmp`)
- RcppArmadillo in `LinkingTo` only (NOT `Imports` — causes NOTE)
- NAMESPACE file must exist before `Rcpp::compileAttributes()` works

## Dependencies

**Imports**: Rcpp, Rfast, dplyr, igraph, tidyr, rlang, stats, utils
**LinkingTo**: Rcpp, RcppArmadillo
**System**: GNU make, C++17, OpenMP (optional)
