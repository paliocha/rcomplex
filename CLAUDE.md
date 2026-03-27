# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

rcomplex is an R package for comparative co-expression network analysis across species. It maps orthologous genes (via ortholog groups / HOGs from OrthoFinder, FastOMA, PLAZA, etc.), builds co-expression networks independently per species, then tests conservation at three levels:

- **Gene / HOG-level**: Hypergeometric tests with q-value correction (`compare_neighborhoods()` + `summarize_comparison()`) and gene-identity permutation with adaptive stopping (`permutation_hog_test()`)
- **Module-level**: Community detection (Leiden / Infomap / SBM) with multi-resolution consensus, cross-species comparison via hypergeometric or Jaccard permutation tests
- **Clique-level**: C++ Bron-Kerbosch clique detection, leave-k-out jackknife stability for trait-exclusive cliques, hub gene identification

Based on [Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106).

## Build & Test

```bash
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'devtools::document()'
R CMD INSTALL .
Rscript -e 'devtools::test()'
R CMD check .                     # expect 1 WARNING from R_ext/Boolean.h
```

## Package Architecture

### R layer
| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()`, `reduce_orthogroups()` |
| `R/network.R` | `compute_network()` — correlation, MR/CLR, density threshold, torch GPU |
| `R/comparison.R` | `compare_neighborhoods()` — pair-level hypergeometric |
| `R/summary.R` | `summarize_comparison()`, `permutation_hog_test()`, torch FE helpers |
| `R/modules.R` | `detect_modules()` (single + consensus), `compare_modules()`, `classify_modules()` |
| `R/cliques.R` | `find_cliques()`, `clique_stability()`, `clique_hubs()` |
| `R/rcomplex-package.R` | Package-level roxygen, namespace imports |

### C++ layer (src/, RcppArmadillo + OpenMP)
| File | Purpose |
|------|---------|
| `src/mutual_rank.cpp` | MR normalization with column-major access |
| `src/clr.cpp` | CLR normalization |
| `src/density_threshold.cpp` | Quantile-based density thresholding |
| `src/neighborhood_comparison.cpp` | Pairwise neighborhood overlap (hypergeometric) |
| `src/hog_permutation.cpp` | HOG permutation engine (bit-vector/flag-vector, Besag & Clifford) |
| `src/fe_permutation.cpp` | GPU-precomputed FE permutation engine |
| `src/module_jaccard_permutation.cpp` | Batched Jaccard permutation engine |
| `src/reduce_orthogroups.cpp` | Ward.D2 paralog merging |
| `src/coclassification.cpp` | Co-classification matrix for consensus modules |
| `src/find_cliques_common.h` | Shared clique primitives (BK, backtracking, trait, Jaccard) |
| `src/find_cliques.cpp` | C++ clique detection wrapper |
| `src/find_cliques_stability.cpp` | Leave-k-out stability engine with OpenMP |
| `src/sample_k_distinct.h` | Shared rejection-sampling utility |

### Tests
| File | Covers |
|------|--------|
| `tests/testthat/test-network.R` | Network construction, MR/CLR, density threshold, torch |
| `tests/testthat/test-comparison.R` | Neighborhood comparison, effect sizes |
| `tests/testthat/test-summary.R` | Pair-level q-value correction |
| `tests/testthat/test-permutation.R` | HOG permutation (correctness, adaptive stopping, torch) |
| `tests/testthat/test-modules.R` | Module detection, comparison, classification, consensus |
| `tests/testthat/test-cliques.R` | Clique detection (igraph + C++ backends) |
| `tests/testthat/test-stability.R` | Leave-k-out jackknife stability |
| `tests/testthat/test-clique-hubs.R` | Hub gene identification |
| `tests/testthat/test-reduce-orthogroups.R` | Paralog reduction |
| `tests/testthat/helper-reference.R` | Pure-R reference implementations |

## Key Design Decisions

### Column-major memory access
Armadillo stores matrices column-major. All hot loops use `colptr()` for sequential reads. ~3x speedup on MR normalization and neighbor list extraction.

### Integer indices in C++, string mapping in R
Homebrew clang ABI issue with `std::unordered_map<std::string, ...>`. All C++ uses integer indices; R wrappers handle string-to-int mapping. The codebase also avoids `std::unordered_map` entirely — uses sorted vectors + binary search instead.

### HOG-level testing uses permutation, not Fisher's method
Fisher's method is anti-conservative for multi-copy HOGs (correlated tests). `permutation_hog_test()` permutes gene identities instead.

### Adaptive consensus module detection
Multi-resolution Leiden sweep + Jeub-inspired adaptive thresholding. Subtracts expected co-classification under random assignment (E_k = sum_m(s_m/N)²), discounting mega-modules at low resolutions.

### Build system
- `Makevars` / `Makevars.win`: C++23, `$(SHLIB_OPENMP_CXXFLAGS)` for portable OpenMP
- RcppArmadillo in `LinkingTo` only (NOT `Imports`)

## Dependencies

**Imports**: Rcpp, Rfast, dplyr, igraph, DiscreteQvalue, qvalue, tidyr, parallel, rlang, stats, utils
**Suggests**: sbm, torch, testthat, lintr, pkgdown
**LinkingTo**: Rcpp, RcppArmadillo
**System**: GNU make, C++23, OpenMP (optional)
