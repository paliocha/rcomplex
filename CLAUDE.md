# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

rcomplex is an R package for comparative co-expression network analysis across species. It maps orthologous genes (via ortholog groups / HOGs from OrthoFinder, FastOMA, PLAZA, etc.), builds co-expression networks independently per species, then tests conservation at three levels:

- **Gene / HOG-level**: Hypergeometric tests with q-value correction (`compare_neighborhoods()` + `summarize_comparison()`), gene-identity permutation with adaptive stopping (`permutation_hog_test()`), batch orchestration (`find_coexpressologs()`, `density_sweep()`), degree-preserving edge-swap null (`coexpressolog_null()`)
- **Module-level**: Community detection (Leiden / Infomap / SBM) with multi-resolution consensus, cross-species comparison via hypergeometric or Jaccard permutation tests, module hubs + conservation
- **Clique-level**: C++ Bron-Kerbosch / Tomita clique detection, leave-k-out jackknife stability for trait-exclusive cliques, threshold sweep, HOG conservation classification

Based on [Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106).

## Build & Test

```bash
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'devtools::document()'
R CMD INSTALL .
Rscript -e 'devtools::test()'
R CMD build . && R CMD check --no-manual rcomplex_0.2.0.tar.gz   # expect "Status: OK"
```

Check the built tarball, not the source directory — `Authors@R` only expands at build time, so `R CMD check .` fails with "Author/Maintainer missing". `--no-manual` avoids needing pdflatex. The historical `R_ext/Boolean.h` warning no longer appears with clang 22.

## Package Architecture

### R layer
| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()`, `reduce_orthogroups()`, `prepare_orthologs()` |
| `R/network.R` | `compute_network()` — correlation, MR/CLR, density threshold, sparse extraction, torch GPU |
| `R/network-sparse.R` | Sparse dispatch: `.net_is_sparse()`, `.net_check()`, `.net_cpp_args()` (store guard), `as_sparse_network()`, `dense_to_dgc()` |
| `R/mr_block.R` | `mr_block()` — exact local MR reconstruction for gene subsets (incl. sub-store entries) |
| `R/comparison.R` | `compare_neighborhoods()`, `comparison_to_edges()`, `find_coexpressologs()` (alias: `run_pairwise_comparisons()`), `density_sweep()`, `get_coexpressed_hogs()` |
| `R/coexpressolog_null.R` | `coexpressolog_null()` — degree-preserving edge-swap null |
| `R/summary.R` | `summarize_comparison()`, `permutation_hog_test()`, `compute_qvalues()` (randomized-p pi0), torch FE helpers |
| `R/modules.R` | `detect_modules()` (single + consensus), `compare_modules()`, `compare_modules_paired()`, `classify_modules()`, `coarsen_modules()`, `identify_module_hubs()`, `classify_hub_conservation()`, `characterize_hubs()` |
| `R/tag_permutation.R` | `tag_permutation()` — trait-specific module recurrence test |
| `R/cliques.R` | `find_cliques()`, `clique_stability()`, `clique_persistence()`, `clique_threshold_sweep()`, `clique_perturbation_test()`, `clique_intensity_test()`, `classify_cliques()` |
| `R/se_methods.R` | `extract_orthologs()`, `build_se()` (internal) — SummarizedExperiment helpers |
| `R/rcomplex-class.R` | S3 `rcomplex` container: constructor, print/summary, `.rcomplex` methods for 9 pipeline functions |
| `R/rcomplex-package.R` | Package-level roxygen, namespace imports |

### C++ layer (src/, RcppArmadillo + OpenMP)
| File | Purpose |
|------|---------|
| `src/mutual_rank.cpp` | MR normalization with column-major access; in-place kernel (`mutual_rank_inplace_cpp`) + cached reference |
| `src/clr.cpp` | CLR normalization |
| `src/density_threshold.cpp` | Quantile-based density thresholding |
| `src/sparse_extract.cpp` | Sparse (dgCMatrix-slot) extraction of the thresholded MR matrix |
| `src/neighbor_lists.h` | Shared neighbour-list construction: dense `arma::mat` or validated dgCMatrix slots |
| `src/neighborhood_comparison.cpp` | Pairwise neighborhood overlap (hypergeometric, self-excluded urn); dense + sparse entry points |
| `src/hog_permutation.cpp` | HOG permutation engine (bit-vector/flag-vector, Besag & Clifford); dense + sparse entry points |
| `src/fe_permutation.cpp` | GPU-precomputed FE permutation engine |
| `src/module_jaccard_permutation.cpp` | Batched Jaccard permutation engine |
| `src/reduce_orthogroups.cpp` | Ward.D2 paralog merging |
| `src/coclassification.cpp` | Co-classification matrix with per-pair null subtraction (Jeub et al. 2018) |
| `src/find_cliques_common.h` | Shared clique primitives (BK/Tomita, backtracking, trait, Jaccard) |
| `src/find_cliques.cpp` | C++ clique detection wrapper |
| `src/find_cliques_stability.cpp` | Leave-k-out stability engine with OpenMP |
| `src/sample_k_distinct.h` | Shared rejection-sampling utility |

### Tests
| File | Covers |
|------|--------|
| `tests/testthat/test-network.R` | Network construction, MR/CLR, density threshold, in-place MR, torch |
| `tests/testthat/test-network-sparse.R` | Sparse network object: dense-vs-sparse equality across all consumers, store guard |
| `tests/testthat/test-equivalence.R` | Equivalence against canonical ComPlEx (seeded fixture, 149 calls) |
| `tests/testthat/test-comparison.R` | Neighborhood comparison, effect sizes, sparse dispatch, pval_combine |
| `tests/testthat/test-summary.R` | Pair-level q-value correction |
| `tests/testthat/test-pi0.R` | Randomized-p pi0 estimation (simulation vs truth, Storey, BH) |
| `tests/testthat/test-permutation.R` | HOG permutation (correctness, adaptive stopping, sparse, torch) |
| `tests/testthat/test-mr-block.R` | `mr_block()` exact reconstruction vs dense network |
| `tests/testthat/test-coexpressolog-null.R` | Edge-swap null (seeded, parallel reproducibility) |
| `tests/testthat/test-coexpressed-hogs.R` | `get_coexpressed_hogs()` cross-species partner queries |
| `tests/testthat/test-modules.R` | Module detection, comparison, classification, consensus |
| `tests/testthat/test-module-hubs.R` | Hub identification, tie-breaking, hub conservation |
| `tests/testthat/test-tag-permutation.R` | Trait-specific module recurrence test |
| `tests/testthat/test-cliques.R` | Clique detection (igraph + C++ backends) |
| `tests/testthat/test-stability.R` | Leave-k-out jackknife stability |
| `tests/testthat/test-threshold-sweep.R` | Threshold sweep structural survival |
| `tests/testthat/test-perturbation.R` | Clique perturbation (noise robustness) |
| `tests/testthat/test-intensity-test.R` | Clique intensity permutation null |
| `tests/testthat/test-classify-cliques.R` | HOG classification waterfall pipeline |
| `tests/testthat/test-reduce-orthogroups.R` | Paralog reduction |
| `tests/testthat/test-se.R` | SummarizedExperiment integration (build_se, extract_orthologs, S4 compute_network) |
| `tests/testthat/test-rcomplex-class.R` | S3 container construction, printing, method pass-through |
| `tests/testthat/helper-reference.R` | Pure-R reference implementations + shared net helpers (`sparse_net()`, ...) |
| `tests/testthat/helper-clique-fixtures.R` | Shared clique test fixtures |

## Key Design Decisions

### Sparse network object (dgCMatrix + store_density)
`compute_network(sparse = TRUE)` (default since 0.2.0) stores a `dgCMatrix` with both triangles of entries at or above the `store_density` quantile (default `max(density, 0.05)`); the analysis `threshold` is still computed from the full dense MR, so dense and sparse results are identical. `.net_cpp_args()` / `.net_check()` in `R/network-sparse.R` are the single choke points: every consumer validates there, and a threshold below `store_threshold` errors (the sparse object cannot represent that density). Sub-store values are reconstructible exactly via `mr_block()`.

### Membership-only consumers
Rewired null networks from `coexpressolog_null()` are binary (`threshold = 1`, `store_threshold = 1`): only membership-based consumers (neighborhood comparison, HOG permutation, adjacency extraction) are valid downstream — nothing that reads edge weights (density_sweep multipliers < 1, perturbation, module weights).

### Self-excluded urn
The anchor gene is never its own neighbour, so it leaves the ortholog-mapped set (k) and the hypergeometric population (N - 1) in `compare_neighborhoods()`, both permutation engines, and the torch FE matrix. O(1/N) p-value shift vs canonical ComPlEx; the reported `*.p.val.con` keeps the canonical `x > 1` gate.

### Randomized-p pi0 (pair level)
Exact hypergeometric p-values pile up at 1 and force Storey's pi0 to 1. `summarize_comparison(pi0_method = "randomized")` (default) estimates pi0 on draws of `p.val.gt + U * p.val.eq` (exactly uniform under H0) and applies it to the exact p-values. Draws use the global RNG — `set.seed()` for reproducible q-values. HOG-level q-values stay `DiscreteQvalue::DQ(method = "Liang")` (Besag–Clifford p-values have discrete support).

### pval_combine default "max"
`comparison_to_edges()` / `summarize_comparison()` / `find_coexpressologs()` / `density_sweep()` combine directional q-values with `pmax` by default: both directions must be significant (reciprocal criterion of Netotea et al. 2014, the `Max.p.val` filter). `"min"` is the permissive either-direction option (pre-0.2.0 behaviour).

### Column-major memory access
Armadillo stores matrices column-major. All hot loops use `colptr()` for sequential reads. ~3x speedup on MR normalization and neighbor list extraction. MR normalization runs in place on the correlation matrix (`mutual_rank_inplace_cpp()`), holding the `compute_network()` peak transient at ~1.5 n^2 doubles.

### Integer indices in C++, string mapping in R
Homebrew clang ABI issue with `std::unordered_map<std::string, ...>`. All C++ uses integer indices; R wrappers handle string-to-int mapping. The codebase also avoids `std::unordered_map` entirely — uses sorted vectors + binary search instead.

### HOG-level testing uses permutation, not Fisher's method
Fisher's method is anti-conservative for multi-copy HOGs (correlated tests). `permutation_hog_test()` permutes gene identities instead.

### Iterative consensus module detection
Multi-resolution Leiden sweep + iterative consensus per Jeub et al. (2018). Per-pair null subtraction: E(i,j) = (1/K) sum_k (s_m(i)/N)(s_m(j)/N), not a scalar mean. Iterates co-classification → Leiden sweep on consensus graph until all resolutions converge (ARI > 0.999).

### Build system
- `Makevars` / `Makevars.win`: C++23, `$(SHLIB_OPENMP_CXXFLAGS)` for portable OpenMP
- RcppArmadillo in `LinkingTo` only (NOT `Imports`)

## Dependencies

**Imports**: methods, Rcpp, Rfast, dplyr, igraph, Matrix (>= 1.5-0), DiscreteQvalue, qvalue, tidyr, parallel, rlang, stats, utils
**Suggests**: DT, knitr, purrr, rmarkdown, S4Vectors, sbm, stringr, SummarizedExperiment, tibble, torch, testthat, lintr, withr, pkgdown
**LinkingTo**: Rcpp, RcppArmadillo
**System**: GNU make, C++23, OpenMP (optional)
