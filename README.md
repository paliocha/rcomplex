# rcomplex

<!-- badges: start -->
[![R-CMD-check](https://github.com/paliocha/rcomplex/actions/workflows/r.yml/badge.svg)](https://github.com/paliocha/rcomplex/actions/workflows/r.yml)
<!-- badges: end -->

Comparative co-expression network analysis across species.

**rcomplex** maps orthologous genes between two species, builds
co-expression networks independently, then tests conservation at two
complementary levels:

- **Gene / HOG level** -- Are individual genes' co-expression
  neighborhoods preserved across species?
- **Module level** -- Are higher-order communities (modules) of
  co-expressed genes conserved, partially conserved, or species-specific?

Based on [Netotea *et al.*
(2014)](https://doi.org/10.1186/1471-2164-15-106), extended with
permutation-based HOG-level testing, module detection via Leiden /
Infomap / Stochastic Block Model, and cross-species module comparison
with Jaccard permutation or hypergeometric tests.

## Installation

```r
devtools::install_github("paliocha/rcomplex")
```

**System requirements:** C++17 compiler, GNU make. OpenMP is optional
but recommended for parallel permutation tests.

## Quick start

```r
library(rcomplex)

# 1. Parse ortholog groups (PLAZA format)
orthologs <- parse_orthologs("orthogroups.txt", "species1", "species2")

# 2. Build co-expression networks
net1 <- compute_network(expr1, norm_method = "MR", density = 0.03)
net2 <- compute_network(expr2, norm_method = "MR", density = 0.03)

# 3. Compare neighborhoods (pair-level hypergeometric tests)
comparison <- compare_neighborhoods(net1, net2, orthologs)
```

From here, two analysis paths are available.

### Gene / HOG-level analysis

```r
# Pair-level FDR with Storey q-values
summary <- summarize_comparison(comparison)

# HOG-level permutation test (recommended for multi-copy gene families)
hog_results <- permutation_hog_test(net1, net2, comparison, n_cores = 4L)
```

### Module-level analysis

```r
# Detect modules in each network
mod1 <- detect_modules(net1, method = "leiden")
mod2 <- detect_modules(net2, method = "leiden")

# Compare modules across species (hypergeometric or Jaccard permutation)
comp <- compare_modules(mod1, mod2, orthologs, method = "jaccard", n_cores = 4L)

# Classify: conserved, partially conserved, or species-specific
classes <- classify_modules(comp)
```

## Main functions

| Function | Purpose |
|----------|---------|
| `parse_orthologs()` | Parse PLAZA ortholog group files |
| `compute_network()` | Correlation + MR/CLR normalization + density threshold |
| `compare_neighborhoods()` | Pair-level hypergeometric neighborhood tests |
| `summarize_comparison()` | Storey q-values and summary statistics |
| `permutation_hog_test()` | Permutation-based HOG-level conservation test |
| `detect_modules()` | Community detection (Leiden / Infomap / SBM) |
| `compare_modules()` | Cross-species module overlap (hypergeometric or Jaccard permutation) |
| `classify_modules()` | Three-tier module conservation classification |
| `find_coexpression_cliques()` | Conserved clique detection across species |

## Statistical methods

### Neighborhood comparison

For each ortholog pair, the co-expression neighborhood in species 1 is
mapped to species 2 through the ortholog table and tested for significant
overlap with the species-2 neighborhood using `phyper()`. The test is
performed in both directions. Effect sizes are fold-enrichments over the
hypergeometric expectation.

### HOG-level permutation test

Fisher's method for combining pair-level p-values is invalid within HOGs
because pairs share network neighborhoods and the ortholog mapping,
violating the independence assumption.

`permutation_hog_test()` instead uses a gene-identity permutation null:
for each HOG, M random species-1 genes and N random species-2 genes are
drawn (matching the HOG's gene counts) and the sum-of-fold-enrichments
statistic is recomputed. The permutation p-value is exact regardless of
the dependency structure.

Besag & Clifford (1991) adaptive stopping terminates early once enough
exceedances are observed, and Liang (2016) discrete q-values handle the
non-uniform null distribution that adaptive stopping produces.

### Module comparison

`compare_modules()` supports two methods:

- **Hypergeometric**: Maps module genes through the ortholog table and
  tests overlap with `phyper()`. Q-values via Storey (2003).
- **Jaccard + permutation**: Computes observed Jaccard index on
  ortholog-mappable genes. The null permutes the ortholog mapping
  (Fisher-Yates shuffle), preserving module and network structure. Batched
  permutation shares one shuffle per iteration across all active pairs.
  Q-values via Liang (2016) discrete method.

`classify_modules()` assigns each module to one of three categories:

| Classification | Criteria |
|----------------|----------|
| Conserved | Best-match q < alpha AND Jaccard >= threshold |
| Partially conserved | Best-match q < alpha AND Jaccard < threshold |
| Species-specific | No significant match in the other species |

## Architecture

### R layer

| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()` |
| `R/network.R` | `compute_network()` -- correlation, MR/CLR, density threshold |
| `R/comparison.R` | `compare_neighborhoods()` -- pair-level hypergeometric |
| `R/summary.R` | `summarize_comparison()`, `permutation_hog_test()`, shared q-value helpers |
| `R/modules.R` | `detect_modules()`, `compare_modules()`, `classify_modules()` |
| `R/cliques.R` | `find_coexpression_cliques()` |

### C++ layer (RcppArmadillo + OpenMP)

| File | Purpose |
|------|---------|
| `src/mutual_rank.cpp` | MR normalization with column-major access |
| `src/clr.cpp` | CLR normalization |
| `src/density_threshold.cpp` | Quantile-based density thresholding |
| `src/neighborhood_comparison.cpp` | Pairwise neighborhood overlap |
| `src/hog_permutation.cpp` | HOG permutation engine (bit-vector / flag-vector intersections) |
| `src/module_jaccard_permutation.cpp` | Batched Jaccard permutation engine |

All C++ functions use integer indices only (string mapping is done in R)
due to Homebrew clang ABI constraints. Network matrices are accessed
column-major for cache-friendly reads on symmetric Armadillo matrices.

## References

- Netotea, S. *et al.* (2014). ComPlEx: conservation and divergence of
  co-expression networks in *A. thaliana*, *Populus* and *O. sativa*.
  *BMC Genomics*, 15, 106.
  [doi:10.1186/1471-2164-15-106](https://doi.org/10.1186/1471-2164-15-106)
- Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
  *Biometrika*, 78(2), 301--304.
  [doi:10.1093/biomet/78.2.301](https://doi.org/10.1093/biomet/78.2.301)
- Storey, J. D. & Tibshirani, R. (2003). Statistical significance for
  genomewide studies. *PNAS*, 100(16), 9440--9445.
  [doi:10.1073/pnas.1530509100](https://doi.org/10.1073/pnas.1530509100)
- Liang, K. (2016). False discovery rate estimation for large-scale
  homogeneous discrete p-values. *Biometrics*, 72(2), 639--648.
  [doi:10.1111/biom.12429](https://doi.org/10.1111/biom.12429)
- Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to
  Leiden: guaranteeing well-connected communities. *Scientific Reports*,
  9, 5233.
  [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
- Rosvall, M. & Bergstrom, C. T. (2008). Maps of random walks on complex
  networks reveal community structure. *PNAS*, 105(4), 1118--1123.
  [doi:10.1073/pnas.0706851105](https://doi.org/10.1073/pnas.0706851105)

## License

MIT
