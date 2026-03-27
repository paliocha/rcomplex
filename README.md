# rcomplex

<!-- badges: start -->
[![R-CMD-check](https://github.com/paliocha/rcomplex/actions/workflows/r.yml/badge.svg)](https://github.com/paliocha/rcomplex/actions/workflows/r.yml)
<!-- badges: end -->

Comparative co-expression network analysis across species.

**rcomplex** maps orthologous genes between species, builds co-expression
networks independently, then tests conservation at three complementary
levels:

- **Gene / HOG level** -- Are individual genes' co-expression
  neighborhoods preserved across species?
- **Module level** -- Are higher-order communities (modules) of
  co-expressed genes conserved, partially conserved, or species-specific?
- **Clique level** -- Which fully connected subsets of ortholog groups
  are conserved, and how robust is their trait exclusivity to species
  removal?

Based on [Netotea *et al.*
(2014)](https://doi.org/10.1186/1471-2164-15-106), extended with
permutation-based HOG-level testing, iterative multi-resolution
consensus modules (Jeub *et al.*, 2018), C++ clique detection with
Bron-Kerbosch / Tomita pivoting, and leave-k-out jackknife stability
analysis.

## Installation

```r
devtools::install_github("paliocha/rcomplex")
```

**System requirements:** C++23 compiler, GNU make. OpenMP is optional
but recommended for parallel permutation and stability tests.

## Quick start

```r
library(rcomplex)

# 1. Parse ortholog groups (tab-delimited, see format below)
orthologs <- parse_orthologs("orthogroups.txt", "species1", "species2")

# 2. Build co-expression networks
net1 <- compute_network(expr1, norm_method = "MR", density = 0.03)
net2 <- compute_network(expr2, norm_method = "MR", density = 0.03)

# 3. Compare neighborhoods (pair-level hypergeometric tests)
comparison <- compare_neighborhoods(net1, net2, orthologs)
```

From here, three analysis paths are available.

### Gene / HOG-level analysis

```r
# Pair-level q-values (Storey & Tibshirani, 2003)
summary <- summarize_comparison(comparison)

# HOG-level permutation test (recommended for multi-copy gene families)
hog_results <- permutation_hog_test(net1, net2, comparison, n_cores = 4L)
```

### Module-level analysis

```r
# Detect modules (single resolution)
mod1 <- detect_modules(net1, method = "leiden", resolution = 1.0)

# Or: iterative multi-resolution consensus (Jeub et al., 2018)
mod1 <- detect_modules(net1, resolution = seq(0.1, 5, by = 0.1), seed = 42)
# mod1$resolution_scan shows n_modules, modularity, ARI at each resolution
# mod1$params$n_consensus_iterations shows how many iterations until convergence

mod2 <- detect_modules(net2, resolution = seq(0.1, 5, by = 0.1), seed = 42)

# Compare modules across species (hypergeometric or Jaccard permutation)
comp <- compare_modules(mod1, mod2, orthologs, method = "jaccard", n_cores = 4L)

# Classify: conserved, partially conserved, or species-specific
classes <- classify_modules(comp)
```

### Clique-level analysis

```r
# Find conserved cliques (C++ Bron-Kerbosch / Tomita pivoting + backtracking)
# Edges come from summarize_comparison() with q-values
cliques <- find_cliques(edges, target_species, min_species = 3L)

# Leave-k-out stability for trait-exclusive cliques
# species_trait: any discrete trait (life habit, climate zone, ploidy, ...)
stability <- clique_stability(
  edges, target_species,
  species_trait = c(SP_A = "annual", SP_B = "annual",
                    SP_C = "perennial", SP_D = "perennial"),
  min_species = 3L, max_k = 3L, n_cores = 4L
)

# Identify hub genes recurring across many trait-exclusive cliques
hubs <- clique_hubs(cliques, target_species,
                    species_trait = c(SP_A = "annual", SP_B = "annual",
                                     SP_C = "perennial", SP_D = "perennial"),
                    min_hogs = 3L)
```

## Main functions

| Function | Purpose |
|----------|---------|
| `parse_orthologs()` | Parse ortholog group files (tab-delimited) |
| `reduce_orthogroups()` | Merge correlated paralogs within HOGs (Ward.D2 clustering) |
| `compute_network()` | Correlation + MR/CLR normalization + density threshold |
| `compare_neighborhoods()` | Pair-level hypergeometric neighborhood tests |
| `summarize_comparison()` | Storey q-values and summary statistics |
| `permutation_hog_test()` | Permutation-based HOG-level conservation test |
| `detect_modules()` | Community detection (Leiden / Infomap / SBM); iterative multi-resolution consensus |
| `compare_modules()` | Cross-species module overlap (hypergeometric or Jaccard permutation) |
| `classify_modules()` | Three-tier module conservation classification |
| `find_cliques()` | C++ clique detection via Bron-Kerbosch with Tomita pivoting |
| `clique_stability()` | Leave-k-out jackknife stability for trait-exclusive cliques |
| `clique_hubs()` | Rank genes by recurrence across trait-exclusive cliques |

## Ortholog file format

`parse_orthologs()` reads a tab-delimited file with at least these
columns:

| Column | Description |
|--------|-------------|
| `species` | Species code for the anchor species |
| `gene_id` | Gene identifier in the anchor species |
| `gene_content` | Semicolon-delimited list of `species_code:gene1,gene2,...` entries |

Each row defines one ortholog group from the anchor species' perspective.
The `gene_content` field lists members from other species in the format
`species_code:gene_id1,gene_id2,...`, separated by semicolons.

This format is produced by PLAZA, OrthoFinder, and FastOMA (with
appropriate reformatting). Any tool that can produce a tab-delimited file
with these three columns will work.

### Paralog reduction (optional)

Multi-copy gene families often contain recent duplicates with nearly
identical expression. `reduce_orthogroups()` merges these within each
HOG using Ward.D2 hierarchical clustering on Pearson correlation
distance (1 - r). Paralogs above `cor_threshold` (default 0.7) are
replaced by their averaged expression. Subfunctionalized copies with
distinct expression programs are preserved as separate clusters.

```r
reduced <- reduce_orthogroups(expr1, orthologs, cor_threshold = 0.7)
net1 <- compute_network(reduced$expr_matrix, norm_method = "MR", density = 0.03)
```

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

### Multi-resolution consensus modules

When `detect_modules()` receives a vector of resolutions, it runs
iterative multi-resolution consensus clustering (Jeub *et al.*, 2018):

1. Run Leiden at each resolution on the original network.
2. Build an N x N co-classification matrix C (fraction of resolutions
   placing each gene pair in the same module).
3. Subtract the per-pair expected co-classification under random
   assignment: E(i,j) = (1/K) sum_k (s_m(i)/N) * (s_m(j)/N), keeping
   only excess signal above chance.
4. Run Leiden at all resolutions on the thresholded consensus graph.
5. Repeat from step 2 until all resolutions yield the same partition
   (ARI > 0.999), or `max_consensus_iter` is reached.

The per-pair null subtraction (step 3) replaces the naive scalar
threshold that over-retains edges from low-resolution mega-modules.
Re-running the full resolution sweep on the consensus graph (step 4)
avoids the resolution limit that afflicts single-resolution Leiden on
dense graphs (Fortunato & Barthelemy, 2007).

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

### Clique detection

`find_cliques()` uses a two-level decomposition to find conserved
cliques across species within each ortholog group:

1. **Species-level**: Bron-Kerbosch with Tomita pivoting (Tomita *et al.*,
   2006) on the species adjacency graph (at most 16 nodes) to enumerate
   all maximal species cliques. The pivot is chosen as the vertex in
   P ∪ X maximising |N(u) ∩ P|, giving worst-case O(3^(n/3)) time.
2. **Gene-level**: Backtracking search assigns one gene per species,
   minimising mean q-value across all C(k, 2) edges. Early pruning
   rejects partial assignments where any required edge is missing.

This avoids the combinatorial explosion of enumerating cliques directly
on the gene-level graph when HOGs contain many paralogs.

### Clique stability

`clique_stability()` performs leave-k-out jackknife analysis to assess
phylogenetic robustness of trait-exclusive cliques. A clique is
*trait-exclusive* if all its species share the same value of a discrete
trait (e.g., life habit, climate zone, ploidy level).

For k = 1, 2, ..., max_k species removed at a time:

1. All edges involving the removed species are dropped
2. Full clique detection re-runs on the reduced edge set
3. Reduced cliques are matched to full-dataset cliques by Jaccard
   similarity of gene assignments (ignoring removed species)
4. Trait exclusivity preservation is checked for each match

A clique's *stability score* at level k is the fraction of C(N, k)
species-removal subsets where its trait exclusivity is preserved. The
*stability class* is the highest k at which the score equals 1.0.

The analysis is parallelised over subsets with OpenMP (`n_cores`
parameter) and uses uint64_t bitmask filtering for species membership
(supports up to 64 species).

## Architecture

### R layer

| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()`, `reduce_orthogroups()` |
| `R/network.R` | `compute_network()` -- correlation, MR/CLR, density threshold |
| `R/comparison.R` | `compare_neighborhoods()` -- pair-level hypergeometric |
| `R/summary.R` | `summarize_comparison()`, `permutation_hog_test()`, shared q-value helpers |
| `R/modules.R` | `detect_modules()`, `compare_modules()`, `classify_modules()` |
| `R/cliques.R` | `find_cliques()`, `clique_stability()`, `clique_hubs()` |

### C++ layer (RcppArmadillo + OpenMP)

| File | Purpose |
|------|---------|
| `src/reduce_orthogroups.cpp` | Ward.D2 paralog merging engine |
| `src/coclassification.cpp` | Co-classification matrix with per-pair null subtraction |
| `src/mutual_rank.cpp` | MR normalization with column-major access |
| `src/clr.cpp` | CLR normalization |
| `src/density_threshold.cpp` | Quantile-based density thresholding |
| `src/neighborhood_comparison.cpp` | Pairwise neighborhood overlap |
| `src/hog_permutation.cpp` | HOG permutation engine (bit-vector / flag-vector intersections) |
| `src/fe_permutation.cpp` | GPU-precomputed FE permutation engine |
| `src/module_jaccard_permutation.cpp` | Batched Jaccard permutation engine |
| `src/find_cliques_common.h` | Shared clique primitives (Bron-Kerbosch / Tomita, backtracking, Jaccard, trait) |
| `src/find_cliques.cpp` | C++ clique detection wrapper |
| `src/find_cliques_stability.cpp` | Leave-k-out stability engine with OpenMP |

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
- Lancichinetti, A. & Fortunato, S. (2012). Consensus clustering in
  complex networks. *Scientific Reports*, 2, 336.
  [doi:10.1038/srep00336](https://doi.org/10.1038/srep00336)
- Jeub, L. G. S., Sporns, O. & Fortunato, S. (2018). Multiresolution
  consensus clustering in networks. *Scientific Reports*, 8, 3259.
  [doi:10.1038/s41598-018-21352-7](https://doi.org/10.1038/s41598-018-21352-7)
- Fortunato, S. & Barthelemy, M. (2007). Resolution limit in community
  detection. *PNAS*, 104(1), 36--41.
  [doi:10.1073/pnas.0605965104](https://doi.org/10.1073/pnas.0605965104)
- Bron, C. & Kerbosch, J. (1973). Algorithm 457: finding all cliques of
  an undirected graph. *Communications of the ACM*, 16(9), 575--577.
  [doi:10.1145/362342.362367](https://doi.org/10.1145/362342.362367)
- Tomita, E., Tanaka, A. & Takahashi, H. (2006). The worst-case time
  complexity for generating all maximal cliques and computational
  experiments. *Theoretical Computer Science*, 363(1), 28--42.
  [doi:10.1016/j.tcs.2006.06.015](https://doi.org/10.1016/j.tcs.2006.06.015)
- Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to
  Leiden: guaranteeing well-connected communities. *Scientific Reports*,
  9, 5233.
  [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
- Rosvall, M. & Bergstrom, C. T. (2008). Maps of random walks on complex
  networks reveal community structure. *PNAS*, 105(4), 1118--1123.
  [doi:10.1073/pnas.0706851105](https://doi.org/10.1073/pnas.0706851105)

## License

MIT
