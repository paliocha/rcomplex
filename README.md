# rcomplex

<!-- badges: start -->
[![R-CMD-check](https://github.com/paliocha/rcomplex/actions/workflows/r.yml/badge.svg)](https://github.com/paliocha/rcomplex/actions/workflows/r.yml)
<!-- badges: end -->

Comparative analysis of plant co-expression networks in R

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

# 2. Build co-expression networks (accepts matrix or SummarizedExperiment)
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

# Convert to edge format for clique analysis (multi-species)
edges_AB <- comparison_to_edges(summary$results, "SP_A", "SP_B")

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

# Identify hub genes within modules (6-tier tie-breaking cascade)
hubs1 <- identify_module_hubs(mod1, net1, orthologs,
                              comparison = summary$results)
hubs2 <- identify_module_hubs(mod2, net2, orthologs,
                              comparison = summary$results)
hubs1[hubs1$is_hub, ]

# Classify hub conservation across traits
trait <- c(SP_A = "annual", SP_B = "perennial")
hub_class <- classify_hub_conservation(
  list(SP_A = hubs1, SP_B = hubs2), trait,
  module_comparisons = list("SP_A.SP_B" = comp)
)
hub_class[hub_class$classification != "non_hub", ]

# Query co-expression partners of a hub HOG across all species
partners <- get_coexpressed_hogs("HOG42", networks, orthologs,
                                  species_trait = trait, min_species = 2L,
                                  edges = edges)
partners[partners$coexpressed_traits == "annual", ]      # annual-only partners
partners[grepl(",", partners$coexpressed_traits), ]       # cross-trait partners
```

### Clique-level analysis

```r
# Species and trait definitions
annual_sp    <- c("BDIS", "HVUL", "BMAX", "VBRO")
perennial_sp <- c("BSYL", "HJUB", "BMED", "FPRA")
all_sp       <- c(annual_sp, perennial_sp)
trait <- setNames(rep(c("annual", "perennial"), each = 4), all_sp)

# Find annual-exclusive cliques (C++ Bron-Kerbosch / Tomita pivoting)
cliques <- find_cliques(edges, annual_sp)

# Leave-k-out stability over the full 8-species universe
# max_k defaults to length(all_sp) - 2 = 6, testing all meaningful depths
stab <- clique_stability(edges, annual_sp, trait,
                         all_species = all_sp,
                         full_cliques = cliques, n_cores = 4L)

# Phylogenetically stable cliques (survive any single species dropout)
k1 <- stab$stability[stab$stability$k == 1, ]
stable_cliques <- cliques[k1$clique_idx[k1$stability_score == 1], ]

# Co-expressolog persistence (robustness to threshold tightening)
persist <- clique_persistence(cliques, annual_sp, networks, edges)
persist[persist$persistence > 2.0, ]  # survive 2x stricter thresholds

# Edge-weight robustness metrics (already in find_cliques output)
cliques$intensity   # Onnela geometric mean of (1 - q)
cliques$coherence   # edge weight homogeneity (1 = all equal)

# Bootstrap perturbation test (noise robustness)
pert <- clique_perturbation_test(cliques, annual_sp, networks, orthologs,
                                  n_boot = 100, noise_sd = 0.1)

# Permutation null for clique intensity
z_test <- clique_intensity_test(cliques, annual_sp, networks, orthologs,
                                 edges = edges, n_perm = 500)
```

## Main functions

| Function | Purpose |
|----------|---------|
| `parse_orthologs()` | Parse ortholog group files (tab-delimited) |
| `reduce_orthogroups()` | Merge correlated paralogs within HOGs (Ward.D2 clustering) |
| `extract_orthologs()` | Derive ortholog pairs from two SummarizedExperiment objects by HOG |
| `compute_network()` | Correlation + MR/CLR normalization + density threshold (S4 generic: matrix or SE) |
| `compare_neighborhoods()` | Pair-level hypergeometric neighborhood tests |
| `summarize_comparison()` | Storey q-values and summary statistics |
| `comparison_to_edges()` | Convert comparison results to edge format for clique analysis |
| `permutation_hog_test()` | Permutation-based HOG-level conservation test |
| `detect_modules()` | Community detection (Leiden / Infomap / SBM); iterative multi-resolution consensus |
| `compare_modules()` | Cross-species module overlap (hypergeometric or Jaccard permutation) |
| `classify_modules()` | Three-tier module conservation classification |
| `identify_module_hubs()` | Within-module hub identification with 6-tier conservation-aware tie-breaking |
| `classify_hub_conservation()` | Hub conservation across traits (conserved / rewired / trait-specific) |
| `get_coexpressed_hogs()` | Query co-expression partners of a candidate HOG across species |
| `find_cliques()` | C++ clique detection via Bron-Kerbosch with Tomita pivoting |
| `clique_stability()` | Leave-k-out jackknife stability for trait-exclusive cliques |
| `clique_persistence()` | Co-expressolog persistence scores (robustness to threshold tightening) |
| `clique_threshold_sweep()` | Structural survival of cliques across stricter density thresholds |
| `clique_perturbation_test()` | Bootstrap noise robustness for clique edge weights |
| `clique_intensity_test()` | Permutation null model for clique intensity (Z-score) |
| `classify_cliques()` | Waterfall HOG classification (complete/partial/differentiated/trait_specific) |

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

### SummarizedExperiment integration

`compute_network()` is an S4 generic that accepts either a numeric matrix
or a `SummarizedExperiment`. When passed an SE, it extracts the specified
assay and proceeds as usual:

```r
net1 <- compute_network(se1, assay = "vst.count", norm_method = "MR", density = 0.03)
```

`extract_orthologs()` derives ortholog pairs by matching HOG identifiers
in `rowData()` of two SE objects, producing the same data frame format as
`parse_orthologs()`:

```r
orthologs <- extract_orthologs(se1, se2, hog_col = "hog")
comparison <- compare_neighborhoods(net1, net2, orthologs)
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

1. (Optional) Test K = 1 null via spectral norm permutation: compare the
   leading eigenvalue of the excess co-classification matrix against a
   null from degree-preserving rewiring. If p > 0.05, return a single
   module.
2. Run Leiden at each resolution on the original network.
3. For each gene pair connected in the original thresholded network,
   compute co-classification C (fraction of resolutions placing the pair
   in the same module) and per-pair expected E(i,j) = (1/K) sum_k
   (s_m(i)/N) * (s_m(j)/N). Memory is O(|E|), not O(N^2).
4. Build a sparse consensus graph from edges with positive excess C - E.
5. Run Leiden at all resolutions on the consensus graph.
6. Repeat from step 3 until all resolutions yield the same partition
   (ARI > 0.999), or `max_consensus_iter` is reached.

The sparse edge-restricted co-classification (step 3) replaces the
dense N x N matrix, reducing memory from ~3.2 GB to ~48 MB for N = 20K
genes at 3% density. Re-running the full resolution sweep on the
consensus graph (step 5) avoids the resolution limit that afflicts
single-resolution Leiden on dense graphs (Fortunato & Barthélemy, 2007).

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

### Hub identification

`identify_module_hubs()` identifies hub genes within each species' modules
and reports all three within-module centrality measures (weighted degree,
betweenness, eigenvector), mean edge weight, and global degree. **Weighted
degree is the recommended primary centrality** for co-expression networks:
it directly measures the total co-expression strength to module neighbours,
which is the standard definition of a hub in the co-expression literature.
Betweenness identifies bridge genes between sub-clusters (a different
property), while eigenvector centrality can be unstable on disconnected
subgraphs.

Hub selection uses a 6-tier biologically informed tie-breaking cascade
when genes share the same primary centrality:

1. Primary centrality (user-selected)
2. Global weighted degree (importance beyond the module)
3. Alternative centrality (betweenness if primary is degree; degree otherwise)
4. Mean within-module edge weight (strong edges vs many weak edges)
5. Per-gene conservation effect size (optional; from `summarize_comparison()`)
6. Per-HOG minimum q-value (optional; lower = more conserved)

`classify_hub_conservation()` then classifies each HOG by whether it acts
as a hub across trait groups: **conserved_hub** (hub in both traits, in
corresponding modules), **rewired_hub** (hub in both traits, in
non-corresponding modules), **trait-specific** (hub in one trait only),
**sporadic**, or **non_hub**.

### Clique detection

`find_cliques()` uses a two-level decomposition to find conserved
cliques across species within each ortholog group:

1. **Species-level**: Bron-Kerbosch with Tomita pivoting (Tomita *et al.*,
   2006) on the species adjacency graph (up to 64 species) to enumerate
   all maximal species cliques. The pivot is chosen as the vertex in
   P ∪ X maximising |N(u) ∩ P|, giving worst-case O(3^(n/3)) time.
2. **Gene-level**: Backtracking search assigns one gene per species,
   minimising a composite cost across all C(k, 2) edges. By default
   cost is the mean q-value; the `cost_weights = c(q = 1, effect = 0)`
   argument lets users blend in effect size to favour paralogs with
   stronger enrichment. Early pruning rejects partial assignments where
   any required edge is missing.

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

### Clique persistence

`clique_persistence()` measures how robust each clique's conservation
signal is to threshold tightening. For each clique edge (species pair),
it identifies co-expressologs -- genes that are co-expression neighbours
of the clique gene in both species -- and computes the ratio of the
weakest co-expressolog edge's MR value to the species' density threshold.

A persistence of 1.0 means the weakest supporting edge is exactly at
threshold (marginal). Values above 1.0 indicate the conservation signal
would survive at stricter density thresholds.

### Clique edge-weight robustness

`find_cliques()` returns three per-clique edge-weight summary statistics
following Onnela *et al.* (2005):

- **Intensity**: geometric mean of (1 - q) across all clique edges.
  Higher values indicate uniformly strong conservation signal.
- **Coherence**: ratio of geometric mean to arithmetic mean of (1 - q),
  equal to 1.0 when all edge weights are identical. Low coherence flags
  cliques with a mix of strong and weak edges.
- **min_effect_size**: minimum fold-enrichment across clique edges,
  identifying the bottleneck enrichment.

`clique_threshold_sweep()` returns a `$persistence` element with
`birth`, `death`, and `persistence` per clique (persistence-diagram
compatible), quantifying the range of density thresholds over which each
clique exists.

`clique_perturbation_test()` assesses noise robustness by adding
Gaussian noise to MR scores, re-running network thresholding and clique
detection, and measuring survival of original cliques across bootstrap
replicates.

`clique_intensity_test()` builds a permutation null for clique intensity
by shuffling ortholog assignments, re-running neighborhood comparison
and clique detection, and reporting a Z-score for each clique's observed
intensity relative to the null distribution.

## Architecture

### R layer

| File | Purpose |
|------|---------|
| `R/orthologs.R` | `parse_orthologs()`, `reduce_orthogroups()` |
| `R/network.R` | `compute_network()` -- S4 generic (matrix / SE), correlation, MR/CLR, density threshold |
| `R/comparison.R` | `compare_neighborhoods()`, `comparison_to_edges()`, `get_coexpressed_hogs()` -- pair-level hypergeometric, edge conversion, HOG co-expression queries |
| `R/summary.R` | `summarize_comparison()`, `permutation_hog_test()`, shared q-value helpers |
| `R/modules.R` | `detect_modules()`, `compare_modules()`, `classify_modules()`, `identify_module_hubs()`, `classify_hub_conservation()` |
| `R/cliques.R` | `find_cliques()`, `clique_stability()`, `clique_persistence()`, `clique_threshold_sweep()`, `clique_perturbation_test()`, `clique_intensity_test()`, `classify_cliques()` |
| `R/se_methods.R` | `extract_orthologs()`, `build_se()` (internal) -- SummarizedExperiment helpers |

### C++ layer (RcppArmadillo + OpenMP)

| File | Purpose |
|------|---------|
| `src/reduce_orthogroups.cpp` | Ward.D2 paralog merging engine |
| `src/coclassification.cpp` | Dense and sparse co-classification with per-pair null subtraction; spectral norm for K=1 test |
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
| `src/sample_k_distinct.h` | Shared rejection-sampling utility for subset generation |

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
- Senbabaoglu, Y. *et al.* (2014). Critical limitations of consensus
  clustering in class discovery. *Scientific Reports*, 4, 6207.
  [doi:10.1038/srep06207](https://doi.org/10.1038/srep06207)
- Fortunato, S. & Barthélemy, M. (2007). Resolution limit in community
  detection. *PNAS*, 104(1), 36--41.
  [doi:10.1073/pnas.0605965104](https://doi.org/10.1073/pnas.0605965104)
- Bron, C. & Kerbosch, J. (1973). Algorithm 457: finding all cliques of
  an undirected graph. *Communications of the ACM*, 16(9), 575--577.
  [doi:10.1145/362342.362367](https://doi.org/10.1145/362342.362367)
- Onnela, J.-P. *et al.* (2005). Intensity and coherence of motifs in
  weighted complex networks. *Physical Review E*, 71, 065103.
  [doi:10.1103/PhysRevE.71.065103](https://doi.org/10.1103/PhysRevE.71.065103)
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
