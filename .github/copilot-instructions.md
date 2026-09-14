# Copilot instructions for rcomplex

## Project scope

`rcomplex` is an R package for comparative co-expression network analysis
across species. It maps genes through ortholog groups (HOGs), builds a
co-expression network per species, and tests conservation at four levels:

- gene/HOG: neighborhood overlap and gene-identity permutation tests;
- module: community detection, paralog-aware projection, and topology
  preservation;
- clique: conserved multi-species structures and robustness analyses;
- trait: all-pairs preservation relabelling and recurrent divergence tests.

The public API is implemented in `R/`. Performance-critical kernels are in
`src/` and use RcppArmadillo plus optional OpenMP. R wrappers own validation,
gene-name mapping, statistical assembly, and dense/sparse dispatch; C++ code
works primarily with integer indices and numeric arrays.

## Build, test, lint, and documentation

Run commands from the repository root.

```bash
# Regenerate RcppExports.R and RcppExports.cpp after changing Rcpp exports
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate NAMESPACE and man/*.Rd from roxygen comments
Rscript -e 'devtools::document()'

# Compile and install the package
R CMD INSTALL .

# Run the full testthat suite
Rscript -e 'devtools::test()'

# Run one test file, e.g. tests/testthat/test-network.R
Rscript -e 'devtools::test(filter = "^network$")'

# Run package lint exactly as CI does
Rscript -e 'lintr::lint_package()'

# Build first, then check the generated source tarball
R CMD build .
R CMD check --no-manual rcomplex_0.3.0.tar.gz
```

Use the generated tarball for `R CMD check`, not the source directory:
`Authors@R` is expanded during build, so `R CMD check .` reports missing
Author/Maintainer fields. `--no-manual` avoids requiring pdflatex.

CI also compiles C++ with `R CMD INSTALL --preclean .` and runs
`testthat::test_local()`. There is no `.lintr`; CI sets
`LINTR_ERROR_ON_LINT=true`, and the default 80-character line limit applies.
`src/Makevars` and `src/Makevars.win` require C++23 and use
`$(SHLIB_OPENMP_CXXFLAGS)` for portable optional OpenMP.

Do not edit `NAMESPACE`, `R/RcppExports.R`, `src/RcppExports.cpp`, or `man/*.Rd`
by hand. Change roxygen/Rcpp declarations and regenerate them.

## Architecture and data flow

The main data flow is:

1. `parse_orthologs()` / `prepare_orthologs()` in `R/orthologs.R` create the
   `Species1`, `Species2`, `hog` mapping used throughout the package.
2. `compute_network()` in `R/network.R` computes correlation, applies MR or
   CLR normalization, derives a density threshold, and normally returns a
   sparse network.
3. `compare_neighborhoods()` in `R/comparison.R` validates the two networks,
   maps gene names to zero-based indices, and calls the dense or sparse C++
   neighborhood kernel. `summarize_comparison()` and
   `permutation_hog_test()` in `R/summary.R` produce pair- and HOG-level
   inference. `find_coexpressologs()` orchestrates species-pair batches.
4. The resulting edge table feeds separate downstream analyses:
   - `R/modules.R`, `R/ortholog_map.R`, and `R/module_preservation.R` detect
     modules, resolve paralogs, and test preservation in both directions;
   - `R/cliques.R` analyzes species-graph cliques and their stability;
   - `R/clique_gene_graph.R` enumerates per-HOG gene/paralog cliques for the
     five-tier conservation taxonomy;
   - `R/preservation_matrix.R` and `R/tag_permutation.R` perform distinct
     trait-level tests.
5. `R/rcomplex-class.R` provides an S3 container that carries named networks,
   orthologs, species metadata, and result slots through the same pipeline.
   `R/se_methods.R` provides the SummarizedExperiment boundary.

Dense and sparse networks share public behavior. Every C++ network consumer
must validate and unpack networks through `.net_check()` and
`.net_cpp_args()` in `R/network-sparse.R`; two-network consumers also use
`.net_pair_sparse()`. Keep string-to-index mapping in R and parallel dense and
sparse entry points in C++.

## Repository-specific conventions and invariants

### Network storage

- Sparse storage is the default. A network is a list containing `network`,
  `threshold`, and metadata; sparse matrices must be `dgCMatrix` objects with
  both triangles, matching non-null dimnames, and no stored diagonal.
- `store_threshold` is a hard representational boundary. An analysis threshold
  below it must error because discarded edges cannot be recovered. Use
  `mr_block()` for exact local reconstruction of sub-store MR values.
- Networks returned by `coexpressolog_null()` are binary rewired networks.
  Only membership-based consumers are valid for them; weighted module,
  perturbation, and looser-density analyses are not.

### Statistical semantics

- Neighborhood tests use a self-excluded hypergeometric urn: the anchor leaves
  both the mapped set and the `N - 1` population. Keep this aligned across the
  R wrappers, C++ kernels, and torch implementation.
- Directional q-values default to `pval_combine = "max"` so both directions
  must be significant. `"min"` is intentionally the permissive alternative.
- Pair-level q-value estimation uses randomized exact p-values; HOG-level
  permutation results use `DiscreteQvalue::DQ(method = "Liang")`.
- Module preservation measures connectivity, not membership. Its call is
  carried by `avg.weight` and `cor.degree`, combines them with `pmax`, and
  applies the configured joint-null calibration before BH adjustment.
  Classification thresholds use `Zsummary_std`, not raw `Zsummary`.
- Permutation p-values have finite resolution. Preserve floor/tie diagnostics
  and rank effects with `Zsummary_std`; use p/q values for significance calls,
  not effect ordering.
- `preservation_matrix_test()` and `tag_permutation()` answer different trait
  questions. The former uses every unordered species pair and reports free and
  phylogeny-blocked relabelling nulls; the latter tests recurrent HOGs across
  designated contrasts.

### Randomness and parallelism

Every exported randomized function accepts `seed = NULL` and enters
`.seed_scope()` from `R/rng.R`:

- `seed = NULL` consumes and advances the caller's ambient RNG stream;
- an explicit seed uses a private stream and restores the caller's state
  exactly, including absence of `.Random.seed`;
- batch wrappers seed once and pass `seed = NULL` to nested calls;
- core-count reproducibility comes from stable per-task seeds, not forked RNG
  propagation.

When adding a randomized entry point, use this contract and add it to the
table-driven coverage in `tests/testthat/helper-rng-contract.R` and
`tests/testthat/test-rng-contract.R`.

### Clique models

Do not interchange the two clique APIs:

- `find_cliques()` builds a species graph and returns one best gene assignment
  per species clique; stability, persistence, threshold, perturbation, and
  intensity analyses consume this shape.
- `gene_clique_graph()` builds a per-HOG `(species, gene)` graph and enumerates
  every maximal paralog combination. Runs at different `alpha_graph` values
  need distinct `id_prefix` values before row-binding. Its classifier needs
  the unfiltered edge table to distinguish absent, untested, non-significant,
  and extendable gaps.

### C++ implementation

- Armadillo matrices are column-major; hot loops use `colptr()` and sequential
  column access.
- Keep gene/string mapping in R. Avoid `std::unordered_map`, especially with
  string keys; existing kernels use integer indices, sorted vectors, and
  binary search for compiler/ABI portability.
- Preserve dense-versus-sparse equivalence whenever changing a network
  consumer. `tests/testthat/helper-reference.R` contains pure-R reference
  implementations used to check kernels.

## Active design direction

Plan a density-integrated strength measure for an individual co-expressolog
edge: an orthologous gene pair whose network neighbourhood is conserved across
two species. This is deliberately narrower than a Marbach/Seidr-style
multi-algorithm crowd. Density settings are dependent views of one inference
pipeline, not independent voters.

The feature should complement, not replace, existing robustness signals:

- pair-level q-values, effect sizes, and Jaccard overlap;
- `density_sweep()` across analysis thresholds;
- clique intensity and coherence from `find_cliques()`;
- `clique_persistence()` and `clique_threshold_sweep()`;
- `clique_perturbation_test()` for MR-score noise.

`clique_stability()` removes species, not expression samples. Density
integration measures sensitivity to network sparsification; species
leave-k-out measures taxon-composition robustness; sample bootstrap or
subsampling would measure expression-data uncertainty. Keep these concepts and
outputs separate.

### Proposed public entry point

Add `coexpressolog_strength()` as an S3 generic with:

- `coexpressolog_strength.default(networks, orthologs, densities,
  reference_density, species_pairs, alpha, ...)`;
- `coexpressolog_strength.rcomplex()` in `R/rcomplex-class.R`, storing the
  result in a new `edge_strength` slot;
- implementation and helpers in a focused `R/coexpressolog-strength.R` file.

Return a list rather than only a summary data frame:

- `strength`: one row per species-pair/coexpressolog edge;
- `profiles`: the long-form per-edge, per-density measurements;
- `reference`: ordinary inferential results at `reference_density`;
- `params`: density grid, reference density, method, and scoring choices.

The edge summary should initially expose `strength`, `median_rank`,
`density_support`, `strictest_density`, `rank_iqr`, and `n_transitions`.
`strength` is a bounded descriptive ranking score, not a probability.

### Density and network entry points

Accept explicit matched densities such as
`c(0.01, 0.02, 0.03, 0.05, 0.075)`, not only threshold multipliers.
Equivalent densities make sparsity comparable across species; a common
threshold multiplier does not.

Add one internal threshold-at-density helper in `R/network-sparse.R` and use it
from both `density_sweep()` and the new feature where practical. Its result
must exactly match `density_threshold_cpp()` for dense networks. For a sparse
network, every requested density must be no greater than `store_density`;
otherwise error and instruct the caller to rebuild with a larger store. Respect
ties and the package's both-triangles storage convention.

At each density, make shallow copies of network objects with the derived
analysis thresholds, then call `compare_neighborhoods()`. Do not recompute
correlations or normalization. Keep every tested ortholog pair, including
zero-overlap and non-significant rows, so ranks use a stable comparison
universe.

### Per-density evidence

Add an internal helper that derives a directional hypergeometric Z-score from
the existing neighborhood counts. For direction 1:

- population `M = nrow(net1$network) - 1`, preserving the self-excluded urn;
- draws `m = Species1.neigh`;
- successes `k = Species1.ortho.neigh`;
- observed `x = Species1.neigh.overlap`;
- expectation `m * k / M`;
- variance from the standard hypergeometric finite-population formula.

Direction 2 uses the corresponding species-2 fields. Define reciprocal
evidence as `pmin(z1, z2)`, consistent with the package's requirement that
both directions support conservation. Undefined zero-variance cases should
remain explicit rather than being silently assigned strong evidence.

Within each species pair and density, convert reciprocal evidence to a
normalized percentile rank in `[0, 1]`, with 1 strongest. Resolve ties
deterministically using, in order, reciprocal Z, minimum directional fold
enrichment, geometric-mean Jaccard, minimum overlap, and stable gene IDs.
Do not use `1 - q.value` as the primary strength: q-values are inferential,
depend on the tested family, and may tie at a permutation floor.

### Across-density summary

Summarize densities within each edge; do not treat densities as independent
votes. The proposed primary score is the area under normalized rank versus
log-density, scaled to `[0, 1]`. Also report:

- median normalized rank and its IQR;
- fraction of densities satisfying the chosen reference call criterion;
- sparsest density supporting that call;
- number of supported/unsupported transitions across the ordered grid.

Support need not be monotonic because neighbourhood sizes and multiple-testing
correction change with density. Retain the full profile rather than reducing
the curve to its strictest supported density.

Run q-value inference at one explicit `reference_density`; do not combine
dependent p-values across the density grid. The first implementation should
use analytical neighborhood comparisons for the curve and reserve HOG
permutation inference for the reference analysis.

### Integration and validation entry points

- Add the generic and methods to roxygen-generated `NAMESPACE`.
- Extend the `rcomplex` constructor, print, and summary methods with the new
  result slot without changing existing pipeline slots.
- Keep current `find_coexpressologs()` edge columns stable. Downstream clique
  code should opt into a strength column only after the edge-level behavior is
  established; do not silently replace q-based clique weights.
- Add focused tests in `tests/testthat/test-coexpressolog-strength.R`.
- Test exact dense/sparse agreement, density/store-boundary errors,
  self-excluded Z calculations against an R reference, deterministic ties,
  species-pair isolation, input-order invariance, non-monotone profiles,
  reference-density equivalence, empty inputs, and S3 pass-through.
- Compare ranking reproducibility under sample subsampling against
  single-density q-value, effect-size, and Jaccard rankings before promoting
  the score as a robustness measure.

Potential multi-method rank aggregation remains a later extension. If added,
first aggregate densities within each inference method and then give each
method one vote; never let a method gain weight merely because it was evaluated
at more densities.
