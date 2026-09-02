# rcomplex 0.2.0

Sparse-network release. The network object is now a thresholded sparse
matrix, the whole pipeline runs without densifying, and three statistical
decisions land: the reciprocal `pval_combine = "max"` default, pi0
estimated from randomized p-values, and the self-excluded hypergeometric
urn.

## Sparse network object

- `compute_network(sparse = TRUE, store_density = NULL)` is the new
  default. The network is returned as a `dgCMatrix` holding both triangles
  of the entries at or above the `store_density` quantile (default
  `max(density, 0.05)`; diagonal absent), extracted in C++
  (`extract_sparse_cpp()`) so the dense n x n matrix never leaves
  `compute_network()`; the analysis `threshold` is unchanged (computed from
  the full dense matrix). `sparse = FALSE` returns the pre-0.2.0 dense
  object. New exported `as_sparse_network(net, store_density)` converts
  dense network objects to the identical sparse representation.
- All consumers (`density_sweep()`, `get_coexpressed_hogs()`,
  `detect_modules()` single + consensus, `identify_module_hubs()`,
  `clique_persistence()`, `clique_threshold_sweep()`,
  `clique_perturbation_test()`, `clique_intensity_test()`, `rcomplex()`)
  validate networks through the shared `.net_check()` (storage class,
  dimnames, no stored diagonal, store guard) and return results identical
  to the dense path (`test-network-sparse.R`) -- except
  `clique_perturbation_test()`, which perturbs only the stored entries of a
  sparse network (different RNG stream; edges discarded below the store can
  never be promoted). SBM module detection densifies with a warning.
  Analyses that would need entries below the stored superset (e.g.
  `density_sweep()` with a loose multiplier) error with a message asking
  for a larger `store_density`. `print()`/`summary()` of `rcomplex`
  objects report per-species storage.
- `compare_neighborhoods()` and `permutation_hog_test()` (C++ and torch
  backends) accept a `dgCMatrix` (both triangles stored) as `net$network`
  and return results identical to the dense path. Neighbour-list
  construction is shared in `src/neighbor_lists.h` (dense `arma::mat` or
  validated `dgCMatrix` slots, OpenMP over columns) with internal entry
  points `compare_neighborhoods_sparse_cpp()` /
  `hog_permutation_test_sparse_cpp()`; R-side dispatch, storage-class check
  (square, identical dimnames) and store-threshold guard
  (`store_threshold`, or `min(@x)` for hand-built nets) live in
  `R/network-sparse.R`; torch adjacency edge extraction is factored into a
  testable `.adj_edges()`. `Matrix (>= 1.5-0)` added to Imports.
- In-place MR: normalization now runs in place on the correlation matrix
  (`mutual_rank_inplace_cpp()`); the R-side clip / `abs()` / `diag<-`
  temporaries are gone, so the peak transient of `compute_network()` drops
  from ~4 n^2 to ~1.5 n^2 doubles. The in-place kernel errors on non-double
  input and on `NaN` correlations (e.g. a constant gene kept via
  `min_var = NULL`), which the previous path ranked silently.
- Memory smokes: network object at n = 5000 genes (density 0.03,
  store_density 0.05) is 191.3 MB dense vs 14.9 MB sparse (1.25M stored
  entries); `Rprofmem()` around `find_coexpressologs()` on two sparse 5k
  networks shows zero R-heap allocations >= n^2 * 8 B (200 MB), largest
  single allocation 5 MB. For `compute_network()` (n = 3000) the only
  R-heap allocation >= n^2 * 8 B is the correlation result (72 MB), vs 7
  such allocations (~360 MB) on the 0.1.x path.
- `get_coexpressed_hogs()`: `mean_weight` is now storage-independent. Per
  candidate copy the mean runs only over partner-HOG genes in that copy's
  own neighborhood (edge values at or above the network threshold, hence
  always stored), instead of the union of neighbors across all copies.
  Previously, for multi-copy candidate HOGs, sub-threshold cross-copy
  entries entered the mean, so a sparse network (implicit zeros below the
  store) silently returned different weights than the dense path. Dense and
  sparse results are now identical; single-copy candidate HOGs are
  unaffected.

## Statistical changes

- Self-excluded urn (D5): in `compare_neighborhoods()`, both
  `permutation_hog_test()` engines (bit-vector / flag-vector) and the torch
  fold-enrichment matrix, the anchor gene now leaves the ortholog-mapped
  set (k) and the hypergeometric population (N - 1); it was never its own
  neighbour, so the overlap x is unchanged. P-values shift by O(1/N)
  relative to canonical ComPlEx (fixture, N = 320: identical 149 calls,
  max |diff| of the combined BH value 2.7e-3, mean relative difference
  7e-3; `test-equivalence.R` tolerance widened to 1e-2 with an absolute
  5e-3 guard and a boundary rule for flipped calls).
  `Species*.ortho.neigh` now excludes the anchor; the reported
  `*.p.val.con` keeps the canonical `x > 1` gate. `Species*.jaccard` also
  shifts for pairs whose anchor gene is itself ortholog-mapped: the union
  now uses the post-exclusion k (intended consequence of the excluded
  anchor). Pure-R reference (`helper-reference.R`) updated to the same urn.
- Pair-level pi0 from randomized p-values (D4): `compare_neighborhoods()`
  adds `Species*.p.val.gt` (P(X > x), ungated) and `Species*.p.val.eq`
  (P(X = x)).
  `summarize_comparison(pi0_method = c("randomized", "storey", "none"),
  B = 20L)` estimates Storey's pi0 on `B` draws of the randomized p-value
  `p.val.gt + U * p.val.eq` (lower tail for `alternative = "less"`), which
  is exactly uniform under H0 (Dickhaus et al. 2012), and applies it to the
  exact p-values; the exact hypergeometric p-values pile up at 1 and drive
  the plain Storey estimate to pi0 = 1. `"storey"` is the pre-0.2.0
  behaviour, `"none"` is BH. The estimate is recorded as `$summary$pi0`
  (`c(sp1 = , sp2 = )`). The default draws from the global RNG:
  `set.seed()` first for reproducible q-values.
  `find_coexpressologs(pi0_method = )` passes it through (analytical
  method); `compare_modules()` keeps `"storey"`. Comparison tables without
  the new columns error under the default (use `pi0_method = "storey"` or
  rerun `compare_neighborhoods()`). New `test-pi0.R` (simulation:
  randomized pi0 within 0.05 of the truth, Storey = 1, `"none"` == BH,
  null-only run has no q < 0.05).
- `pval_combine = c("max", "min")` (D2) on `comparison_to_edges()`,
  `summarize_comparison()`, `find_coexpressologs()` and `density_sweep()`:
  by default a co-expressolog is now called only when BOTH directional
  tests are significant -- the reciprocal criterion of Netotea et al.
  (2014), the `Max.p.val` filter of the original ComPlEx. Default results
  therefore contain FEWER conserved edges than 0.1.x; pass
  `pval_combine = "min"` to restore the old either-direction behaviour.
  `density_sweep()` forwards both new arguments to `find_coexpressologs()`
  (unknown values now error instead of vanishing into `...`);
  `clique_threshold_sweep()` is pinned to `pi0_method = "storey"`
  (deterministic q-values across multipliers, matching
  `compare_modules()`). The canonical fixture calls reproduce under the
  defaults: `find_coexpressologs(..., pi0_method = "none")`
  (`test-equivalence.R`).
- Review fixes (roborev 156): the randomized-p pi0 draw now runs over
  ALL comparison rows, before the `filter_zero` filter -- the randomized
  p-value is uniform under H0 only unconditionally, and conditioning on
  overlap > 0 truncated the null and deflated pi0 (anti-conservative
  q-values whenever the expected overlap is small); pi0 therefore refers
  to the full ortholog-pair set. `clique_perturbation_test()` and
  `clique_intensity_test()` gain `pval_combine` (default `"max"`) and
  `pi0_method` (default `"storey"`, deterministic, matching
  `clique_threshold_sweep()`) and forward them to every internal
  `find_coexpressologs()` rerun, so survival rates and the intensity
  null are computed under the same edge-calling criterion that built
  the baseline cliques/edges. `as_sparse_network()` now fails fast when
  the requested `store_density` would give a store threshold above the
  network's analysis threshold (the store would silently drop analysis
  edges), naming `store_density` and both thresholds in the error.
  `coexpressolog_null()`: the serial path (`n_cores = 1` or Windows)
  now saves `.Random.seed` and restores it on exit, so the caller's
  ambient RNG stream continues where the observed run left it, exactly
  as under `mclapply()`; the ineffective `OMP_NUM_THREADS` save/set/
  restore around the workers was removed; and a rewired null run that
  leaves a species pair with no overlap > 0 rows now records 0 for that
  pair under the built-in conserved-count statistic instead of aborting
  (a user-supplied statistic missing a name still errors).

## New functions

- `mr_block(x, genes, net)`: exact local reconstruction of mutual-rank
  values for a gene subset (e.g. a module heatmap), including entries a
  sparse network discarded below `store_threshold`. MR_ij only needs the
  correlations of genes i and j ranked over all n network genes, so the
  k x k block is computed from a k x n correlation slice using the
  network's stored parameters (`cor_method`, `abs_cor`,
  `mr_log_transform`); matches the dense `compute_network(sparse = FALSE)`
  block to 1e-8 for raw MR, log MR, `abs_cor` and Spearman
  (`test-mr-block.R`).
- `coexpressolog_null()`: degree-preserving edge-swap null for
  co-expressolog statistics. Each sparse network is binarised at its
  analysis threshold (not the stored superset), rewired with
  `igraph::keeping_degseq()` (`swap_factor * ecount` swaps; the graph stays
  simple), and `find_coexpressologs()` reruns on the rewired networks with
  exactly the observed run's arguments (`...`). The default statistic
  counts conserved calls per species pair plus `"total"`; any
  `function(edges) -> named numeric` can replace it. Returns a data frame
  (observed, null mean/sd/max, fold, empirical p-value
  `(1 + sum(null >= observed)) / (n_perm + 1)`) with the full
  `n_perm x k` null matrix as `attr(, "null")`. Permutation `b` seeds its
  worker with `seed + b`, so results are identical for any `n_cores`
  (`parallel::mclapply()` on Unix; serial on Windows). Rewired networks are unweighted (`threshold = 1`,
  `store_threshold = 1`), so only membership-based consumers are valid
  downstream. Requires sparse networks (`as_sparse_network()`).

## Validation and documentation

- Equivalence fixture and test against canonical ComPlEx
  (`tests/testthat/test-equivalence.R`): a seeded synthetic dataset with
  expected calls from natstreet/ComPlEx_python (itself validated against
  Hvidsten's `RComPlEx.Rmd`); identical neighborhood overlaps and
  co-expressolog calls (149 pairs) under the ortholog-restricted gene
  universe.
- Roxygen: `@section Gene universe` on `compute_network()`,
  `compare_neighborhoods()` and `find_coexpressologs()` (all-genes universe
  vs the ortholog-restricted canonical implementations);
  `compare_neighborhoods()` `@details` documents the `x > 1` gate, the
  effect size for all x, and the self-excluded urn.
- Internal: the unreachable invalid-index branch of the neighborhood
  comparison now encodes p = 1 as `p.val.gt = 0`, `p.val.eq = 1`, keeping
  the randomized-p decomposition valid in the degenerate case.
- Internal testing hook: `options(rcomplex.force_flag_vector = TRUE)`
  forces the HOG permutation engines into the flag-vector intersection
  mode (the > 100K-gene path), which is now oracle-covered: identical to
  the seeded bit-vector run, to `reference_T_obs()` under the
  self-excluded urn, and dense vs sparse.

# rcomplex 0.1.0

Initial development series (`DESCRIPTION` stayed at 0.1.0 throughout; no
release was tagged). Features accumulated under that version number:

- Networks: `compute_network()` -- Pearson/Spearman correlation (Rfast),
  Mutual Rank (Obayashi log-transformed or raw) or CLR normalization in
  C++ with column-major access, quantile-based density thresholding,
  low-variance gene filter; S4 generic accepting a matrix or a
  `SummarizedExperiment`.
- Gene / HOG level: `compare_neighborhoods()` bidirectional hypergeometric
  neighborhood tests with fold-enrichment effect sizes;
  `summarize_comparison()` with Storey q-values;
  `permutation_hog_test()` gene-identity permutation with Besag & Clifford
  (1991) adaptive stopping (bit-vector / flag-vector C++ engines) and
  Liang (2016) discrete q-values -- replacing Fisher's method, which is
  anti-conservative for multi-copy HOGs; `comparison_to_edges()`.
- torch GPU backend for the permutation test: fold-enrichment matrix
  precomputed on GPU via GEMM, then a C++ permutation loop with table
  lookups; CUDA / MPS / CPU devices with float32 fallback where float64 is
  unsupported.
- Batch orchestration: `find_coexpressologs()` (alias
  `run_pairwise_comparisons()`) across all species pairs;
  `density_sweep()` re-running the pipeline over threshold multipliers.
- Module level: `detect_modules()` (Leiden / Infomap / SBM) with iterative
  multi-resolution consensus clustering (Jeub et al. 2018) -- sparse
  edge-restricted co-classification with per-pair null subtraction and a
  spectral-norm K = 1 permutation test; `compare_modules()`
  (hypergeometric or batched Jaccard permutation),
  `compare_modules_paired()`, `classify_modules()`, `coarsen_modules()`.
- Module hubs: `identify_module_hubs()` (weighted degree, betweenness and
  eigenvector centrality; 6-tier conservation-aware tie-breaking),
  `classify_hub_conservation()` (conserved / rewired / trait-specific /
  sporadic hubs) and `characterize_hubs()`; `tag_permutation()`
  trait-recurrence permutation test; `get_coexpressed_hogs()`
  cross-species co-expression partner queries.
- Clique level: `find_cliques()` C++ two-level decomposition
  (Bron-Kerbosch with Tomita pivoting on the species graph, up to 64
  species, plus backtracking gene assignment with composite q/effect
  cost); `clique_stability()` leave-k-out jackknife with OpenMP and uint64
  bitmask filtering; `clique_persistence()`; `clique_threshold_sweep()`
  with birth/death persistence; `clique_perturbation_test()` bootstrap
  noise robustness; `clique_intensity_test()` permutation null with
  Phipson & Smyth (2010) p-values; Onnela intensity/coherence statistics;
  `classify_cliques()` waterfall classification.
- Orthologs: `parse_orthologs()` (PLAZA / OrthoFinder / FastOMA formats),
  `reduce_orthogroups()` Ward.D2 paralog merging (C++),
  `prepare_orthologs()`.
- SummarizedExperiment integration: `extract_orthologs()` and an internal
  `build_se()`.
- `rcomplex` S3 container class threading networks, orthologs and results
  through the pipeline functions, with print/summary methods.
- Vignette: eight-species Pooideae annual/perennial tutorial with a
  bundled dataset.
