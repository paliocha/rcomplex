# rcomplex 0.3.0

Module-preservation release. Module-level conservation is now a
topology-preservation test, not a gene-overlap test. Overlap called a module
conserved when its membership survived, even if the wiring was gone --
connectivity, not membership, is what selection acts on (Mähler et al.
2017). The overlap hypergeometric was also anti-conservative on multi-copy
HOGs: one HOG with three paralogs contributed three correlated draws to the
same urn. `module_preservation()` replaces it with the permutation test
NetRep computes when only an adjacency matrix is available (Ritchie et al.
2016), and `resolve_ortholog_map()` reduces each multi-copy HOG to one
counterpart per gene before any module label is projected.

## Reproducibility

- `detect_modules(seed = )` now reaches the parallel workers. Under the
  default `RNGkind` a forked `mclapply` child deletes `.Random.seed` and
  re-seeds from clock and PID, so `set.seed()` in the parent never reached
  `igraph::cluster_leiden()`: the same call at `n_cores > 1` could return
  different partitions, and every downstream result inherited that. Each
  parallel task now derives its stream from the seed and its own task index,
  making the result identical at any core count. Reproducibility holds per
  machine and per igraph build, and depends on `RNGkind` as well as `seed`.

- The K = 1 test's early-stopping grid was batched on `n_cores`, so
  `n_perm_completed` and its p-value were functions of the machine (0.040 at
  eight cores versus 0.0476 at one). It now batches on the significance grid
  the serial path already used. Serial results are unchanged; parallel
  results move to match them. **Anyone who ran with `n_cores > 1` should
  expect different K = 1 numbers, and a flipped verdict rewrites the whole
  partition.**

## Breaking changes

- `tag_permutation()` no longer permutes trait labels across all species.
  The design it is applied to is a randomised block --- the phylogenetic
  pair is the block, its two species are the two levels of the trait, and
  the pairs are the replication --- so an unconditional shuffle draws
  labellings in which a pair carries the same trait on both sides. Such a
  pair contributes nothing to the statistic, so the old null mixed the
  observed design with designs having fewer contributing pairs. For four
  pairs, 77% of its support came from those other designs, which is where
  its inflated variance came from (null mean 178 with sd 178) and made it
  anti-conservative: measured type I error 17.2% at a nominal 5%. The null
  now swaps the two labels *within* each pair, so every draw keeps the
  observed design.

  **This changes published numbers, and it moves them the wrong way.** On
  the eight-species Pooideae set the annual-side result goes from
  p = 0.088 to the exact p = 6/16 = 0.375; the "2.7x enrichment" was the
  observed statistic measured against a null mean 77% composed of
  degenerate designs, and against the correct null the ratio is about
  1.04x. Anyone who reported a `tag_permutation()` p-value must recompute
  it.

  The new null is enumerated exactly whenever `2^k <= n_perm` for `k`
  swappable pairs, so `n_perm` is ignored in that case: with four pairs
  the label space has 16 points and `n_perm = 10000` claimed a resolution
  it never had. The return value gains `p_min`, `exact` and `n_swappable`.
  Because the observed labelling is always one of the points,
  **`p_min = 2^-k`, and a design with fewer than five swappable pairs
  cannot produce p < 0.05 for any signal whatsoever** --- at four pairs
  the floor is 0.0625. The function warns when `p_min > 0.05` so a
  non-significant result is not read as evidence of absence when
  significance was unreachable by construction. Note also that running the
  test for two complementary trait values reads two entries of the *same*
  null distribution; they are not independent tests.

  `pairs` must now be disjoint (each species in at most one pair), since a
  within-pair swap would otherwise change another pair's labels.

  Two things measured while fixing this, neither of them changed in code.
  The retired null's type I error grows with the number of pairs, because
  the fraction of its support that reproduces the observed design falls
  as `C(k, a) C(k - a, a) 2^D / C(2k, k)`: 0.16 at four pairs, 0.69 at
  six, and **1.000 at eight or more** --- it rejected on every replicate
  with no signal present. And the "2.7x enrichment" is what zero signal
  produces against that null: with nothing planted the ratio of observed
  to null mean measures 2.80-3.00, so the reported 2.68 sits slightly
  *below* the no-signal expectation.

  Separately, `min_recurrence` does not scale with the number of pairs.
  A HOG reaches two of `k` sides by chance with probability about 0.03 at
  `k = 4` but 0.19 at `k = 10`, so with the default of 2 the statistic
  saturates on chance recurrence and simulated power becomes non-monotone
  in `k`, collapsing by `k = 10`; `max(2, round(k / 2))` restores it. The
  default is unchanged for continuity and the trade-off is now documented
  on the argument --- but anyone acting on the "add pairs" advice above
  should raise it at the same time.

- Removed, with no deprecation shim: `compare_modules()` (with its
  `compare_modules_hypergeometric()` and `compare_modules_jaccard()`
  engines, `best_match_direction()` and `compute_best_matches()`),
  `classify_modules()`, `coarsen_modules()` and `compare_modules_paired()`
  (generic, default and `rcomplex` methods), together with
  `src/module_jaccard_permutation.cpp`. Module-level conservation now runs
  through `module_preservation()` / `preservation_paired()`; the separate
  "which module corresponds to which" question is `module_correspondence()`.
- Classification vocabulary: conserved / partially_conserved /
  species_specific, with the species held in `sp1` / `sp2`, becomes
  conserved / moderate / diverged / untested, with the real species names in
  a `species` column and new `reference` / `test` columns naming the two
  sides of the direction the row describes. `untested` is a `q.value` of
  `NA` -- `cor.degree` undefined because intramodular connectivity is
  constant, so nothing was measured. It is deliberately not folded into
  `diverged`, which would be a positive claim of divergence.
- `classify_hub_conservation(module_comparisons = )` now takes a named list
  of `module_correspondence()` results keyed by the ALPHABETICALLY SORTED
  species pair, and validates that each element carries a `pairs` data frame
  with `module_sp1`, `module_sp2`, `jaccard` and `q.value`. Passing
  something else (e.g. `preservation_paired()$raw`) errors instead of
  leaving every HOG at `NA`, which was indistinguishable from supplying no
  comparison at all.
- `tag_permutation()` retargeted at `preservation_paired()$classification`:
  it takes each contrast's two sides from `reference` / `test` and counts
  only `"diverged"` rows. Preservation reports only modules with at least
  `min_module_size` mapped genes, so the HOG pool is smaller than the
  retired gene-overlap engine's.
- `rcomplex` container: the `module_comparisons` slot is replaced by
  `preservation` (a `preservation_paired()` result) and `correspondence` (a
  `module_correspondence()` result per alphabetically sorted species pair,
  built from the resolved map `preservation_paired()` already computed for
  that orientation). `classify_hub_conservation.rcomplex()` reads the
  latter; `print()` / `summary()` report preservation.
- Fixed: `compare_modules_paired()` keyed its `raw` list unsorted
  (`paste(sp1, sp2)`) while `classify_hub_conservation()` looks up the
  ALPHABETICALLY SORTED pair, so any `phylo_pairs` row with `sp1 > sp2`
  silently missed every lookup and degraded all of its HOGs to
  `multi_trait_hub`. `preservation_paired()` runs the sorted direction as
  one of its two, and the container keys correspondence by the sorted pair.
- Fixed: `classify_preservation()` crashed on a zero-row preservation table
  (the scalar `species` / `pair_name` clashed with the zero-length columns);
  it now returns a zero-row classification.
- `module_correspondence()` computes `jaccard` on the one-to-one projected
  map, so the values run systematically HIGHER than the retired engine's and
  the unchanged `jaccard_threshold = 0.1` of `classify_hub_conservation()`
  is now slightly more permissive. Its q-values default to randomized-pi0,
  which draws from the global RNG, so the `conserved_hub` / `rewired_hub`
  verdict is seed-dependent: `set.seed()` first.
- Preservation is directional -- whether A's modules survive in B is a
  different question from the reverse -- and `preservation_paired()` always
  runs both, so its output is roughly twice the size of the old
  per-contrast output; a contrast listed twice in either orientation is now
  rejected rather than silently overwritten. There is no coarsening and no
  `matched_scale`: preservation never partitions the test species, so a
  module-count ratio is meaningless.

## New functions

- `module_preservation(modules_ref, net_ref, net_test, ...)`: permutation
  test for whether a reference species' modules keep their topology in a
  test network. Two statistics carry the call -- the pair NetRep computes
  when only an adjacency matrix is available (Ritchie et al. 2016):
  `avg.weight` (`sum(kIM) / (m^2 - m)`, the module density) and
  `cor.degree` (Pearson correlation of intramodular connectivity between
  the reference and test networks, i.e. whether hub identity is conserved).
  `meanClusterCoeff` and `meanMAR` are reported as diagnostics only and
  take NO part in the call: a hard-thresholded MR network leaves the
  surviving edge weights nearly constant (max/min ratio about 1.04 at
  density 0.03), so both lose their dynamic range -- including them in a
  median-of-three collapsed a density signal of Z = 124 to Z = 5.3 and
  misclassified a perfectly preserved module. The null shuffles gene
  identities with the edges held constant, each module taking a contiguous
  block of the shuffled genes of its own size, and only ortholog-mappable
  test-species genes enter the shuffle (the NetRep overlap null model).
  One-sided `p = (exceedances + 1) / (n_perm + 1)` per statistic, combined
  across the two with `pmax`, so a module is preserved only when BOTH are
  significant -- the same reciprocal criterion as `pval_combine = "max"`
  elsewhere in the package. `n_perm` defaults to 10000 and sets the p-value
  floor at `1 / (n_perm + 1)`; at 1000 permutations every strongly
  preserved module ties. The combined p-value is recalibrated against the
  permutation joint null of the two statistics (`p.calibrated`) before
  Benjamini-Hochberg: raw `pmax` is valid but was measured roughly 400x
  conservative, with 0.10 the smallest q-value it could emit. The rejection
  region is unchanged, and `calibrate = "none"` recovers the raw result.
  `Zsummary = (Z_avg.weight + Z_cor.degree) / 2` is reported alongside for
  continuity with the WGCNA literature (Langfelder et al. 2011); on
  adjacency-only inputs it is the GWENA `z_summary()` formula reduced to
  the statistics available. Dense and sparse networks share the kernel
  (`src/module_preservation.cpp`, OpenMP over permutations, each iteration
  drawing the same permutation whatever `n_cores` is).
- `resolve_ortholog_map(orthologs, genes1, genes2, ...)`: reduces
  multi-copy HOGs to one counterpart per gene. Cliques first -- a clique
  fixes one gene per species simultaneously, so its copy choices are
  globally consistent across every species at once -- then mutual-best
  coexpressologs (`gene1`'s highest-ranked partner must also rank `gene1`
  highest), then whatever is left, which `module_preservation()` resolves
  by majority vote with ties dropped. Resolution may only choose WHICH
  paralog copy carries a label, never which genes are mappable: filtering
  the mappable set on coexpressolog evidence would select the tested genes
  on the statistic being tested, so every candidate pair whose species-2
  gene no resolved pair claims is kept and the mappable species-2 gene set
  is identical to the one `orthologs` implies.
  `module_preservation(sensitivity = TRUE)` re-runs under a naive map --
  built from `orthologs` alone, no resolution -- and reports both results
  side by side with the attribute recording whether the two runs covered
  the identical test-species genes; a mismatch means the resolution layer
  is filtering rather than choosing. Note the *tested* set can still differ:
  resolving a copy rescues genes whose candidate labels would otherwise tie
  in the majority vote, so the guard reports `same_projected_set` alongside,
  and the circularity itself is measured by the `p_copy` columns -- a null
  over random copy choices that holds the projected gene set fixed.
- `classify_preservation(pres, alpha = 0.05, z_conserved = 10)`: conserved
  (`q.value < alpha` and `Zsummary >= z_conserved`), moderate
  (`q.value < alpha`, weaker), diverged (`q.value >= alpha`) and untested
  (`q.value` is `NA`), carrying the species and pair labels through and
  warning with the module names whenever anything is untested.
- `module_correspondence(modules_ref, modules_test, map)`: module-pair
  cross-tabulation over the resolved map with a hypergeometric
  excess-overlap test, returning `module_sp1`, `module_sp2`, `size_sp1`,
  `size_sp2`, `overlap`, `jaccard`, `p.value` and `q.value` -- the columns
  `classify_hub_conservation()` expects. The test is no longer
  anti-conservative: the resolved map contributes one draw per test-species
  gene instead of one per paralog.
- `preservation_paired(modules, networks, orthologs, pairs, group = )`:
  `module_preservation()` over both directions of every contrast, returning
  `classification` (one row per module per direction), `summary` (counts
  per contrast and direction) and `raw` keyed `"<reference>.<test>"`.
  `"untested"` modules are counted under `"untested"` rather than
  attributed to the reference species' trait group -- nothing was measured,
  so attributing them would overstate the evidence -- and are carried as
  their own level rather than `NA`, which `aggregate()` would drop.
  `rcomplex` method included.

## Validation and documentation

- New `tests/testthat/test-module-preservation.R` (44 blocks): the kernel
  statistics and the observed per-module statistics against a pure-R
  reference, dense vs sparse exactly equal, results independent of
  `n_cores`, permutation p-values uniform and `Zsummary` centred near zero
  under the null, shared module structure detected as preserved and an
  unrelated test network not, a module with constant intramodular
  connectivity reported untested rather than diverged, the `sensitivity`
  re-run (including a copy choice that changes the result and a naive map
  covering different genes), `module_correspondence()` output sanity, and
  `preservation_paired()` covering every module in both directions.
- New `tests/testthat/test-ortholog-map.R` (22 blocks): the resolution
  waterfall and its precedence (a coexpressolog cannot re-claim a gene the
  clique layer resolved), the preserved-gene-set invariant, exactly one
  species-1 partner per resolved species-2 gene, `rank_by` handling
  (`"q.value"` errors on the HOG-level constant q-values that
  `find_coexpressologs(method = "permutation")` broadcasts), and input
  validation.
- `tests/testthat/test-modules.R` reduced to `detect_modules()` coverage
  (single-resolution and consensus); everything it held for the retired
  comparison engine is gone with it.

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
