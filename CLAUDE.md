# CLAUDE.md

Guidance for Claude Code in this repository. Keep sentences short.

## Overview

rcomplex compares gene co-expression networks across species. It maps
orthologous genes (hogs), builds one network per species, and tests
conservation ([Netotea *et al.*,
2014](https://doi.org/10.1186/1471-2164-15-106); clique taxonomy from
[Rodriguez *et al.*, 2026](https://doi.org/10.1038/s41467-026-75624-2)).

Start with the driver
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md).
It takes `expr` (or `networks`) and `orthologs`, and returns edges,
cliques and a classification.
[`print()`](https://rdrr.io/r/base/print.html) shows `r_threshold` per
species: the weakest correlation that passed the density threshold.
[`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
writes the tables.

- Inputs:
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md)
  reads the long table (`species`, `gene`, `hog`).
  [`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md)
  imports a network.
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
  imports modules.
- `clades`: named list of species vectors. It may nest, never cross.
  [`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md)
  derives it. A clique’s home clade is the smallest clade holding all
  its species.
- `partition`: per-sample factor for per-level correlation. `block`: a
  designed sample factor (time point, tree, zone). `sign`: `"positive"`
  or `"negative"` correlation.
- Edge `score` is `-log2(p)` in bits. `evalue` is the expected count of
  pairs this strong by chance.

The four levels: - **Gene / hog**: hypergeometric q-values
([`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
[`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md)),
hog permutation, rank specificity, edge-swap null. - **Module**:
[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md),
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md),
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md),
[`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md),
hubs. - **Clique**:
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
(species graph),
[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
(gene graph),
[`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md),
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md). -
**Trait**:
[`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
on an all-pairs matrix.

## Build & Test

``` bash
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'devtools::document()'
R CMD INSTALL .
NOT_CRAN=true Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint_package()'
Rscript dev/sharpen-prose.R
R CMD build . && R CMD check --no-manual rcomplex_0.4.1.tar.gz   # expect "Status: OK"
```

Check the built tarball, not the directory: `Authors@R` expands at build
time. CI runs
[`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html)
with `LINTR_ERROR_ON_LINT` and there is no `.lintr` file, so keep lines
at or under 80 characters. Only one `rcomplex` installs per library: use
`R CMD INSTALL -l <lib> .` and `R_LIBS=<lib>`.

## Surface budget

Guard: `tests/testthat/test-surface.R`. - Exports: 32. - Named formals
per export: 8 or fewer. Exceptions:
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md),
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md),
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
[`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md),
the matrix method of
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md). -
Each `R/*.R` file: 1,500 lines or fewer. `README.md`: 150 or fewer.

## Prose rules (STE-lite)

Applies to everything a user reads: roxygen, README, quickstart,
messages. - `@description`: 3 sentences at most. Rationale, literature
and caveats go to `@details`. - Sentences under 20 words. One action
each. Active voice. - One term per concept: species, gene, hog, edge,
clique, module, block, partition, clade. Keep identifiers such as
`lineage_specific`, `trait_specific`, `trait_groups`. - Every error
names the argument, says what it got, says what it needs. - Link
exported functions as `[fn()]`. Write internals as `\code{fn()}`. -
`Rscript dev/sharpen-prose.R` checks `R/*.R`, README and quickstart.

## Key Design Decisions

- **Sparse network.** `compute_network(sparse = TRUE)` stores a
  `dgCMatrix` at the `store_density` quantile; `threshold` still comes
  from the dense MR. `.net_check()` in `R/network-sparse.R` is the one
  choke point; a threshold below `store_threshold` errors.
  [`mr_block()`](https://paliocha.github.io/rcomplex/reference/mr_block.md)
  rebuilds sub-store values.
- **Blockwise build (opt-in).** `block_size = b` never forms the n x n
  matrix (`src/network_block.cpp`). It keeps top-fraction candidates per
  column and widens the fraction when it is not exact. The result
  matches the dense build except at floating-point near-ties. Default
  stays dense.
- **Membership-only consumers.** Null networks from
  [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md)
  are binary. Nothing that reads edge weights is valid on them.
- **Self-excluded urn.** The anchor gene leaves the mapped set and the
  population in
  [`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)
  and both permutation engines.
- **One RNG contract.** Every function that draws random numbers takes
  `seed = NULL` and uses `.seed_scope()` (`R/rng.R`). The caller’s
  stream advances by what the function drew, and by nothing else. Forks
  use `.task_seed()`. `test-rng-contract.R` fails on a missing entry
  point.
- **Randomized-p pi0.** Exact hypergeometric p-values pile up at 1.
  `pi0_method = "randomized"` (default) estimates pi0 on
  `p.val.gt + U * p.val.eq`. Hog q-values use
  [`DiscreteQvalue::DQ`](https://rdrr.io/pkg/DiscreteQvalue/man/DQ.html).
- **Specificity score (opt-in).** `method = "rank"` uses
  [`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md)
  (`src/specificity.cpp`) and calibrates against
  [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  partners. Gate rank edges at `min_power = 0.9`. `"hypergeometric"`
  stays the default for ComPlEx equivalence.
- **`pval_combine = "max"`.** Both directions must be significant
  (Netotea *et al.* 2014). `"min"` is the either-direction option.
- **Column-major, integer indices.** Hot loops use `colptr()`. C++ uses
  integer indices and sorted vectors, never `std::unordered_map` with
  string keys (Homebrew clang ABI). R maps strings to integers.
- **Hog permutation, not Fisher’s method.** Fisher is anti-conservative
  for multi-copy hogs.
  [`permutation_hog_test()`](https://paliocha.github.io/rcomplex/reference/permutation_hog_test.md)
  permutes gene identities.
- **Modules: connectivity, not membership.**
  [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
  permutes gene identities with edges fixed and tests `avg.weight` and
  `cor.degree` with `pmax`. The `pmax` p-value is recalibrated against
  the joint null (`p_calibrated`).
  [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.md)
  cuts on `Zsummary_std`.
  [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md)
  picks the paralog copy.
- **Two clique backends.**
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  returns one gene assignment per species clique; stability and sweeps
  need that. A species clique with no shared copy set falls back to its
  largest subsets that have one. No copy cap; ties go to effect size,
  then gene id.
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  returns every maximal clique, one per paralog combination. Twin
  paralogs (same species, same partners) collapse first, so a clique
  carries `multiplicity` and the classifier splits classes by status. A
  hog over `max_cliques` (default 1e5) maximal cliques or copy
  combinations is reported as species blocks: one best copy per species;
  `capped` lists those hogs. A block is one combination to the
  classifier.
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  runs a tier waterfall with thresholds from `length(species)`; `clades`
  and home clades set the specificity and divergence tiers. It takes the
  unfiltered edge table, because a failed test is negative evidence only
  if it had `power`.
  [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  runs three gene graphs (at `alpha`, at 0.9, and every tested pair of
  the unclaimed hogs) and classifies them together; `run` names the
  graph.
- **Trait level.**
  [`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
  reads `Zsummary_std` from every contrast of
  [`all_species_pairs()`](https://paliocha.github.io/rcomplex/reference/all_species_pairs.md).
  It reports `p_free` and `p_blocked`; their gap shows phylogenetic
  confounding.
- **Joint modules are lined-up species modules.**
  [`joint_modules()`](https://paliocha.github.io/rcomplex/reference/joint_modules.md)
  never re-clusters genes. It joins per-species modules by hog overlap
  against a species-pair null (probe: literature note sections 10-12).
  Family-steered clustering lost on planted, wood and leaf.
- **P-values saturate.** A relabelling p-value cannot go below the tail
  of a finite label space. Read `p_attainable` and
  [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md).
  Rank on `Zsummary_std`; keep p and q for the call.
- **Consensus modules.**
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  is Leiden only, with an iterative multi-resolution consensus (Jeub *et
  al.* 2018). Default objective: modularity. Few samples give modules
  that do not replicate.
- **Build.** C++23, `$(SHLIB_OPENMP_CXXFLAGS)`, RcppArmadillo in
  `LinkingTo` only, `-DARMA_DONT_USE_OPENMP` (Armadillo’s OpenMP
  deadlocks forked workers).

## Design notes

`dev/design-notes/` holds plans and research. The operative plan is
`sharpen-plan.md`.

## Known CI/tooling gotchas

- **Forked workers must not enter OpenMP with more than one thread.**
  Linux libgomp is not fork-safe; macOS libomp survives, so it never
  reproduces locally. Kernels guard pragmas with `if(n_cores > 1)` and
  workers pass `n_cores = 1`.
- **`mclapply()` deadlocks under `covr`** (r-lib/covr#322).
  `.can_fork()` in `R/rng.R` returns `FALSE` when `R_COVR == "true"`.
  Every `mclapply()` site must call `.can_fork()` and pass its results
  through `.check_fork_results()`.
- Workflows carry `timeout-minutes: 30`, so a hang fails fast.
