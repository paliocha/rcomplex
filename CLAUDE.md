# CLAUDE.md

Guidance for Claude Code in this repository. Keep sentences short.

## Overview

rcomplex compares gene co-expression networks across species. It maps
orthologous genes (hogs), builds one network per species, and tests
conservation ([Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106);
clique taxonomy from [Rodriguez *et al.*, 2026](https://doi.org/10.1038/s41467-026-75624-2)).

Start with the driver `rcomplex()`. It takes `expr` (or `networks`) and
`orthologs`, and returns edges, cliques and a classification. `print()`
shows `r_threshold` per species: the weakest correlation that passed the
density threshold. `write_rcomplex()` writes the tables.

- Inputs: `read_orthologs()` reads the long table (`species`, `gene`,
  `hog`). `as_network()` imports a network. `as_modules()` imports modules.
- `clades`: named list of species vectors. It may nest, never cross.
  `clades_from_tree()` derives it. A clique's home clade is the smallest
  clade holding all its species.
- `partition`: per-sample factor for per-level correlation. `block`: a
  designed sample factor (time point, tree, zone). `sign`: `"positive"`
  or `"negative"` correlation.
- Edge `score` is `-log2(p)` in bits. `evalue` is the expected count of
  pairs this strong by chance.

The four levels:
- **Gene / hog**: hypergeometric q-values (`find_coexpressologs()`,
  `density_sweep()`), hog permutation, rank specificity, edge-swap null.
- **Module**: `detect_modules()`, `module_preservation()`,
  `preservation_paired()`, `module_correspondence()`, hubs.
- **Clique**: `find_cliques()` (species graph), `gene_clique_graph()`
  (gene graph), `classify_cliques()`, `classify_gene_cliques()`.
- **Trait**: `preservation_matrix_test()` on an all-pairs matrix.

## Build & Test

```bash
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'devtools::document()'
R CMD INSTALL .
NOT_CRAN=true Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint_package()'
Rscript dev/sharpen-prose.R
R CMD build . && R CMD check --no-manual rcomplex_0.4.0.tar.gz   # expect "Status: OK"
```

Check the built tarball, not the directory: `Authors@R` expands at build
time. CI runs `lintr::lint_package()` with `LINTR_ERROR_ON_LINT` and there
is no `.lintr` file, so keep lines at or under 80 characters. Only one
`rcomplex` installs per library: use `R CMD INSTALL -l <lib> .` and
`R_LIBS=<lib>`.

## Surface budget

Guard: `tests/testthat/test-surface.R`.
- Exports: 31.
- Named formals per export: 8 or fewer. Exceptions: `rcomplex()`,
  `module_preservation()`, `find_coexpressologs()`, `density_sweep()`, the
  matrix method of `compute_network()`.
- Each `R/*.R` file: 1,500 lines or fewer. `README.md`: 150 or fewer.

## Prose rules (STE-lite)

Applies to everything a user reads: roxygen, README, quickstart, messages.
- `@description`: 3 sentences at most. Rationale, literature and caveats
  go to `@details`.
- Sentences under 20 words. One action each. Active voice.
- One term per concept: species, gene, hog, edge, clique, module, block,
  partition, clade. Keep identifiers such as `lineage_specific`,
  `trait_specific`, `trait_groups`.
- Every error names the argument, says what it got, says what it needs.
- Link exported functions as `[fn()]`. Write internals as `\code{fn()}`.
- `Rscript dev/sharpen-prose.R` checks `R/*.R`, README and quickstart.

## Key Design Decisions

- **Sparse network.** `compute_network(sparse = TRUE)` stores a
  `dgCMatrix` at the `store_density` quantile; `threshold` still comes
  from the dense MR. `.net_check()` in `R/network-sparse.R` is the one
  choke point; a threshold below `store_threshold` errors. `mr_block()`
  rebuilds sub-store values.
- **Blockwise build (opt-in).** `block_size = b` never forms the n x n
  matrix (`src/network_block.cpp`). It keeps top-fraction candidates per
  column and widens the fraction when it is not exact. The result matches
  the dense build except at floating-point near-ties. Default stays dense.
- **Membership-only consumers.** Null networks from `coexpressolog_null()`
  are binary. Nothing that reads edge weights is valid on them.
- **Self-excluded urn.** The anchor gene leaves the mapped set and the
  population in `compare_neighborhoods()` and both permutation engines.
- **One RNG contract.** Every function that draws random numbers takes
  `seed = NULL` and uses `.seed_scope()` (`R/rng.R`). The caller's stream
  advances by what the function drew, and by nothing else. Forks use
  `.task_seed()`. `test-rng-contract.R` fails on a missing entry point.
- **Randomized-p pi0.** Exact hypergeometric p-values pile up at 1.
  `pi0_method = "randomized"` (default) estimates pi0 on
  `p.val.gt + U * p.val.eq`. Hog q-values use `DiscreteQvalue::DQ`.
- **Specificity score (opt-in).** `method = "rank"` uses
  `compare_specificity()` (`src/specificity.cpp`) and calibrates against
  `null_network()` partners. Gate rank edges at `min_power = 0.9`.
  `"hypergeometric"` stays the default for ComPlEx equivalence.
- **`pval_combine = "max"`.** Both directions must be significant
  (Netotea *et al.* 2014). `"min"` is the either-direction option.
- **Column-major, integer indices.** Hot loops use `colptr()`. C++ uses
  integer indices and sorted vectors, never `std::unordered_map` with
  string keys (Homebrew clang ABI). R maps strings to integers.
- **Hog permutation, not Fisher's method.** Fisher is anti-conservative
  for multi-copy hogs. `permutation_hog_test()` permutes gene identities.
- **Modules: connectivity, not membership.** `module_preservation()`
  permutes gene identities with edges fixed and tests `avg.weight` and
  `cor.degree` with `pmax`. The `pmax` p-value is recalibrated against the
  joint null (`p_calibrated`). `classify_preservation()` cuts on
  `Zsummary_std`. `resolve_ortholog_map()` picks the paralog copy.
- **Two clique backends.** `find_cliques()` returns one gene assignment
  per species clique; stability and sweeps need that. `gene_clique_graph()`
  returns every maximal clique, one per paralog combination.
  `classify_gene_cliques()` runs a tier waterfall with thresholds from
  `length(species)`; `clades` and home clades set the specificity and
  divergence tiers. It takes the unfiltered edge table, because a failed
  test is negative evidence only if it had `power`.
- **Trait level.** `preservation_matrix_test()` reads `Zsummary_std` from
  every contrast of `all_species_pairs()`. It reports `p_free` and
  `p_blocked`; their gap shows phylogenetic confounding.
- **P-values saturate.** A relabelling p-value cannot go below the tail
  of a finite label space. Read `p_attainable` and `pvalue_resolution()`.
  Rank on `Zsummary_std`; keep p and q for the call.
- **Consensus modules.** `detect_modules()` is Leiden only, with an
  iterative multi-resolution consensus (Jeub *et al.* 2018). Default
  objective: modularity. Few samples give modules that do not replicate.
- **Build.** C++23, `$(SHLIB_OPENMP_CXXFLAGS)`, RcppArmadillo in
  `LinkingTo` only, `-DARMA_DONT_USE_OPENMP` (Armadillo's OpenMP
  deadlocks forked workers).

## Design notes

`dev/design-notes/` holds plans and research. The operative plan is
`sharpen-plan.md`.

## Known CI/tooling gotchas

- **Forked workers must not enter OpenMP with more than one thread.**
  Linux libgomp is not fork-safe; macOS libomp survives, so it never
  reproduces locally. Kernels guard pragmas with `if(n_cores > 1)` and
  workers pass `n_cores = 1`.
- **`mclapply()` deadlocks under `covr`** (r-lib/covr#322). `.can_fork()`
  in `R/rng.R` returns `FALSE` when `R_COVR == "true"`. Every
  `mclapply()` site must call `.can_fork()` and pass its results through
  `.check_fork_results()`.
- Workflows carry `timeout-minutes: 30`, so a hang fails fast.
