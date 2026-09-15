# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

rcomplex is an R package for comparative co-expression network analysis across species. It maps orthologous genes (via ortholog groups / HOGs from OrthoFinder, FastOMA, PLAZA, etc.), builds co-expression networks independently per species, then tests conservation at four levels:

- **Gene / HOG-level**: Hypergeometric tests with q-value correction (`compare_neighborhoods()` + `summarize_comparison()`), gene-identity permutation with adaptive stopping (`permutation_hog_test()`), batch orchestration (`find_coexpressologs()`, `density_sweep()`), degree-preserving edge-swap null (`coexpressolog_null()`)
- **Module-level**: Community detection (Leiden / Infomap / SBM) with multi-resolution consensus, connectivity preservation testing (`module_preservation()` — permutation null on `avg.weight` + `cor.degree`) over paralog-resolved ortholog maps (`resolve_ortholog_map()`), module correspondence (`module_correspondence()`), module hubs + conservation
- **Clique-level**: two backends. The *species* graph (`find_cliques()` — C++ Bron-Kerbosch / Tomita plus a backtracker that picks one best gene assignment per species clique) carries leave-k-out jackknife stability for trait-exclusive cliques, persistence, threshold sweep, perturbation and intensity tests, and `classify_cliques()`. The per-HOG *gene* graph (`gene_clique_graph()`) enumerates every maximal clique, one per paralog combination, and feeds the five-tier taxonomy of `classify_gene_cliques()`
- **Trait-level**: `preservation_matrix_test()` — relabelling null over an all-pairs preservation matrix (`all_species_pairs()` + `preservation_paired()`), reported under a free and a phylogeny-blocked null; `tag_permutation()` — recurrence of HOGs in diverged modules across designated contrasts. `pvalue_resolution()` reports how much resolution a finite label space leaves

Based on [Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106). The gene-graph clique taxonomy follows [Rodriguez *et al.*, 2026](https://doi.org/10.1038/s41467-026-75624-2).

## Repository status (as of 2026-09-14)

`main` is green: `R CMD check` is `Status: OK`, `devtools::test()` passes
3036/0/7/10 (pass/fail/warn/skip; warns and skips are pre-existing and
environment-gated on missing `sbm`/`torch`), `lintr::lint_package()` reports
no lints, and all four `main` CI workflows (lint, R-CMD-check x2, test-coverage,
pkgdown) pass. There is one open draft PR, `#4`
(`experiment/clique-module-deployment`, "ZDS interpretation"): a large (133
files, 61 commits) speculative branch that predates today's fixes, has
diverged from `main`, and fails CI (including the pre-fix `covr` hang below).
It has been deliberately set aside — do not merge, rebase, or otherwise act
on it without being asked.

## Build & Test

```bash
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'devtools::document()'
R CMD INSTALL .
Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint_package()'
R CMD build . && R CMD check --no-manual rcomplex_0.3.0.tar.gz   # expect "Status: OK"
```

Check the built tarball, not the source directory — `Authors@R` only expands at build time, so `R CMD check .` fails with "Author/Maintainer missing". `--no-manual` avoids needing pdflatex. The historical `R_ext/Boolean.h` warning no longer appears with clang 22. CI runs `lintr::lint_package()` with `LINTR_ERROR_ON_LINT` and there is no `.lintr` file, so any lint fails the build — keep lines at or under 80 characters.

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
| `R/modules.R` | `detect_modules()` (single + consensus), `identify_module_hubs()`, `classify_hub_conservation()`, `characterize_hubs()` |
| `R/ortholog_map.R` | `resolve_ortholog_map()` — reduce multi-copy HOGs toward one counterpart per gene (cliques, then mutual-best coexpressologs); anything left is carried as `unresolved` and settled by majority vote inside `module_preservation()` / `module_correspondence()` |
| `R/module_preservation.R` | `module_preservation()`, `classify_preservation()`, `module_correspondence()`, `preservation_paired()` |
| `R/preservation_matrix.R` | `all_species_pairs()` — every `choose(n, 2)` contrast as a `pairs` table; `preservation_matrix_test()` — trait relabelling null over the all-pairs preservation matrix (free + within-block, enumerated below `enum_max` and sampled above it) |
| `R/pvalue_saturation.R` | `pvalue_resolution()` — distinct values, ties at the minimum and at 1, and permutation-floor status (`permutation-limited` / `evidence-limited` / below floor) of a p- or q-vector |
| `R/tag_permutation.R` | `tag_permutation()` — trait-specific module recurrence test |
| `R/tag_blocks.R` | Internal exchangeability blocks for `tag_permutation()`: `.tp_blocks()` (connected components of the species/contrast graph), `.tp_block_labellings()`, `.tp_expected()` |
| `R/coexpressolog-strength.R` | `coexpressolog_strength()` — density-integrated robustness score for individual coexpressolog edges (re-examines existing networks at several matched densities, no recomputed correlation/normalization); `suggest_reference_density()` — WGCNA-style scale-free-fit diagnostic for picking a reference density |
| `R/cliques.R` | `find_cliques()`, `clique_stability()`, `clique_persistence()`, `clique_threshold_sweep()`, `clique_perturbation_test()`, `clique_intensity_test()`, `classify_cliques()` |
| `R/clique_gene_graph.R` | `gene_clique_graph()` — maximal cliques of the per-HOG (species, gene) graph, every paralog combination reported; `classify_gene_cliques()` — five-tier conservation taxonomy with thresholds derived from `length(species)` |
| `R/se_methods.R` | `extract_orthologs()`, `build_se()` (internal) — SummarizedExperiment helpers |
| `R/rcomplex-class.R` | S3 `rcomplex` container: constructor, print/summary, and a `.rcomplex` method for every pipeline generic registered in `NAMESPACE` |
| `R/rcomplex-package.R` | Package-level roxygen, namespace imports |
| `R/rng.R` | `.seed_scope()` / `.seed_restore()` — the one RNG seeding contract every seeded entry point routes through; `.can_fork()` — gates every `mclapply()` call site on whether forking is safe (see the RNG section below for why) |

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
| `src/module_preservation.cpp` | Module preservation permutation engine (`avg.weight`, `cor.degree`, plus diagnostics) and per-gene intramodular statistics; dense + sparse entry points |
| `src/rewire_degseq.cpp` | Degree-preserving edge-swap kernel for `coexpressolog_null()`: igraph's `keeping_degseq` trial on a bit-matrix adjacency, drawing from R's RNG |
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
| `tests/testthat/test-modules.R` | Module detection and consensus |
| `tests/testthat/test-module-determinism.R` | Consensus determinism on a deliberately ambiguous network: bit-reproducible across core counts, ambient RNG stream unchanged, K = 1 null resolution, consensus fixed point |
| `tests/testthat/test-ortholog-map.R` | Paralog resolution layers, precedence, mappable-set invariant |
| `tests/testthat/test-module-preservation.R` | Preservation kernel vs R reference, null calibration, classification, correspondence, paired directions |
| `tests/testthat/test-preservation-matrix.R` | `all_species_pairs()`, all-pairs relabelling test: label-space sizes and floors, binary difference vs multilevel dispersion, within-block exclusion, saturation reporting |
| `tests/testthat/test-pvalue-saturation.R` | `pvalue_resolution()` counts, floor status, off-grid detection, `suggested_n_perm`, print method |
| `tests/testthat/test-module-hubs.R` | Hub identification, tie-breaking, hub conservation |
| `tests/testthat/test-tag-permutation.R` | Trait-specific module recurrence test |
| `tests/testthat/test-cliques.R` | Clique detection (igraph + C++ backends) |
| `tests/testthat/test-clique-gene-graph.R` | Gene-graph cliques and the five-tier taxonomy: tier waterfall, `missing_reason` states, derived thresholds, duplicate-row collapse, q-floor diagnostics |
| `tests/testthat/test-stability.R` | Leave-k-out jackknife stability |
| `tests/testthat/test-threshold-sweep.R` | Threshold sweep structural survival |
| `tests/testthat/test-perturbation.R` | Clique perturbation (noise robustness) |
| `tests/testthat/test-intensity-test.R` | Clique intensity permutation null |
| `tests/testthat/test-classify-cliques.R` | HOG classification waterfall pipeline |
| `tests/testthat/test-reduce-orthogroups.R` | Paralog reduction |
| `tests/testthat/test-se.R` | SummarizedExperiment integration (build_se, extract_orthologs, S4 compute_network) |
| `tests/testthat/test-rcomplex-class.R` | S3 container construction, printing, method pass-through |
| `tests/testthat/test-rng-contract.R` | The RNG seeding contract, table-driven over every seeded entry point |
| `tests/testthat/helper-reference.R` | Pure-R reference implementations + shared net helpers (`sparse_net()`, ...) |
| `tests/testthat/helper-rng-contract.R` | Fixtures and the seeded-entry-point table for the RNG contract test |
| `tests/testthat/helper-clique-fixtures.R` | Shared clique test fixtures |

## Key Design Decisions

### Sparse network object (dgCMatrix + store_density)
`compute_network(sparse = TRUE)` (default since 0.2.0) stores a `dgCMatrix` with both triangles of entries at or above the `store_density` quantile (default `max(density, 0.05)`); the analysis `threshold` is still computed from the full dense MR, so dense and sparse results are identical. `.net_cpp_args()` / `.net_check()` in `R/network-sparse.R` are the single choke points: every consumer validates there, and a threshold below `store_threshold` errors (the sparse object cannot represent that density). Sub-store values are reconstructible exactly via `mr_block()`.

### Membership-only consumers
Rewired null networks from `coexpressolog_null()` are binary (`threshold = 1`, `store_threshold = 1`): only membership-based consumers (neighborhood comparison, HOG permutation, adjacency extraction) are valid downstream — nothing that reads edge weights (density_sweep multipliers < 1, perturbation, module weights).

### Self-excluded urn
The anchor gene is never its own neighbour, so it leaves the ortholog-mapped set (k) and the hypergeometric population (N - 1) in `compare_neighborhoods()`, both permutation engines, and the torch FE matrix. O(1/N) p-value shift vs canonical ComPlEx; the reported `*.p.val.con` keeps the canonical `x > 1` gate.

### One RNG seeding contract
Every exported function that consumes randomness takes `seed = NULL` and routes it through `.seed_scope()` (`R/rng.R`). The rule: **the caller's stream advances by exactly what the function drew from it, and by nothing else.** `seed = NULL` draws from the ambient stream and leaves it advanced (as `sample()` does); a seed draws from a private stream and restores the caller's on exit, byte for byte, including removing `.Random.seed` when it did not exist before. This replaced three coexisting contracts in 0.3.0 — pinning the exit state at `set.seed(seed)` (which handed a downstream unseeded draw a stream decided by the *upstream* function's seed), restoring the ambient stream, and a bare `set.seed(seed)` with no restore. Batch wrappers (`find_coexpressologs()`, `density_sweep()`, `preservation_paired()`) seed once and pass `seed = NULL` down, so their inner calls draw in sequence instead of all reusing one set of uniforms. `mclapply()` forks never propagate `.Random.seed` to the parent, so core-count reproducibility comes from `.task_seed()`, not from the ambient stream. `tests/testthat/test-rng-contract.R` enforces this table-driven over every seeded entry point and fails if a new one appears without joining the table.

### Randomized-p pi0 (pair level)
Exact hypergeometric p-values pile up at 1 and force Storey's pi0 to 1. `summarize_comparison(pi0_method = "randomized")` (default) estimates pi0 on draws of `p.val.gt + U * p.val.eq` (exactly uniform under H0) and applies it to the exact p-values. `summarize_comparison()`, `find_coexpressologs()` and `density_sweep()` take `seed` (default `NULL`) to pin those draws under the contract above. HOG-level q-values stay `DiscreteQvalue::DQ(method = "Liang")` (Besag–Clifford p-values have discrete support).

### pval_combine default "max"
`comparison_to_edges()` / `summarize_comparison()` / `find_coexpressologs()` / `density_sweep()` combine directional q-values with `pmax` by default: both directions must be significant (reciprocal criterion of Netotea et al. 2014, the `Max.p.val` filter). `"min"` is the permissive either-direction option (pre-0.2.0 behaviour).

### Column-major memory access
Armadillo stores matrices column-major. All hot loops use `colptr()` for sequential reads. ~3x speedup on MR normalization and neighbor list extraction. MR normalization runs in place on the correlation matrix (`mutual_rank_inplace_cpp()`), holding the `compute_network()` peak transient at ~1.5 n^2 doubles.

### Integer indices in C++, string mapping in R
Homebrew clang ABI issue with `std::unordered_map<std::string, ...>`. All C++ uses integer indices; R wrappers handle string-to-int mapping. The codebase also avoids `std::unordered_map` entirely — uses sorted vectors + binary search instead.

### HOG-level testing uses permutation, not Fisher's method
Fisher's method is anti-conservative for multi-copy HOGs (correlated tests). `permutation_hog_test()` permutes gene identities instead.

### Module-level testing measures connectivity, not membership
Gene overlap called a module conserved whenever its membership survived, even when the wiring inside it was gone. `module_preservation()` permutes gene identities with the edges held fixed and tests `avg.weight` (module density) and `cor.degree` (Pearson correlation of intramodular connectivity), combined with `pmax` so both must be significant. `pmax` is valid for that intersection-union null but was measured ~400x conservative, so the combined p-value is first recalibrated against the permutation joint null of the two statistics (`p.calibrated`, argument `calibrate`) and q-values are Benjamini-Hochberg on that; `classify_preservation()` reads its cut point against `Zsummary_std`, which has unit null variance, not raw `Zsummary`. `resolve_ortholog_map()` may only choose *which* paralog copy carries a label, never which genes are **mappable** — but it can still change which genes are **tested**, by rescuing genes whose candidate labels would otherwise tie in the majority vote. The circularity defence therefore rests on the `p_copy` columns of `sensitivity`, a null over random copy choices holding the projected gene set fixed, not on the mappable-set invariant.

### Two clique backends, because they answer different questions
`find_cliques()` cliques the *species* graph — nodes are species, edges are species pairs joined by a co-expressolog — and the C++ backtracker returns **one** best gene assignment per species clique (fewest missing edges first, then the lowest composite cost `cost_weights["q"] * mean_q - cost_weights["effect"] * mean_effect`, whose `c(q = 1, effect = 0)` default is mean-q only). That single-object-per-clique shape is what `clique_stability()`, `clique_persistence()`, `clique_threshold_sweep()` and `clique_perturbation_test()` need: each asks what happens to *a* clique under resampling. It cannot say which paralog combinations form cliques, because it reports one. `gene_clique_graph()` cliques the per-HOG *gene* graph instead — nodes are (species, gene) pairs, edges are the co-expressolog calls — and enumerates every maximal clique, so a multi-copy HOG contributes each combination separately. That is what the published taxonomy needs. Because co-expressolog edges are cross-species, a clique normally carries at most one gene per species; `n_species` records this rather than assuming it, so a caller who supplies within-species edges sees `n_species < n_members` instead of silently wrong species-pair arithmetic. `classify_gene_cliques()` then runs a five-tier waterfall (`complete_conserved`, `lineage_specific`, `partial_significant`, `partial_present`, `differentiated`, else `unclassified`; `underpowered` takes the place of `lineage_specific` or `differentiated`), first match wins, with every threshold derived from `length(species)` — e.g. `partial_significant` tolerates `S - 2` non-significant edges, the largest tolerance under which no member can be isolated. A clique maximal on the strict graph need not be maximal on a looser one, so the intended use is to combine `gene_clique_graph()` runs at several `alpha_graph` values under distinct `id_prefix` values. The tiers separate absent evidence from negative evidence: `missing_reason` is `absent` / `untested` / `tested_ns` / `underpowered` / `extendable`, which is why `classify_gene_cliques()` takes the **unfiltered** edge table — a gap tier has to see the rows that were tested and failed in order to refuse them. A failed test is negative evidence only if it could have succeeded: hypergeometric power rises with degree, so `comparison_to_edges()` writes a per-edge `power` (probability of a call had a fraction `f0` of the neighbourhood been conserved; `NA` on the permutation path) and both classifiers read it through `min_power`. The rule is that a specificity or divergence call must survive treating every underpowered edge as possibly conserved; one that does not is `underpowered`. `NA` power, or no `power` column, is the old behaviour.

### Trait-level testing is all-pairs preservation, not HOG recurrence
Two trait-level tests with different units. `tag_permutation()` asks whether HOGs recur in *diverged* modules across designated contrasts; its null swaps the two trait labels within each contrast, so the label space is the product over connected components of the species/contrast graph (`R/tag_blocks.R`) — `2^k` for a disjoint pairing, the special case rather than the assumption — and `2^-4 = 0.0625` means a design with fewer than five independent contrast groups cannot reach p < 0.05 however strong the signal. `preservation_matrix_test()` drops the recurrence requirement: it reads `Zsummary_std` from every `choose(n, 2)` contrast of an all-pairs `preservation_paired()` run, and the unit is the module-direction, not the HOG. The statistic is `mean(Zsummary_std | concordant) - mean(Zsummary_std | discordant)` for two trait levels (one-sided upward, positive meaning discordant pairs are the less preserved), and for more levels the row-weighted between-class SD of the class means, with concordant rows pooled so the binary case is exactly the difference. Running *all* pairs is what buys the resolution: eight species split 4/4 give `choose(8, 4) = 70` free labellings against the 16 a within-genus null over four genera of two can reach. Both nulls are reported — `p_free` over all species, `p_blocked` within `block` — because their agreement is the phylogenetic-control diagnostic: a free p-value far below the blocked one bought its resolution by breaking phylogenetic control. Within-block pairs are excluded by default, since in a paired design every within-genus pair is trait-discordant and every concordant pair is between-genus, so trait status and phylogenetic distance are perfectly confounded and the confound runs *against* the hypothesis; the exclusion is by block membership, which no relabelling changes, so the null stays valid. The statistic never reads `p.value` or `q.value`.

### Permutation p-values saturate by design
A relabelling p-value cannot go below the tail of the maximum of a **finite** label space, and that space is a property of the species set, not of how long the analysis runs. Relabellings that leave every pair's concordance alone reproduce the statistic exactly — swapping the two trait values everywhere is one — so a binary design always carries at least a twofold tie at the top: `p_attainable` is 2/70 = 0.029 free on eight species split 4/4 and 2/16 = 0.125 inside four blocks of two, and no amount of sampling lowers either. `preservation_matrix_test()` reports `p_min`, `p_attainable`, `n_tied_max` and `exact` per null and warns when the binding floor exceeds 0.05; a null whose space exceeds `enum_max` is sampled and then reports its floor over draws (`1 / (n_perm + 1)`), which can sit far below the floor the design really imposes — only `exact = TRUE` bounds the design. Downstream q-values collapse the same way: on the eight-species Pooideae run, 511 module-directions took 172 distinct q-values with 35 tied at the floor of 0.00071 and 20 at exactly 1, and among those 35 `Zsummary_std` spanned 6.4 to 66.7 — a tenfold range of effect the q-value cannot see. `preservation_paired()` also corrects BH *within* each contrast, which under-corrects an all-pairs matrix by roughly the number of contrasts, and a global BH could not put any module below `n_tests / (n_perm + 1)` (0.26 for 511 tests at `n_perm = 2000`); pass `n_perm_pres` and `$saturation$q_floor_global` reports that number. `pvalue_resolution()` exists to report this resolution wherever a p-value is reported. Rank on `Zsummary_std`, which is standardised to unit null variance; reserve p and q for the significance call.

### Iterative consensus module detection
Multi-resolution Leiden sweep + iterative consensus per Jeub et al. (2018). Per-pair null subtraction: E(i,j) = (1/K) sum_k (s_m(i)/N)(s_m(j)/N), not a scalar mean. Iterates co-classification → Leiden sweep on consensus graph until all resolutions converge (ARI > 0.999).

### Build system
- `Makevars` / `Makevars.win`: C++23, `$(SHLIB_OPENMP_CXXFLAGS)` for portable OpenMP
- RcppArmadillo in `LinkingTo` only (NOT `Imports`)
- `-DARMA_DONT_USE_OPENMP` in both Makevars: Armadillo's *internal* OpenMP deadlocks `mclapply()` workers forked after the parent used OpenMP (see Known CI/tooling gotchas). Package kernels parallelise with their own `#pragma omp`, which the flag does not touch

## Dependencies

**Imports**: methods, Rcpp, Rfast, collapse, data.table, igraph, kit, Matrix (>= 1.5-0), DiscreteQvalue, qvalue, parallel, rlang, stats, utils
**Suggests**: DT, dplyr, knitr, purrr, rmarkdown, S4Vectors, sbm, stringr, SummarizedExperiment, tibble, torch, testthat, lintr, withr, pkgdown
**LinkingTo**: Rcpp, RcppArmadillo
**System**: GNU make, C++23, OpenMP (optional)

## Design notes

`dev/design-notes/` holds research/handoff documents for work that has been
investigated but not yet implemented. `network-sparsification-plan.md`
(2026-09-15) is the operative plan for automatic network sparsification:
gated work packages T1-T6, where T1-T3 must reach a go/no-go decision before
any package code changes. `mdl-engine.md` is the earlier MDL literature
review it supersedes, kept as historical research rather than a
specification. The directory deliberately lives under `dev/`, not `docs/`:
`docs/` is pkgdown's GitHub Pages output directory (see `.Rbuildignore`'s
`^docs$`), and `pkgdown::build_site_github_pages()` refuses to touch a
non-empty `docs/` that it did not itself build.

`prepare_data/` is gitignored and exists only in the maintainer's checkout.
`prepare_data/convert_to_se.R` is the one-time script that converted the
Pooideae long-format `vst_hog.RDS` into the per-species SummarizedExperiment
objects in `prepare_data/data/`. It is not part of the package.

## Known CI/tooling gotchas

- **Forked workers must not enter OpenMP with more than one thread.**
  On Linux, libgomp is not fork-safe: once the parent has run an OpenMP
  region with `n_cores > 1`, a `mclapply()` worker that enters another
  OpenMP region blocks for good (macOS LLVM libomp survives, so it never
  reproduces locally). Package kernels guard their pragmas with
  `if(n_cores > 1)` and workers pass `n_cores = 1`; Armadillo's own
  internal OpenMP is disabled with `ARMA_DONT_USE_OPENMP` for the same
  reason (it hung `detect_modules()`'s K = 1 workers inside
  `arma::eigs_sym()`). Only C/C++ CI (`cpp-check.yml`, plain
  `testthat::test_local()`) reaches the `n_cores = 3` determinism test --
  R-CMD-check's `--as-cran` makes it skip via `_R_CHECK_LIMIT_CORES_` and
  covr disables forking -- so the `n_cores = 2` regression test next to it
  is the one R-CMD-check runs. Both workflows carry `timeout-minutes: 30`.
- **`mclapply()` forking deadlocks under `covr` coverage instrumentation.**
  `covr::package_coverage()` compiles the package with gcov instrumentation;
  forking under gcov is a documented deadlock (a forked child can inherit a
  coverage-counter file lock the parent held at fork time and never release
  it — r-lib/covr#322), not a bug in the forked code itself. `.can_fork()`
  in `R/rng.R` returns `FALSE` whenever `Sys.getenv("R_COVR") == "true"`
  (the same signal `covr::in_covr()` uses), and every `mclapply()` call site
  (`coexpressolog_null()`, and `detect_modules()`'s two sweeps plus its
  batched significance loop in `R/modules.R`) is gated on it. A new fork
  call site must route through `.can_fork()` too, or it will reintroduce
  the hang under `.github/workflows/test-coverage.yml` (which also carries
  a `timeout-minutes: 30` safety net now, so a regression fails fast
  instead of burning hours).
