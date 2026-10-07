# Sharpen plan: MDO I + II, BLAST-shaped UX

2026-10-07. Baseline `main` b245ac0 (0.3.2.9000). Target release 0.4.0.

Goal. Hvidsten proposal MDO I (multi-species, phylogeny-aware) and II
(nulls, effect sizes) are already covered by rcomplex. What is missing is
the BLAST/Clustal shape: one verb, data in, table out, defaults good.
Package is 57 exports, 18.7k R lines, 1.5k-line tutorial. Most of it is
research scaffolding from the 2026-09/10 Orion benchmarks. Cut it, keep
the math, add one driver and nested clades.

Framing. rcomplex is an ortholog-anchored *local* GCN aligner. Pairwise
alignment = co-expressologs (the BLAST hit); multiple alignment = cliques
across N species (the Clustal column); modules = aligned blocks. The
aligner must accept any GCN, from any expression design, not only the
Pooideae annual/perennial or the EVOTREE wood data. Three consequences
(Martin, 2026-10-07): input generalisation (WP10), harmonisation of
heterogeneous expression designs before alignment (WP11, TEA-GCN as the
reference), and a controlled prose style for everything a user reads
(WP8, Karpathy's ASD-STE100 note, section 7).

Regime rule (Martin, 2026-10-07: most people cannot run giant
experiments). rcomplex works the same way at n = 6 and n = 6,000
samples: same verbs, same output columns, no default that assumes a
sample count. What changes with n is power and resolution, and the
package says so instead of hiding it: `print()` shows n and the
correlation each network needed to pass its density threshold, the
`power` column already scales with degree, `pvalue_resolution()`
already reports the floor, and an opt-in null pass reports the false
calls at the user's alpha on the user's n (WP12). Calibration must
hold at every n; power may not, and is reported. Harmonisation has one
tool per regime (WP11): `block` for a designed experiment with
replicates (needs two samples per level), `partition` for a compendium
(needs `min_partition_n` per level). The maintainer's own data is the
small-n case (Pooideae, 20 samples per species, five time points by
four replicates), so small n is the default test regime, not an edge
case.

Rules for every work package: ponytail ladder (deletion over addition, no
new abstraction, no new object slots), karpathy guidelines (surgical
edits, state assumptions), lines <= 80, lint clean, no Claude attribution
lines in commits or PRs. Worktree agents commit before finishing.

## 1. Baseline (measured, `dev/sharpen-census.R`)

| Metric | Now | Target |
|---|---|---|
| Exports | 57 | <= 30 |
| R lines | 18,704 | <= 12,000 |
| C++ lines | 5,763 | <= 5,000 |
| Test lines | 23,209 | <= 15,000 |
| Largest R file | cliques.R 2,719 | <= 1,500 any file |
| Max named formals | module_preservation 19 | <= 8 any export |
| README | 511 lines | <= 150 |
| Tutorial | 1,472 lines | quickstart <= 150, walkthrough <= 800 |
| Species column names | `Species1` (64), `species1` (20), `sp1` (14) | one |
| Stat column names | `q.value` (16) vs `effect_size` (11) | one |
| `rcomplex()` container | 621 lines, 15 S3 methods, 0 uses in tutorial | driver |

Export census (uses by location; `internal` = calls from other package
code, `readme`/`tutorial`/`methods` = mentions, `tests` = calls):

```
fn                          internal readme tutorial methods tests
split_layers                       1      0        0       0     7
run_pairwise_comparisons           0      1        0       0     8
rcomplex                           1      1        0       1    44
mr_block                           2      1        0       2    16
summarize_specificity              6      1        0       2    20
as_sparse_network                  9      1        0       2    23
summarize_comparison              25      1        0       3    55
coexpressolog_strength             0      2        0       1     1
bicm_species_z                     1      2        0       0     1
conservation_lattice               1      2        0       0     3
as_preservation_matrix             2      2        0       0     2
extract_orthologs                  2      2        0       3     5
all_species_pairs                  3      1        1       2    14
suggest_reference_density          4      2        0       1     0
coexpressolog_null                 5      2        0       3    31
conservation_pattern_table         6      2        0       0     5
compare_specificity                9      2        0       2    32
comparison_to_edges               15      1        1       2    31
resolve_ortholog_map              16      2        0       3    58
as_modules                        19      2        0       1    56
clique_intensity_test              0      1        2       3    34
module_auroc_reciprocal            0      3        0       0     3
module_replication                 0      3        0       0     6
clique_perturbation_test           0      1        2       3    25
recurrence_modules                 1      3        0       0     5
clique_persistence                 0      1        2       3    15
characterize_hubs                  3      1        2       1    11
subspace_preservation              3      3        0       0     8
permutation_hog_test               4      2        1       3    52
classify_preservation              5      3        0       3    20
parse_orthologs                    5      3        0       3     0
recurrence_graph                   6      3        0       0     5
identify_module_hubs               8      2        1       4    38
compare_neighborhoods             28      2        1       3    77
get_coexpressed_hogs               1      2        2       2    20
density_sweep                      3      3        1       2    29
tag_permutation                    4      2        2       2    65
clique_stability                   5      2        2       3    42
module_correspondence             15      2        2       3    14
detect_modules                    28      3        1       5    90
reduce_orthogroups                 1      2        3       3    15
prepare_orthologs                  2      3        2       1    13
clique_threshold_sweep             3      1        4       2    20
pvalue_resolution                  8      3        2       4    40
null_network                      13      5        0       2    22
module_auroc                      17      5        0       0    12
compute_network                   56      4        1      12   137
classify_hub_conservation          6      2        4       4    29
module_preservation               26      5        1       7    71
classify_cliques                   4      2        5       2    32
preservation_paired               19      3        4       3    28
preservation_matrix_test          12      5        3       4    71
classify_gene_cliques              9      5        4       4    85
gene_clique_graph                 10      6        3       5   104
find_cliques                      29      2        7       7    71
find_coexpressologs               43      6        3       4    83
```

## 2. Export triage

**Tier A, stays exported (29).** The user path plus the building blocks a
power user needs by name. `split_layers` is the harmonisation primitive
(WP11), not scaffolding. `read_orthologs` and `as_network` are new
(WP10); `parse_orthologs` folds into `read_orthologs`.

```
read_orthologs  prepare_orthologs  reduce_orthogroups
split_layers  compute_network  as_network  null_network  coexpressolog_null
find_coexpressologs  density_sweep
gene_clique_graph  classify_gene_cliques
find_cliques  clique_stability  clique_threshold_sweep  classify_cliques
as_modules  detect_modules  resolve_ortholog_map  module_preservation
classify_preservation  module_correspondence  preservation_paired
preservation_matrix_test  identify_module_hubs  classify_hub_conservation
get_coexpressed_hogs  pvalue_resolution  rcomplex (driver, WP6)
```

**Tier B, demote to internal (11).** Called by Tier A, never needed by
name. Keep the function, drop `@export`, add `@keywords internal`.
Tests call them as `rcomplex:::fn()`.

```
compare_neighborhoods  summarize_comparison  comparison_to_edges
permutation_hog_test  compare_specificity  summarize_specificity
run_pairwise_comparisons  all_species_pairs  mr_block  as_sparse_network
extract_orthologs
```

**Tier C, remove from package (17).** Research probes from the Orion
benchmark (PRs #36-#43, design note module-engine-borrow-scope 10.x) and
clique robustness tests whose nulls are structurally unavailable
(memory: clique metric ceiling). None on the user path. Source moves
verbatim to `dev/probes/<file>.R` with a two-line header (was
`rcomplex::fn` until 0.4.0, see git history); tests, Rd, README and
tutorial sections deleted. Git keeps everything.

```
module_auroc  module_auroc_reciprocal  module_replication
subspace_preservation  as_preservation_matrix
recurrence_graph  recurrence_modules
conservation_pattern_table  conservation_lattice  bicm_species_z
coexpressolog_strength  suggest_reference_density
clique_persistence  clique_perturbation_test  clique_intensity_test
tag_permutation  characterize_hubs
```

Estimated removal: ~5,900 R lines (module_auroc.R 452, subspace.R 264,
recurrence_graph.R 446, clique_patterns.R 338, coexpressolog-strength.R
783, tag_permutation.R 980, tag_blocks.R 140, cliques.R lines 893-2309
~1,400, modules.R characterize_hubs ~170, rcomplex-class.R methods for
Tier C ~150), ~4,600 test lines, src/module_auroc.cpp 117,
src/subspace.cpp 37, src/fe_permutation.cpp if only strength uses it
(check).

`tag_permutation` is the one judgement call: it is one of the two trait
tests in CLAUDE.md, but `preservation_matrix_test()` supersedes it on
resolution (2^-4 floor vs 2/70) and the tutorial already labels it
"secondary". Gate G1 decides.

## 3. Work packages

Integration branch `refactor/sharpen` off `main`. Each WP is one PR into
it. One PR `refactor/sharpen` -> `main` after WP9, version 0.4.0.
Acceptance commands run from the package root after
`Rscript -e 'devtools::document()' && R CMD INSTALL .`.

Order: WP0 -> {WP1 || WP2 || WP3} -> WP4 -> WP5 -> {WP7 || WP10 || WP11}
-> {WP13 || WP14} -> WP6 -> WP12 -> WP8 -> WP9 = 0.4.0. Then, gated
and not blocking the release: WP15 || WP16.

### WP0 Branch + surface snapshot (serial, first)

- Files: `tests/testthat/test-surface.R` (new), `dev/sharpen-census.R`.
- Do: `git checkout -b refactor/sharpen main`. Test snapshots sorted
  `getNamespaceExports("rcomplex")` and, per export, the count of named
  formals (excluding `...`). `expect_snapshot()` so every later removal
  shows in the diff.
- Accept: `Rscript -e 'testthat::test_local(filter = "surface")'` passes;
  `_snaps/surface.md` lists 57 exports.
- Deps: none.

### WP1 Demote Tier B (parallel with WP2, WP3)

- Files: `R/comparison.R`, `R/specificity.R`, `R/summary.R`,
  `R/mr_block.R`, `R/network-sparse.R`, `R/se_methods.R`,
  `R/preservation_matrix.R` (`all_species_pairs`),
  `R/rcomplex-package.R` (index), `NAMESPACE`, `man/`, `README.md`
  function index rows, tests that call Tier B by name.
- Do: remove `@export`, add `@keywords internal`. Do not rename, do not
  move. `rcomplex:::` in tests. Delete README index rows.
- Accept: `R CMD build . && R CMD check --no-manual rcomplex_*.tar.gz`
  OK; `lintr::lint_package()` clean; surface snapshot updated to 45.
- Deps: WP0.

### WP2 Remove Tier C, non-modules files (parallel with WP1, WP3)

- Files: delete `R/module_auroc.R`, `R/subspace.R`,
  `R/recurrence_graph.R`, `R/clique_patterns.R`,
  `R/coexpressolog-strength.R`, `R/tag_permutation.R`, `R/tag_blocks.R`,
  `src/module_auroc.cpp`, `src/subspace.cpp`; cut `clique_persistence`,
  `clique_perturbation_test`, `clique_intensity_test` and their helpers
  from `R/cliques.R` (`clique_threshold_sweep` and `clique_stability`
  stay; `find_cliques_stability.cpp` is `clique_stability`'s kernel,
  keep it); cut the matching `.rcomplex` methods from
  `R/rcomplex-class.R`; delete `tests/testthat/test-module-auroc*.R`,
  `test-subspace.R`, `test-recurrence-graph.R`, `test-tag-permutation.R`,
  `test-intensity-test.R`, `test-perturbation.R`, `test-restricted-nulls.R`
  (check: keep any case that covers Tier A); strip Tier C rows from
  `helper-rng-contract.R` table; README and tutorial sections
  ("Secondary: do the same orthogroups recur", "Clique robustness
  diagnostics"); methods.Rmd paragraphs; CLAUDE.md lines naming them.
- Do: `git mv` R sources to `dev/probes/` first, then delete from
  package. `Rcpp::compileAttributes()` after C++ removal.
- Accept: check OK; `grep -rE 'module_auroc|subspace_preservation|
  recurrence_graph|clique_intensity|clique_perturbation|clique_persistence|
  tag_permutation|coexpressolog_strength' R/ src/ tests/ README.md
  vignettes/` returns nothing; `wc -l R/*.R` total < 13,500; surface
  snapshot updated.
- Deps: WP0, gate G1. Do not touch `R/modules.R` (WP3 owns it).

### WP3 detect_modules diet (parallel with WP1, WP2)

- Files: `R/modules.R`, `src/coclassification.cpp`
  (`sparse_excess_spectral_norm_cpp` only), `R/rng.R`
  (`.blas_fork_safe()` if the K = 1 fork was its only user), `DESCRIPTION`
  (drop `sbm` from Suggests), `tests/testthat/test-modules.R`,
  `test-module-determinism.R`, `test-fork-safety.R`, tutorial section
  "Testing for community structure (K = 1 null)", CLAUDE.md gotcha
  paragraphs that only describe the K = 1 fork site.
- Do: keep `method = "leiden"` with the multi-resolution consensus; delete
  `infomap`, `sbm`, `test_community_structure()` and its C++ kernel,
  `characterize_hubs()` and its `.rcomplex` method. `detect_modules()`
  loses the `method` argument.
- Accept: check OK; `R/modules.R` < 900 lines; the `n_cores = 2`
  determinism test still passes; surface snapshot updated.
- Deps: WP0, gate G2.

### WP4 One name style (serial after WP1-3 merge)

- Files: every `R/*.R`, `src/*.cpp` that builds a data.frame, tests,
  README, vignettes, `inst/extdata/orthologs_small.txt` header.
- Do: columns `species1 species2 gene1 gene2 hog p_value q_value
  effect_size power`; keep `Zsummary_std` (WGCNA term). Arguments
  already `n_cores seed alpha`; check `n_perm` vs `n_perm_pres` and keep
  both only where two nulls really exist. No deprecation shims (pre-1.0,
  0.4.0 is the break).
- Accept: `grep -rE '"(Species[12]|sp[12]|q\.value|p\.val[a-z.]*)"' R/
  src/ tests/` empty; check OK.
- Deps: WP1-3 merged, gate G3.

### WP5 Argument diet (serial after WP4)

- Files: `R/module_preservation.R`, `R/cliques.R`
  (`clique_threshold_sweep` 14), `R/modules.R`, `R/comparison.R`,
  `R/preservation_matrix.R`, `R/coexpressolog_null.R`,
  `R/ortholog_map.R` (`resolve_ortholog_map` 9), tests.
- Do: for each formal of a Tier A export that neither README, tutorial,
  nor the gitignored root vignette (`prepare_data/vignettes/
  root-workflow.Rmd`, Martin greps) ever sets, delete it and hardcode
  the default. No `control = list()`. Keep formals CLAUDE.md names as
  design decisions (`calibrate`, `pval_combine`, `rho0`, `min_power`,
  `swap_factor`, `store_density`, `block_size`).
- Accept: surface test: no export has > 8 named formals; check OK.
- Deps: WP4.

### WP13 Scores and E-values (parallel with WP14)

- Files: `R/comparison.R` (`comparison_to_edges()`), `R/specificity.R`
  if the rank path builds its frame elsewhere, `R/module_preservation.R`
  (one line), tests, `man/`.
- Why: BLAST reports a bit score (evidence, comparable across searches)
  and an E-value (expected chance hits in a database this size). A
  user reads E < 1e-3 without a statistics course. rcomplex reports
  `p_value`, `q_value`, `effect_size`, `power`; the ingredients are
  there, the two BLAST numbers are not.
- Do: two columns on every edge table, same definitions on every path
  (hypergeometric, rank, permutation), computed on the combined
  p-value (`pval_combine`, so `max` by default):
  `score = -log2(p_value)` in bits, the evidence against chance;
  `evalue = n_tests * p_value`, where `n_tests` is the number of
  ortholog pairs tested in that species-pair comparison (the rows of
  the direction with more tests; both directions share the pair set
  up to genes absent from a network). This is BLAST's own relation,
  `E = N * 2^(-S)`, so a user who knows BLAST knows these. `effect_size`
  stays as the "percent identity" of the hit: magnitude, not evidence.
  `score` is not `log2(effect_size)` on purpose: a fold enrichment of
  8 on a 3-gene neighbourhood is no evidence, and the hypergeometric p
  already weighs magnitude against neighbourhood size. Permutation
  p-values floor at `1 / (n_perm + 1)`, so `score` floors with them;
  `pvalue_resolution()` already says so. `module_preservation()` gains
  the same `evalue` on its calibrated p with `n_tests` = modules
  tested. No new argument anywhere; `n_tests` is a column, not a
  parameter. `print.rcomplex()` (WP6) and the quickstart show
  `score` and `evalue` first, then `q_value`, `effect_size`, `power`.
  Cliques: both clique tables (`gene_clique_graph()`, `find_cliques()`)
  gain `score` = the sum of the member edges' `score`. Log-odds add,
  so this is the NetworkBLAST log-likelihood-ratio of a conserved
  subnetwork and the BLAST bit score of a longer alignment: it grows
  with clique size on purpose. `mean_q` stays for the classifiers. No
  clique `evalue`: the clique null is not analytic; `rcomplex(null =
  TRUE)` (WP12) counts null-network cliques at or above each score
  instead and `summary()` reports that count.
- Accept: `evalue == n_tests * p_value` and `score == -log2(p_value)`
  on the fixture edges for all three methods; `sum(evalue < 1)` on a
  `null_network()` comparison is at most about 1 (one chance hit
  expected); clique `score` equals the sum over its `n_edges` member
  rows on the clique fixture; check OK.
- Deps: WP5.

### WP14 Signed networks: anticorrelation (parallel with WP13)

- Files: `R/network.R`, `R/network-sparse.R` (`.net_cpp_args()`),
  `src/mutual_rank.cpp`, `src/network_block.cpp`, `R/null_network.R`
  and `R/comparison.R` only where `abs_cor` is forwarded, tests,
  `man/`, `vignettes/articles/methods.Rmd` (one paragraph).
- Why: today a negative correlation ranks last and never enters a
  network (`abs_cor = FALSE`, the default), or loses its sign
  (`abs_cor = TRUE`). Neither can report conserved anticorrelation
  (the same repressive partners in both species) or a sign flip (the
  same partners, positive in one species and negative in the other),
  which is rewiring with a mechanism and is what the proposal's
  "regulatory rewiring" means at the co-expression level.
- Do: `compute_network(sign = c("positive", "negative", "unsigned"))`
  replaces `abs_cor`. The kernels already apply `abs_cor` in one line
  per column (`col[r] = abs_cor ? fabs(v) : v` in `mutual_rank.cpp`
  and `network_block.cpp`); it becomes a three-way: negate for
  `negative`, `fabs` for `unsigned`, identity for `positive`. MR,
  density threshold, sparse store, blockwise validity rule and every
  consumer are unchanged, because MR ranks whatever it is given. `sign`
  joins `params` and `.net_cpp_args()` so `null_network()` and
  `density_sweep()` rebuild with the same sign. Nothing else changes:
  `find_coexpressologs(list(A = net_A_pos, B = net_B_neg))` already
  tests whether A's positive partners are B's negative partners,
  because the test reads membership only. The driver (WP6) takes the
  same `sign` argument plus `"both"`: builds positive and negative
  networks per species, runs the four comparisons, and the edge table
  gains a `sign` column with values `+` (positive in both), `-`
  (negative in both), `flip` (positive in one, negative in the other,
  either way). Gene-graph cliques keep the column and the classifiers
  ignore it; `sign` is a filter for the user, not a tier.
- Caveats to write into `@details`, not into code: negative
  correlations are rarer and weaker in RNA-seq, a top-3 % negative
  network exists at any n, and WP12's `r_threshold` line is what tells
  the user whether it holds anything (r >= -0.3 at n = 20 is noise).
  No default changes: `sign = "positive"`.
- Accept: `compute_network(x, sign = "negative")$network` equals
  `compute_network(-x_cor_proxy)`'s where the proxy flips half the
  genes (`x[flip, ] <- -x[flip, ]` makes those pairs anticorrelated:
  the negative network of `x` must contain exactly the flipped pairs
  the positive network of the proxy contains); `sign = "unsigned"`
  reproduces the old `abs_cor = TRUE` network bit for bit on the
  fixture; blockwise and dense agree for every sign; `grep abs_cor R/
  src/ tests/` empty; check OK.
- Deps: WP11 (same files; WP11 merges first).

### WP6 Driver (serial, after WP7, WP10, WP11, WP13, WP14)

- Files: `R/rcomplex-class.R` (rewrite, target < 300 lines),
  `tests/testthat/test-rcomplex-class.R` (rewrite), `man/`.
- Do: `rcomplex()` becomes the BLAST verb.

  ```r
  rcomplex(expr = NULL, orthologs, networks = NULL, block = NULL,
           clades = NULL, density = 0.03,
           sign = c("positive", "negative", "unsigned", "both"),
           method = c("hypergeometric", "rank"), alpha = 0.1,
           modules = FALSE, null = FALSE, n_cores = 1L, seed = NULL)
  ```

  Twelve formals: over the WP5 budget of eight by design, the driver
  is the one place the knobs meet; the surface test (WP9) exempts it
  by name.

  `expr`: named list (species) of matrices or SummarizedExperiments.
  `networks`: named list of network objects instead, from
  `compute_network()` or `as_network()` (WP10); exactly one of `expr`
  and `networks`. `orthologs`: long data.frame (`species gene hog`) or
  a file path for `read_orthologs()`. `block`: named list of per-sample
  factors, one per species in `expr`; each species goes through
  `split_layers()` and its wiring layer is what `compute_network()`
  sees (WP11); `null_network(block = )` for any null the driver runs.
  Runs `compute_network()` per species, `find_coexpressologs()` over
  all pairs, `gene_clique_graph()` + `classify_gene_cliques()`; with
  `modules = TRUE` adds `detect_modules()` + `preservation_paired()`
  over all pairs. Returns a plain list of class `rcomplex`: `networks`,
  `edges`, `cliques`, `classification`, optional `modules`,
  `preservation`, plus `call`. Methods: `print()` (species, genes,
  samples, edges, tier counts), `summary()` = classification table,
  `as.data.frame()` = classification, `write_rcomplex(x, dir)` writes
  one TSV per element. Delete every existing `.rcomplex` S3 method; the
  driver calls the defaults. Seed once, pass `seed = NULL` down (RNG
  contract); add `rcomplex` to the rng-contract table.
- Accept: `test-rcomplex-class.R` runs the driver on `inst/extdata` in
  < 30 s and checks tier counts, once from `expr` and once from
  `networks = lapply(expr, compute_network)` with identical edges;
  `print()` snapshot; rng-contract test passes; check OK.
- Deps: WP7, WP10, WP11, WP13, WP14, gate G4.

### WP7 Nested clades, MDO I (parallel with WP10, WP11)

- Files: `R/clique_gene_graph.R`, `R/cliques.R` (`classify_cliques`,
  `clique_stability`), `R/modules.R` (`classify_hub_conservation`),
  `R/preservation_matrix.R` (`block`), `R/orthologs.R` or new
  `R/clades.R` (< 60 lines), tests, `DESCRIPTION` (Suggests `ape`).
- Do: `species_trait` (flat named vector) becomes `clades`: a named list
  of species vectors, nested allowed, must be laminar (every pair nested
  or disjoint; validate once in `.check_clades()`). A flat trait is
  `split(names(trait), trait)`. `lineage_specific`, `trait_specific`,
  `differentiated` and the `underpowered` rule are evaluated against
  the smallest clade that contains the clique's species; output gains a
  `clade` column (name of that clade). `clades_from_tree(phy, min_size
  = 2L)` turns an `ape::phylo` into that list via `ape::prop.part()`;
  no other tree code. `preservation_matrix_test(block = )` accepts the
  same list (first level only). Rename, no alias.
- Accept: new test with a 3-level nested clade list on the clique
  fixtures (a clique lineage-specific to an inner clade is reported with
  that clade, not the outer); every existing trait test passes after
  `clades = split(...)`; `clades_from_tree()` round-trips a 4-tip tree;
  check OK with and without `ape` installed.
- Deps: WP5, gate G5.

### WP10 Generalise inputs (parallel with WP7, WP11)

- Files: `R/orthologs.R`, `R/network-sparse.R` or new `R/as_network.R`
  (< 120 lines), `inst/extdata/` (one OrthoFinder `N0.tsv` fixture of
  ~20 rows), tests, `man/`.
- Do: today `parse_orthologs()` reads PLAZA only, two species at a time,
  and nothing imports a network built elsewhere. Two additions.
  (a) `read_orthologs(file, species = NULL, format = c("auto",
  "orthofinder", "plaza", "long"))` returns the long table `species
  gene hog` for every species in the file (or the `species` subset).
  Formats: OrthoFinder `N0.tsv` / `Orthogroups.tsv` (one column per
  species, comma-separated genes), PLAZA (today's parser), long
  (already the output shape; passthrough with column check). `auto`
  reads the header: `HOG`/`Orthogroup` column means OrthoFinder,
  `gene_content` means PLAZA, `species gene hog` means long. FastOMA
  and eggNOG users convert to long themselves; one documented sentence,
  no parser. `parse_orthologs()` is deleted; `prepare_orthologs()`
  takes the long table where it took SE rowData.
  (b) `as_network(x, density = 0.03, genes = NULL)`: `x` is a symmetric
  numeric matrix (dense or `dgCMatrix`, genes in dimnames) or an edge
  list `gene1 gene2 weight`. Builds the same network object
  `compute_network()` returns (`network` as `dgCMatrix`, `threshold`
  at the `density` quantile of the stored weights, `store_threshold`,
  `params`), so every consumer works unchanged through `.net_check()`.
  Weights are taken as given; a TEA-GCN `zScore(Co-exp_Str_MR)` or a
  WGCNA adjacency both qualify. No MR recomputation: that is
  `compute_network()`'s job on expression.
- Accept: `read_orthologs()` on the N0 fixture and on
  `orthologs_small.txt` give the same long table for the two shared
  species; `as_network()` on `compute_network(x)$network` round-trips to
  identical `find_coexpressologs()` edges; `find_cliques()` and
  `module_preservation()` run on `as_network()` input; check OK.
- Deps: WP5.

### WP11 Harmonise heterogeneous designs (parallel with WP7, WP10)

- Files: `R/split_layers.R`, `R/network.R`, `R/null_network.R`, tests,
  `vignettes/articles/methods.Rmd` (one section).
- Why: species rarely share a design. Pooled leaf + root networks are
  tissue networks (memory: tissue confound), and a 10-tissue compendium
  against a 3-tissue one aligns tissue coverage, not regulation. TEA-GCN
  (Lim et al. 2026, Nat Commun 17:5906) solves the public-compendium
  version: k-means partitions of samples in PCA space, correlation per
  partition, rectified average (negatives to zero, then mean) across
  partitions, then MR and a global z-score so networks compare across
  species. rcomplex already has the small-n version: `split_layers()`
  projects out a designed block and the wiring layer is correlation
  within blocks pooled over them; MR + density threshold already make
  networks comparable (TEA-GCN's z-score is the same idea on raw
  weights). Three steps close the gap without reimplementing TEA-GCN.
- Do:
  (a) `split_layers()` keeps its name and gains nothing; it is the
  documented answer to "my species have different conditions" and the
  driver's `block` argument (WP6). Its `@description` becomes three STE
  sentences; the Breschi/Cote/Parsana rationale moves to methods.Rmd.
  (b) `compute_network(partition = NULL)`: a per-sample factor. When
  given, correlation is computed within each level and combined by
  rectified average before MR (dense path only; `block_size` with
  `partition` errors with one sentence). Levels with fewer than
  `min_partition_n = 5` samples are dropped with a message. This is
  TEA-GCN's partition aggregation on a user-supplied partition; rcomplex
  does not cluster samples, the user brings the factor (tissue, study,
  or a k-means they ran). One kernel call per level, existing
  correlation code, ~40 lines.
  (c) `harmonise_blocks(block)`: no. The driver checks that every
  species' `block` levels are a subset of the union and messages which
  species lack which level; nothing is dropped automatically. Ten
  lines inside WP6, no new export.
- Accept: `compute_network(x, partition = f)` on a 2-level factor equals
  the hand-computed `pmax(cor_1, 0) + pmax(cor_2, 0)) / 2` passed
  through MR; `split_layers()` tests unchanged; `null_network(block =)`
  still matches the wiring layer; check OK.
- Deps: WP5.

### WP12 Sample-size honesty (serial after WP6)

- Files: `R/rcomplex-class.R` (driver, `print`, `summary`),
  `tests/testthat/test-rcomplex-class.R`, `tests/testthat/
  helper-preservation.R` (`pres_expr()` already takes `n_samp`).
- Do: three small things, no new export.
  (a) `print.rcomplex()` adds one line per species: `n_genes`,
  `n_samples`, `density`, and `r_threshold`, the smallest correlation
  that passed the density threshold (read it from the network object;
  if `compute_network()` does not keep it, store it in `params`, one
  field). A user sees "20 samples, r >= 0.71" next to "600 samples,
  r >= 0.18" and knows what each network can resolve.
  (b) `rcomplex(null = FALSE)`: when `TRUE`, the driver runs
  `null_network(block = block)` per species and `find_coexpressologs()`
  on the null networks with the same settings, and `summary()` reports
  `calls_null` beside `calls` per species pair, plus their ratio as the
  empirical false-call rate at `alpha`. Opt-in because it doubles the
  run time; `summary()` prints `null: not run` otherwise. Same seed
  contract as the rest of the driver.
  (c) Test across regimes: one generator (`pres_expr()`), three sample
  sizes (6, 20, 200), same genes, same orthologs. The driver runs at
  every n with identical output columns; with `null = TRUE` the
  false-call rate is at most `alpha` at every n; power (calls among
  planted conserved pairs) is non-decreasing in n. Keep it under 20 s:
  200 genes, `density = 0.05`.
- Accept: the regime test passes; `print()` snapshot shows the per
  species line; check OK.
- Deps: WP6.

### WP8 Docs (serial after WP12)

- Files: `README.md`, `vignettes/quickstart.Rmd` (new),
  `vignettes/rcomplex-tutorial.Rmd` -> `vignettes/articles/walkthrough.Rmd`,
  `vignettes/articles/methods.Rmd`, `_pkgdown.yml` (new), `CLAUDE.md`,
  every `@description` and `@param` in `R/`, every `stop()` /
  `warning()` / `message()` string.
- Do: README <= 150 lines: one-paragraph pitch, a mermaid pipeline
  diagram (GitHub renders it, no tooling), install, the driver on
  `inst/extdata` in <= 12 lines, output columns, three pitfalls (few
  samples, density, p-value floor), citation. Function index moves to
  `_pkgdown.yml` reference groups: Run, Inputs, Harmonise, Networks,
  Co-expressologs, Cliques, Modules, Traits, Nulls and diagnostics.
  Quickstart <= 150 lines, driver only. Walkthrough = current tutorial
  minus Tier C and K = 1 sections, <= 800 lines, uses the building
  blocks. methods.Rmd: drop Tier C paragraphs, gain the WP11 section.
  CLAUDE.md: rewrite the overview to the new surface, drop Tier C design
  notes, add the surface budget and the prose rules.
  Prose rules (section 7, STE-lite) apply to every file in this WP:
  `@description` <= 3 sentences; rationale and literature go to
  `@details` or methods.Rmd; one term per concept (the WP4 column names
  are the vocabulary: species, gene, hog, edge, clique, module, block,
  clade); every error message names the argument, says what it got,
  says what it needs, one action per sentence.
- Accept: `wc -l README.md` <= 150; both vignettes knit;
  `pkgdown::build_site()` OK; every Tier A export appears in exactly one
  reference group; `Rscript dev/sharpen-prose.R` (WP8 writes it: flags
  `@description` blocks over 3 sentences or any sentence over 25 words
  in `R/*.R`, `README.md`, `vignettes/quickstart.Rmd`) reports zero.
- Deps: WP12.

### WP9 Surface guard (last)

- Files: `tests/testthat/test-surface.R`.
- Do: replace the snapshot with budgets: exports <= 30; named formals
  <= 8 per export except `rcomplex()`; `R/*.R` <= 1,500 lines each;
  `README.md` <= 150 lines. Snapshot of the export list stays so
  additions are explicit.
- Accept: test passes on `refactor/sharpen`; break any budget locally
  and watch it fail.
- Deps: WP8.

### WP15 Joint modules across species (gated G9, after 0.4.0)

- Files: new `R/joint_modules.R` (< 150 lines), tests, `man/`,
  `dev/probes/joint_modules_pooideae.R` (the Orion probe script).
- Why: rcomplex detects modules per species, then tests preservation.
  Two frameworks find modules that span species in one step. The
  multilayer view (Mucha et al. 2010): species are layers, orthology
  couples layers, one modularity over the supra-graph. The hypergraph
  view (Kaminski et al. 2019; Chodrow, Veldt & Benson 2021): each HOG
  is one hyperedge over all its copies, and hypergraph modularity with
  a partial affinity rule lets a HOG spread over modules, so the
  partition itself says which copies carry the conserved role. The
  two agree on one-copy HOGs and differ on multi-copy HOGs: pairwise
  coupling pulls every copy pair with the same force and big HOGs
  dominate; one hyperedge per HOG with a size-aware null does not.
  Both are a Leiden run on one graph, so one function probes both.
- Do: `joint_modules(nets, orthologs, coupling = c("star", "pairwise"),
  weight = NULL, objective = c("modularity", "CPM"), resolution =
  NULL, seed = NULL)`. Build one igraph: nodes `(species, gene)` for
  every gene in every network; intra-species edges from each sparse
  store, divided by their maximum as `detect_modules()` does for CPM;
  inter-species links by `coupling`: `star` adds one node per HOG and
  one edge of weight `weight` from it to each copy (hypergraph star
  expansion, the 50-line approximation of hypergraph modularity);
  `pairwise` adds an edge of weight `weight` between every ortholog
  copy pair across species (multislice coupling). `weight = NULL` is
  the median intra-species edge weight. One `igraph::cluster_leiden()`
  on it, HOG nodes dropped from the result. Return: per-species
  partitions in `as_modules()` shape plus a HOG table: `hog`,
  `n_copies`, `n_modules` (how many modules the HOG spreads over),
  `module_main`, `copies_main` (the copies in it), so a HOG split
  across modules is readable as a subfunctionalisation candidate. No
  consensus sweep inside; the caller sweeps `weight` and feeds
  `detect_modules()`'s consensus if wanted.
- Accept: on `pres_fixture()` (planted modules, two species) both
  couplings recover the planted modules with ARI > 0.9 per species; a
  HOG with one planted non-conserved copy puts that copy outside
  `module_main`; `pairwise` on a 5 x 5 HOG fixture merges the copies
  where `star` does not (the documented difference); rng-contract
  table gains `joint_modules`; check OK. The Orion probe (not a test)
  is split-half and cross-species replication on leaf and wood at
  three weights, against the existing two-step path.
- Deps: WP3, WP7, WP10 (long ortholog table). Blocks nothing.

### WP16 Edge gain and loss on the species tree (gated G10, after 0.4.0)

- Files: `R/clades.R` (WP7 owns it; this WP adds one function) or new
  `R/edge_history.R` (< 120 lines), tests, `DESCRIPTION` (Suggests
  `phangorn`), `man/`.
- Why: `clades` (WP7) is flat. The proposal's "phylogeny-aware
  inference" means states along a tree. A clique's membership over
  species is a binary character; parsimony on the species tree
  reconstructs where co-expression was gained and lost. This is the
  formal version of `lineage_specific`: a clade-restricted clique is
  one gain on that clade's stem, or one loss outside it, and the tree
  says which is cheaper.
- Do: `edge_history(cliques, tree, species)`: for each clique row of
  `classify_gene_cliques()` or `classify_cliques()`, build the
  character from `species_present` (member 1, tested-and-absent 0,
  `untested` / `underpowered` as `?`), run Fitch parsimony with
  ambiguity (`phangorn::ancestral.pars`, `type = "ACCTRAN"`), and
  return per clique: `state_root`, `n_gain`, `n_loss`, `branches`
  (node labels where the state changes), `n_mpr` (equally
  parsimonious reconstructions; ties are reported, not hidden). A
  second table per branch: `n_gain`, `n_loss`, so a branch with many
  losses names a lineage where co-expression diverged. `tree` is an
  `ape::phylo`, tips named by species, same object `clades_from_tree()`
  reads. No ML, no Dollo: Fitch is the one rule, `n_mpr` is the honest
  uncertainty.
- Accept: 4-tip balanced tree, clique present in one clade of two
  gives `n_gain = 1` on that stem and `n_loss = 0` under ACCTRAN;
  present in three of four tips gives one loss on the missing tip's
  branch; a `?` tip never counts as a change; check OK with and
  without `phangorn`.
- Deps: WP7. Blocks nothing.

## 4. Gates (Martin decides, before the WP starts)

| Gate | WP | Question | Default if silent |
|---|---|---|---|
| G1 | WP2 | Remove all 17 Tier C functions? Any live Orion script or the root vignette sourcing one? `tag_permutation` in or out? | remove all, `dev/probes/` keeps source |
| G2 | WP3 | Drop Infomap, SBM, K = 1 test from `detect_modules()`? | drop |
| G3 | WP4 | snake_case columns, no deprecation shims, 0.4.0 break? | yes |
| G4 | WP6 | Driver returns plain list; container S3 methods deleted? | yes |
| G5 | WP7 | `species_trait` renamed to `clades` with no alias? | yes |
| G6 | WP10 | `parse_orthologs()` deleted in favour of `read_orthologs()`? Long table `species gene hog` as the one ortholog shape? | yes |
| G7 | WP13 | `score = -log2(p)` bits and `evalue = n_tests * p` as the two headline columns, `effect_size` kept as magnitude? | yes |
| G8 | WP14 | `abs_cor` replaced by `sign`; driver `sign = "both"` runs four comparisons and tags `+ - flip`? | yes |
| G9 | WP15 | Build the joint-module probe (both couplings)? Adopt into the package only if the Orion probe beats `detect_modules()` + `module_preservation()` on replication | build the probe |
| G10 | WP16 | Edge gain/loss on a species tree, `phangorn` in Suggests? | yes |

Out of scope, on purpose: rank-vs-hypergeometric default (needs Orion
validation, design note 11.15); Bioconductor conventions (MDO V);
plotting and enrichment (MDO IV); regulatory layers (MDO III); PR #4;
sample clustering for partitions (TEA-GCN's k-means: the user brings
the factor, or builds the network in TEA-GCN and imports it with
`as_network()`); global network alignment (IsoRank-style topology
matching without orthologs): rcomplex is ortholog-anchored by design.

## 5. Agents

Budget rule (Martin, 2026-10-07: Fable quota is short this week). Fable
is the parent only: it answers gates, launches, reads reports, merges.
No WP runs on Fable. Model routing follows the Claude Code subagent
docs ("route high-volume tasks to cheaper models"): census on haiku,
mechanical cut / rename / docs / review on sonnet, new code on opus.

Five project subagents in `.claude/agents/` (checked in; `.claude` is
Rbuildignored). Per the docs: short `description` with the trigger,
everything else in the body; `tools` allowlisted; `skills` preloaded so
agents do not spend turns discovering them; `permissionMode:
acceptEdits` for writers; `maxTurns` as a ceiling, not a target. One
definition per kind of work; the WP number travels in the launch
prompt, the plan file is the source of truth.

| File | model | tools | WPs |
|---|---|---|---|
| `sharpen-census.md` | haiku | Read, Grep, Glob, Bash | pre-WP2 and after every merge: `dev/sharpen-census.R`, grep `prepare_data/` for Tier C names |
| `sharpen-cutter.md` | sonnet | Read, Edit, Write, Grep, Glob, Bash | WP0, WP1, WP2, WP3, WP4, WP5, WP9: work defined by a list, no design |
| `sharpen-builder.md` | opus | Read, Edit, Write, Grep, Glob, Bash | WP6, WP7, WP10-WP16: new functions with a stated signature |
| `sharpen-docs.md` | sonnet | Read, Edit, Write, Grep, Glob, Bash | WP8 |
| `sharpen-reviewer.md` | sonnet | Read, Grep, Glob, Bash | every PR before merge, read-only |

Launch prompt shape (caveman, under 15 lines):

```
WP<n> of dev/design-notes/sharpen-plan.md section 3. Gate G<x>: <answer>.
Branch: worktree off refactor/sharpen. Read the WP, then the files it
lists, then start. Stop and report if a file outside the list must
change.
```

Writers run with `isolation: "worktree"` and must commit before
finishing (memory: worktree commits). Parallel sets go in one Agent
call. Worktree agents that share a file must not run together: WP2 and
WP3 both touch `R/rcomplex-class.R` only through method deletion for
disjoint functions, so WP2 deletes the Tier C methods and WP3 leaves
`rcomplex-class.R` alone except `characterize_hubs`. WP7, WP10, WP11
are disjoint by file (`clique_*`/`cliques.R`/`preservation_matrix.R`;
`orthologs.R`/`as_network.R`; `split_layers.R`/`network.R`/
`null_network.R`). WP13 (`comparison.R`) and WP14 (`network.R`,
kernels) are disjoint; WP14 waits for WP11's `network.R` merge.

Fable cost per WP: one launch, one report read, one reviewer report
read, one merge. Reports are capped at 30 lines in the agent bodies so
the parent context stays small.

## 6. Expected end state

29 exports, ~12k R lines, ~5.3k C++, ~15k test lines, README 150
lines, a 10-line quickstart that is the whole BLAST analogue:

```r
library(rcomplex)
expr <- list(BDIS = bdis_vst, HVUL = hvul_vst)          # genes x samples
block <- list(BDIS = bdis_meta$tissue, HVUL = hvul_meta$tissue)
res <- rcomplex(expr, "N0.tsv", block = block,
                clades = list(annual = "BDIS", perennial = "HVUL"))
res
#> rcomplex: 2 species, 1,842 hogs, 3 tiers
#>   BDIS  18,201 genes  20 samples  density 0.03  r >= 0.71
#>   HVUL  19,877 genes  20 samples  density 0.03  r >= 0.69
#>   edges 4,112   cliques 1,203   complete_conserved 611 ...
summary(res)       # classification table; with null = TRUE, calls_null
head(res$edges)    # gene1 gene2 hog sign score evalue q_value effect_size power
write_rcomplex(res, "out/")
```

or, from networks built elsewhere:

```r
nets <- list(ATH = as_network(tea_gcn_ath), OSA = as_network(tea_gcn_osa))
res <- rcomplex(networks = nets, orthologs = "N0.tsv")
```

MDO I and II then read: multi-species (any N, nested clades, species
tree in), any ortholog source, any network source, designed-block
harmonisation and partition aggregation for heterogeneous designs,
network-aware nulls (edge-swap, shuffled expression within block, rank
test), effect sizes and power on every edge, p-value resolution on every
test. Nothing in the proposal's P1 that rcomplex does not do except the
rank-test default, which waits on Orion.

## 7. Prose: Karpathy's ASD-STE100 note, audited

Karpathy (X, 2026-10-02): ask the model to write in ASD-STE100,
Simplified Technical English, the aerospace maintenance standard; "80 %
STE" is enough. Rules that matter: sentences under 20 words, one action
per sentence, active voice, one term per concept. His other tips
(diagrams, HTML pages, 3b1b-style video) are about explaining model
output, not package text.

Audit against rcomplex. The package prose fails three of the four
rules. `@description` blocks run to 40-line essays with literature
(`split_layers`, `compute_network`); CLAUDE.md sentences average well
over 20 words; the same concept has several names (`Species1`,
`species1`, `sp1`; `trait`, `lineage`, `clade`; `block`, `layer`,
`partition`). Error messages are mostly fine already (they name the
argument).

Decision. Adopt the four rules as STE-lite for everything a user reads:
`@description`, `@param`, `@return`, README, quickstart, messages.
Not the 900-word dictionary (the domain vocabulary is outside it) and
not `@details`, methods.Rmd, or design notes, which are allowed to
argue. The vocabulary is the WP4 column names plus `block` (a designed
sample factor), `partition` (a factor for per-level correlation),
`clade` (a species set). `layer` is retired from user text except as
`split_layers()`'s own output names. Enforcement: `dev/sharpen-prose.R`
in WP8's acceptance, and the reviewer's checklist. Diagrams: one
mermaid pipeline figure in README; no HTML, no video. For agents: no
change, caveman is already shorter than STE and agents read the plan,
not the prose.
