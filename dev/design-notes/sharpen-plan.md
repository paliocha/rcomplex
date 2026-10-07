# Sharpen plan: MDO I + II, BLAST-shaped UX

2026-10-07. Baseline `main` b245ac0 (0.3.2.9000). Target release 0.4.0.

Goal. Hvidsten proposal MDO I (multi-species, phylogeny-aware) and II
(nulls, effect sizes) are already covered by rcomplex. What is missing is
the BLAST/Clustal shape: one verb, data in, table out, defaults good.
Package is 57 exports, 18.7k R lines, 1.5k-line tutorial. Most of it is
research scaffolding from the 2026-09/10 Orion benchmarks. Cut it, keep
the math, add one driver and nested clades.

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

**Tier A, stays exported (27).** The user path plus the building blocks a
power user needs by name.

```
compute_network  parse_orthologs  prepare_orthologs  reduce_orthogroups
find_coexpressologs  density_sweep  null_network  coexpressolog_null
gene_clique_graph  classify_gene_cliques
find_cliques  clique_stability  clique_threshold_sweep  classify_cliques
as_modules  detect_modules  resolve_ortholog_map  module_preservation
classify_preservation  module_correspondence  preservation_paired
preservation_matrix_test  identify_module_hubs  classify_hub_conservation
get_coexpressed_hogs  pvalue_resolution  rcomplex (driver, WP6)
```

**Tier B, demote to internal (12).** Called by Tier A, never needed by
name. Keep the function, drop `@export`, add `@keywords internal`.
Tests call them as `rcomplex:::fn()`.

```
compare_neighborhoods  summarize_comparison  comparison_to_edges
permutation_hog_test  compare_specificity  summarize_specificity
run_pairwise_comparisons  all_species_pairs  mr_block  as_sparse_network
extract_orthologs  split_layers
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

Order: WP0 -> {WP1 || WP2 || WP3} -> WP4 -> WP5 -> {WP6 || WP7} -> WP8
-> WP9.

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
  `R/split_layers.R`, `R/preservation_matrix.R` (`all_species_pairs`),
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

### WP6 Driver (parallel with WP7)

- Files: `R/rcomplex-class.R` (rewrite, target < 250 lines),
  `tests/testthat/test-rcomplex-class.R` (rewrite), `man/`.
- Do: `rcomplex()` becomes the BLAST verb.

  ```r
  rcomplex(expr, orthologs, clades = NULL, density = 0.03,
           method = c("hypergeometric", "rank"), alpha = 0.1,
           modules = FALSE, n_cores = 1L, seed = NULL)
  ```

  `expr`: named list (species) of matrices or SummarizedExperiments.
  `orthologs`: data.frame or file path (via `parse_orthologs()`).
  Runs `compute_network()` per species, `find_coexpressologs()` over all
  pairs, `gene_clique_graph()` + `classify_gene_cliques()`; with
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
  < 30 s and checks tier counts; `print()` snapshot; rng-contract test
  passes; check OK.
- Deps: WP5, gate G4.

### WP7 Nested clades, MDO I (parallel with WP6)

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

### WP8 Docs (serial after WP6 + WP7)

- Files: `README.md`, `vignettes/quickstart.Rmd` (new),
  `vignettes/rcomplex-tutorial.Rmd` -> `vignettes/articles/walkthrough.Rmd`,
  `vignettes/articles/methods.Rmd`, `_pkgdown.yml` (new), `CLAUDE.md`.
- Do: README <= 150 lines: one-paragraph pitch, install, the driver on
  `inst/extdata` in <= 12 lines, output columns, three pitfalls (few
  samples, density, p-value floor), citation. Function index moves to
  `_pkgdown.yml` reference groups: Run, Networks, Co-expressologs,
  Cliques, Modules, Traits, Nulls and diagnostics. Quickstart <= 150
  lines, driver only. Walkthrough = current tutorial minus Tier C and
  K = 1 sections, <= 800 lines, uses the building blocks. methods.Rmd:
  drop Tier C paragraphs. CLAUDE.md: rewrite the overview to the new
  surface, drop Tier C design notes, add the surface budget.
- Accept: `wc -l README.md` <= 150; both vignettes knit;
  `pkgdown::build_site()` OK; every Tier A export appears in exactly one
  reference group.
- Deps: WP6, WP7.

### WP9 Surface guard (last)

- Files: `tests/testthat/test-surface.R`.
- Do: replace the snapshot with budgets: exports <= 30; named formals
  <= 8 per export; `R/*.R` <= 1,500 lines each; `README.md` <= 150
  lines. Snapshot of the export list stays so additions are explicit.
- Accept: test passes on `refactor/sharpen`; break any budget locally
  and watch it fail.
- Deps: WP8.

## 4. Gates (Martin decides, before the WP starts)

| Gate | WP | Question | Default if silent |
|---|---|---|---|
| G1 | WP2 | Remove all 17 Tier C functions? Any live Orion script or the root vignette sourcing one? `tag_permutation` in or out? | remove all, `dev/probes/` keeps source |
| G2 | WP3 | Drop Infomap, SBM, K = 1 test from `detect_modules()`? | drop |
| G3 | WP4 | snake_case columns, no deprecation shims, 0.4.0 break? | yes |
| G4 | WP6 | Driver returns plain list; container S3 methods deleted? | yes |
| G5 | WP7 | `species_trait` renamed to `clades` with no alias? | yes |

Out of scope, on purpose: rank-vs-hypergeometric default (needs Orion
validation, design note 11.15); Bioconductor conventions (MDO V);
plotting and enrichment (MDO IV); regulatory layers (MDO III); PR #4.

## 5. Agents

One agent per WP, `general-purpose` in a worktree unless noted. Launch
prompts in caveman. Every launch prompt starts with this block:

```
Rules. Invoke andrej-karpathy-skills:karpathy-guidelines first. Ponytail
ladder: delete before edit, edit before add, no new abstraction, no new
object slots, no refactor off the WP path. Lines <= 80. Touch only the
files the WP lists; if another file must change, stop and say why.
Before commit: Rscript -e 'Rcpp::compileAttributes()' (if src/ changed),
Rscript -e 'devtools::document()', R CMD INSTALL ., Rscript -e
'devtools::test()', Rscript -e 'lintr::lint_package()', then the WP
acceptance command; paste its output in the report. Commit on the
worktree branch before finishing or the work is lost. No Co-Authored-By,
no Claude-Session, no "Generated with" anywhere. Report in caveman:
files touched, lines +/-, test count before/after, acceptance output,
anything you did not do.
```

| Agent | Type | WP | Prompt adds |
|---|---|---|---|
| census | cavecrew-investigator | pre-WP2, re-run after each merge | run `Rscript dev/sharpen-census.R`, report table; for each Tier C name grep `prepare_data/` (gitignored, maintainer checkout) and report hits |
| snapshot | general-purpose, worktree | WP0 | file + command from WP0 |
| demoter | general-purpose, worktree | WP1 | Tier B list verbatim; "no renames" |
| pruner | general-purpose, worktree | WP2 | Tier C list verbatim; `git mv` to `dev/probes/` first; "do not touch R/modules.R" |
| modules | general-purpose, worktree | WP3 | WP3 list; "keep leiden + consensus only" |
| renamer | general-purpose, serial on branch | WP4 | the column list; `grep` acceptance |
| dieter | general-purpose, serial on branch | WP5 | the keep-list of formals; per export paste the formals it dropped and the grep showing nobody set them |
| driver | general-purpose, worktree | WP6 | signature verbatim; "plain list, delete container methods" |
| clades | general-purpose, worktree | WP7 | laminar rule; `ape::prop.part`; test cases |
| docs | general-purpose, worktree | WP8 | line budgets; reference groups |
| guard | cavecrew-builder | WP9 | budget numbers |
| reviewer | cavecrew-reviewer | every PR | review the diff for: new abstraction, new argument, new export, file outside WP list, column name outside the style, attribution line. One line per finding |

Parallel sets run as one Agent call with several tool uses. Worktree
agents that share a file must not run together: WP2 and WP3 both touch
`R/rcomplex-class.R` only through method deletion for disjoint
functions, so WP2 deletes the Tier C methods and WP3 leaves
`rcomplex-class.R` alone except `characterize_hubs`.

## 6. Expected end state

27 exports, ~11-12k R lines, ~5.2k C++, ~15k test lines, README 150
lines, a 10-line quickstart that is the whole BLAST analogue:

```r
library(rcomplex)
expr <- list(BDIS = read_expr("bdis.tsv"), HVUL = read_expr("hvul.tsv"))
res <- rcomplex(expr, "orthogroups.tsv",
                clades = list(annual = "BDIS", perennial = "HVUL"))
res
summary(res)
write_rcomplex(res, "out/")
```

MDO I and II then read: multi-species (any N, nested clades, species
tree in), network-aware nulls (edge-swap, shuffled expression, rank
test), effect sizes and power on every edge, p-value resolution on every
test. Nothing in the proposal's P1 that rcomplex does not do except the
rank-test default, which waits on Orion.
