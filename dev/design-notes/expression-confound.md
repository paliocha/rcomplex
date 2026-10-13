# Expression level and the conservation calls (WP-d)

2026-10-13. Branch `null/orthologs-strata`, base `81c48aa` (main after
WP-c2). Data: the EVOTREE wood run of `prepare_data/evotree/e5_rcomplex_0.4.R`
(6 species, density 0.03, alpha 0.1; angiosperm = Asp Birch Cher,
conifer = Nor Scots Lodge, pinus = Scots Lodge). Probe scripts
`e6_matched_go.R` and `e7_abundance_null.R` sit beside it; the data are
not in the repository. Expression is VST.

## 0. Verdict

- Nothing in the package conditions on expression level. The
  hypergeometric test conditions on the two neighbourhood sizes and the
  urn; MR ranks within gene; `power` is combinatorics on `m`, `k`, `N`
  and `rho0`. A low-expressed gene has a full-size neighbourhood at a
  density threshold, is tested with power, fails, and lands in
  `trait_specific` or `unclassified` as negative evidence.
- The tiers track expression. Odds of `complete_conserved` rise 1.91x
  per VST unit (logistic on 16,941 tested hogs, p 1e-132): P(core) is
  0.8 % at VST 4 and 10 % at VST 8. The core's GO enrichment shrinks
  20-40 % in fold against an expression-matched background and loses
  20-54 % of its significant terms; the wall, cytoskeleton and cell
  cycle terms survive.
- The hypergeometric null is not the right reference for a call.
  Random non-ortholog gene pairs in the real networks reach q < 0.1 in
  9-16 % of cases (module structure), twice as often when both genes
  are in the top expression deciles (18-24 %) as when one is low
  (7-14 %). The driver's shuffled-expression null gives 0 calls in
  every species pair and so measures none of this.
- Decision: `rcomplex(null = TRUE)` shuffles the hog labels among genes
  of the same expression stratum, per species, and repeats the whole
  pipeline on that table: edges, cliques, tiers. `calls_null` then
  counts what abundance and network structure alone produce, and the
  tier table gets a null column. No new argument. Section 4 has the
  spec, section 3 the alternatives.

## 1. What the pipeline controls

| quantity | conditioned on | blind to |
|---|---|---|
| hypergeometric p (`compare_neighborhoods()`) | `m`, `k`, `N - 1` (degree) | expression level, module structure |
| MR normalisation | rank within gene | level |
| `power` (`R/comparison.R:302`) | `m`, `k`, `N`, `rho0` | level: a noisy neighbourhood is full size |
| `underpowered` tier | too few mapped partners | too little signal |
| edge-swap null (`coexpressolog_null()`) | degree sequence | level, modules |
| `null_network()` (driver `null = TRUE`) | nothing: shuffles samples per gene | everything: 0 calls on wood (15/15 pairs) and on leaf |

## 2. Measurements

### 2.1 Tiers against expression (e5 4d, e6)

Per hog: mean VST over genes and species; CV over sections.

| tier | hogs | median VST | median CV |
|---|---|---|---|
| complete_conserved | 938 | 7.16 | 0.152 |
| unclassified | 2,702 | 6.25 | 0.122 |
| specific: angiosperm | 1,479 | 6.14 | 0.158 |
| partial_present | 247 | 6.14 | 0.200 |
| specific: conifer | 2,926 | 5.85 | 0.210 |
| underpowered | 511 | 5.57 | 0.124 |
| no clique | 8,138 | 4.82 | 0.175 |

Kruskal-Wallis p < 1e-300 on both. Logistic fits over the 16,941
tested hogs, `tier ~ mean + cv`:

| set | OR per VST unit | p | OR per CV unit | P at VST 4 | P at VST 8 |
|---|---|---|---|---|---|
| complete_conserved | 1.91 | 1.7e-132 | 22.4 | 0.008 | 0.101 |
| partial_present | 1.53 | 6.7e-17 | 15.2 | 0.003 | 0.017 |
| specific_angiosperm | 1.33 | 5.4e-40 | 3.2 | 0.035 | 0.103 |
| specific_conifer | 1.25 | 6.5e-15 | 3.0 | 0.022 | 0.051 |

Conditional on the mean, a higher CV raises the odds: flat genes have
no co-expression to conserve. The marginal "core has low CV" is a level
effect.

### 2.2 GO enrichment against a matched background (e6)

Universe: 7,367 tested hogs with Potra22 BP terms, 1,741 terms with 10
or more hogs. Unmatched: one hypergeometric. Matched: strata = mean VST
deciles over the universe; one hypergeometric per stratum conditional
on the set's count in it; exact p from the convolution of the stratum
pmfs (the exact test behind Cochran-Mantel-Haenszel). Sensitivity:
mean quintile x CV tertile.

| set | hogs | median VST (universe 5.94) | terms q < 0.1: unmatched / mean / mean x CV | top-15 median fold: unmatched / mean / mean x CV |
|---|---|---|---|---|
| complete_conserved | 722 | 7.23 | 110 / 92 / 51 | 2.33 / 1.89 / 1.41 |
| specific_angiosperm | 1,073 | 6.31 | 20 / 15 / 6 | 1.46 / 1.42 / 1.36 |
| partial_present | 159 | 6.66 | 5 / 8 / 9 | 3.63 / 3.92 / 3.84 |
| specific_conifer | 591 | 6.08 | 0 / 0 / 0 | 3.23 / 3.30 / 3.16 |

Generic terms shrink most (cellular component organization 1.48 to
1.31; vesicle-mediated transport 2.45 to 1.73). Cell wall biogenesis
(2.41 to 2.04, q 1e-4), microtubule-based process (2.33 to 2.02), cell
cycle (1.88 to 1.89) and endocytosis (7.1 to 3.2) stay. partial_present
gains: its DNA replication and break repair terms sit above the tier's
expression, so the unmatched test under-called them. Tables:
`go_matched_summary.tsv`, `go_matched_top15.tsv`, `go_matched_full_*.tsv`.

### 2.3 Random gene pairs through the test (e7)

30,000 random non-ortholog pairs per species pair, 300 per cell of
(decile of gene 1, decile of gene 2), appended to the real ortholog
pairs and run through `compare_neighborhoods()`, so the urn and the
mapping are the real ones. `pmax` is the larger of the two direction
p-values (the `pval_combine = "max"` rule). q by BH over real and random
rows together. Classes: low = deciles 1-3, high = 8-10.

| pair | random pairs q < 0.1 | high-high | low-low | low-high | real pairs q < 0.1 | real high-high | real low-high |
|---|---|---|---|---|---|---|---|
| Asp-Birch | 0.155 | 0.238 | 0.153 | 0.126 | 0.399 | 0.592 | 0.247 |
| Nor-Scots | 0.155 | 0.237 | 0.173 | 0.137 | 0.408 | 0.626 | 0.287 |
| Asp-Nor | 0.087 | 0.177 | 0.075 | 0.073 | 0.236 | 0.448 | 0.172 |

Mean fold enrichment of the random pairs rises from 0.93-1.00 (one gene
low) to 1.10-1.17 in the (10, 10) cell, with 1.06-1.24 in the (1, 1)
cell: genes of similar abundance share neighbours. Neighbourhoods have
about 700 genes, so a fold of 1.15 is many standard deviations. Nominal
p < 0.05 holds in 10-17 % of random pairs (5 % expected), p < 0.001 in
5-10 %.

Two findings. First, the uniform hypergeometric calls about one random
pair in seven at q < 0.1: two genes in one large module overlap by
construction. The real orthologs are called 2.5 times as often in every
class, so the signal is real, but the shuffled-expression null (0
calls) is not its reference. Second, abundance doubles the random-pair
call rate from the low-high to the high-high cells. The gene-identity
permutation null (`find_coexpressologs(method = "permutation")`, H0:
random genes from the same networks) is the right null for both; it
draws uniformly over genes and so lacks the stratification.

## 3. Alternatives

- A. Stratified draws in the permutation engines. `sample_k_distinct()`
  at `src/hog_permutation.cpp:407-408` and `src/fe_permutation.cpp:152-153`
  draws each hog gene's replacement from its expression stratum. Gives
  a per-hog p-value "conserved beyond abundance". Three engines (dense,
  sparse, torch), new int vectors through `permutation_hog_test()`,
  and the driver does not reach `method = "permutation"`. Later, if
  per-hog conditional p-values are wanted.
- B. Shuffle the hog labels within expression strata and rerun the
  pipeline (this WP). One `find_coexpressologs()` more per run (about
  5 s per species pair on wood), the same q-value machinery, a null
  tier table for free. A rate, not a per-hog p; one shuffle is enough
  (30,000 pairs put the rate within 1 %).
- C. Report only: per-hog expression on the result and a matched
  background helper. B needs the same primitive (`gene_mean`), so C's
  table can follow from it later without a new export.

## 4. Spec: WP-d expression-stratified ortholog null (about 150 lines)

- Files: `R/rcomplex.R` (driver 92-155, `.driver_orthologs` 160-173,
  `.driver_nulls` 199-204, `summary.rcomplex()` 295-319,
  `write_rcomplex()` roxygen 329-332), `R/network.R` (the three return
  lists at 398, 458, 480), `tests/testthat/test-rcomplex.R` (139, 165,
  198, 252), `tests/testthat/helper-rng-contract.R:230`, `README.md:74`,
  `vignettes/quickstart.Rmd` (56, 85, 91), `NEWS.md`, `man/`.
- Signature: unchanged. `null = TRUE` changes meaning; `rcomplex()` no
  longer needs `expr` for it (`method = "rank"` and `block` still do;
  fix the condition at line 109 and its message).
- `compute_network()` stores `gene_mean`: a named numeric, `rowMeans()`
  of the matrix it correlates, over the genes it keeps. The driver
  overwrites it from `expr[[s]]` (assay 1 of a SummarizedExperiment)
  after the build, as it does `params$n_samples` at line 121, because
  with `block` the wired matrix is centred. `as_network()` imports have
  no `gene_mean`.
- Shuffle, per species `s`, on the long ortholog table before
  `prepare_orthologs()` (split `.driver_orthologs()` so the filtered
  long table is kept): rows whose gene is in
  `rownames(networks[[s]]$network)`; stratum = `cut()` of `gene_mean`
  at its deciles over those genes, `include.lowest = TRUE`, breaks made
  `unique()`; permute the `gene` column within each stratum with
  `sample()`. Hog sizes per species and each gene's stratum are
  preserved; orthology is broken. Without `gene_mean`: one stratum and
  `message("null: <s> has no gene_mean, so the shuffle ignores
  expression")`. The draw runs inside the driver's `.seed_scope(seed)`.
- Null run: `edges_null <- .driver_edges(networks, xs, ortho_null, ...)`
  on the real networks; `g0 <- gene_clique_graph(edges_null, ...)` and
  `classification_null <- classify_gene_cliques(g0, edges_null, ...)`
  with the same `min_size`, `alpha`, `clades`, `min_power` as the real
  run. Store `edges_null` and `classification_null` only. For
  `method = "rank"`, build the rank calibration networks once and pass
  them to both `.driver_edges()` calls. `.driver_nulls()` stays for
  that.
- `summary.rcomplex()`: the tier table gains `cliques_null` (null
  cliques per tier, 0 where none; `NA` when `null = FALSE`). The
  `null` table keeps its columns; `calls_null` now counts the shuffled
  table's calls. `write_rcomplex()` writes `classification_null.tsv`
  by its existing loop; list it in the roxygen.
- Prose (STE-lite): `@param null`: "If `TRUE`, repeat the run with hog
  labels shuffled among genes of the same expression level. The null
  calls and tiers estimate what abundance and network structure alone
  produce." Replace every "shuffled networks" in the files above.
  `null_network()` keeps its own documentation.
- Tests: (1) on a fixture with a `gene_mean` gradient, every shuffled
  row keeps its gene's decile and every (hog, species) count is
  unchanged; (2) `edges_null` has the columns of `edges`;
  `classification_null` has the columns of `classification`; `summary()`
  returns `cliques_null`; (3) `seed` reproduces `edges_null`; (4)
  `rcomplex(networks = nets, null = TRUE)` runs, with the message when
  a network lacks `gene_mean`; (5) `test-rng-contract.R` passes with no
  new entry point (the shuffle is inside `rcomplex()`).

### Acceptance

- `R CMD check` OK on the built tarball, `lintr::lint_package()` clean,
  `Rscript dev/sharpen-prose.R` clean, surface test unchanged.
- Wood: `R_LIBS=<lib> Rscript prepare_data/evotree/e5_rcomplex_0.4.R`
  (local script, absolute data paths; skip if the data directory is
  missing and say so). Report `summary(res)$null` and the tier table.
  Expected: `calls_null / calls` about 0.3-0.4 in every species pair
  (it was 0.0), and a non-empty null tier table. Runtime within 2 min
  of the 3.8 min run.

### NEWS (development version)

- `rcomplex(null = TRUE)` now shuffles hog labels among genes of the
  same expression level and reruns edges, cliques and tiers on the
  shuffled table. `calls_null` counts the calls that abundance and
  network structure alone produce; the tier table gains `cliques_null`.
  The shuffled-expression null gave 0 calls on every data set and is
  gone from the driver; `null_network()` stays. `null = TRUE` no longer
  needs `expr`.
- `compute_network()` stores `gene_mean`.

## 5. After WP-d

- Re-read the wood report's "0 null calls" line; rerun and update it.
- The matched-background GO test (section 2.2) stays an analysis
  recipe; a package helper waits for a second use.
- Alternative A if per-hog conditional p-values are wanted.
