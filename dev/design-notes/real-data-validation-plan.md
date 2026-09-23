# Real-data validation of PRs #13 and #14 (and #7)

Date: 2026-09-15. Status: **plan only, nothing run.** Audience: a fresh
session (and Martin) picking this up once #13 and #14 are pushed. Everything
below was written from the state of the repo, the PRs and the Orion notes on
that date; verify each path and SHA before relying on it.

## 0. What is being validated

| PR | Issue | Branch (head at writing) | Change |
|---|---|---|---|
| #13 | #11 | `fix/clique-intensity-range` (`36c7054`) | `find_cliques()` `intensity` / `coherence` weight each clique edge by its **Jaccard-index percentile** among all tested pairs of its species pair (ranked before the `edge_type` filter), instead of `1 - q.value`. Warnings `rcomplex_prefiltered_edges` and `rcomplex_missing_jaccard`; `NA` intensity when any clique edge lacks a finite Jaccard. |
| #14 | #12 | `fix/degree-power-classification` (`4bbb57d`), **stacked on #13** | Per-edge `power` column (detection probability under a reference conserved fraction `f0`) in `comparison_to_edges()` / `find_coexpressologs()` (analytical path only) / `summarize_comparison()`; new class `"underpowered"` in `classify_cliques()` and `classify_gene_cliques()` (`min_power = 0.8`). |
| #7 (merged) | — | `main` | C++ degree-preserving rewire kernel in `coexpressolog_null()`. Production edge-swap null on Orion still outstanding. Optional section 6. |

Session evidence these checks should be compared against (all simulation,
none real data; details in issue comments #11 and #12):

- Pooideae clique means (nf-rcomplex March TSVs): `1 - mean_q` spanned
  0.958-0.997; Martin observed intensity 0.992-0.998, coherence
  0.9995-1.000 across all classes.
- Degree-controlled simulation (true conservation equal at every degree):
  effect-size percentile Spearman with degree **-0.97**; Jaccard **0.18**
  (0.12 among conserved edges under weak signal).
- Weak-signal simulation: call rate by degree decile 0.49 -> 1.00; median
  power at `f0` 0.52 -> 0.99; strict unreachability flagged 0/451 ns edges.
- Gene-graph fixture, all power 0.1: HOG5 differentiated -> underpowered,
  HOG4 unclassified -> partial_present; nothing else moved.

## 1. Critical caveat: the production pipeline cannot validate either PR as-is

Verified in `~/Documents/R/nf-rcomplex/bin/` (local copy):

- `run_evotree_cliques.R::build_edge_table()` (~line 120) and
  `run_clique_analysis.R` (~line 198) build edge tables **by hand** from
  **HOG-level permutation q-values** (`permutation_hog_test()`, Liang), with
  columns `gene1, gene2, species1, species2, hog, q.value, effect_size` only:
  **no `jaccard`, no `power`.**
- `find_differentiated_cliques()` subsets to `q.value < relaxed_q` and sets
  `type = "conserved"` on every remaining row (~line 182) before
  `find_cliques()`; `run_clique_analysis.R` passes a filtered
  `parallel_edges` (~line 318).

Consequences if the pipeline is simply re-run with the PRs installed:
#13 gives `NA` intensity everywhere plus `rcomplex_missing_jaccard` and
`rcomplex_prefiltered_edges` warnings; #14 changes no classification
(no `power` column; the permutation path would write `power = NA` anyway).

**So the validation builds its own analytical-path edge tables with
`find_coexpressologs()`** on the same networks, and treats the production
permutation-path outputs as a separate, unchanged reference. Adapting the
pipeline scripts (carry `jaccard`, stop pre-filtering, decide how power maps
onto HOG-level q) is a follow-up decision for Martin, not part of this plan.

## 2. Hypotheses and pass criteria

### #13 (intensity / coherence)
- **H13.1 range.** Real-clique intensity is no longer confined to a ~0.01
  band: report min / IQR / max per class and tissue. Pass: IQR spans a
  substantial part of (0, 1] (e.g. > 0.2) in every class with >= 30 cliques.
- **H13.2 degree neutrality.** Spearman(intensity, clique mean member degree)
  per tissue. Pass: |rho| <= ~0.25. Control computed ad hoc on the same
  cliques: the effect-size-percentile weight should show a strongly negative
  rho (simulation -0.97); if it does not, the simulation's premise does not
  hold on real data and #13's rationale must be revisited.
- **H13.3 membership unchanged.** Same cliques (hog + genes), `mean_q`,
  `n_edges`, `min_effect_size` as `main` on identical edge tables. Pass: identical.
- **H13.4 warnings.** On analytical-path unfiltered tables: no
  `rcomplex_prefiltered_edges`, no `rcomplex_missing_jaccard`. Record whether
  the production hand-built tables trigger both (expected yes, section 1).
- **H13.5 intensity test at real scale.** `clique_intensity_test()` on a
  stratified clique subset: share with `n_matched > 0`, z / p distribution,
  runtime. The null only matches by chance on small fixtures; real scale is
  where it must work. Pass: a meaningful share of cliques gets `n_matched`
  >= 20 at the chosen `n_perm`.

### #14 (power / underpowered)
- **H14.1 urn recovery.** No `.urn_size()` mismatch warning on any of the
  2 x 28 species pairs; recovered `N - 1` equals network gene count - 1.
- **H14.2 f0.** Per species pair and direction; plausible (simulation 0.167),
  no `NA` on pairs with calls.
- **H14.3 power by degree.** Share of non-conserved edges with
  `power < 0.8` by member-degree decile. Pass: strong monotone decrease with
  degree (the mechanism #12 claims).
- **H14.4 classification transitions.** Old = same edges with the `power`
  column dropped (tested to reproduce the pre-#14 behaviour); new = with
  `power`. Transition matrices per classifier and tissue. **Hard check:**
  only lineage_specific / differentiated / trait_specific -> underpowered and
  -> partial_present (gene graph) moves; any other direction is a bug.
- **H14.5 degree confound (issue #12 acceptance).** Member degree by class
  (Kruskal-Wallis; logistic or multinomial model class ~ log degree + trait
  group). Pass: the degree gap between trait_specific / differentiated and
  complete / partial shrinks after reclassification; report the before/after
  effect size, not just p.
- **H14.6 sensitivity.** `min_power` in {0.5, 0.8, 0.9}; `f0` in
  {NULL, 0.1, 0.2, 0.3}. Report class counts; flag any result whose sign
  depends on the choice.
- **H14.7 life-history balance.** Per-species degree distributions by trait
  group at density 0.03; if they differ, the underpowered share differs by
  group and the trait-level conclusions need that caveat.

## 3. Environment (Orion; from memory notes, verify first)

- SSH aliases `orion` (login, commands killed after 5 min: submit, don't
  run) and `orion-filemanager` (transfers). Use the `ohpcc-nmbu` skill.
- `module load R/4.4.2` (Apptainer; its g++ 13.2 handles C++23). Container
  prints "WARNING: ignoring R_HOME" to stdout: use `tail -1` on
  `$(Rscript ...)`.
- Production lib `/mnt/users/martpali/R/library-4.4`; pipeline
  `/mnt/project/FjellheimLab/martpali/AnnualPerennial/rcomplex-fastoma`;
  scratch `/mnt/project/FjellheimLab/martpali/rcomplex-testrun/`.
- Partitions `orion` / `GPU` / `OOD`. Prefer `.tsv(.gz)` outputs for tables,
  `.rds` for networks.
- **Do not overwrite the production lib.** Install each version into its own
  lib, e.g. `.../R/lib-main-<sha>` and `.../R/lib-pr14-<sha>`, and set
  `R_LIBS` per job. #14's head contains #13, so one PR lib covers both.

## 4. Design choices

- **One network build, reused by every version.** Networks do not depend on
  #13 / #14. Compute or reuse per species x tissue `.rds`.
- **Old vs new inside one run where the code allows it:** #14 old behaviour =
  drop `power` (tested path). #13 old intensity = recompute `1 - q` weights ad
  hoc in the analysis script (documented formula), plus a `main`-lib run on
  one tissue for H13.3 membership identity.
- **Parameters match the pipeline** (`nextflow.config`): 8 species
  `BDIS, HVUL, BMAX, VBRO` (annual) / `BSYL, HJUB, BMED2BMAX, FPRA`
  (perennial), tissues `leaf, root`, density 0.03, Pearson + MR,
  `min_var = 0.1`, `q_conserved = 0.10`, congener pairs
  `BDIS=BSYL, HVUL=HJUB, BMAX=BMED2BMAX, VBRO=FPRA`, FastOMA `RootHOGs.tsv`,
  seed 42. PANN / PSUP excluded.
- **Baseline is `main` at the same inputs, never the March TSVs:** 0.2.0 /
  0.3.0 changed pi0, `pval_combine` default and the urn, so March numbers
  differ for reasons unrelated to these PRs.
- **Memory:** up to ~53k genes (HJUB before filtering). Dense MR transient
  ~1.5 n^2 doubles ~ 34 GB at 53k; the pipeline requested 128 GB / 16 CPU
  for network builds. Keep that.

## 5. Steps

- **S0 prepare (login node, short).** `git fetch`; build tarballs for `main`
  and #14's head (`R CMD build`); install to the two libs via a short job;
  record SHAs and `sessionInfo()` into the output directory.
- **S1 networks.** Check whether per-species x tissue networks at density
  0.03 already exist under the pipeline results (and whether they are sparse
  `compute_network()` objects from >= 0.2.0). Reuse if so; otherwise an array
  job of 16 (`compute_network(..., density = 0.03, sparse = TRUE)`), 128 GB,
  16 CPU.
- **S2 analytical edges (#14 lib).** Per tissue:
  `find_coexpressologs(networks, orthologs, method = "analytical",
  alpha = 0.10, pval_combine = "max", pi0_method = "randomized", seed = 42,
  out_file = "edges_<tissue>.tsv.gz")` over all 28 pairs. Capture warnings
  (H14.1). One tissue also with the `main` lib for H13.3.
- **S3 cliques + intensity (#13).** `find_cliques()` on the unfiltered table
  for the full 8 species and each trait group (min_species as in the
  pipeline). Add ad hoc `intensity_1mq` and `intensity_effpct`; add clique
  mean member degree from the networks. Tables for H13.1-H13.4.
- **S4 classification (#14).** `classify_cliques(edges, species,
  species_trait = life_cycle)` and
  `classify_gene_cliques(gene_clique_graph(edges), edges, species,
  lineage = <trait groups>)`, each with and without `power`, plus the
  sensitivity grid. Transition matrices (H14.4, H14.6).
- **S5 degree diagnostics (#12).** Degree per gene per network; per class and
  per trait group (H14.3, H14.5, H14.7).
- **S6 intensity test (#13, heavy, optional).** Stratified subset (e.g.
  <= 50 cliques per class), `n_perm = 5` timing run first, then `n_perm` sized
  from it (H13.5). Every permutation reruns all pair comparisons.
- **S7 edge-swap null (#7, optional).** See section 6.
- **S8 report.** A markdown report with the tables above and the pass / fail
  of each hypothesis; post short summaries on PR #13 and #14 before they are
  merged; update `project_clique_metric_ceiling` memory.

## 6. Optional: #7 edge-swap null in production

- Script `meristems/paper/analysis/_source/61_coexpressolog_edgeswap_null.R`
  (reads cached `networks/<sp>_network.rds` and `edges.tsv.gz`; env
  `PERM=100 SWAP=10 ALPHA=0.1`), runscript `meristems/slurm_logs/coexp_null.slurm`.
- Previous run: job 1213069, 12 CPU / 64 G; igraph rewiring ~100% of cost,
  forks > 22 min each, cgroup OOM kills. Benchmarks in
  `/mnt/project/FjellheimLab/martpali/rcomplex-testrun/bench_edgeswap/`.
- Rerun with a `main` lib; compare wall time per permutation and peak RSS
  (`sacct -j <id> --format=JobID,Elapsed,MaxRSS`) against 1213069. Seeded
  null values change numerically by design; compare distributions, not
  numbers. Requesting physical cores (`n_cores` <= physical) avoids the
  hyperthread contention seen before.

## 7. Clean-session start checklist

1. Read this file, memory `project_clique_metric_ceiling.md` and
   `project_ci_state.md`, the PR #13 / #14 descriptions, and the comments on
   issues #11 / #12.
2. `gh pr view 13 --json headRefOid,state`; `gh pr view 14 --json
   headRefOid,state,baseRefName`. If #13 has merged, #14 targets `main`.
3. Verify on Orion: lib paths, pipeline results directory, whether networks
   exist (S1), free scratch space.
4. Confirm with Martin: whether S6 and section 6 are in scope, and the
   `min_species` / `lineage` definitions to use in S3-S4.

## 8. Session pickup, 2026-09-17 (state verified, decisions taken)

Verified this session; supersedes the "head at writing" SHAs in section 0.

- **#13 is merged.** `main` is `7a57766` ("Merge pull request #13"), which
  contains `36c7054`. Section 0's `fix/clique-intensity-range` row is history.
- **#14** head `4bbb57d`, base already retargeted to `main`, `MERGEABLE` /
  `CLEAN`, but **no CI had ever run on it** (it was stacked on #13 while the
  workflows only fire for PRs into `main`). Marked ready this session so lint,
  the two R-CMD-checks, test-coverage and pkgdown run, and Copilot reviews,
  **before** anything is submitted on Orion.
- **Orion, read-only probe.** SSH aliases work. `R/4.4.2` is available
  (modules go up to `R/4.6.0`). `/mnt/users/martpali/R/library-4.4` holds
  **rcomplex 0.2.0** (built 2026-09-02) and is the only lib with the package,
  so S0 must build both new libs. `/mnt/project/FjellheimLab` has 117 TB free.
- **S1 inputs exist.** 16 `results/networks/<tissue>__<species>_network.rds`
  (1.1-5.7 GB each) plus a duplicate set under `results/modules/<tissue>/`,
  all dated **2026-03-29** — i.e. built *before* 0.2.0, so their sparseness and
  object layout must not be assumed. Job `1317571` inspects
  `root__BDIS_network.rds` for class, `dgCMatrix` vs dense, dim and stored
  density; if they are pre-0.2.0 dense objects, S1 rebuilds rather than reuses.
  `FastOMA-HOGs/RootHOGs.tsv` (30 MB, 2026-03-29) is in place.
- **Working directory** for everything below:
  `/mnt/project/FjellheimLab/martpali/rcomplex-testrun/validation-2026-09-17/`.

Martin's answers to the section 7 checklist:

1. **Both optional sections are in scope**: S6 (`clique_intensity_test()` at
   real scale, timing run first) and section 6 (the #7 edge-swap null rerun).
2. **CI and review come before the Orion run** — the lib is built from the
   reviewed #14 head, not from `4bbb57d` as-is.
3. **S3-S4 parameters: all 8 species, plus one run per trait group**
   (4 annual `BDIS, HVUL, BMAX, VBRO`; 4 perennial `BSYL, HJUB, BMED2BMAX,
   FPRA`), `min_species = 3` within a group. The pipeline's own
   `min_lineages = 3` / `min_within_pairs = 2` / `min_sp_offset = 1` stay with
   the production permutation-path reference, which this run does not touch.

### S1 resolved, 2026-09-17: reuse the March networks

Job `1318139` rebuilt `root__BDIS` from `vst_hog.RDS` with the pre-#13
baseline lib and compared it against the March object:

| Quantity | March object | Rebuilt today |
|---|---|---|
| genes into `compute_network()` | 23,661 | 19,960 |
| `n_genes` / `n_removed` | 12,240 / 11,421 | 10,293 / 9,667 |
| threshold | 11,829.5 | 9,946.8 |
| MR values, 2000 common genes | Spearman **0.99987** | |
| thresholded edge sets | Jaccard **0.946**, 96.95% of old edges kept | |

`compute_network()` therefore reproduces: every one of the 10,293 rebuilt
genes is among the old 12,240, and the MR values agree to three decimal
places of rank correlation. The whole difference enters **upstream**, in
`reduce_orthogroups()`: the March run merged ~1,000 paralogs where today's
merges 4,709. The March command line (recovered from
`work/*/*/.command.sh`) is `--cor_threshold 0.7`, the same value used in
the rebuild, so this is a version difference in the Ward.D2 merge, not a
parameter difference.

**Decision: reuse the 16 March networks as the common input.** Both arms of
every comparison read the same networks, so the reduction difference cancels;
rebuilding would cost 16 x 128 GB jobs and would *change* the gene set
relative to every earlier result. S1 needs no array job. All 16 are
structurally valid for today's code: plain dense matrices with `dimnames`,
`list(network, threshold, n_genes, n_removed, params)`, which is exactly what
`compute_network(sparse = FALSE)` still returns, and `.net_check()`'s dense
branch accepts them (the store guard applies to `dgCMatrix` only).

Gene lists and thresholded degrees for all 16 are extracted once into
`validation-2026-09-17/netinfo/` for S5.

### #14 review round, 2026-09-17

CI on `6716f0d` is green on all four workflows (lint, both R-CMD-checks,
test-coverage). Copilot raised three inline findings:

1. **`pmax(..., na.rm = TRUE)` returns `-Inf` when both directions are `NA`**
   -- *rejected, not reproducible.* R returns `NA` when every input is `NA`;
   only `max(numeric(0))` gives `-Inf`. Checked directly:
   `pmax(NA_real_, NA_real_, na.rm = TRUE)` is `NA`, and the vector form
   keeps `NA` in the all-`NA` position.
2. **The `underpowered` check sat in the `else` of "no significant test"**
   -- *real, fixed.* A species with one significant edge that still cannot
   join the clique is kept out by its failures, so the check now runs on the
   failed rows after the extendable test.
3. **`filter_zero`** -- *real, fixed per Martin's decision*: the analytical
   batch now passes `filter_zero = FALSE`, exposed as an argument on
   `find_coexpressologs()` and `density_sweep()`. Zero-overlap tested pairs
   reach the edge table with their `power` instead of arriving at the
   classifiers as `absent`/`untested`. Every analytical q-value rises, since
   the multiple-testing set grows; the ComPlEx_python equivalence tests pin
   `filter_zero = TRUE` to keep correcting over the canonical set.

This enlarges the S2 edge tables by the number of zero-overlap ortholog
pairs, which is why the per-pair counts are measured before S2 is sized.

### Decisions and measurements, 2026-09-17 (second half)

**`coexpressolog_null()` keeps `filter_zero = TRUE`.** It forwards `...` to
`find_coexpressologs()`, so the new `filter_zero = FALSE` default reached it
and would have added every zero-overlap ortholog pair to *every* permutation.
Its statistic reads called edges only, so those rows are pure cost. It now
takes `filter_zero = TRUE` as its own default and passes it to both the
observed and the permuted runs. Section 6's edge-swap null is therefore
unaffected by the `filter_zero` change.

This surfaced as a test failure rather than by inspection: the validator
behind `"statistic is missing"` only fires when a rewired permutation yields
an *empty* edge table, which retention prevented. The control mattered --
the same test passes on pristine `6716f0d`, which is what identified the
change as the cause rather than test-order flakiness.

**S2 sizing (measured, job 1318392).** Ortholog pairs to be tested per
tissue, over all 28 species pairs: **leaf 471,918**, **root 419,434** --
roughly 15-20k per pair. Retaining zero-overlap rows is therefore cheap.

**S2 dry run (job 1318599, baseline lib, `root` BDIS-HVUL).** 12,870 ortholog
pairs in, 12,669 edges out, 6,071 conserved at `alpha = 0.10`, **68 s wall**
including the 220 MB `vst_hog.RDS` load. The 201-row gap is the old
`filter_zero = TRUE` default dropping zero-overlap pairs (1.6% of tested
pairs); under the fixed #14 head the same pair should emit all 12,870, which
is the check to run first.

Accordingly `s2.slurm` is right-sized from 250 GB / 12 h / 16 CPU to
**64 GB / 2 h / 8 CPU**.

**Ortholog source for S2.** `FastOMA-HOGs/RootHOGs.tsv` keys on *protein*
IDs (`Bmax_egapxtmp_028053-P2`) which do not join to the networks' gene IDs
(`Bradi1g45090.v3.2`). S2 therefore takes HOG membership from
`vst_hog.RDS` (`HOG`, `abbrev`, `gene_id`) and builds per-pair tables with
the pipeline's own `prepare_orthologs_from_hog()` (`max_paralogs = 10`,
expressed-gene filtered), exactly as `run_networks.R` does.

**Open regression.** `R CMD check` with vignettes built fails on my branch in
`rcomplex-tutorial.Rmd` with `map must be a data frame from
resolve_ortholog_map()`; the same vignette builds clean on pristine
`6716f0d`. `R/module_preservation.R` and `R/ortholog_map.R` are byte-identical
between the two, and `.map_coexpressolog_layer()` already filters
`type == "conserved"`, so paralog resolution never sees the retained `ns`
rows. Under diagnosis; the `--no-build-vignettes` shortcut is what hid it, so
the documented workflow (build *with* vignettes) is the one to run.

### S2 complete, 2026-09-17: analytical edges for both tissues

Library `lib-pr14-1c37c7c` (PR #14 head after the review fixes), verified on
load: `.jaccard_percentile` and `.edge_power` present, `filter_zero` default
`FALSE`.

**Single-pair verification first** (job 1320706), against the same pair the
baseline dry run used:

| Run | ortho pairs | edges | conserved | power NA |
|---|---|---|---|---|
| baseline lib (pre-#13) | 12,870 | 12,669 | 6,071 | n/a (no column) |
| `lib-pr14-1c37c7c` | 12,870 | **12,870** | 6,066 | **0** |

The 201-row gap was exactly the zero-overlap pairs, as inferred from the dry
run and now measured. Every retained row carries a computed power, so those
low-degree failures reach the classifiers as *underpowered* rather than as
`absent` / `untested` gaps -- the artefact #12 exists to remove, now shown on
real data rather than simulation. Conserved moved 6,071 -> 6,066: five calls,
the expected small rise in q-values under the larger multiple-testing set.

**Full runs** (jobs 1320735 root, 1320736 leaf), 28 species pairs each:

| Tissue | pairs | ortho pairs | edges | conserved | power NA | compute |
|---|---|---|---|---|---|---|
| root | 28 | 419,434 | 419,434 | 164,386 | 0 | 0.38 min |
| leaf | 28 | 471,918 | 471,918 | 153,938 | 0 | 0.42 min |

`edges == ortho_pairs` on **all 56 pairs** (zero mismatching rows), and
`power_na = 0` throughout: retention holds everywhere, not just on the pair
spot-checked. The totals match the S2 sizing job exactly (471,918 / 419,434),
so that estimate is confirmed rather than merely plausible. Wall time is ~3
min per tissue including the network loads; the 64 GB / 2 h request is ample.

Outputs: `s2_root/` (18 MB) and `s2_leaf/` (20 MB), one
`edges_<tissue>_<sp1>_<sp2>.tsv.gz` per pair plus `s2_summary_<tissue>.tsv`.

These tables are the input to S3 (cliques + intensity) and S4
(classification). They are the **unfiltered** analytical-path tables the
gap tiers need: `classify_gene_cliques()` must see the rows that were tested
and failed in order to refuse them.

### S3 complete, 2026-09-18: cliques + intensity, and a finding against #13

Jobs 1320969 (root), 1321152 (leaf), plus diagnostics 1321159 and 1333019.
`find_cliques()` on the unfiltered analytical edge tables, three runs per
tissue: all 8 species and each trait group, `min_species = 3`.

| tissue | run | cliques | intensity IQR | coherence IQR | IQR `1 - q` | IQR effect pct |
|---|---|---|---|---|---|---|
| root | all8 | 14,665 | 0.136 | 0.007 | 0.012 | 0.132 |
| root | annual | 4,431 | 0.160 | 0.006 | 0.006 | 0.164 |
| root | perennial | 3,571 | 0.117 | 0.004 | 0.008 | 0.127 |
| leaf | all8 | 13,289 | 0.105 | 0.004 | 0.013 | 0.109 |
| leaf | annual | 4,572 | 0.136 | 0.003 | 0.004 | 0.141 |
| leaf | perennial | 3,306 | 0.103 | 0.002 | 0.006 | 0.112 |

**H13.1 range: partial.** The ceiling is gone -- intensity spans 0.35-0.999
where `1 - q` spans 0.91-1.000, an IQR roughly 11x wider -- but 0.10-0.16 is
below the plan's `IQR > 0.2` bar. Recorded as partial, not passed.

**Coherence was not fixed.** IQR 0.002-0.007, median 0.997-0.999: still
ceiling-bound under Jaccard weights, exactly as it was under `1 - q`. #11 was
framed as fixing intensity *and* coherence; only intensity moved.

**H13.2 degree neutrality: fails, and inverts the simulation.** Spearman of
intensity against clique mean member degree:

| run | Jaccard pct | `1 - q` | effect pct |
|---|---|---|---|
| root all8 | **+0.581** | +0.295 | -0.064 |
| root annual | +0.521 | +0.329 | +0.054 |
| root perennial | +0.488 | +0.253 | -0.145 |
| leaf all8 | **+0.530** | +0.301 | +0.127 |
| leaf annual | +0.634 | +0.510 | +0.364 |
| leaf perennial | +0.461 | +0.275 | +0.101 |

The simulation predicted Jaccard ~+0.18 and effect-size percentile -0.97. On
real data the ordering is reversed: the weight #13 moved *to* is the most
degree-dependent, the weight it moved *from* is the near-neutral one. Leaf
replicates root, so it is not a tissue quirk.

**Four controls rule out species composition**, in increasing strictness:

1. `rho(mean_degree, n_species)` = 0.076.
2. Within `n_species` strata, rho stays 0.637 / 0.606 / 0.605 / 0.559 / 0.483
   / 0.323 for n = 3..8.
3. Within *identical* species composition (219 compositions per tissue,
   strata n >= 100): root median rho 0.643 over 34 strata (0.323-0.846),
   leaf median 0.593 over 39 strata (0.294-0.754).
4. Within a single species pair at edge level, the Jaccard index itself
   tracks `sqrt(deg1 * deg2)`: rho +0.30 to +0.54 over all 8 pairs tested,
   against +0.04 to +0.34 for effect size.

The mechanism is mechanical, not biological: for random neighbourhoods of
size `k`, `m` from `N` genes, `E[Jaccard] ~ km / (N(k + m) - km)`, which rises
with degree. At fixed density degree scales with network size, so Jaccard's
*null* expectation drifts upward with degree before any conservation signal
enters, while the hypergeometric effect size is observed/expected and so is
degree-normalised by construction. Ranking within a species pair does not
help, because the degree variation is within the pair too.

Connectivity being what selection acts on (Mahler et al. 2017) argues for
measuring conserved connectivity well, not for a statistic whose null rises
with degree: a hub and a low-degree gene with equally conserved
neighbourhoods should not score differently because one neighbourhood is
larger. If hub conservation is the claim, it should be measured and reported,
not inherited from the metric.

Filed as **#15**; the shipped weighting is left in place so #14 can finish.

**Also: `power` correlates negatively with degree on real data** (rho -0.31
to -0.62 over the same 8 pairs), where #12's simulation had detection
probability rising with degree (0.52 -> 0.99). This does not affect #14's
correctness but changes how its classification output should be read, and
needs its own check before S4 results are interpreted.

### S4 root, 2026-09-18: classification moves only where #14 says it may

Job 1333360, `f0 = NULL`, root, 419,434 edges (all carrying power, mean 0.708).

**HOG level** (`classify_cliques`, 16,405 HOGs). Counts move only in the
sanctioned direction, and only one transition type occurs at any threshold:

| `min_power` | moved | transition |
|---|---|---|
| 0.5 | 182 | `trait_specific` -> `underpowered` |
| 0.8 | 228 | `trait_specific` -> `underpowered` |
| 0.9 | 243 | `trait_specific` -> `underpowered` |

`complete` (359), `partial` (7,764) and `unclassified` (8,011) are identical
across all four arms. `trait_specific` falls 271 -> 89 -> 43 -> 28.

**Gene-graph level** (`classify_gene_cliques`, 53,844 cliques from 276,537
graph rows). Again one transition type, `unclassified -> partial_present`:

| arm | complete_conserved | partial_present | unclassified | underpowered |
|---|---|---|---|---|
| no power | 7,910 | 1,941 | 43,993 | 0 |
| 0.5 | 7,910 | 2,192 | 43,742 | 0 |
| 0.8 | 7,910 | 2,541 | 43,393 | 0 |
| 0.9 | 7,910 | 3,152 | 42,773 | 9 |

**H14.4 passes.** The hard check was that only `lineage_specific` /
`differentiated` / `trait_specific` may become `underpowered`, and only
`partial_present` may rise. Observed exactly that: `trait_specific ->
underpowered` at HOG level, `unclassified -> partial_present` in the gene
graph, nothing else, no backward moves, `complete_conserved` fixed at 7,910
in every arm.

**But read the counts against #16.** Power falls with degree on this data, so
the 243 HOGs withdrawn at `min_power = 0.9` are the high-degree ones, not the
low-degree failures #12 set out to protect. The mechanism is doing what the
code says; what it selects is not what the issue text describes.

`underpowered` barely fires in the gene graph (9 cliques, and only at 0.9)
because most edges have high power. Combined with the f0 grid saturating at
`f0 >= 0.2` (#16), the informative range for this parameter on real data is
narrow.

### H13.1 retired, 2026-09-21: the weight is fixed, the range criterion is not meaningful

`feat/onnela-maxent-weight` replaces the Jaccard-percentile Onnela weight
with an ensemble connection probability: association strength
(`effect_size`) mapped through `p = z w / (1 + z w)`, with `z` fitted per
species pair by maximum entropy (Garlaschelli, Ahnert, Fink & Caldarelli
2013). Validated on the root tissue, same 14,665 cliques as job 1320969:

| weight | intensity IQR | intensity rho(degree) | coherence IQR | coherence rho |
|---|---|---|---|---|
| Jaccard percentile (#13) | 0.136 | **+0.581** | 0.0070 | **+0.371** |
| ensemble probability | 0.0821 | **-0.077** | 0.0034 | **-0.063** |

`NA intensity: 0 of 14,665` -- the degenerate-bracket fallback never fires
on real pairs, which carry 15-20k tested edges each.

**H13.2 passes** (|rho| <= 0.25 on both statistics, both tissues' cliques
built from the same edge tables). **H13.1 is retired rather than failed.**
Its `IQR > 0.2` bar was written when intensity was a percentile, which is
uniform by construction and so has an IQR of 0.5 before any aggregation.
A probability-valued weight has no such guarantee: `p = zw/(1 + zw)`
concentrates toward the middle of (0, 1), so a narrow IQR is a property of
the map, not evidence of a dead metric. Ranking cliques on intensity stays
valid within a comparable set; the number is not meant to spread across
the unit interval.

Martin's decision: ship the MaxEnt weight, drop the IQR criterion, and
document the narrow range rather than chase it.

**Still open:** coherence has never discriminated on this data under any
weight (IQR 0.002-0.007), except under an `AS / max(AS)` scaling that gave
IQR 0.044 but uses the raw maximum the 2013 paper supersedes. Whether
coherence earns its place at all is a separate question from #15.

## 9. Execution plan, 2026-09-21: phases 0-5

Everything downstream derives from the networks, and Martin's `log = TRUE`
decision changes them, so the order below is forced: land code, rebuild
networks, rebuild derived layers, then diagnose.

### Phase 0 -- land the open PRs

- **#17** (`feat/onnela-maxent-weight`) merges when its two R-CMD-checks
  finish; lint and coverage already pass. Closes **#15**. Delete branch.
- **#16** stays open: a finding about `power`, not a defect in merged code.
  #14's merge comment already tells readers to interpret through it.
- **#4** stays untouched (CLAUDE.md).

### Phase 1 -- rebuild the 16 networks with `mr_log_transform = TRUE`

Raw MR is unbounded with a surviving-weight `max/min` of 1.035 at density
0.03: after thresholding the weights carry almost no information, so
"weighted" Leiden and `avg.weight` are weighted in name only. `S = 1 -
log(MR)/log(n)` is bounded in [0, 1] and measured `max/min = 2.464`
(range 0.375-0.925) on `root__BDIS`.

- Array job over 8 species x 2 tissues, 128 GB, 16 CPU, from
  `vst_hog.RDS` via the pipeline's `prepare_expression_matrix()` and
  `reduce_orthogroups(cor_threshold = 0.7)`.
- **Gate:** edge-set Jaccard against the March networks ~0.92 and
  surviving-weight `max/min` ~2.46, both measured in the S1 spot-check.
  A miss stops the phase rather than propagating into S2.

Note the 8% edge turnover: the log transform flips the rank direction, and
`sqrt(r_i r_j)` on descending ranks is not a monotone function of the
ascending version, so the selected edge set genuinely differs. Every
number in sections 5-8 above was computed on the old networks.

### Phase 2 -- rebuild the derived layers

Build one lib from merged `main` (carrying #17's weight and #14's power
column), then:

- **S2** edges, both tissues, 28 pairs each, ~3 min/tissue.
  Gate: `edges == ortho_pairs` on all 56 pairs, `power_na = 0`.
- **S3** cliques + intensity, ~10 min/tissue.
  Gate: intensity rho with mean member degree stays near 0 -- confirms the
  #17 fix survives log-weighted networks.
- **S4** classification, ~7 min/tissue.
  Gate: transitions only `trait_specific -> underpowered` and
  `unclassified -> partial_present`.

S3 and S4 are re-runs, but they are cheap and they keep every number in
the final report from one configuration rather than two.

### Phase 3 -- S5, degree diagnostics (~15 min)

Member degree by class and by trait group; Kruskal-Wallis plus
`class ~ log(degree) + trait`. H14.5 asks whether the degree gap between
trait-specific and conserved classes shrinks after reclassification;
H14.7 whether the trait groups differ in degree distribution at all.
Given #16, expect the opposite of #12's prediction and report it that way.

### Phase 4 -- S6, the intensity test, which is also the dynamic-range answer

Intensity's narrow IQR has three compounding causes: clique edges are
pre-selected to `q < alpha` (restricted range), the MaxEnt map pulls
weights toward `p = 0.5` by construction (entropy is maximal there), and
a geometric mean over `E` edges shrinks spread like `1/sqrt(E)`.

Spread is not the goal; discrimination is. The principled quantity is the
clique's intensity standardised against a **size- and
composition-matched** null: `z = (observed - null_mean) / null_sd`. That
is unbounded, has real spread, and reads as "how unusual is this clique".

- `n_perm = 5` timing run on a stratified subset first, then size the real
  run from it.
- Report `z` and `p` per clique; rank on `z`, not on raw intensity.
- Also report `mean(log p)` (log-intensity), which does not saturate.
- **Gate:** a meaningful share of cliques reaching `n_matched >= 20`. If
  the null rarely rebuilds a matching clique, that is a finding about the
  test, not about the data.

### Phase 5 -- report and close

- S8 report against every H13/H14 criterion, H13.1 recorded as retired
  with the S6 `z` distribution standing in for it.
- Open a follow-up on **coherence**: IQR 0.002-0.007 under every weight
  tested. If its `z` does not discriminate either, propose removing it
  rather than shipping a constant.
- Decide **#16** with S5 and S6 evidence in hand.

---

## 10. Outcome, 2026-09-21/22 — what landed, and two retractions

Phases 0-4 ran. This section is the handoff: what is true now, what was
withdrawn, and what is still open.

### Merged to main

| PR | What | Commit on main |
|----|------|----------------|
| #18 | within-HOG shuffle for `clique_intensity_test()` | `df11047` |
| #19 | `null_model = "matched_edges"` (matched edge-set null) | `bdf045f` |
| #21 | exact hypergeometric urn carried out of `compare_neighborhoods()` | `9e3f134` |
| #23 | matched-edge pools matched on clique size | `b191d83` |
| #24 | `underpowered` demoted from a classification to a flag column (breaking: `classification` no longer takes `"underpowered"`) | `0a54a27` |
| #26 | `coherence` retired from `find_cliques()` (#20) | `3a99584` |
| #27 | `power` computed against fold enrichment; `rho0` replaces `f0` (#22, #16) | `de1c88d` |
| #28 | `gap` and `n_edges` columns from `clique_intensity_test()` (#25) | `fc3bd4a` |

### Two results were retracted. Read this before trusting older sections.

**1. The `within_hog` real-data figure is withdrawn.** Section 9 Phase 4
anticipated `within_hog` rescuing the null. The reported "13 of 13 cliques
got a usable z" was produced by a **bug**: the shuffle grouped by `hog`
alone, and on a stacked all-pairs ortholog table that moves genes across
species boundaries. `compare_neighborhoods()` then drops those rows by
network membership — silently. Measured: 10 of 30 rows (33%) on a stacked
3-species table carried a wrong-species gene. Fixed in #18 by partitioning
on `(hog, species1, species2)`. With the fix, `within_hog` is the
**identity** for a single-copy HOG, so its null has no spread at all.

**2. A "re-measurement" that confirmed nothing.** Re-running the
power/degree correlation against the new library appeared to give
bit-identical numbers, which read as "the urn does not matter". It was
invalid: the script reads `power` from stored `s2_*_log/edges_*.tsv.gz`
tables written hours earlier by a *different* library, and those tables
carry no urn columns. It compared a vector with itself. The valid test is
`urn_isolate.R`, which recomputes both paths in one process.

### The intensity null: the permutation approach is structurally unavailable

Clique intensity is computed **conditional on the clique existing**, and
clique existence depends on the ortholog mapping being right. Any null
that perturbs that mapping destroys what it measures. Measured on a
stacked 3-species fixture:

| null | rows changed | species kept | outcome |
|------|--------------|--------------|---------|
| `global` | 0.96 | **0.55** | 0/30 matched — ~45% of rows dropped by the membership filter |
| `within_hog` | **0.00** | 1.00 | identity for single-copy HOGs, `null_sd = 0` |
| species-preserving relabel (prototype) | 1.00 | 1.00 | **0/30 permutations produced any clique at all** |

There is no setting in between: "enough randomisation to have a null" and
"little enough that a clique survives" are the same knob. Note also that
`global`'s failure was never "HOGs destroyed" — the `hog` column is never
touched — it is the same silent species-mismatch drop, worse.

**Replacement: `null_model = "matched_edges"`.** Hold the mapping fixed;
replace each clique edge with one drawn from that species pair's own pool.
Nothing is re-clustered, so no clique has to be rebuilt.

### S6 under the matched null, pair-only pools (Phase 4 answered)

**Correction, 2026-09-22.** These tables were produced by
`lib-matched-fbe547b`, which is the **#19-only** commit (15:23 on 09-21,
no `match_clique_size`), although it had been recorded as the "#19+#23
line". Everything in this subsection and the next therefore describes
the **pair-only** null; #23's size-matched pools first ran on real data
on 2026-09-22 — see "#23 on real data" below. The mislabel was caught by
`s7_coherence.R`'s reproduction guard, not by anyone reading a log.

Root `s6_root_matched/`, leaf `s6_leaf_matched/`, `n_perm = 2000`, all
cliques (no `MAXCL` subset needed).

- usable `z` for **14597/14597** root and **13190/13190** leaf
- ran in **5.8 / 5.5 min** at **MaxRSS 820 MB / 924 MB** — against a 12 h,
  340 GB job that returned all-`NA`. Future runs can request ~4 GB.
- **degree bias absorbed**: `cor(z, mean_degree)` = **-0.041** root,
  **+0.116** leaf. This was the stated risk from #15 and it did not
  materialise; no degree stratification needed.

### The clique-size effect, and why #23 exists

Median `z` ran **-0.51** (3 species) to **+6.32** (8 species). Two
compounding causes, both measured:

1. `null_sd` shrinks as `1/sqrt(E)` — `null_sd * sqrt(E)` constant to
   within 1% (0.0820-0.0829 root).
2. **Selection, the larger term.** A clique exists only because all its
   edges passed together, so larger cliques are built from stronger edges.
   Median edge weight by the clique size it belongs to: 0.615, 0.641,
   0.663, 0.687, 0.720, **0.752** for sizes 3-8, against a pair-only pool
   median of **0.638** and 0.609 for edges in no clique.

So ranking cliques on raw `z` was largely ranking them on size. #23 pools
on `(species pair, clique size)`; a **membership**, not an edge, is the
unit, because 3904 of 8212 HOGs produce more than one clique.

**Still report the gap (`observed - null_mean`) separately** — it runs
-0.024 to +0.098 (root) and is the substantive descriptive finding.
Folding it into a single score hides it.

### #23 on real data (2026-09-22): the size gradient is gone

`lib-main-0a54a27` (main after #24), `s6_root_matched23/`,
`s6_leaf_matched23/`, `n_perm = 2000`, `seed = 42`, submitted through the
ohpcc-nmbu wrappers (jobs 1348920 / 1348921; 5.9 / 5.8 min, MaxRSS 828 /
929 MB).

| | pair-only (pre-#23) | size-matched (#23) |
|---|---|---|
| median `z` by size 3→8, root | -0.51, -0.16, 0.53, 1.63, 3.47, 6.32 | -0.16, -0.28, -0.28, -0.28, -0.23, 0.21 |
| median `z` by size 3→8, leaf | -0.75, -0.29, 0.46, 1.96, 5.10, 6.82 | -0.19, -0.19, -0.29, -0.12, 0.31, 0.49 |
| `cor(z, mean_degree)` root / leaf | -0.041 / +0.116 | -0.096 / -0.044 |
| usable `z` | all | all (14597 / 13190) |

- Within a size class the ranking is unchanged: Spearman(old `z`, new
  `z`) ≥ 0.989 at every size. Across sizes it is 0.82 root / 0.76 leaf,
  and the top 100 by `z` went from 99-100% size 7-8 to a spread over
  every size (root 1/10/17/33/28/11 for sizes 3-8; 40 of the old top 100
  survive).
- **What remains is a spread effect, not a location bias.** The share of
  |z| > 1.96 rises with size, 0.12 → 0.67 root and 0.10 → 0.54 leaf, and
  it does so symmetrically (root: z > 1.96 goes 0.09 → 0.35, z < -1.96
  goes 0.03 → 0.32). Cause: `null_sd` still shrinks as 1/sqrt(E)
  (0.041 → 0.014), but the IQR of the observed gap does not (0.067 →
  0.080), because a clique's edges all belong to one HOG and share its
  conservation level — they are not E independent draws. So `z` across
  sizes conflates effect with evidence: a size-8 clique reaches |z| >
  1.96 on a gap a size-3 clique cannot. Rank on `z` *within* size, or on
  the gap, which is reported for exactly this reason. Not acted on.

### S7: coherence under the matched null (#20, option 2 measured and rejected)

`s7_coherence.R` rebuilds the (species pair, clique size) pools with the
package's own helpers, draws the same null as `clique_intensity_test()`
— bit-identical: `cor(z_int, s6 z) = 1`, `max |null_mean diff| = 0` —
and scores coherence (`gm / mean`) on the same draws. Tables
`s7_root_matched/`, `s7_leaf_matched/`.

| | root | leaf |
|---|---|---|
| raw coherence IQR | 0.00317 | 0.00312 |
| coherence `z` IQR (range) | 0.91 (-6.8..3.8) | 0.97 (-6.6..3.3) |
| median coherence `z` | +0.79 | +0.75 |
| share z_coh > 1.96 / < -1.96 | 0.065 / 0.014 | 0.078 / 0.025 |
| share \|z_int\| > 1.96, same cliques | 0.243 | 0.229 |
| Spearman(z_coh, clique size) | 0.43 | 0.41 |
| Spearman(z_coh, z_int) overall | 0.05 | -0.01 |
| Spearman(z_coh, mean_degree) | -0.001 | +0.065 |

Standardising gives coherence a range, but the range is clique size:

- median z_coh / sqrt(E) is 0.34-0.38 at every size (root). The gap
  `observed - null_mean` is a near-constant +0.0026 — positive for 82% /
  79% of cliques, i.e. a clique's edges are slightly more homogeneous than
  random same-pair, same-size edges, a HOG-level conservation level — and
  `null_sd` shrinks as 1/sqrt(E), so z_coh ≈ 0.35·sqrt(E) + noise.
- within a size class, Spearman(z_coh, raw coherence) is 0.98-0.999: the
  standardised statistic is the raw ratio re-scaled, and the raw ratio is
  the one with no range. Within size, z_coh is independent of z_int for
  small cliques (-0.04 at size 3) and largely z_int for large ones (+0.62
  at size 8).

**Recommendation: retire `coherence` (option 1).** The matched null was
the cheapest rescue and it does not discriminate; option 3 (a
dispersion-preserving weight) would be a second weight scale beside the
ensemble probability, and nothing downstream reads coherence. Not acted
on.

### S8: which alternative hypothesis for power (#22) — measured, and a retraction

`s8_alt_power.R`, all 28 pairs × both directions × both tissues,
recomputed in one process from `networks_log` (jobs 1349036 / 1349037,
9-10 min, 6-7 GB, `lib-main-0a54a27`). Guard: the script's recomputation
of the current power equals the package's `power` column (`max diff =
0`) before either alternative is read. Tables `s8_alt_power/s8_*.tsv`.

**Retraction 3.** The mechanism posted on #22 earlier that day — that
the current alternative `f0 * min(k, m)` sits *below* the null mean
`k * m / (N - 1)` for hubs — is false: median **0.0%** of tested rows
(max 3.0% root, 0.3% leaf). It was a plausible derivation posted before
it was measured. #16's mechanism stands: `x*` outruns `f0 * n`, and the
called fraction `x / m` rises with degree (+0.39 root, +0.55 leaf), so
one median fraction under-serves hubs and over-serves low-degree genes.

Medians over 56 pair-directions (reference values `f0` 0.06-0.13,
`rho0` 2.1-3.2, `e0` 0.05-0.08):

| alternative | rho(power, degree) root / leaf | underpowered all / lowest decile / highest decile, root | leaf |
|---|---|---|---|
| fraction (current) | -0.30 / -0.28 | 0.110 / 0.041 / 0.297 | 0.044 / 0.004 / 0.119 |
| fold `rho0 * k * m / N` | +0.998 / +0.998 | 0.079 / 0.785 / 0.000 | 0.050 / 0.495 / 0.000 |
| excess `E0 + e0 * (n - E0)` | +0.48 / +0.48 | 0.005 / 0.021 / 0.000 | 0.001 / 0.002 / 0.000 |

The choice decides who `underpowered` flags: hubs, low-degree genes
(almost deterministically under fold), or nobody. The degree-flatness
test among called pairs (`x / m` +0.39 / +0.55, fold -0.12 / +0.16,
excess +0.21 / +0.40) cannot decide it, because the called set is
selected through `x*` in opposite directions for fraction and fold, so
each carries a selection term of the sign observed. It is a modelling
choice, as #22 said. Recommendation posted: fold, on consistency with
`effect_size` and the intensity weights; excess would empty the flag
and argues for retiring it instead. Decision open.

### S9: the shipped `power` after #27 (closes #16)

`s9_power_check.R` on `lib-main-de1c88d` (main after #27), all 28 pairs
per tissue, job 1349050 (20 min, 6 GB). Guard: the package's `power`
column equals the fold formula recomputed from the comparison frame,
`max |diff| = 0` on every pair.

| | fraction (#16, `lib-matched-fbe547b`) | fold (`lib-main-de1c88d`) |
|---|---|---|
| Spearman(`power`, degree), root | -0.658..-0.187, 28/28 negative | +0.863..+0.990, median +0.95, 28/28 positive |
| Spearman(`power`, degree), leaf | -0.668..-0.253, 28/28 negative | +0.950..+0.987, median +0.98, 28/28 positive |
| underpowered, lowest decile, root / leaf | 0.04 / 0.00 | 0.81 / 0.62 |
| underpowered, highest decile, root / leaf | 0.30 / 0.12 | 0.00 / 0.00 |
| underpowered, all, root / leaf | 0.11 / 0.04 | 0.10 / 0.07 |

`rho0` per direction 2.1-2.9. #16 closed with this table.

### #16 settled on mechanism, open on meaning (historical; closed by S9)

The power/degree inversion **reproduces on current code**: negative on
**56 of 56** species pairs (root -0.658..-0.187, leaf -0.668..-0.253,
median ~-0.36), while `jaccard`/degree is positive (+0.38..+0.60).

It is **intrinsic, not a reconstruction artifact**. `urn_isolate.R` on 6
root pairs: exact vs reconstructed urn give `max |power difference| =
0.000e+00`, zero `NA` either way. So #21 is a robustness fix (it removes a
failure mode that *can* fire) rather than a numerical correction.

What remains open is what `underpowered` should *mean*: it flags **hub**
genes, the opposite of the low-degree failures #12 introduced it for.
That turns on whether the reference conserved fraction is a fraction or a
count — **#22**.

### Orion artifacts

Under `validation-2026-09-17/`:

- libs: `lib-main-de1c88d` (main after #27, the current one),
  `lib-main-0a54a27` (main after #24: #23 pools + #24 flag), `lib-matched-fbe547b` (**#19 only**, despite its earlier
  "#19+#23" label), `lib-urn-57bb90d` (#21). `lib-null-a997ce7` is
  **pre-fix** — it has the buggy hog-only shuffle; do not reuse.
- scripts: `s6_matched.R` + `.slurm` (matched null; needs no networks),
  `s7_coherence.R` (coherence under the matched null; recomputes the
  package's intensity null in-process and refuses to continue unless it
  matches the stored s6 table), `install_main.sh <sha>` (installs
  `pkg/main-<sha>/` and asserts `match_clique_size`),
  `power_degree_remeasure.R` + `.slurm` (**reads stored `power`** — see
  retraction 2), `urn_isolate.R` + `.slurm` (recomputes both urn paths),
  `install_matched.slurm`, `install_urn.slurm`. From 2026-09-22 new
  scripts live in the laptop's gitignored
  `prepare_data/validation-2026-09-17/` and go up through the ohpcc-nmbu
  wrappers (`hpc.env` at the repo root, gitignored); their logs are in
  `../slurm_logs/<jobid>.out`, not `logs/`.
- results: `s9_power/` (`s9_power_check.R`: shipped power vs degree
  after #27), `s8_alt_power/` (`s8_alt_power.R`: power under three
  alternatives, per pair-direction), `s6_root_matched23/`,
  `s6_leaf_matched23/` (#23 pools),
  `s7_root_matched/`, `s7_leaf_matched/` (coherence + intensity null
  moments and `z`), `s6_root_matched/`, `s6_leaf_matched/` (pair-only),
  `power_degree_root.tsv`, `power_degree_leaf.tsv`, `urn_isolate_root.tsv`.

### Still open

- **#22** — fraction vs count; decides what `underpowered` means.
- **#20** — closed: `coherence` retired in #26 (`3a99584`) after the
  S7 measurement above.
- **#25** — closed by #28: `gap` and `n_edges` columns, docs say to
  compare `z` within a size and `gap` across sizes. Design-effect and
  empirical standardisation were rejected: the between-HOG spread they
  would absorb is the signal.
- **#22, #16** — closed: Martin chose fold (biologically: `rho0` is
  how many times more often than chance a pair's partners are shared);
  #27 shipped it and S9 confirmed the shipped column.
- **Two clique classifiers** — Martin (2026-09-22): having both
  `classify_cliques()` (species graph) and `classify_gene_cliques()`
  (gene graph) is awkward. The `classify_species_cliques()` rename is
  shelved; the direction is consolidation, not clearer names for two.
  Not started.

## 11. Second dataset, 2026-09-22: the improved package on the EVOTREE wood data

Rodriguez et al. 2026 (Nat Commun, doi 10.1038/s41467-026-75624-2) —
the paper the gene-graph taxonomy follows. Six species (aspen, birch,
cherry; Norway spruce, Scots pine, lodgepole pine), 65-106 wood
cryosection samples each, HOGs from OrthoFinder N1. Their input data and
outputs are on Zenodo (10.5281/zenodo.21025760), unpacked on Orion at
`rcomplex-testrun/evotree/`; scripts in `prepare_data/evotree/`
(`e1_reproduce.R`, `e2_taxonomy.R`), run with `lib-main-f176fb9`.
Their published tiers are the **baseline**; the question was what the
improvements change.

**rcomplex reproduces their test exactly.** Same 444,213 tested pairs,
overlap counts identical on all 15 pairs, all 196,183 of their calls
recovered plus 77 borderline extras from the self-excluded urn, q-values
inside their 3-significant-digit rounding. Every difference below is
therefore due to the improvements, not the network or the test.

**Pooideae findings replicate, and are stronger here.**

| | Pooideae (root / leaf) | wood |
|---|---|---|
| `rho0` per direction | 2.1-3.2 | 2.1-3.1 |
| Spearman(`power`, degree) | +0.95 / +0.98 median | +0.98 median, +0.93..+1.00 |
| underpowered: all / lowest decile / highest | 0.10 / 0.81 / 0.00 (root) | 0.20 / 1.00 / 0.00 |
| median z by clique size (matched null) | −0.28..+0.21 (sizes 3-8) | −0.10, −0.16, −0.14, +0.29 (sizes 3-6) |
| IQR of gap by size | 0.067-0.080 | 0.060-0.088 |
| IQR of z by size | 1.6 → 5.6 | 2.0 → 4.1 |
| share |z| > 1.96 by size | 0.12 → 0.67 | 0.17 → 0.53 |

So #23 (no size gradient in z) and #25 (spread grows with size while
the gap does not) both hold on a second dataset, and the `rho0` default
of "about 2.5" is not a Pooideae number.

**Gene-graph taxonomy, HOG level, best tier per HOG (16,892 HOGs).**

| tier | theirs | baseline (power-blind) | with power | note |
|---|---|---|---|---|
| complete | 940 | 941 (940 shared) | 941 | exact |
| partial_significant | 538 | 537 (537) | 537 | exact |
| lineage_specific dicot / conifer | 1003 / 1975 | 1004 / 1977 (all shared) | same | exact |
| differentiated | 293 | 369 (293 shared) | 354 (282) + 11 → `underpowered` | needs `alpha_graph = Inf` |
| partial_present | 620 | 211 (211) | 226 (226) | 333 of theirs are `unclassified`: sixth species tested and failed |

The only tier the package refuses is partial_present when the missing
species was *tested* and not called (`missing_reason = tested_ns`);
that is the absent-versus-negative-evidence rule, and power rescues 15
of those (13 + 2). Power gating moved 114 HOGs in all: 86 unclassified
and 13 differentiated to `underpowered`, 13 unclassified and 2
differentiated to partial_present. Trap: `gene_clique_graph(alpha_graph
= 1)` drops q = 1 rows (documented; `Inf` is the unfiltered graph) —
with 1 the differentiated tier found 44 of 293.

**Species-graph classifier (`classify_cliques()`, lineage as trait):**
of 416 `differentiated` HOGs **319** carry the `underpowered` flag; of
7,817 `trait_specific`, 648. Three quarters of the divergence calls on
wood rest on a cross-lineage edge whose gene could not have been called
at typical enrichment. That is the flag doing what #12 introduced it
for, on data the package was not tuned on.

Complete cliques on the strict graph: 50,750 vs their 50,697.

**E3: the 333 refused partial-present HOGs are negative evidence.** In
all 8,944 of their five-member cliques the sixth species was tested
against the members (`missing_reason = tested_ns`); the maximum power
among those tests is 1.0 in three quarters of cliques (minimum 0.33),
the smallest q is 0.82 at the median, and `min_power` rescues 12 HOGs at
0.8, 17 at 0.9, 45 at 0.99. The rule stands. The missing species is
aspen in 3,391 of those cliques, birch and cherry about 2,200 each, the
conifers 187-604; aspen's data come from a different study.

## 12. Pooideae under 0.3.0 with power, 2026-09-22 (P2)

Analytical-path edge tables rebuilt with `lib-main-f176fb9`
(`s2_edges_log.R` → `s2_<tissue>_v030/`, alpha 0.1, `pval_combine =
"max"`, randomized pi0), then both classifiers with and without the
`power` column (`p2_classify_v030.R` → `s4_<tissue>_v030/`). Power is
the fold-enrichment power of #27. Mean edge power 0.94 root / 0.94 leaf;
9.4% / 7.8% of edges below 0.8.

**Species graph (`classify_cliques()`, trait = annual/perennial,
`min_species = 3`).** No HOG is `differentiated` in either tissue; the
divergence signal of the study is entirely `trait_specific`.

| | root | leaf |
|---|---|---|
| complete / partial / trait_specific / unclassified (power-blind) | 346 / 7745 / 297 / 8017 | 440 / 6609 / 384 / 8502 |
| trait_specific flagged `underpowered`, `min_power = 0.5` | 17 (6%) | 0 |
| `min_power = 0.8` (default) | **130 (44%)** | **118 (31%)** |
| `min_power = 0.9` | 192 (65%) | 209 (54%) |

So at the default gate roughly a third to a half of the trait-specific
HOGs rest on a cross-trait edge whose gene could not have been called at
typical enrichment; the flag qualifies the call, it does not remove it,
but any downstream claim about those HOGs needs the flag beside it.

**Gene graph (`classify_gene_cliques()`, 0.1 / 0.9 / Inf, lineage =
trait; the unfiltered `Inf` graph was dropped, see below):** 372,618
cliques over 11,774 HOGs (root), 249,340 over 11,233 (leaf).

| best tier per HOG | root | leaf |
|---|---|---|
| complete_conserved | 346 | 440 |
| lineage_specific | 0 | 4 |
| partial_significant | 483 | 324 |
| partial_present | 228 | 329 |
| differentiated | 0 | 0 |
| unclassified | 10,717 | 10,136 |

Power gating moves **nothing** at `min_power` 0.5 or 0.8 and two HOGs
(unclassified → partial_present) at 0.9, in either tissue. The tiers
power can gate, `lineage_specific` and `differentiated`, are empty here:
`lineage_specific` needs every outside species absent or untested, and
in Pooideae the other trait group is present and tested. So the trait
signal of this study is visible only to `classify_cliques()`'s
`trait_specific` (297 / 384 HOGs, 31-44% flagged), and the gene graph
has no tier for "one trait group conserved, the other present but not
co-conserved". That is the decisive input for the classifier
consolidation note.

**With the `trait_specific` tier (PR #29, `lib-main-a305ac4`):** 93 root /
97 leaf HOGs, against the species graph's 297 / 384 at `min_species = 3`
and 129 / 116 at `min_species = 4` (overlap 39 / 35). 226 / 307 of the
species-graph calls use three of the group's four species; 45 / 33 are
outranked by `partial_significant`; 45 / 44 remain unclassified. The tier's HOGs are 70 annual + 23
perennial (root) and 81 + 16 (leaf). Of the species graph's extra calls
(282 / 364), the fourth group member is absent in 96 / 122, tested and
rejected in 88 / 157, and in 90 / 76 the group triangle sits inside a
larger mixed clique with one species of the other group. See the
consolidation note, step 0.

Trap: `gene_clique_graph(alpha_graph = Inf)` on all Pooideae HOGs ran
for an hour on one core at 39 GB without finishing (eight species, up to
ten paralogs each); it took 0.9 min on wood. Restrict it to leftover
HOGs, or skip it when no HOG is differentiated.

## 13. Rotation gene-set pilot, 2026-09-23 (Dørum et al. 2009, limma::roast)

Question: is the rotation test a useful layer for rcomplex? Pilot on
HOG:0022829, an annual-specific root clique (Bradi1g51520 with its
barley, B. mexicanum and V. bromoides orthologs). Per species and
tissue, the focal gene's top-100 co-expressed genes as the set, a
linear trend in log(day) over the five time points (four replicates
each) as the design, `roast()` with 9,999 rotations, a random set of
100 as control, and the single-gene trend beside it. Script and table:
`prepare_data/validation-2026-09-17/rotation_pilot/`.

- The four annual clique genes are age-stable in root: single-gene
  |t| 0.5-1.3 (p 0.22-0.59), set p 0.22-0.68, active proportion
  0-0.34. The HOG's perennial copies (8 in B. sylvaticum, 14 in the
  tetraploid H. jubatum, 8 in the B. mediterraneum assembly, 2 in
  F. pratensis) are often strongly age-responsive: |t| up to 5.3
  (Brasyl.7G013700, down) and 7.7 (HJUBATUM 4H_1G00353160, up), with
  set p down to 1e-4. So the conserved annual programme is not an
  age-tracking one, while the expanded perennial family is.
- Caveat that limits the set p-values: random 100-gene sets reach
  p < 0.05 in 19 of 61 tests, because much of the root transcriptome
  tracks age. `roast()` is self-contained (null = no association), so
  set p-values are inflated for any set here; a competitive test or a
  random-set reference is needed before a set p means "this set more
  than others". The single-gene contrast above does not depend on it.
- The neighbourhood is defined on the same samples it is tested on,
  which is circular; a cross-species conserved neighbourhood would be
  the honest set.
- Rotation itself buys nothing on this design: 20 samples per tissue
  give ample permutations. It would on EVOTREE's three trees per
  species, which is where a sample-level clique test belongs.

**With `camera()` (competitive) and `fry()` added** (same sets, same
design; table `rot_test_HOG0022829_camera.tsv`): `fry` reproduces
`roast` throughout. `camera` on the random control sets gives 0 of 59
below 0.05 against roast's 22, so the competitive null is the
calibrated one here. On the focal sets `camera` is extreme for the
age-responsive perennial copies (p down to 1e-57) and, unlike the
self-contained tests, also significant for two annual clique
neighbourhoods (BDIS 7.6e-5, VBRO 2.7e-7, both "down") although their
single-gene |t| is about 1.2 and roast/fry give p 0.2-0.3. That is the
competitive question showing: the annual sets are *flatter* than a
transcriptome that drifts with age, so they stand out from the
background without moving themselves. Read together: the annual
programme is age-stable in absolute terms and unusually so relative to
the genome; the perennial family is age-responsive on both readings.

Not a package feature; kept as the reference run for that decision.
