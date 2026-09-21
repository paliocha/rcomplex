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
