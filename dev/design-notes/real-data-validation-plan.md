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
