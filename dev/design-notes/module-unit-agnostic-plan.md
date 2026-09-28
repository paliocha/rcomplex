# Plan: decouple module testing from module detection (branch `feature/module-unit-agnostic`)

2026-09-28. Martin asked whether the module engine is worth keeping. The
answer was to split it: keep the testing half, retire per-species community
detection as the *source* of modules. This plan covers the first two steps.
The anchored-regulon engine (the replacement unit) is out of scope until it
passes a cross-fitted validation (see "Later").

## Why

Evidence (design note `module-engine-redesign.md`, sections 11.1-11.14, and
the 2026-09-28 pathway-convergence audit):

- `detect_modules()` modules are seed-stable (ARI 0.97) but do not replicate:
  split-half ARI 0.06-0.15 on n = 20 leaf; 8-11 Leiden modules on shuffled
  expression. The best consensus / multilayer engine capped at ARI ~0.25.
  Cause: at n = 20 the correlation graph is a random geometric graph under
  the null, so module structure exists on noise.
- The K = 1 test compares against degree-preserving rewiring, which destroys
  that geometry; on this argument it should reject K = 1 on shuffled
  expression. **Not yet measured** -- WP2 measures it before anything changes.
- The testing half (`module_preservation()`, `preservation_paired()`,
  `module_correspondence()`, `preservation_matrix_test()`,
  `tag_permutation()`) is calibrated and separates discovery (reference) from
  test (other species), but only accepts `detect_modules()` output.

What the testing half reads (main, 67d256f): `module_preservation()`
(R/module_preservation.R:245-348) and `module_correspondence()` (:1182-1204)
read only `$modules` (named vector gene -> label, NA = unassigned) and
`$module_genes` (named list label -> genes). `identify_module_hubs()`
(R/modules.R:975) also needs `$graph`.

## Outcome

1. `as_modules()` turns any partition -- anchored regulons, root cores,
   curated pathways -- into the object the testing half takes, so it can be
   tested without `detect_modules()`.
2. `detect_modules()` is marked superseded: documented limits, a once-per-
   session notice, and `test_k1 = FALSE` by default **if** WP2 shows the K = 1
   test rejects on shuffled data. Nothing is removed; nf-rcomplex keeps working.

No new dependency. No change to any statistic.

## Work packages

Branch: `feature/module-unit-agnostic` from `main` (done). WP1 and WP2 are
independent (disjoint files) and may run in parallel; WP3 needs WP2's
number; WP4 needs WP1 and WP3. Each WP ends green on its acceptance command
and with one commit (no Co-Authored-By trailer, last line the
`Claude-Session:` line).

### WP1: `as_modules()` (new file `R/as_modules.R`, new test file)

```r
as_modules(x, genes = NULL, min_size = 1L)
#  x: named vector (gene -> module label; NA = unassigned), or a named list
#     (label -> character vector of genes), or a detect_modules() result
#     (returned unchanged, so the call is idempotent).
#  genes: optional universe; genes not in x get NA (unassigned).
#  min_size: modules smaller than this become unassigned.
#  returns list(modules = <named chr, labels as character>,
#               module_genes = <named list>, n_modules, method = "external",
#               params = list(min_size = min_size))
```

- A list input whose sets overlap errors, naming up to five shared genes and
  their modules: the testing half projects one label per gene, so it needs a
  partition. (Overlapping regulons: the caller splits them into disjoint
  batches. Supporting overlap is out of scope.)
- Empty labels, duplicated gene names in a vector, and non-character genes
  error with the offending value.
- Labels are stored as character (detect_modules uses integers; the testing
  half already `as.character()`s them, R/module_preservation.R:317).
- Update the two "must be output from detect_modules()" messages in
  R/module_preservation.R (:247, :1183) to "must be a module assignment from
  detect_modules() or as_modules()". The hubs check (R/modules.R:977) keeps
  its message: hubs need `$graph`, which `as_modules()` does not build.
- Roxygen with one example: a named list of two gene sets fed to
  `module_preservation()`.

Acceptance (`devtools::test(filter = "as-modules")` plus the existing
preservation tests unchanged):
- `module_preservation(as_modules(m$modules), ...)` is `identical()` to
  `module_preservation(m, ...)` for a `detect_modules()` result `m` on the
  existing preservation fixture (same seed); the same for
  `module_correspondence()` and `preservation_paired()`.
- list input and vector input describing the same partition give identical
  objects; `as_modules(as_modules(x))` is identical to `as_modules(x)`.
- overlap, duplicate-gene, empty-label errors fire with the named values.
- `min_size` turns small modules into NA and drops them from `module_genes`.

### WP2: measure the K = 1 test on null data (script, no package change)

`dev/bench/k1_null_check.R`: for (a) seeded `rnorm` expression, 20 samples,
and (b) Pooideae leaf per species with sample labels shuffled within each
gene (the shuffled-expression null used everywhere else; skip if
`prepare_data/data/` is absent), 2,000 top-variance genes, `compute_network()`
at the package defaults, then `detect_modules(resolution = seq(0.5, 2, 0.5),
test_k1 = TRUE, n_perm_k1 = 100)`; 10 seeds for (a), all 8 species for (b).
Report per run `has_structure`, `p_value`, `n_modules`; write a TSV to
`dev/bench/k1_null_check.tsv` and print the rejection rate.

Acceptance: the script runs end to end locally in under 30 min and prints the
rejection rate for (a) and (b). **Decision rule for WP3:** if the rejection
rate on null data exceeds 0.2 (alpha is 0.05), flip the default; otherwise
keep `test_k1 = TRUE` and only document.

### WP3: soft-deprecate detection (R/modules.R, docs)

- `detect_modules()` roxygen: a "Status: superseded" paragraph at the top of
  `@description`, in plain words: per-species modules at small sample sizes
  are reproducible across seeds but not across independent samples (numbers
  above); use `as_modules()` to test gene sets from elsewhere. Same paragraph,
  shortened, on `identify_module_hubs()` / `classify_hub_conservation()` /
  `characterize_hubs()` (hubs are defined on those modules).
- One notice per session from `detect_modules()`:
  `rlang::inform(<one sentence + pointer to ?as_modules>, .frequency = "once",
  .frequency_id = "rcomplex_detect_modules")`. `inform`, not `warn`: the
  function still works and must not fail `R CMD check` examples or
  `expect_no_warning()` tests.
- If WP2 triggered: `test_k1 = FALSE` default in `detect_modules.default`,
  `detect_modules.rcomplex` and `detect_modules_consensus`; the `test_k1`
  param doc states the measured null rejection rate. Tests that relied on the
  default set it explicitly (grep: tests passing no `test_k1` to a multi-
  resolution call).
- `NEWS.md` bullet under `# rcomplex 0.3.1`; CLAUDE.md: one paragraph in Key
  Design Decisions ("Detection is superseded; testing takes any partition")
  and the Project Overview module bullet trimmed to match. README: one line
  pointing module users to `as_modules()` (keep it one line; PR #32 rewrites
  the README).

Acceptance: `devtools::test()` green; `lintr::lint_package()` clean;
`R CMD check` on the built tarball "Status: OK"; the notice appears once in a
fresh session calling `detect_modules()` twice (`testthat::expect_message`
then `expect_no_message`, with `rlang::reset_message_verbosity()` between
test files if needed).

### WP4: PR

Push, open the PR (no merge). Body: the evidence table, WP2 numbers, the
default decision, and a **nf-rcomplex compatibility** section:
- `bin/run_species_modules.R:200` passes `test_k1` explicitly -> unaffected
  by the default change; it will print the notice once per process.
- `bin/run_module_comparison.R:171-382` calls `compare_modules()` /
  `classify_modules()`, which 0.3.0 removed: that script is already broken
  against rcomplex >= 0.3.0 (pre-existing, not caused here; flag for the
  nf-rcomplex side).

## Out of scope (and why)

- Removing `detect_modules()` or the hubs: nf-rcomplex depends on them;
  removal is a later release once the pipeline has switched.
- A shuffled-expression K = 1 null: it needs the expression matrix, which
  `detect_modules()` does not take; not worth building for a superseded
  function.
- Overlapping modules in the testing half.
- Porting the anchored-regulon engine (below).

## Later (not this branch)

Anchored regulons become the package's module unit only after a cross-fitted
validation: build programs on half of each species' samples (Pooideae: two
replicates per time point; wood: half the trees) and measure preservation and
trait contrasts on the other half. The 2026-09-28 audit showed same-data
readings are inflated in scope species (attachment +0.25-0.38 selection vs
+0.05-0.10 held-out). With `as_modules()` in place, that validation can use
`module_preservation()` directly.

## Verification

```bash
Rscript -e 'devtools::document()'
R CMD INSTALL .
Rscript -e 'devtools::test(filter = "as-modules|module")'
Rscript dev/bench/k1_null_check.R
Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint_package()'
R CMD build . && R CMD check --no-manual rcomplex_0.3.1.tar.gz
```
