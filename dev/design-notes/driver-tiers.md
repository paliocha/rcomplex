# Driver tiers, ortholog reader, release (WP-g to WP-j)

2026-10-14. Branch `docs/driver-tiers`, base `16f2879` (main after
WP-f). Follows `clique-block.md`. Martin: plan items 1 to 5 of the
open list; the fourteen leaf families were a side quest.

## 0. Verdict and order

- WP-g. `rcomplex()` builds one gene graph, at `alpha`. Two of the
  published tiers need cliques from looser graphs, so the driver never
  reports them. It gets the three-pass workflow that
  `prepare_data/evotree/e2_taxonomy.R` ran by hand. Biology first:
  this changes what the tables say.
- WP-i. `read_orthologs()` misses the EVOTREE HOG file over one header
  word. One line, one fixture.
- Release. Tag `v0.4.0` on `f693b33`; `main` becomes 0.5.0. After
  WP-g and WP-i, so the release carries them.
- Report. Rerun the wood report (`e5_*.R`) on the released code and
  republish the artifact. A run, not a work package.
- WP-j. EVOTREE's "best unique cliques" thinning. Gated on Martin.

Each WP is one branch off `main`, one builder, one reviewer, one PR.
No Claude attribution anywhere. Rules of `CLAUDE.md` apply: 80-column
lines, STE-lite prose, surface budget (31 exports, 8 formals).

## 1. WP-g: three graphs in the driver

### What the published workflow did

EVOTREE (`CliquesWorkflow/`, local copy
`~/.cache/rcomplex/evotree_code/`) runs three clique passes per
orthogroup, each on its own edge set:

1. `CompleteCliques.R`: pairs with `Max.p.Val < 0.1`.
2. `PartialCliques.R`: pairs with `Max.p.Val < 0.9`, orthogroups with a
   complete clique excluded.
3. `DifferentiatedCliques.R`: every tested pair, only orthogroups with
   no complete or partial claim (`leftoverOG`).

A clique maximal on one graph need not be maximal on another, so the
tiers that read the loose graphs need cliques from them.
`classify_gene_cliques()` was built for this: it takes a list of
`gene_clique_graph()` runs, offsets their ids, adds a `run` column and
takes `hog_class` as the earliest tier any run reaches. `e2_taxonomy.R`
lines 49-85 are the reference: strict at 0.1, loose at 0.9, then
all-tested (`alpha_graph = Inf`) on hogs whose best tier so far is not
`complete_conserved`, `partial_significant` or `partial_present`.

### What the driver does today

`R/rcomplex.R` lines 151-160: one `gene_clique_graph()` at
`alpha_graph = alpha`, then `classify_gene_cliques()` with its default
`alpha_graph = 0.9`. Nothing in `cliques` can be `partial_significant`
(needs a clique edge in `[alpha, 0.9)`) or `differentiated` (needs
every pair of the home tested, rejected pairs included). On wood the
hand run gave 293 `differentiated` and 537 `partial_significant` hogs;
the driver gives 0 and 0, and those hogs sit in `unclassified`.

### Spec

`cliques(e)` in `rcomplex()` becomes three passes:

1. `strict <- gene_clique_graph(e, min_size, alpha_graph = alpha)`.
2. `loose <- gene_clique_graph(e, min_size, alpha_graph = 0.9)`.
3. Classify `list(strict, loose)`. Hogs whose `hog_class` is
   `complete_conserved`, `partial_significant` or `partial_present`
   are claimed. `tested <- gene_clique_graph(e[!claimed, ], min_size,
   alpha_graph = Inf)`.
4. Classify `list(strict, loose, tested)` with `alpha_call = alpha`,
   `alpha_graph = 0.9`, `clades`, `min_power` as today. That is
   `res$classification`; it carries `run`.

No new driver argument. 0.9 is the classifier's default and the
published value. `min_size` stays `min(3L, length(sp))`. The `null =
TRUE` path reuses the same closure, so the null gets the same three
passes.

Result tables. `res$cliques`, `res$members`, `res$twins` and
`res$capped` are the three runs stacked, each with a leading `run`
column (1, 2, 3) and `alpha_graph` column (`alpha`, 0.9, `Inf`).
`clique_id` stays the run's own id, as `classification$clique_id`
does, so `(run, clique_id)` joins the tables. `write_rcomplex()` needs
no change; the columns flow through. `print()` gains one line naming
the three graphs and their clique counts. `summary()` is unchanged.

Roxygen of `rcomplex()`: `@details` gets one paragraph on the three
passes and which tiers each one feeds; `@return` documents `run` and
`alpha_graph`. Quickstart: the tier table it shows must come from a
run that has all tiers; update the prose where it says which tiers
appear.

Cost. The budget bounds every hog (`max_cliques = 1e5`), and pass 3
runs on leftover hogs only. The builder measures wood and full leaf
(`~/.cache/rcomplex/E_full.rds`, 2.66 M rows; wood edges from
`prepare_data/evotree/`) and reports wall time, peak memory, clique
count per pass and `nrow(capped)` per pass. If full leaf exceeds
twice the current driver time (8 min) or 32 GB, stop and report before
changing the spec.

### Tests

1. Fixture with one hog whose six genes form a complete clique at
   `q < 0.1` in three species and a fourth species joined at
   `q = 0.5`: classification has a run-2 clique
   `partial_significant`; run 1 has only the three-species clique.
2. Fixture with two clades, each a complete clique within, every cross
   pair tested at `q = 1`: a run-3 clique `differentiated`;
   `hog_class` is `differentiated`; that hog has no run-1 or run-2
   clique over all species.
3. A hog claimed in pass 3's exclusion (complete clique) has no run-3
   rows.
4. `(run, clique_id)` of `classification` joins `cliques` and
   `members` with no misses; `write_rcomplex()` and a reread
   round-trip the `run` column.
5. `rcomplex(null = TRUE)` on the package example returns
   `classification_null` with a `run` column.
6. Wood: the driver's `hog_class` counts equal a hand three-pass run
   with the same `edges`, `species`, `clades` and thresholds, tier by
   tier. The builder reports the counts next to the WP-f baseline
   (complete 938 at species level; gene graph complete 940,
   `partial_significant` 537, `differentiated` 293 in the e2 run).

### Acceptance

`R CMD check` OK on the tarball, all tests pass, lint and
`dev/sharpen-prose.R` clean. NEWS bullet states the behaviour change
with the wood counts before and after. CLAUDE.md "Two clique
backends" gains one sentence: the driver runs three graphs. Reviewer
checks that `capped` of pass 3 is not polluted by hogs that pass 1
resolved (the exclusion is what keeps it clean).

## 2. WP-i: OrthoFinder header variants in `read_orthologs()`

### Today

`R/orthologs.R` line 54: an OrthoFinder file is one whose header has
`HOG` or `Orthogroup`. The EVOTREE file
(`Orthogroups_20240823_clean_N1.tsv.gz`) has `OrthoGroup` (capital G),
then `OG`, `Gene Tree Parent Clade`, then one column per species with
genes separated by `", "`. `read_orthologs()` stops with "Cannot read
the file format from the header". `e1_reproduce.R` lines 60-80 reshape
it by hand.

### Spec

Recognise the hog column case-insensitively: `HOG`, `Orthogroup`,
`OrthoGroup`, and `N<k>.HOG`-style names, in both the format sniff and
`.orthofinder_long()`. Keep the skip list (`OG`, `Gene Tree Parent
Clade`). `", "` already works through `trimws()`. Gene id suffixes
(`.p1`, `.1.p1`) are the user's mapping, not the reader's; say so in
`@details` in one sentence.

### Tests

1. Fixture `inst/extdata/N1_variant.tsv` (three hogs, three species,
   `OrthoGroup` header, `", "` separators, `NA` cells): `read_orthologs()`
   returns the long table; row count equals the number of genes.
2. `format = "auto"` picks `orthofinder` for it; `format = "long"` on
   it errors naming the missing columns.
3. The six wood species of the EVOTREE file read with
   `species = c("Populus_tremula", ...)` give the hog set that
   `e1_reproduce.R` uses (builder reports counts, no committed data).

### Acceptance

Check OK, tests pass, prose clean. NEWS bullet. One sentence in
`@details`.

## 3. Release: tag 0.4.0, bump to 0.5.0

0.4.0 was released on 2026-10-12 as `f693b33` (merge of #65) and never
tagged; the last tag is `v0.3.2`. Since then the clique tables changed
shape (`cliques` + `members`), the null was redefined (hog shuffle
within expression deciles), the copy cap went and blocks arrived.
Those break 0.4.0 scripts, so the next release is 0.5.0.

Steps, after WP-g and WP-i merge:

1. `git tag -a v0.4.0 f693b33 -m "rcomplex 0.4.0"` and push the tag.
2. On a branch: `DESCRIPTION` Version 0.5.0; `NEWS.md` header
   `# rcomplex 0.5.0` over the development bullets, with a two-line
   lead like the 0.4.0 one; `CLAUDE.md` line 46 tarball name.
3. Check OK, PR, merge, tag `v0.5.0`.

Martin decides the number; 0.5.0 is the recommendation.

## 4. Report refresh

The wood report (https://claude.ai/artifact/J62cYfRu5UdTfMJzVjtgPH,
version 2) shows WP-d numbers: 51,775 cliques, no blocks, no
`differentiated` or `partial_significant`. After the 0.5.0 tag, rerun
`prepare_data/evotree/e5_wood_report_data.R` against the released
library and republish as version 3 (same link). Numbers to carry:
125,790 cliques (WP-f), 347 species-level hogs moved by the subset
fallback, the tier counts from WP-g test 6, and the `calls_null /
calls` ratio from WP-d. Networks and edges are cached under
`prepare_data/evotree/data/` and `prepare_data/cache/`; scratchpad
copies are gone.

## 5. WP-j: best unique cliques (gated)

EVOTREE's `bestUniqueCliques.R` walks each orthogroup's cliques by
mean FDR and keeps a clique when none of its genes is in a kept
clique: 70,458 cliques became 2,145 unique. We report every clique.

If Martin wants that view: a logical column `representative` on
`classification` (not a new export; the surface is at 31). Per hog,
order rows by `mean_effect_size` descending, then `mean_q` ascending,
then `clique_id`; a row is representative when no member gene of it
belongs to an earlier representative row of the same hog. Blocks
(`capped`) count as one gene per species. `summary()` gains a
`representative` count per tier. Test: a hog with two overlapping and
one disjoint clique gives two representatives, the better of the pair
first. Not started; Martin's call.

## 6. Left open

- Twin paralogs adjacent within a species, and copy resolution from
  network evidence (expression only).
- WP15a (star-expansion baseline) unblocked; WP15b held; WP16/17 wait
  on Orion probes. See `sharpen-plan.md`.
- Three unmerged worktrees under `.claude/worktrees/` (README draft,
  blockwise kernel, benchmark script) and one Copilot worktree. Not
  touched.
