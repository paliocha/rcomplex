# Shape of the gene-clique table (issue #66, memory half)

2026-10-08. Branch `docs/clique-table-shape`, base `5f31d6f` (0.4.0 with
the #66 loop fix). Measurements on a 64 GB, 10-core Mac, R 4.6.1,
igraph 2.3.3. Data: the Pooideae leaf edge tables at 8k, 10k, 12k
genes per species and full size (8 species, 18.8k-41.3k genes,
2,664,155 rows). Clades: annual = BDIS HVUL BMAX VBRO, perennial =
BSYL HJUB BMED FPRA, as in the pre-merge driver run.

## 0. Verdict

- The member-row table is 99.6 bytes per row. 56 % of it is the
  per-clique statistics repeated on every member row, 20 % is the
  character `clique_id`. At full size it is about 20 GB.
- The shape is not the only memory problem. Two others are larger.
  The builder's peak comes from materialising the cliques of three
  HOGs (86 % of all cliques). The classifier needs 11 GB and 511 s
  for the 1.45 M cliques at 12k. At full size, extrapolated, that is
  2.8 h and far more than 64 GB. The driver cannot finish at full
  size under any table shape until the classifier changes.
- The clique explosion repeats information. At 12k the 1.45 M cliques
  fall into 7,413 (HOG, species set) groups. Only 2 of those groups
  hold more than one tier. At full size 28.3 M cliques fall into
  41,866 groups (677x).
- Recommendation: shape B (a clique table plus an integer member
  table), built HOG by HOG in disjoint pieces. Add a classifier that
  works on the sparse gene x clique incidence. Report D (tiers by
  cliques, HOGs and (HOG, species set) groups) in `summary()`, computed
  rather than stored. No new argument. Two WPs, about 600 changed
  lines in total. Section 5 has the spec.

## 1. Where the bytes are

Status quo at 12k: 1,452,683 cliques, 10,113,744 member rows (6.96
per clique). Sizes from `object.size()`, which counts each distinct
string once.

| column | MB | B/row | distinct values |
|---|---|---|---|
| clique_id (chr) | 187.5 | 19.4 | 1,452,683 |
| hog (chr) | 77.4 | 8.0 | 3,742 |
| species (chr) | 77.2 | 8.0 | 8 |
| gene (chr) | 78.8 | 8.2 | 20,171 |
| n_members, n_species, n_edges | 3 x 38.6 | 3 x 4 | 6 |
| mean_q, max_q, mean_effect_size, score | 4 x 77.2 | 4 x 8 | up to 1.45 M |
| mean_q_floor | 77.2 | 8.0 | 1 |
| n_cliques_at_q_floor | 38.6 | 4.0 | 1 |
| **total** | **961.0** | **99.6** | |

- Repeated per-clique statistics: 540 MB, 56 %. Stored once per
  clique they take 78 MB. The repetition is the `rep(..., times =
  nm_v)` at `R/clique_gene_graph.R:441-451` and the floor columns at
  117-118.
- `clique_id`: 187.5 MB, 20 %. 77 MB of it is one pointer per row,
  110 MB is 1.45 M distinct strings of about 76 bytes each
  (`paste0(id_prefix, h, "_", seq_len(n_k))`, line 390). A global
  integer id costs 39 MB.
- `species` as a factor: 38.6 MB, half of the character column.
- `gene` is only pointers: 20,171 distinct strings, already shared
  with the edge table. A character column costs 8 bytes per row, not
  the length of the id.

Full size: 28,343,526 cliques, 215,914,224 rows (21.35x the rows,
19.51x the cliques of 12k).

| part | 12k MB | full, extrapolated |
|---|---|---|
| repeated statistics | 540 | 11.5 GB |
| clique_id pointers + strings | 77 + 110 | 1.6 + 2.2 GB |
| hog, species, gene pointers | 233 | 5.0 GB |
| **table** | **961** | **about 20 GB** |

The observed peak was 26.5 GB in 335 s. So the table is about three
quarters of the peak. The rest is the per-HOG transient (section 3.1).

Classification table: 288 MB at 12k, 208 B per clique, 121.5 MB of
it `clique_id`. About 5.9 GB at full size, 3.6 GB with an integer id.

Runtime and peak memory of the status quo (`m1.R`):

| genes / species | cliques | `gene_clique_graph()` | `classify_gene_cliques()` |
|---|---|---|---|
| 8k | 41,329 | 1.5 s | 18 s |
| 10k | 267,992 | 3.9 s | 98 s |
| 12k | 1,452,683 | 16.1 s | 511 s, 11.1 GB process peak |
| full | 28,343,526 | 335 s, 26.5 GB | about 2.8 h, not run |

The classifier costs about 350 microseconds per clique: an R
`lapply()` over cliques (lines 863-900) with a `combn()`, pasted pair
keys and an 18-element list per clique.

## 2. Consumers of the member-row shape

| consumer | reads | B breaks? | edit |
|---|---|---|---|
| `classify_gene_cliques()` (`R/clique_gene_graph.R:732-953`) | members per clique: species set, member keys for the pair lookup (829-883) and for `.gcg_missing_reason()` (489-536); `n_members` for the merge check (52-81) | yes | rewrite on the incidence (WP-b) |
| `rcomplex()` (`R/rcomplex.R:121-129`) | passes `cliques` through | no | 3 lines: store `cliques` and `members` |
| `print.rcomplex()` (223-256), `summary.rcomplex()` (261-293) | only `classification` | no | `summary()` gains a `hogs` column |
| `as.data.frame.rcomplex()` (295) | `classification` | no | none |
| `write_rcomplex()` (309-321) | every top-level data frame (`Filter(is.data.frame, ...)`, 312) | no | none: `members.tsv` appears by itself |
| `find_cliques()`, `classify_cliques()`, `clique_threshold_sweep()`, `clique_stability()` | the species-graph table from `find_cliques()` | untouched | none |
| `get_coexpressed_hogs()` (`R/comparison.R:1146`) | networks, not cliques | untouched | none |
| `tests/testthat/test-clique-gene-graph.R` (1,468 lines, 75 builder calls, 63 classifier calls) | mostly classifier output; member rows at 105-106, 128, 563, 663, 673, 1143; hand-built member tables at 573, 593, 695, 742, 782 | a few | about 15 tests, about 80 lines, if the classifier keeps accepting a member-row table |
| `tests/testthat/test-scores.R:117-124` | loops over member rows of each clique | yes | 6 lines |
| walkthrough `gene-cliques` chunk | `table(gene_class$classification)` | no | none |
| quickstart (`vignettes/quickstart.Rmd:119`) | names `cliques.tsv` | no | add `members.tsv` |

## 3. Candidate shapes

### 3.1 Builder memory is mostly not the table

Clique counts per HOG are extreme. At full size 10,443 HOGs have
cliques. The median HOG has 3. Three HOGs have 11.0 M, 8.6 M and 4.7 M
(86 %), the top five 91 %. HOG:0016636 holds 10.96 M cliques with all
8 species: one species pattern. At 12k one HOG (HOG:0025026, 7
species, HVUL absent, 6-15 copies before the 10-copy cap) holds 1.40 M
of the 1.45 M cliques. All of them are `partial_present`.

One `igraph::max_cliques()` call on HOG:0016636 alone peaks at 4.6 GB
(`m8.R`): R holds 11 M small integer vectors and igraph its own copy.
`igraph::count_max_cliques()` over all HOGs takes 45 s and 1.0 GB at
full size (`m4.R`, `MODE=COUNT`). So enumeration itself is cheap. What
costs is materialising it in R:

| prototype (full size) | time | peak RSS | result |
|---|---|---|---|
| status quo `gene_clique_graph()` | 335 s | 26.5 GB | about 20 GB table |
| count only, C (`m4` COUNT) | 45 s | 1.0 GB | counts |
| D aggregates, one `max_cliques()` per HOG (`m4` D) | 253 s | 23.5 GB | 9 MB |
| B, one `max_cliques()` per HOG (`m4` B) | 271 s | 26.1 GB | 4.9 GB |
| B, blocks of 2e5 cliques, one call per HOG (`m7`) | 114 s | 16.6 GB | 4.8 GB |
| B, blocks, one call per start vertex in HOGs of 40+ nodes (`m7`) | 129 s | 11.8 GB | 4.8 GB |
| same, loop only, before assembly | | 6.9 GB | |

D's 9 MB result did not lower the peak: any shape needs the
piecewise build, `max_cliques(g, subset = v)` for every vertex id `v`.
On all 167 12k HOGs with 20+ nodes the pieces are disjoint and their
union is the full set (`m9.R`). Caution: `subset` over only the
non-capped vertices lost 3,906 cliques in 24 HOGs. On HOG:0016636 the
largest piece is 1.02 M cliques, peak 0.86 GB, 6.6 s against 8.1 s.

The remaining 4.9 GB between the loop and the peak is the final
`unlist()` while the per-HOG lists are alive. A count pass (45 s)
lets the builder preallocate the columns. That should bring the peak
near 8 GB, but it is not measured.

`m7.R` also swaps `match()` on pair keys (377-380) for a dense
per-HOG edge-index matrix, and `split()` + `vapply()` for `collapse`.
At 12k it takes 7.1 s against 16.1 s, with all seven statistics
`all.equal()` to the status quo for all 1,452,683 cliques.

### 3.2 Sizes of each shape

| shape | 12k MB | full | TSV at 12k |
|---|---|---|---|
| status quo | 961 | about 20 GB | 1,584 MB |
| A: same rows, int id, factor species, stats dropped | 233 | about 5.0 GB, stats lost | |
| B: cliques + members (int id, factor species, chr gene) | 78 + 156 | 1.5 + 3.3 GB, measured | 170 + 359 MB |
| B, members as int gene index + gene table | 78 + 116 | 1.5 + 2.5 GB | |
| B as `ngCMatrix` genes x cliques | 44 | 0.93 GB, measured | needs flattening |
| C: one row per clique, list-columns | 1,765 | about 35 GB | needs flattening |
| W: one row per clique, one gene column per species | 190 | about 3.7 GB | 412 MB |
| D: (HOG, species set) aggregates | 1.7 | 9 MB, measured | trivial |

- **A** is B with the statistics thrown away; not a separate option.
- **B** is the normalised member table. TSV is two files,
  `cliques.tsv` and `members.tsv`; `write_rcomplex()` writes both
  unchanged if the driver stores them as two top-level data frames.
  At full size that is about 3.3 + 7.7 GB of TSV against about 34 GB
  now.
- **C** is the worst: 1.45 M small vectors cost 1.67 GB (a 48-byte
  header each) against 156 MB as member rows, and TSV needs
  flattening. Rejected.
- **W** is `find_cliques()`'s shape. It is competitive (64 B per
  clique for the gene pointers at S = 8) and gives one TSV. It cannot
  hold the case `gene_clique_graph()` documents: a within-species edge
  gives `n_species < n_members`, two genes of one species. Its columns
  depend on the species set, so it is not an incidence. Second choice.
- **D** is in 3.3.
- **E**, a per-HOG clique cap or a species-scaled copy cap, would cut
  the three giant HOGs. It changes results, and silently: the cap
  would pick which paralog combinations exist. The 10-copy cap already
  does this, with a message (lines 332-356, 412-418). Section 3.1
  shows the explosion costs 45 s and 1 GB to enumerate once the
  materialisation is fixed, so a cap buys nothing that B does not.
  Keep it as the fallback only if Orion shows a HOG past the 64 GB
  node even under B.

### 3.3 D: does the explosion carry information?

The tier can depend on which paralog copies form the clique. With
`alpha_graph == alpha_call`, as the driver calls it
(`R/rcomplex.R:121-128`), every clique edge is significant, so `n_sig
= n_pairs` and the home clade follow from the species set alone. But
`missing_reason` depends on the members (`.gcg_missing_reason()`,
lines 499-534). Whether a species is `untested`, `tested_ns`,
`underpowered` or `extendable` depends on the rows between its genes
and *these* members. `lineage_specific`, `trait_specific`,
`underpowered` and `partial_present` read those reasons (1090-1107).

Measured (`m3.R`):

| | 8k | 12k | full |
|---|---|---|---|
| cliques | 41,329 | 1,452,683 | 28,343,526 |
| HOGs | 1,410 | 3,742 | 10,443 |
| (HOG, species set) | 2,290 | 7,413 | 41,866 |
| (HOG, species set, tier) | 2,290 | 7,415 | not classifiable now |
| groups with more than one tier | 0 | 2 (163 cliques) | |
| groups with more than one `missing_reason` | 42 | 156 | |
| collapse, cliques per (HOG, set) | 18x | 196x | 677x |

Tiers at 12k by unit:

| tier | cliques | (HOG, set) | HOGs |
|---|---|---|---|
| complete_conserved | 10,660 | 4 | 4 |
| partial_present | 1,401,514 | 19 | 19 |
| trait_specific | 460 | 65 | 65 |
| underpowered | 19 | 8 | 8 |
| unclassified | 40,030 | 7,319 | 3,721 |

The two mixed groups are `trait_specific/unclassified` and
`trait_specific/underpowered`. So the explosion mostly repeats
information. A clique count of 1.4 M `partial_present` is one HOG.
But it does not *only* repeat it: a tier at the pattern level is not
exact.

The published workflow made the same reduction. Rodriguez *et al.*
(2026) report "70,458 cliques from 2098 orthogroups, of which 2145
cliques were non-overlapping (unique)", and count their lineage
tiers in unique cliques too (PubMed, PMC13503837,
[doi:10.1038/s41467-026-75624-2](https://doi.org/10.1038/s41467-026-75624-2)).
Their headline unit is the orthogroup. rcomplex already reports it as
`hog_class` (lines 947-951).

Can D be computed without enumerating? Partly.

- `n_cliques` per pattern: only by enumeration. igraph's count gives
  a per-HOG total. A per-pattern count needs a C callback igraph's R
  API does not offer. Enumeration in C is cheap (45 s), so this does
  not matter.
- Constant per pattern: `n_members`, `n_species`, and `n_edges =
  choose(|P|, 2)` for cross-species edges.
- Extremes per pattern (best `mean_q`, top `score`, top
  `mean_effect_size`, `max_q`): a max-weight clique search per
  pattern. Feasible by branch and bound at S <= 8, c <= 10, but it
  needs a C++ kernel. Not prototyped.
- Means over cliques, and the tier (per-clique `missing_reason`):
  enumeration.

So D is the right report and the wrong store. Computing it from B and
the classification is a grouped count, a few seconds at full size.

### 3.4 B as an incidence, and the hypergraph

B's member table is the COO form `(clique_id, gene)` of a genes x
cliques incidence; `Matrix::sparseMatrix()` makes it an `ngCMatrix`
(44 MB at 12k, 932 MB at full: 4-byte `i` per nonzero, 4-byte `p` per
clique, no `x`). The status-quo table is the same COO incidence plus
repeated columns. B is the same data, normalised, not a new
capability.

That incidence is what a hypergraph engine reads (node x hyperedge).
But it is not WP15's hypergraph (`dev/design-notes/sharpen-plan.md`,
WP15; `hypergraph-community-literature.md` section 7). There the
hyperedge is the HOG, over all its copies, with one nonzero per gene,
taken from the ortholog table. Cliques as hyperedges would be
wrong. 86 % would come from three HOGs, as near-duplicates. Hypergraph
modularity sums over hyperedges, so a HOG with c^S cliques would
weigh c^S. If a clique hypergraph is ever wanted, its hyperedges are
D's patterns, not the cliques. WP15 needs nothing from this note.

The incidence does pay off in the classifier. Every `missing_reason`
is a function of sparse products of the incidence `M` with edge
matrices `T` (tested), `S` (significant), `F` (failed) and `U`
(failed, power below `min_power`), summed per species. `absent`: no
node of the species in the HOG. `untested`: `T M` zero on the
species. `extendable`: a gene with `(S M)[g, c] = m_c`.
`underpowered`: `F M` equals `U M` on the species and is positive.

`m6.R` computes all reasons this way. At 8k: 41,329 of 41,329
identical to `classify_gene_cliques()` in 1.1 s, against 18 s. At 12k:
1,452,683 of 1,452,683 identical in 15.8 s, against 511 s (32x),
3.5 GB process peak, which includes loading the 961 MB status-quo
table for the comparison. Pair counts follow the same way: `n_sig =
colSums(M * (S M)) / 2`, and likewise `n_present` and the q sums.
`max_q` and the within/cross split need the blockwise pair expansion
the builder already does.

## 4. Recommendation

B stored, built piecewise; a vectorised classifier on the incidence;
D computed in `summary()`. No enumeration on request: enumeration in C
costs 45 s, materialising it in R was the cost. `find_cliques()` stays:
it is one row per clique already and has no paralog explosion.

Driver at full size after both WPs: cliques 1.5 + members 3.3 +
classification about 3.6 + edges 0.3 GB, about 8.7 GB held. Builder
peak 11.8 GB measured, near 8 GB with preallocation (not measured).

## 5. Spec

### WP-a Clique table shape (about 250 lines)

- Files: `R/clique_gene_graph.R` (builder 232-461, `.gcg_empty`,
  `.gcg_graph_attrs`), `R/rcomplex.R`, `tests/testthat/test-scores.R`,
  `tests/testthat/test-clique-gene-graph.R`, `vignettes/quickstart.Rmd`,
  `NEWS.md`, `man/`.
- Signature: `gene_clique_graph(edges, min_size = 3L, alpha_graph =
  0.1, ...)`. `id_prefix` goes: integer ids cannot carry it. No new
  argument.
- Return: a list of class `gene_cliques` with `cliques` (one row per
  clique: `clique_id` int 1..K, `hog`, `n_members`, `n_species`,
  `n_edges`, `mean_q`, `max_q`, `mean_effect_size`, `score`,
  `mean_q_floor`, `n_cliques_at_q_floor`) and `members` (`clique_id`
  int, `species` factor over the edge table's species, `gene` chr).
  Attributes as now on `cliques`.
- Build: per HOG as now (graph, collapse, cap). Dense node x node
  edge-index matrix for pair lookup. `max_cliques(subset = v)` for
  every vertex id in HOGs of 40 or more kept nodes, else one call.
  Statistics in blocks of 2e5 cliques with `collapse`. Optional:
  `count_max_cliques()` pass to preallocate.
- Driver: `res$cliques <- g$cliques; res$members <- g$members`.
  `write_rcomplex()` unchanged; it writes `members.tsv`.
- Combining runs (several `alpha_graph`): pass a list of results to
  `classify_gene_cliques()`, which numbers them and adds a `run`
  column. That replaces `id_prefix` and the three merge signatures
  of `.gcg_check_clique_ids()` (lines 52-81) become impossible by
  construction.
- Tests: `test-scores.R:117-124` loops over `members`; builder tests
  read `g$members` where they read member rows; id-prefix tests
  rewritten as list-of-runs tests.

### WP-b Classifier on the incidence (about 350 lines)

- Files: `R/clique_gene_graph.R` (`classify_gene_cliques.default`
  732-953, `.gcg_missing_reason` 489-536, `.gcg_classify_one`
  1020-1176), tests, `NEWS.md`.
- Input: a `gene_cliques` result, a list of them, or (kept, about 10
  lines) any member-row table with `clique_id`, `hog`, `species`,
  `gene`, which is split into the two tables. Hand-built test tables
  keep working.
- Do: node index over the edge table; `T`, `S`, `F`, `U` as symmetric
  `dgCMatrix`; incidence `M`; reasons and pair counts by sparse
  products in blocks of 2e5 cliques (`m6.R`); `max_q` and within/cross
  counts by blockwise pair expansion; home clade per species pattern
  (few patterns), then the tier waterfall as vector `ifelse()` in the
  order of lines 1144-1164; the `trait_specific` core pass and
  `hog_class` as now.
- Return: unchanged columns, `clique_id` integer (plus `run` for a
  list of runs).
- `summary.rcomplex()`: the tier table gains `hogs` (distinct HOGs per
  tier) and `patterns` (distinct HOG x species sets per tier). This is
  D.

### Acceptance

- Identity: on the 8k and 12k edge tables, cliques matched by (HOG,
  sorted members) give `all.equal()` statistics, and every
  classification column except the `clique_id` type is identical to
  `5f31d6f`. Same on `pres_fixture()`, the quickstart and the
  walkthrough data.
- Full leaf (Orion or this Mac): `gene_clique_graph()` <= 150 s and
  peak RSS <= 12 GB (stretch: <= 8 GB with preallocation); stored
  `cliques` + `members` <= 5 GB; `classify_gene_cliques()` <= 15 min
  and peak RSS <= 16 GB; whole `rcomplex(null = TRUE)` finishes.
- 12k: classifier <= 30 s (from 511 s).
- `R CMD check` OK, lint clean, rng-contract table unchanged (no
  randomness here).

### NEWS (0.4.0, which is already not backward compatible)

- `gene_clique_graph()` returns two tables: `cliques`, one row per
  clique with an integer `clique_id`, and `members`, one row per
  clique member. The per-clique statistics are no longer repeated on
  every member row. `id_prefix` is gone: pass several runs to
  `classify_gene_cliques()` as a list. On the full 8-species leaf
  data the result is 4.8 GB instead of about 20 GB, built in about
  2 min instead of 6.
- `classify_gene_cliques()` is about 30 times faster, same tiers.
- `rcomplex()` stores `members` beside `cliques` (`members.tsv`).
  `summary()` counts each tier in cliques, HOGs and species patterns.

## Method

Scripts and logs in the session scratchpad, `shape/`. Package installed
from this worktree into `shape/lib`.

- `m1.R`: status quo builder and classifier at 8k/10k/12k, time, peak.
- `m2.R`: `object.size()` per column; shapes A, B, C, W; TSV sizes.
- `m3.R`: D collapse, tier mixing per (HOG, species set), HOG sizes.
- `m4.R`: `MODE=COUNT`, `D`, `B` prototypes at 12k, 15k, full.
- `m6.R`: `missing_reason` by sparse products against the classifier.
- `m7.R`: piecewise B builder; `CHECK=1` compares statistics,
  `STOP=1` measures the loop alone.
- `m8.R`, `m9.R`: `max_cliques(subset = )` pieces: memory on
  HOG:0016636, completeness on every 12k HOG with 20+ nodes.

Edge tables: `premerge/E_*.rds` (`premerge/E.R`, 0.4.0 defaults).
Full-size status quo (335 s, 26.5 GB): `fix66/big.log`. Paper quote
via PubMed (PMC13503837).
