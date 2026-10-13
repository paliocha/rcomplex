# Large families: species blocks and no copy cap (WP-f)

2026-10-13. Branch `shape/clique-block`, base `cb400ba` (main after
WP-e). Follows `clique-budget.md` section 4. Martin: report a hog over
the clique budget as a block; rework `find_cliques()` accordingly.

## 0. Verdict

- A hog over `max_cliques` is conserved as a block: its copies are so
  densely wired across species that which copy pairs with which cannot
  be resolved. WP-e kept its strongest wiring classes and listed the
  rest in `dropped`. WP-f reports it at the species level instead: one
  row per group of species wired together, with one representative copy
  per species, flagged `capped`. One rule, not two; `dropped` goes.
- `find_cliques()` keeps ten copies per species with an order-dependent
  tie-break (the #79 mechanism) and stops its best-copy search after
  2 million steps without saying so. The cap goes; the search order
  becomes deterministic; the result says when the search was cut.
- Neither published workflow capped copies. See section 1.

## 1. What the published workflows did

ComPlEx (Netotea et al. 2014, doi 10.1186/1471-2164-15-106). OrthoMCL
groups; every ortholog pair tested on its own; a gene is conserved when
at least one of its orthologs has a conserved neighbourhood. Genes had
2 to 4 predicted orthologs on average (up to 6 in Populus). No cap, no
cliques. Multi-copy families were the result, not a nuisance: about
45 % of genes with several orthologs had one conserved and one diverged
copy; in 55 to 60 % of cases the most sequence-similar copy was not the
one with conserved regulation.

EVOTREE (Rodriguez et al. 2026, doi 10.1038/s41467-026-75624-2; code
Zenodo 10.5281/zenodo.21025760, `CliquesWorkflow/`, local copy
`~/.cache/rcomplex/evotree_code/`). `CompleteCliques.R` builds one gene
graph per orthogroup from pairs with `Max.p.Val < 0.1` and calls
`igraph::max_cliques(net, min = 3)`; `PartialCliques.R` and
`DifferentiatedCliques.R` do the same on their own edge sets. No cap on
copies for HOGs. In OG mode (not reported in the paper) four giant OGs
are dropped by hand (`discard_list`). Redundancy is handled afterwards:
`bestUniqueCliques.R` walks each orthogroup's cliques by mean FDR and
keeps those sharing no gene with an earlier pick (70,458 six-species
cliques, 2,145 unique). So EVOTREE equals `gene_clique_graph(alpha_graph
= 0.1, max_cliques = Inf)`; the wood data never needed a bound (largest
family 10,496 cliques).

Our ten-copy rule came from neither. It entered with the C++ rewrite of
the clique finder (`rcomplex-lold/.../find_cliques.cpp`, "replaces the
previous igraph::max_cliques() approach which was exponentially slow"),
a speed shortcut, and `gene_clique_graph()` inherited it until WP-e.

## 2. What `find_cliques()` does today

`src/find_cliques_common.h`. Per hog: collect genes per species and the
species adjacency (lines 256-304); keep the `max_genes_per_sp = 10` most
connected genes per species, ties by entry order (351-400); Bron-
Kerbosch on the species graph (402-477); for each species clique a
backtracking search for the one gene per species with the lowest mean q
(`assign_bt`, 123-207), with `max_iterations = 2000000`; when the budget
is hit the search returns the best assignment so far and nothing records
it. Gene order inside a species follows first appearance in `edges`
(`encode_clique_edges()`, `R/cliques.R:14`), so among equal-cost
assignments the first row wins: order dependence again. Callers:
`R/cliques.R:453-461` (`10L, 0L`), `R/clique_stability.R:225-232`
(`10L`), wrapper signatures in `src/find_cliques.cpp:42-57` and
`src/find_cliques_stability.cpp`.

## 3. Spec

### A. `gene_clique_graph()`: species blocks

- Files: `R/clique_gene_graph.R` (`.gcg_budget()` 89-118 goes;
  `.gcg_no_dropped()` 64; the budget branch 547-575; the `dropped`
  assembly 595-618; roxygen 292-300, 336-338, 375-378),
  `R/clique_gene_classify.R` (`capped` source at `.gcg_members()` 775;
  roxygen 173-176), `R/rcomplex.R` (`res$dropped`, `print()`,
  `write_rcomplex()` roxygen), tests, `NEWS.md`, `man/`.
- Signature unchanged: `gene_clique_graph(edges, min_size = 3L,
  alpha_graph = 0.1, max_cliques = 1e5, ...)`. `Inf` still disables.
- Per hog, after the twin collapse and `count_max_cliques()`: when
  `n_full > max_cliques`, do not enumerate gene cliques. Instead:
  1. Species cliques of the hog: species s and t are adjacent when any
     collapsed edge joins a class of s and a class of t. Take the
     maximal species cliques of at least `min_size` species.
  2. For each species clique, one best assignment of one class per
     species, all mutually adjacent, with the lowest mean q over its
     edges; ties by `effect_size` descending, then `score` descending,
     then representative gene id in radix order. Reuse the level-2
     search of `find_cliques_cpp()` on the collapsed edge list
     (representatives as genes, collapsed mean q and effect). No new
     C++ beyond section B.
  3. Each assignment becomes one clique row through
     `.gcg_clique_stats()`, so `n_members`, `n_species`, `n_edges`,
     `mean_q`, `max_q`, `mean_effect_size`, `score` and `multiplicity`
     (the product of the chosen classes' `n_copies`) have their usual
     meaning, and `members` lists the representatives with `n_copies`.
  4. Record one row per such hog in a table `capped`: `hog`, `n_genes`,
     `n_classes`, `n_cliques_full`, `n_blocks` (rows emitted) and
     `search_complete` (`FALSE` when any block search hit its step
     budget).
- Return: `cliques`, `members`, `twins`, `capped`. `dropped` goes.
  Message when `capped` has rows: "`gene_clique_graph`: N hogs exceeded
  `max_cliques` = B maximal cliques. Each is reported as species blocks;
  see `$capped`."
- `classify_gene_cliques()`: `capped` column as now, sourced from
  `capped$hog`; roxygen: "`TRUE` when the clique is a species block of a
  hog over the clique budget: one representative copy per species."
- Driver: `res$capped <- g$capped`; `write_rcomplex()` writes
  `capped.tsv`; `print()` line "N hogs over the clique budget, reported
  as blocks" only when N > 0.
- Prose (STE-lite). `@param max_cliques`: "Largest number of maximal
  cliques one hog may produce. A hog above it is reported as species
  blocks: one representative copy per species. Default 1e5; `Inf` keeps
  every clique." `@details`: replace the paragraph at 292-300 with
  three or four sentences: paralogs multiply cliques; a copy cap does
  not bound them (the 76-class example); over the budget the hog is
  conserved as a block and which copy pairs with which is unresolved;
  so one row per species clique with the best copies. Mention that
  EVOTREE enumerated every clique and that the wood data never reach
  the budget.

### B. `find_cliques()` and `clique_stability()`: no copy cap

- Files: `src/find_cliques_common.h`, `src/find_cliques.cpp`,
  `src/find_cliques_stability.cpp`, `R/cliques.R`,
  `R/clique_stability.R`, `R/RcppExports.R` (regenerate), tests, NEWS.
- Remove `max_genes_per_sp` and the cap block (351-400) from
  `find_cliques_for_hog()` and both exported wrappers; drop the `10L`
  at `R/cliques.R:459` and `R/clique_stability.R:232`.
- Deterministic search order: `encode_clique_edges()` assigns gene
  integer ids in radix order of `paste(species, gene)`, and
  `find_cliques_for_hog()` sorts each `sp_genes[s]` ascending before
  the search. Among equal-cost assignments the lowest ids win, whatever
  the row order.
- `CliqueResult` gains `bool complete` (`iterations <= max_iterations`
  after `assign_bt`). `find_cliques()` returns it as a logical column
  `exhaustive`: `TRUE` when the best-copy search finished within its
  step budget. The stability wrapper ignores it. `max_iterations` stays
  2,000,000 and hard-coded (no new formal); the C++ wrapper may expose
  it as an argument with that default so a test can hit it.
- Roxygen `find_cliques()` `@details`: two sentences on the search
  budget and `exhaustive`; remove any mention of ten copies.

### Tests

1. Gene graph, complete tripartite fixture (5 classes per species, 125
   cliques): default keeps 125 and `capped` has zero rows;
   `max_cliques = 50` gives exactly one row (the block of all three
   species) whose members are the lowest-mean-q mutually adjacent
   triple, `capped` has one row with `n_cliques_full = 125`,
   `n_blocks = 1`, `search_complete = TRUE`; the same call on shuffled
   rows gives an identical result.
2. A fixture whose species graph has two maximal cliques (ABC and ABD,
   C and D not adjacent) gives two block rows under a budget of 1.
3. `max_cliques = Inf` equals a run that never hits the budget; twin
   collapse still happens before the count (12 twins, 125 cliques,
   budget 1e3, not reduced).
4. `classify_gene_cliques()` carries `capped`; `rcomplex()` stores
   `capped`; the message string.
5. `find_cliques()`: a fixture with 12 copies in one species keeps all
   twelve in play (the old cap dropped two); the result is identical
   under shuffled rows; a fixture with ties among copies picks the
   lowest gene id; `exhaustive` is `TRUE` there and `FALSE` on a dense
   fixture through the C++ wrapper with a tiny `max_iterations`.
6. `clique_stability()` tests still pass without the cap argument.
7. Existing fixtures (`fixtures/gcg_unreduced.rds`) unchanged.

### Acceptance

- `R CMD check` OK on the built tarball, lint clean,
  `dev/sharpen-prose.R` clean, surface test green (31 exports;
  `gene_clique_graph()` 4 named formals; `find_cliques()` at most 8).
- Wood (`prepare_data/evotree/e5_rcomplex_0.4.R` copied to scratch with
  a new output directory; skip if the data are absent): `capped` has
  zero rows; `find_cliques()` on `res$edges` before and after the cap:
  rows, `classify_cliques()` tier counts, hogs whose tier changed, and
  the count of `exhaustive == FALSE`; both functions identical under
  `res$edges[sample(nrow(res$edges)), ]`.
- Full leaf (`gene_clique_graph(readRDS("~/.cache/rcomplex/E_full.rds"))`):
  at most 90 s and 6 GB peak RSS; `capped` names the fourteen families
  of `clique-budget.md` section 2 with their `n_cliques_full`; report
  `n_blocks` and `search_complete` per family and the total clique
  count.

### NEWS (development version)

Edit the WP-e bullet in place (same development version): a hog over
`max_cliques` is reported as species blocks, one row per group of
species wired together with the best copy per species; `capped` lists
those hogs; `dropped` is gone. New bullet: `find_cliques()` and
`clique_stability()` drop the ten-copies-per-species cap; the search
order is deterministic; `find_cliques()` reports `exhaustive`.

## 4. Left open

- A block row's `multiplicity` is the expansion count of that one
  representative clique, not of the family. `n_cliques_full` in
  `capped` holds the family's count.
- EVOTREE's "best unique cliques" (greedy non-overlapping pick) is a
  reporting step, not a clique definition. Whether rcomplex offers it
  is a later question.
