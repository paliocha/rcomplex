# Large families: replace the 10-copy cap by a clique budget (WP-e)

2026-10-13. Branch `shape/clique-budget`, base `65cd20d` (main after
WP-d). Issue #79. Probes `dev/probes/copy-cap.R` (per-hog bounds and
exact counts) and `dev/probes/copy-cap-exact.R` (exact counts, capped
against uncapped). Data: wood (6 species, `e5` run), Pooideae leaf at
12k genes per species (`prepare_data/cache/E_12000.rds`) and the full
leaf edge table (8 species, 2,664,155 rows, rebuilt into
`~/.cache/rcomplex/E_full.rds`; networks 4.4 min, edges 6.0 min).

## 0. Verdict

- The cap is a tie-break bug (#79) and a crude rule, but some bound is
  needed. On wood nothing needs it: the largest family gives 675
  maximal cliques uncapped. On the full leaf data 14 families give more
  than 100,000 cliques each after twin collapse, the largest 6.7
  million, and together the 153 big families give 22.9 million against
  8.0 million with a 10-class cap. That is 6 times today's whole table
  (3.56 million) and does not fit a 64 GB machine in the classifier.
- A cap on genes per species does not bound the output. A family with
  76 wiring classes (fewer than 10 per species) gives 2.41 million
  cliques uncapped and 2.09 million capped: a near-complete 8-partite
  graph is the problem, not the copy number.
- Counting is cheap. `igraph::count_max_cliques()` on the worst family
  (618 classes, 21,355 edges, 6.7 million cliques) takes 1 s. So the
  bound can be measured per hog before anything is enumerated.
- Decision: collapse twins, count, and only when a hog exceeds a clique
  budget drop its weakest wiring classes until it fits, by evidence,
  never by position. Record what was dropped and flag the hog's
  cliques. The 10-copy cap goes. Section 3 has the spec.

## 1. What the cap does today

`R/clique_gene_graph.R:423-446`. Per hog, before the twin collapse, a
species with more than 10 genes keeps the 10 of highest degree over the
called edges; ties keep the gene that entered first. Dropped genes
become isolated vertices. The rule comes from `find_cliques()`
(`R/cliques.R:459`) and the published workflow.

Row-order dependence (#79): the WP-d routing fix re-sorted the edge
table; identical edges, different order, 51,703 against 51,775
cliques on wood; hog tiers unchanged. The degree tie at the cut is the
mechanism: 73 of the 125 capped wood hogs, 122 of 214 on the full leaf
data, have one.

## 2. Measurements

Called edges at q < 0.1. "Classes" are twin classes (same species, same
partner set). The bound is the product over species of class counts
plus one, minus one (every species in or out); it is loose by 4 to 10
orders of magnitude and useless as a budget. Counts are exact
(`count_max_cliques(min = 3)` on the twin-collapsed graph).

| data | hogs with edges | capped now | tie at cut | > 10 classes after collapse | uncapped cliques over capped hogs | largest family |
|---|---|---|---|---|---|---|
| wood | 14,906 | 125 | 73 | 78 | 14,768 | 675 (72 classes) |
| leaf 12k | 6,492 | 15 | 11 | 12 | 16,476 | 1,692 (231 classes) |
| leaf full | 15,101 | 214 | 122 | 191 | 22,922,249 + 10,525 | 6,681,084 (618 classes, 21,355 edges) |

Full leaf, the 153 families the bound could not resolve: uncapped
22.9 M, with the 10 most connected classes per species 8.0 M; the
package today, capping genes before collapse, gives 3.56 M for the
whole table. Fourteen families exceed 100,000 uncapped (6.68 M, 2.82 M,
2.41 M, 2.21 M, 1.87 M, 1.78 M, 1.37 M, 1.09 M, 719 k, 507 k, 439 k,
180 k, 175 k, 153 k); the rest of the 153 sum to about 0.4 M. Two of
the fourteen have fewer than 100 classes (76 and 95) and keep 2.1 M and
1.7 M cliques under the gene cap. Counting took at most 1 s per family.

Under a budget of 100,000 cliques per hog the full leaf table would
hold about 1.4 M from the fourteen plus about 0.5 M from everything
else: smaller than today, with no family thinned that did not need it,
and 2,057 genes on wood (133 hogs) no longer dropped at all.

## 3. Spec: WP-e clique budget (about 180 lines)

- Files: `R/clique_gene_graph.R` (cap 423-446, collapse 448-478,
  enumeration 480-499, message 502-508, return lists), `R/rcomplex.R`
  (`res$dropped`, `print()`), `R/clique_gene_classify.R` (`capped`
  column), `tests/testthat/test-clique-gene-graph.R`,
  `tests/testthat/test-rcomplex.R`, `NEWS.md`, `man/`, `CLAUDE.md`
  ("Two clique backends": replace the 10-copy sentence), README or
  quickstart only where they mention the cap.
- Signature: `gene_clique_graph(edges, min_size = 3L, alpha_graph =
  0.1, max_cliques = 1e5, ...)`. `rcomplex()` keeps its formals and
  passes nothing; the default applies. `Inf` disables the budget.
- Per hog, after the duplicate-row collapse and the twin collapse (the
  twin collapse moves ahead of any reduction):
  1. `n_full <- igraph::count_max_cliques(g, min = min_size)`.
  2. If `n_full <= max_cliques`, enumerate as now.
  3. Else rank the twin classes by evidence: the sum of the collapsed
     edge `score` over the class's edges, descending; ties by the
     representative gene id in the radix order the collapse already
     uses. The top class of every species is always kept. Binary search
     the largest number of kept classes whose subgraph has at most
     `max_cliques` maximal cliques (recount each step; at most about 10
     counts). Drop the other classes with all their genes; they become
     isolated vertices as the capped genes do today.
  4. Record one row per dropped gene (twins included): `hog`,
     `species`, `gene`, `n_cliques_full`, `n_classes`, `n_classes_kept`.
- Return: `cliques`, `members`, `twins` as now, plus `dropped` (the
  table above; zero rows when nothing was dropped). Message when rows
  exist: "`gene_clique_graph`: N ortholog groups exceeded `max_cliques`
  = B maximal cliques; M genes dropped, see `$dropped`." The 10-copy
  cap, its counters and its message go.
- `classify_gene_cliques()`: logical column `capped` on every
  classification row, `TRUE` when the clique's hog has rows in
  `dropped`. A list of runs: `TRUE` if any run dropped from that hog.
- Driver: `res$dropped <- g$dropped`; `write_rcomplex()` writes it by
  its loop (list `dropped.tsv` in the roxygen); `print()` adds "N hogs
  over the clique budget" to the cliques line only when N > 0.
- Prose (STE-lite): `@param max_cliques`: "Largest number of maximal
  cliques one hog may produce. A hog above it keeps its wiring classes
  with the strongest edges until it fits. Default 1e5; `Inf` keeps
  every clique." `@details` paragraph: why a budget and not a copy cap
  (two sentences, the 76-class family as the example).
- Tests:
  1. Complete tripartite fixture, 5 classes per species, 125 cliques:
     default keeps all 125 and `dropped` has zero rows; `max_cliques =
     50` keeps a deterministic class set (the strongest by score), the
     kept count is at most 50, `dropped$n_cliques_full == 125`, and the
     same call on the fixture with its rows shuffled (`seed`ed
     `sample()`) gives an identical result (the #79 test).
  2. The top class of every species survives even when the budget
     would drop it by score.
  3. `max_cliques = Inf` equals a run on a fixture that never hits the
     budget.
  4. Twin collapse happens before the budget: a fixture with 12 twins
     in one species and 125 cliques is not reduced under `max_cliques =
     1e3` (today's cap would drop 2 genes).
  5. `classify_gene_cliques()` carries `capped`; `rcomplex()` stores
     `dropped`; `summary()`/`print()` unchanged except the added line.
  6. The existing cap test(s) in `test-clique-gene-graph.R` are
     rewritten to the budget; identity fixtures
     (`fixtures/gcg_unreduced.rds`) still pass, since no fixture hits
     the budget.
- No new RNG entry point. Lines at or under 80 characters. Surface:
  31 exports; `gene_clique_graph()` has 4 named formals.

### Acceptance

- `R CMD check` OK on the built tarball, lint clean, `dev/sharpen-prose.R`
  clean, surface test green.
- Wood (`prepare_data/evotree/e5_rcomplex_0.4.R` copied to scratch with
  a new output directory; skip if the data are absent): `dropped` has
  zero rows; clique count independent of edge row order (run
  `gene_clique_graph()` on `res$edges` and on `res$edges[sample(nrow(
  res$edges)), ]`, compare `cliques` after sorting by hog and members).
- Full leaf: `gene_clique_graph(readRDS("~/.cache/rcomplex/E_full.rds"))`
  finishes in at most 90 s and 6 GB peak RSS (`/usr/bin/time -l`), gives
  at most 3.6 M cliques, and `dropped` names the fourteen families
  above with their `n_cliques_full`. Report the count of hogs and genes
  dropped and the five largest `n_cliques_full`.

### NEWS (development version)

- `gene_clique_graph()` replaces the 10-copies-per-species cap by a
  clique budget, `max_cliques = 1e5` per hog. Twins collapse first; a
  hog over the budget keeps the wiring classes with the strongest
  edges, chosen by evidence and never by row order (#79). The result
  gains `dropped`, one row per dropped gene with the hog's full clique
  count; `classify_gene_cliques()` flags those hogs' cliques with
  `capped`; `rcomplex()` stores `dropped`. On the wood data nothing is
  dropped any more (2,057 genes were); on the 8-species leaf data
  fourteen families are reduced instead of 214.

## 4. Left open

- `find_cliques()` (species graph) keeps its own copy cap at
  `R/cliques.R:459`; it picks one gene per species and does not
  explode. Align later if wanted.
- A family over the budget is conserved as a block. Whether to report
  such hogs at the species level only (one best copy per species) is a
  design question for after the first runs with `capped`.
