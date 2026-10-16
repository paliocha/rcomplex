# Structify toy: planted context changes on a species tree

Probe, 2026-10-16. `dev/probes/structify.R` (generator), `_check.R`,
`_run.R` (run, score). Outputs: `prepare_data/cache/structify/`.

## Generator

Structify-Net (Cazabet, Citraro & Rossetti 2023,
doi:10.24072/pcjournal.335) ranks every gene pair by a cost `d` and
draws edges from a rank probability `P(r)`. Four rank rules are ported
from the zoo: `blocks` (same block 0, else 1; 20 blocks), `spatial`
(Euclidean on x, y), `core` (the zoo's product
`|x_u - .5| |x_v - .5| |x_u - x_v|`, not `x_u + x_v`) and `nested`
(the zoo's `x_u + x_v`, not `max`). A gene with no structure has
`d = Inf`. 10 % of genes start with none.

Mixture, not the Bezier curve: `P(r) = (1 - eps) [r <= m] + eps m / L`.
The ends match the paper (eps 0: the m best pairs; eps 1: Erdos-Renyi).
The draw keeps exactly m edges with the right marginals:
`K1 ~ Binom(m, p1)` from the top m ranks, `m - K1` from the rest.
Edge weight `(L - r + 1) / L`.

Root attributes pass down the tree `(((A,B)ab,(C,D)cd)abcd,
((E,F)ef,(G,H)gh)efgh)`. Events on a branch apply to its whole clade:
gain (none -> fresh), loss (-> none), switch (-> fresh), dup (twin
copies `.1`, `.2`) and dup_switch (copy `.2` takes a fresh attribute).
Design per run: 1,500 hogs, density 0.03, 60 of each of the five kinds
on non-root branches, plus 30 dups at the root. One event per hog.

## Checks (all pass)

- eps 0: the edge set equals the m lowest-ranked pairs.
- eps 0.5: exactly m edges every draw; top-rank share 0.516 against
  0.515 expected.
- blocks, eps 0, 20 x 75 genes: all 33,728 edges within blocks, 0.608
  of the 55,500 within-block pairs. Ties break at random per species.
- A tip with no events has the root attributes. A loss on `abcd` leaves
  the gene with none in A-D and structured in E-H (degree 0 vs 10-12).
- dup_switch on `ab`: two copies in A and B, `.1` kept, `.2` moved.

## Grid

4 structures x eps {0, 0.2, 0.5, 0.8} x seeds 1-3. `rcomplex()` with
`networks`, `clades_from_tree()`, alpha 0.1, `null = FALSE`; then
`edge_history()`. Per hog, the run-1 clique with most members. Wall
4-6 s per run, plus 2-3 s for the shuffled-ortholog null.

Means over 3 seeds. loss, gain, switch = branch recall (a root tie
counts: a loss on `abcd` reads as a gain on `efgh`); gain only on
`abcd`/`efgh` (rc: root children, the only stems whose complement is a
clade). prec = branch precision. fe = changes per hog without an event.
none = share `complete_conserved` (ceiling 0.92: 90 hogs stay none).

| structure | eps | cliques | loss | gain | switch | prec | fe   | none |
|-----------|-----|--------:|-----:|-----:|-------:|-----:|-----:|-----:|
| blocks    | 0   | 1377 | 0.00 | 0.00 | 0.93 | 1.00 | 0.00 | 0.92 |
| blocks    | 0.2 | 1396 | 0.97 | 1.00 | 0.94 | 0.99 | 0.00 | 0.92 |
| blocks    | 0.5 | 1567 | 1.00 | 0.97 | 0.90 | 0.84 | 0.01 | 0.92 |
| blocks    | 0.8 | 2625 | 0.00 | 0.00 | 0.02 | 0.01 | 0.32 | 0.00 |
| spatial   | 0   | 1376 | 0.00 | 0.00 | 0.82 | 1.00 | 0.00 | 0.92 |
| spatial   | 0.2 | 1401 | 0.99 | 1.00 | 0.80 | 0.98 | 0.00 | 0.92 |
| spatial   | 0.5 | 1565 | 0.99 | 1.00 | 0.83 | 0.80 | 0.01 | 0.92 |
| spatial   | 0.8 | 7493 | 0.15 | 0.36 | 0.24 | 0.03 | 1.21 | 0.00 |
| core      | 0   | 1325 | 0.00 | 0.00 | 0.24 | 0.95 | 0.00 | 0.87 |
| core      | 0.2 | 1718 | 0.67 | 0.95 | 0.39 | 0.28 | 0.16 | 0.66 |
| core      | 0.5 | 2904 | 0.63 | 0.84 | 0.39 | 0.08 | 0.78 | 0.40 |
| core      | 0.8 | 1620 | 0.03 | 0.06 | 0.04 | 0.03 | 1.22 | 0.00 |
| nested    | 0   |  374 | 0.00 | 0.00 | 0.00 |  NA  | 0.00 | 0.25 |
| nested    | 0.2 |  389 | 0.17 | 0.18 | 0.14 | 0.92 | 0.00 | 0.25 |
| nested    | 0.5 |  423 | 0.19 | 0.20 | 0.19 | 0.81 | 0.01 | 0.23 |
| nested    | 0.8 | 1481 | 0.01 | 0.00 | 0.02 | 0.01 | 1.05 | 0.00 |

Root-child tiers and paralogs. gain, loss = share `trait_specific` or
`lineage_specific`; switch = `differentiated`;
copy = dup_switch members in the dup clade that use the kept copy `.1`;
twin = dup cliques with `n_copies` 2 in every dup species; split = dup
hogs with more than one run-1 clique. dup, d_sw = `complete_conserved`.

| structure | eps | gain | loss | switch | dup  | d_sw | copy | twin | split |
|-----------|-----|-----:|-----:|-------:|-----:|-----:|-----:|-----:|------:|
| blocks    | 0   | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |  0 |
| blocks    | 0.2 | 0.80 | 0.84 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |  0 |
| blocks    | 0.5 | 0.89 | 0.93 | 0.93 | 0.99 | 0.99 | 1.00 | 1.00 |  1 |
| blocks    | 0.8 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.91 | 0.00 | 69 |
| spatial   | 0.2 | 0.81 | 0.84 | 0.82 | 1.00 | 1.00 | 1.00 | 1.00 |  0 |
| spatial   | 0.5 | 0.90 | 0.95 | 0.82 | 1.00 | 0.98 | 1.00 | 0.99 |  1 |
| core      | 0.2 | 0.49 | 0.32 | 0.25 | 0.77 | 0.73 | 0.94 | 0.76 | 23 |
| core      | 0.5 | 0.44 | 0.32 | 0.26 | 0.49 | 0.40 | 0.85 | 0.44 | 49 |
| nested    | 0.2 | 0.09 | 0.25 | 0.00 | 0.29 | 0.22 | 0.74 | 0.97 |  1 |

All runs: 0 capped hogs, 0 null cliques, `root_tied` 0-15 %.

## Verdict

- **Losses need a little noise to be seen.** At eps 0 a lost gene has
  no edges, its tests are `underpowered`, and `edge_history()` reads it
  as unknown, so loss recall is 0. From eps 0.2 a lost gene keeps random
  edges, fails the test with power, and recall is 0.97-1.00 on blocks
  and spatial. The same holds for gains on `abcd`/`efgh`.
- **Gains are seen only on 4-species stems.** A gain on a tip or a
  2-species clade cannot make a 3-gene clique. Overall gain recall is
  0.18-0.27; on root-child stems it is 0.95-1.00.
- **Switches are found at eps 0** (0.93 blocks, 0.82 spatial), because
  the moved gene keeps edges. On root-child stems they read as
  `differentiated` (1.00 blocks).
- **Paralogs work.** Perfect twins collapse into one clique of
  multiplicity 2 (twin 1.00 to eps 0.5). The conserved copy is the one
  chosen after a dup_switch (copy 1.00). No hog hits the clique cap.
- **Break point: between eps 0.5 and 0.8.** At 0.8 the 8-species
  cliques fall apart, `complete_conserved` drops to 0, twins split, and
  false events reach 0.3-1.2 per hog. The null stays at 0 cliques, so
  noise does not make cliques from nothing.
- **Blocks and spatial are easy.** Precision 0.98-1.00 at eps 0.2.
- **Core is hard.** Hub genes share most partners, so a switch barely
  changes the neighbourhood (switch recall 0.24-0.39). Precision falls
  to 0.28 at eps 0.2 and false events rise already there.
- **Nested is hardest.** With `x_u + x_v` and 3 % density only genes
  with x below about 0.25 get structured edges, so three quarters of
  hogs never form a clique (tier_none 0.25). What is found is precise
  (0.92 at eps 0.2), but recall stays under 0.2.
