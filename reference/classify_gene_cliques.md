# Classify gene-graph cliques into conservation tiers

Applies the five-tier taxonomy of Rodriguez et al. (2026), plus a
`trait_specific` tier, to the cliques from
[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md).
Every threshold follows from the number of species supplied.

## Usage

``` r
classify_gene_cliques(cliques, ...)

# Default S3 method
classify_gene_cliques(
  cliques,
  edges,
  species,
  clades = NULL,
  alpha_call = 0.1,
  alpha_graph = 0.9,
  min_power = 0.8,
  ...
)
```

## Arguments

- cliques:

  A
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  result, a list of them, or a table of member rows with `clique_id`,
  `hog`, `species` and `gene` columns. Combine runs at several
  `alpha_graph` values to expose every tier: a clique complete at
  `alpha_call` need not be maximal on a looser graph. Pass the runs as a
  list. The result then gains a `run` column, the position of the run in
  the list.

- ...:

  Additional arguments passed to the default method.

- edges:

  The full, unfiltered co-expressolog table. It must not be pre-filtered
  on `q_value`: the gap tier needs to see rows that were tested and
  failed in order to refuse them. An optional `power` column (from
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md))
  enables the `underpowered` tier; without it, or where it is `NA`, the
  classification is unchanged.

- species:

  Character vector of every species in the analysis. Every species
  appearing in `cliques` must be listed; a stranger would be counted
  into the clique's species total while also being reported as missing.

- clades:

  Optional named list of species vectors, one per clade. Clades may nest
  but must not cross. A species in no clade forms its own clade. A flat
  vector `x` of species to clade is `split(names(x), x)`. Required for
  `lineage_specific`, `trait_specific` and `differentiated`. When
  `NULL`, these tiers are skipped.

- alpha_call:

  Significance threshold for calling a species pair conserved (default
  0.1).

- alpha_graph:

  Loose threshold defining the `partial_significant` graph (default
  0.9).

- min_power:

  Detection power below which a non-significant pair is read as
  uninformative rather than as evidence against conservation (default
  0.8). Only used when `edges` has `power`. For rank-test edges
  (`find_coexpressologs(method = "rank")`) use 0.9: under the default
  reference rank `p0` their power is at least 0.5 by construction (0
  when a direction of the species pair has no call under
  `pval_combine = "max"`, or neither has under `"min"`), and it
  overstates detection.

## Value

A data frame with one row per clique:

- run:

  Only for a list of runs: the run of the clique

- clique_id, hog:

  Clique identity. The id is the run's own `clique_id`.

- classification:

  Tier, or `"unclassified"`

- hog_class:

  Earliest tier reached by any clique of that HOG, reproducing the
  HOG-level precedence of the published scripts without discarding the
  per-clique detail

- n_members, n_species:

  Clique size

- n_pairs, n_present, n_sig:

  Member pairs, pairs with a row in `edges`, and pairs significant at
  `alpha_call`

- n_sig_within, n_sig_cross:

  Significant pairs within one child clade of the clique's home, and
  across two (`NA` without `clades`)

- n_underpowered_cross:

  Cross-clade pairs present in `edges`, not significant, and with
  `power` below `min_power` (`NA` without `clades`)

- clade:

  Name of the smallest clade that holds every species of the clique.
  `NA` when no clade holds them all, or without `clades`.

- n_missing, missing_species, missing_reason:

  Species not in the clique, and why: `absent` (no row in the HOG),
  `untested` (in the HOG but never compared to a member), `tested_ns`
  (compared and not significant), `underpowered` (compared and not
  significant, every such test with `power` below `min_power`),
  `extendable` (significant against every member, so the clique came
  from a looser graph). Comma-separated and positionally aligned

- mean_q, max_q:

  Recomputed from `edges` over present pairs

- mean_effect_size:

  Present when `edges` carries `effect_size`; the ranking key to prefer,
  since `mean_q` saturates

- mean_q_floor, n_cliques_at_q_floor:

  Smallest `mean_q` over the classified cliques and how many are tied
  there – the resolution `mean_q` actually has as a ranking key

The two floor diagnostics are also attributes, alongside
`n_species_total`, `n_pairs_total`, `n_within_pairs`, `n_cross_pairs`,
`cross_max`, `alpha_call`, `alpha_graph` and `max_gap`. All are set even
when the result is empty.

## Details

No threshold is tied to the six species and fifteen species-pairs of the
original workflow.

The home clade of a clique is the smallest clade that holds all its
species. With no such clade, the home is the whole species set. The
child clades of the home are the largest clades inside it. A species of
the home in no child clade forms its own child clade. The
`lineage_specific`, `trait_specific` and `differentiated` tiers, and the
within and cross counts, read the home and its child clades. So a clique
can be specific to, or differentiated inside, a nested clade.

With `S = length(species)` species, `P = choose(S, 2)` species-pairs,
top-level clade sizes `n_l`, `W = sum(choose(n_l, 2))` within-clade
pairs and `X = P - W` cross-clade pairs, the tiers are:

- complete_conserved:

  All `S` species present and all `choose(S, 2)` pairs significant at
  `alpha_call`.

- lineage_specific:

  A complete clique over one entire clade, with no species outside it
  testable against the clique. A species that *was* compared against
  every member and came back non-significant is evidence of a boundary
  rather than a gap, so it blocks this tier: that clique is
  `trait_specific`, or with both clades conserved within a candidate for
  `differentiated` on the cliques of the unfiltered graph.

- partial_significant:

  All `S` species present, every clique edge below `alpha_graph`, and at
  least `choose(S - 1, 2) + 1` pairs significant at `alpha_call`.

- partial_present:

  `S - g` species present for `1 <= g <= max_gap`, all
  `choose(S - g, 2)` pairs significant, and every absent species an
  annotation gap or `underpowered` rather than a rejected test.

- differentiated:

  Every species of the home present, every pair tested, at least one
  child clade of two or more species, every such child clade fully
  significant within itself, and at most `choose(H - 1, 2) - W_H`
  cross-clade pairs significant. Here `H` is the home size and `W_H` its
  within-child pairs. The child clades are disjoint, so a nested pair is
  never `differentiated`. With the whole species set as home, the bound
  is `cross_max`.

- trait_specific:

  A complete clique over one entire clade, at least one outside species
  compared against every member and rejected at adequate power
  (`tested_ns`), no outside species `underpowered`, and no complete
  clique of a disjoint clade in the same HOG. Where `lineage_specific`
  reads the other clade's absence as a gap, this reads its presence and
  rejection as a boundary: one clade conserved, the other present but
  not co-conserved. Two clades each conserved within and rejected across
  is `differentiated`, scored on the unfiltered graph's clique, and
  their one-clade cliques stay `unclassified`. This tier is rcomplex's
  addition to the published five.

- underpowered:

  A clique that would be `lineage_specific`, `trait_specific` or
  `differentiated` but for tests that could not have succeeded, read
  from a `power` column in `edges` (see
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)).
  It takes the place of the two specificity tiers when at least one
  outside species is `underpowered` – reading that species as conserved
  would extend the clique, so neither call survives – and of
  `differentiated` when `n_sig_cross + n_underpowered_cross` exceeds its
  bound: a specificity or divergence call must survive treating every
  underpowered pair as possibly significant. A low-degree gene cannot
  reach the call whatever its conservation, so without this its missing
  edges read as a clade boundary.

The waterfall is evaluated in that order, `underpowered` at the position
of the tier it replaces, and the first match wins.

`choose(S - 1, 2) + 1` equals `choose(S, 2) - (S - 2)`, so the
`partial_significant` tolerance is `S - 2` non-significant edges. That
is the largest tolerance under which no member can be isolated: cutting
a member loose needs all `S - 1` of its edges removed, so with `S - 2`
removals every member still keeps a significant edge.

`differentiated` demands that every member pair actually has a row in
`edges`. Without that, a clique whose cross-clade pairs were never
tested would be scored as diverged on absent evidence – the same
conflation `partial_present` refuses through `missing_reason`.

The classifier works on the sparse gene x clique incidence, in blocks of
100,000 cliques. Sparse products with the edge rows give each missing
species its state.

## References

Rodriguez E, Birkeland S, Chapple ED, et al. (2026). Comparative
regulomics of wood formation across dicot and conifer trees. *Nature
Communications* 17(1).
[doi:10.1038/s41467-026-75624-2](https://doi.org/10.1038/s41467-026-75624-2)

## See also

[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)

## Examples

``` r
edges <- data.frame(
  gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
  species1 = c("SP_A", "SP_A", "SP_B"),
  species2 = c("SP_B", "SP_C", "SP_C"),
  hog = "HOG1", q_value = c(0.01, 0.02, 0.03)
)
cl <- gene_clique_graph(edges)
classify_gene_cliques(cl, edges, c("SP_A", "SP_B", "SP_C"))
#>   clique_id  hog     classification clade n_members n_species n_pairs n_present
#> 1         1 HOG1 complete_conserved  <NA>         3         3       3         3
#>   n_sig n_sig_within n_sig_cross n_underpowered_cross n_missing missing_species
#> 1     3           NA          NA                   NA         0                
#>   missing_reason mean_q max_q          hog_class mean_q_floor
#> 1                  0.02  0.03 complete_conserved         0.02
#>   n_cliques_at_q_floor
#> 1                    1
```
