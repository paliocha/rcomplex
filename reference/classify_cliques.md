# Classify HOGs by clique conservation pattern

Convenience wrapper that runs
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
internally (once for all species, once per top-level clade) and applies
a sequential waterfall classification. For fine-grained control over the
clique detection parameters per step, call
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
directly.

## Usage

``` r
classify_cliques(edges, ...)

# Default S3 method
classify_cliques(
  edges,
  target_species,
  clades,
  min_species = 2L,
  stability = NULL,
  min_power = 0.8,
  ...
)
```

## Arguments

- edges:

  Data frame with columns `gene1`, `gene2`, `species1`, `species2`,
  `hog`, `q_value`, `effect_size`, and `type`. Must contain ALL edges
  (conserved + ns + diverged), not pre-filtered, because the
  differentiated check needs to verify absence of cross-group conserved
  edges. An optional `power` column (from
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md))
  enables the `"underpowered"` class; without it, or where it is `NA`,
  the classification is unchanged.

- ...:

  Additional arguments passed to the default method.

- target_species:

  Character vector of all species.

- clades:

  Named list of species vectors, one per clade. Clades may nest but must
  not cross. A species in no clade forms its own clade. A flat trait is
  `split(names(trait), trait)`. A tree gives
  `[clades_from_tree()](phy)`.

- min_species:

  Minimum species for a partial or within-group clique (default 2).

- stability:

  Optional output of
  [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md).

- min_power:

  Detection power below which a non-conserved edge is read as
  uninformative rather than as evidence against conservation (default
  0.8). Only used when `edges` has `power`. For rank-test edges
  (`find_coexpressologs(method = "rank")`) use 0.9: under the default
  reference rank `p0` their power is at least 0.5 by construction (0
  when a direction of the species pair has no call under
  `pval_combine = "max"`, or neither has under `"min"`), and it
  overstates detection.

## Value

A data frame with one row per HOG:

- hog:

  HOG identifier

- classification:

  One of `"complete"`, `"partial"`, `"differentiated"`,
  `"trait_specific"`, `"unclassified"`

- underpowered:

  `TRUE` where a `differentiated` or `trait_specific` call rests on an
  edge whose `power` is below `min_power`. The classification is kept;
  this qualifies it. `FALSE` for every other tier and whenever `edges`
  carries no `power` column.

- n_species:

  Species count in the best clique (NA for unclassified)

- best_mean_q:

  Mean q-value of the best clique (NA for unclassified)

- trait_groups:

  Comma-separated top-level clades that hold the cliques (NA for
  complete/partial/unclassified)

- clade:

  For `trait_specific`, the home clade of the cliques. For
  `differentiated`, the smallest clade that holds all home clades.
  Otherwise, the smallest clade that holds the reported clique. `NA`
  when no clade holds them, and for unclassified HOGs.

- stability_class:

  From stability results (NA if not provided)

- robust:

  Logical: stability class is at least 0 (NA if no stability results
  provided)

## Details

Each HOG is assigned to exactly one category; earlier categories take
precedence.

The pipeline:

1.  **complete**: all target species form a clique (all `C(N,2)` edges
    conserved).

2.  **partial**: a clique with `min_species` to `N-1` species lies in no
    clade.

3.  **differentiated**: cliques sit in two or more disjoint clades, and
    no conserved edge joins two top-level clades.

4.  **trait_specific**: all cliques sit in one clade.

5.  **unclassified**: none of the above.

The home clade of a clique is the smallest clade that holds all its
species. A HOG keeps the home clades that no other of its home clades
holds. Clades are laminar, so these are disjoint. Nested homes count as
one clade, so they are not `differentiated`. Cliques are searched within
each top-level clade. A species in no clade forms its own clade.

## Underpowered calls

`differentiated` and `trait_specific` both rest on edges that were *not*
conserved. An edge too weak to have been called is no evidence either
way, so such a call has to survive reading that edge as conserved. Where
it does not, the row is flagged `underpowered = TRUE` and keeps its
classification – the call is reported and qualified, not replaced, so
nothing is lost and a caller can filter on it. Needs a `power` column in
`edges`; without one, or where it is `NA`, the flag is `FALSE`.

What the flag marks: power is computed against the typical fold
enrichment of a called pair (`rho0`, about 2.5 on the Pooideae data),
and at that enrichment a low-degree gene's shared partners amount to
less than one, which no count test can call. So `underpowered` marks
genes whose neighbourhoods are too small for ordinary conservation to be
visible – on the Pooideae data most of the lowest degree decile and none
of the highest.

## Choosing between the two clique classifiers

This function and
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
answer different questions; neither is deprecated in favour of the
other.

`classify_cliques()` works on the per-orthogroup *species* graph. It
asks which species are joined by conserved co-expression and whether
that pattern respects the trait split, and it returns one row per HOG.
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
commits to one best gene assignment per species clique, so a multi-copy
HOG still gets a single answer, and that answer is what
[`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md)
and
[`clique_threshold_sweep()`](https://paliocha.github.io/rcomplex/reference/clique_threshold_sweep.md)
consume – the `stability_class` / `robust` columns exist only on this
side.

[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
works on the *gene* graph built by
[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md).
It asks which individual gene copies are mutually conserved, so one HOG
can yield several overlapping cliques and the answer names paralogs
rather than species. It applies the taxonomy of Rodriguez et al. (2026),
plus `trait_specific`, with two explicit tolerance tiers
(`partial_significant` for weak wiring, `partial_present` for a missing
gene).

Reach for
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
when which copy sits in the conserved core matters, or when the
published taxonomy is what has to be reported. Reach for
`classify_cliques()` for a one-row-per-HOG trait summary wired into the
stability machinery.

## References

Rodriguez E, Birkeland S, Chapple ED, et al. (2026). Comparative
regulomics of wood formation across dicot and conifer trees. *Nature
Communications* 17(1).
[doi:10.1038/s41467-026-75624-2](https://doi.org/10.1038/s41467-026-75624-2)

## See also

[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
for the species-graph backend this wraps;
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
and
[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
for the copy-level alternative described above.

## Examples

``` r
if (FALSE) { # \dontrun{
clades <- list(annual = c("SP_A", "SP_B"), perennial = c("SP_C", "SP_D"))
result <- classify_cliques(edges, target_species, clades)
table(result$classification)
} # }
```
