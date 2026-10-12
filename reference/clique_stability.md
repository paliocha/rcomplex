# Leave-k-out jackknife structural stability for cliques

Tests how robust each clique is to species removal. Removes k = 1, 2,
..., `max_k` species at a time and re-runs clique detection. Reports
whether matching gene assignments survive (Jaccard similarity).

## Usage

``` r
clique_stability(edges, ...)

# Default S3 method
clique_stability(
  edges,
  target_species,
  clades = NULL,
  all_species = target_species,
  full_cliques = NULL,
  max_k = length(all_species) - 2L,
  n_cores = 1L,
  ...
)
```

## Arguments

- edges:

  Data frame (same format as
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)).

- ...:

  Additional arguments passed to the default method.

- target_species:

  Character vector of species that define clique membership.

- clades:

  Optional named list of species vectors, one per clade. Clades may nest
  but must not cross. A species in no clade forms its own clade.
  Annotations use the top-level clades. If provided, the output gains
  clade annotations (`traits`, `sole_rep`). If `NULL` (default), clade
  columns are `NA`.

- all_species:

  Character vector of ALL species in the analysis universe (default:
  `target_species`). Leave-k-out subsets are drawn from `all_species`.
  Must be a superset of `target_species`. When larger than
  `target_species`, removing a non-target species tests whether the
  clique signal is robust to changes in the broader phylogenetic
  context.

- full_cliques:

  Output of
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md),
  or `NULL` to compute internally (default).

- max_k:

  Maximum number of species to leave out (default:
  `length(all_species) - 2`, leaving at least 2 species).

- n_cores:

  Number of OpenMP threads (default 1).

## Value

A list with components:

- stability:

  Data frame with columns: `clique_idx` (1-based), `hog`, `k`,
  `n_subsets`, `n_stable`, `stability_score`, `species_present`,
  `traits`, `sole_rep`. One row per (clique, k) pair.

- clique_disruption:

  Data frame with columns: `species`, `n_cliques_disrupted` (k=1 only),
  and optionally `trait_value` if `clades` is provided. One row per
  species in `all_species`.

- stability_class:

  Integer vector (length = number of cliques): highest k at which each
  clique is structurally stable across ALL subsets (0 = unstable at k=1)

- novel_cliques:

  Integer: total count of novel cliques across all subsets

## Details

All cliques are tested, whatever their clade composition. Clade
annotations are added after the test if `clades` is provided.

### How it works

For each combination of k species removed from `all_species` (k = 1,
..., max_k):

1.  All edges involving the removed species are dropped

2.  Clique detection re-runs among remaining target species

3.  Reduced cliques are matched to full-dataset cliques by Jaccard
    similarity of gene assignments (considering only non-removed
    species)

A clique is *testable* in a subset if at least 2 of its species remain
active. It is *stable at level k* if its gene assignments are preserved
across ALL C(N, k) subsets where it is testable.

Stability is purely structural — trait annotations are orthogonal and
added post-hoc from `clades` if provided.

### sole_rep column

The `sole_rep` column is `TRUE` if this clique has a single trait value
and that trait value has only one target species representative. Only
populated when `clades` is provided.

### full_cliques parameter

When `full_cliques = NULL` (default), cliques are computed internally
via
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md).
You can also precompute them:

    fc <- find_cliques(edges, target_species)
    stab <- clique_stability(edges, target_species,
                             full_cliques = fc)

## Examples

``` r
if (FALSE) { # \dontrun{
# Structural stability for all cliques
cliques <- find_cliques(edges, all_sp, min_species = 2L)
stab <- clique_stability(edges, all_sp,
  full_cliques = cliques
)

# With trait annotation
trait <- setNames(rep(c("annual", "perennial"), each = 4), all_sp)
stab <- clique_stability(edges, all_sp,
  clades = split(names(trait), trait),
  full_cliques = cliques
)

# Cliques surviving any single species dropout
k1 <- stab$stability[stab$stability$k == 1, ]
stable_cliques <- cliques[k1$clique_idx[k1$stability_score == 1], ]

# Multi-level: stability_class >= 2 survives any pair of dropouts
deeply_stable <- which(stab$stability_class >= 2)
} # }
```
