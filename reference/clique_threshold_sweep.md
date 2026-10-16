# Structural survival of cliques across stricter density thresholds

Convenience wrapper that re-runs the full comparison-to-clique pipeline
([`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)
-\>
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
-\>
[`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)
-\>
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md))
at progressively stricter thresholds. For custom threshold logic, call
the individual functions directly.

## Usage

``` r
clique_threshold_sweep(
  cliques,
  target_species,
  networks,
  orthologs,
  multipliers = c(1.5, 2, 3, 5, 10),
  min_species = length(target_species),
  n_cores = 1L,
  seed = NULL
)
```

## Arguments

- cliques:

  Baseline output of
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md).

- target_species:

  Character vector of species abbreviations.

- networks:

  Named list of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  outputs, keyed by species abbreviation. Each element must have
  `$network` (named numeric matrix) and `$threshold` (scalar).

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog`.

- multipliers:

  Numeric vector of threshold multipliers (each \> 1). Default
  `c(1.5, 2, 3, 5, 10)`.

- min_species:

  Minimum species per clique (default `length(target_species)`).

- n_cores:

  Cores for
  [`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)
  (default 1).

- seed:

  Integer seed for the pi0 draws, or `NULL` (default) to draw from the
  global RNG. A seed makes the sweep reproducible.

## Value

A list with components:

- survival:

  Data frame with one row per (baseline clique, multiplier), including
  `multiplier = 1.0` baseline rows where all cliques trivially survive:
  `clique_idx` (1-based), `hog`, `multiplier`, `survived` (logical),
  `jaccard`, `n_species_orig`, `n_species_new`.

- sweep_cliques:

  Named list of
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  outputs keyed by multiplier.

- sweep_edges:

  Named list of combined edge data frames keyed by multiplier.

- persistence:

  Data frame with one row per baseline clique: `clique_idx`, `hog`,
  `birth` (lowest multiplier where the clique exists; 1.0 for all
  baseline cliques), `death` (first multiplier above birth where the
  clique is lost; `NA` if it survives all tested multipliers),
  `persistence` (`death - birth`; `NA` if death is `NA`).

## Details

No permutations are involved. The sweep uses the same q-values as
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
with randomized-p pi0.

## Examples

``` r
if (FALSE) { # \dontrun{
sweep <- clique_threshold_sweep(cliques, target_species, networks,
  orthologs,
  multipliers = c(1.5, 2, 5)
)
# Survival curve
sapply(
  sort(unique(sweep$survival$multiplier)),
  function(m) mean(sweep$survival$survived[sweep$survival$multiplier == m])
)
} # }
```
