# Sweep density thresholds for robustness assessment

Reruns the pairwise comparison at several threshold multipliers. Shows
how sensitive the conservation results are to the density cutoff.

## Usage

``` r
density_sweep(networks, ...)

# Default S3 method
density_sweep(
  networks,
  orthologs,
  multipliers = seq(0.95, 1.05, by = 0.01),
  method = c("permutation", "hypergeometric", "rank"),
  n_cores = 1L,
  use_torch = FALSE,
  pval_combine = c("max", "min"),
  seed = NULL,
  rho0 = NULL,
  null_networks = NULL,
  ...
)
```

## Arguments

- networks:

  Named list of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  outputs.

- ...:

  Additional arguments passed to the default method.

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog`.

- multipliers:

  Numeric vector of threshold multipliers (default
  `seq(0.95, 1.05, by = 0.01)`).

- method:

  Comparison method passed to
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  `"permutation"` (default), `"hypergeometric"` or `"rank"`.

- n_cores:

  Number of threads (default 1).

- use_torch:

  Logical; GPU acceleration for permutation method (default `FALSE`).

- pval_combine:

  Passed to
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  `"max"` (default) requires both directions to be significant (the
  reciprocal criterion of Netotea et al. (2014)); `"min"` calls a pair
  when either direction is significant.

- seed:

  Integer seed for the sweep's random draws, or `NULL` (default) to draw
  from the global RNG. Applied once here, so the whole sweep is one
  reproducible unit: the multipliers are visited in order and each
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  call continues the same stream, which keeps a level's draws
  independent of its neighbours' while making the sweep as a whole
  repeatable. Because nothing is drawn before the first multiplier, a
  sweep seeded with `s` reproduces `find_coexpressologs(seed = s)`
  exactly at its first level. See
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  for what a seed does and does not pin under each `method`.

- rho0:

  Passed to
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  reference fold enrichment for the analytical `power` column.

- null_networks:

  Passed to
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  (rank method only); every null's threshold is scaled by the same
  multiplier as its species' network.

## Value

A data frame with columns `multiplier`, `eff_density`, `n_significant`,
`edges` (list-column of edge data frames), and `species_densities`
(list-column of named numeric vectors giving the effective density per
species at each threshold level).

## Details

For each multiplier, every network's threshold is scaled by
`threshold * multiplier`, the effective density is recorded, and
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
is called on the modified networks.

## See also

[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
`classify_cliques`

## Examples

``` r
if (FALSE) { # \dontrun{
sweep <- density_sweep(networks, orthologs,
  method = "permutation", use_torch = TRUE
)
} # }
```
