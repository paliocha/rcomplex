# Blockwise sparse mutual-rank network

Blockwise sparse mutual-rank network

## Usage

``` r
mr_block_network_cpp(
  zt,
  log_transform,
  negate,
  density,
  store_density,
  block_size,
  n_cores
)
```

## Arguments

- zt:

  Standardised expression (samples x genes) such that `crossprod(zt)` is
  the correlation matrix.

- log_transform, negate:

  As in
  [`mutual_rank_inplace_cpp()`](https://paliocha.github.io/rcomplex/reference/mutual_rank_inplace_cpp.md).

- density, store_density:

  Analysis and store densities.

- block_size:

  Columns per correlation block.

- n_cores:

  Number of OpenMP threads.

## Value

List with dgCMatrix slots `i`, `p`, `x`, the `threshold` and
`store_threshold`, the number of joined candidate pairs `n_candidates`,
the rank `fraction` that sufficed and the `start_fraction` the build
began with.
