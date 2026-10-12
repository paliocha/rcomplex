# Neighbourhood specificity of ortholog pairs (sparse networks)

Takes the `dgCMatrix` slots of both networks (both triangles stored, no
diagonal, symmetric pattern) as
[`compare_neighborhoods_sparse_cpp()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods_sparse_cpp.md)
does. `pair_*` are the tested pairs, `ortho_*` the full ortholog
mapping; all indices 0-based.

## Usage

``` r
specificity_sparse_cpp(
  p1,
  i1,
  x1,
  thr1,
  p2,
  i2,
  x2,
  thr2,
  pair_sp1_idx,
  pair_sp2_idx,
  ortho_sp1_idx,
  ortho_sp2_idx,
  do_12,
  do_21,
  n_cores,
  grid_frac
)
```

## Arguments

- p1:

  `@p` slot of net1 (column pointers, length n1 + 1)

- i1:

  `@i` slot of net1 (0-based row indices)

- x1:

  `@x` slot of net1 (stored values)

- thr1:

  Co-expression threshold for species 1

- p2:

  `@p` slot of net2

- i2:

  `@i` slot of net2

- x2:

  `@x` slot of net2

- thr2:

  Co-expression threshold for species 2

- pair_sp1_idx:

  0-based index into net1 for each ortholog pair

- pair_sp2_idx:

  0-based index into net2 for each ortholog pair

- ortho_sp1_idx:

  0-based net1 indices for full ortholog table

- ortho_sp2_idx:

  0-based net2 indices for full ortholog table

- do_12, do_21:

  Compute direction 1 -\> 2 / 2 -\> 1.

- n_cores:

  Number of OpenMP threads (default: 1)

- grid_frac:

  Ascending fractions f in (0, 1\]; for each, the AUROC grid holds the
  ceil(f \* n_b)-th largest candidate AUROC of the anchor. Empty skips
  the grid.

## Value

List with, per requested direction, `species1.neigh`, `species1.mapped`,
`species1.auroc`, `species1.p_value`, `species1.jaccard`,
`species1.n.cand` (one element per pair) and the matrix
`species1.auroc.grid` (pairs x fractions), and the `species2.*` set.
