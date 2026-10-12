# Compare co-expression neighborhoods across species (sparse networks)

Same as
[`compare_neighborhoods_cpp()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods_cpp.md)
but takes the slots of a `dgCMatrix` (column-compressed, both triangles
stored) for each network instead of a dense matrix. Column j lists the
neighbours of gene j.

## Usage

``` r
compare_neighborhoods_sparse_cpp(
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
  n_cores = 1L
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

- n_cores:

  Number of OpenMP threads (default: 1)

## Value

DataFrame with comparison results for each ortholog pair
