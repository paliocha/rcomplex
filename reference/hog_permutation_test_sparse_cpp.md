# Permutation-based HOG-level conservation test (sparse networks)

Same as
[`hog_permutation_test_cpp()`](https://paliocha.github.io/rcomplex/reference/hog_permutation_test_cpp.md)
but takes the slots of a `dgCMatrix` (column-compressed, both triangles
stored) for each network instead of a dense matrix. Column j lists the
neighbours of gene j.

## Usage

``` r
hog_permutation_test_sparse_cpp(
  p1,
  i1,
  x1,
  thr1,
  p2,
  i2,
  x2,
  thr2,
  ortho_sp1_idx,
  ortho_sp2_idx,
  hog_sp1_list,
  hog_sp2_list,
  test_greater,
  min_exceedances,
  max_permutations,
  n_cores,
  force_flag_mode = FALSE
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

- ortho_sp1_idx:

  0-based net1 indices for full ortholog table

- ortho_sp2_idx:

  0-based net2 indices for full ortholog table

- hog_sp1_list:

  List of integer vectors: unique 0-based sp1 indices per HOG

- hog_sp2_list:

  List of integer vectors: unique 0-based sp2 indices per HOG

- test_greater:

  If TRUE, test conservation (T \>= T_obs); if FALSE, divergence

- min_exceedances:

  Besag-Clifford stopping parameter (default 50)

- max_permutations:

  Maximum permutations per HOG (default 10000)

- n_cores:

  Number of OpenMP threads (default 1)

- force_flag_mode:

  If TRUE, force the flag-vector intersection mode regardless of network
  size (internal testing hook; default FALSE)

## Value

DataFrame with T_obs, n_perm, n_exceed, p_value per HOG
