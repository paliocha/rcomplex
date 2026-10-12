# Permutation-based HOG-level conservation test

Tests each HOG for co-expression conservation using a gene-identity
permutation null with adaptive stopping (Besag & Clifford, 1991).

## Usage

``` r
hog_permutation_test_cpp(
  net1,
  net2,
  thr1,
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

- net1:

  Co-expression network matrix for species 1 (n1 x n1)

- net2:

  Co-expression network matrix for species 2 (n2 x n2)

- thr1:

  Co-expression threshold for species 1

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
