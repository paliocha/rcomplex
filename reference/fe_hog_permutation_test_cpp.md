# Permutation test using precomputed fold-enrichment matrix

Permutation test using precomputed fold-enrichment matrix

## Usage

``` r
fe_hog_permutation_test_cpp(
  combined,
  hog_sp1_list,
  hog_sp2_list,
  test_greater,
  min_exceedances,
  max_permutations,
  n_cores
)
```

## Arguments

- combined:

  Precomputed combined FE matrix (n1 x n2), where combined\[a, b\] =
  FE1(b-\>a) + FE2(a-\>b).

- hog_sp1_list:

  List of 0-based sp1 gene indices per HOG

- hog_sp2_list:

  List of 0-based sp2 gene indices per HOG

- test_greater:

  If TRUE, test conservation (T \>= T_obs)

- min_exceedances:

  Besag-Clifford stopping parameter

- max_permutations:

  Maximum permutations per HOG

- n_cores:

  Number of OpenMP threads

## Value

DataFrame with T_obs, n_perm, n_exceed, p_value per HOG
