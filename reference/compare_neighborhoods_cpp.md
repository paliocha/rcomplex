# Compare co-expression neighborhoods across species (integer-indexed)

For each ortholog pair, tests the overlap of co-expression neighborhoods
in both directions (sp1-\>sp2 and sp2-\>sp1) using hypergeometric tests.
All gene identifiers are 0-based integer indices (string mapping done in
R).

## Usage

``` r
compare_neighborhoods_cpp(
  net1,
  net2,
  thr1,
  thr2,
  pair_sp1_idx,
  pair_sp2_idx,
  ortho_sp1_idx,
  ortho_sp2_idx,
  n_cores = 1L
)
```

## Arguments

- net1:

  Co-expression network for species 1 (n1 x n1 matrix)

- net2:

  Co-expression network for species 2 (n2 x n2 matrix)

- thr1:

  Co-expression threshold for species 1

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

DataFrame with comparison results for each ortholog pair. The
hypergeometric urn excludes the anchor gene (population n - 1, anchor
dropped from the ortholog-mapped set); `*.p_value_gt` / `*.p_value_eq`
are the ungated upper tail P(X \> x) and point mass P(X = x).
