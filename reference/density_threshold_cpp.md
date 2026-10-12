# Compute density threshold from a symmetric matrix

Extracts the upper triangle (excluding diagonal) and finds the value at
the given density percentile using O(n) partial sort.

## Usage

``` r
density_threshold_cpp(mat, density)
```

## Arguments

- mat:

  Symmetric matrix (n x n)

- density:

  Fraction of top edges to keep (e.g., 0.03 for 3%)

## Value

The threshold value such that `density` fraction of upper triangle
values are \>= the threshold
