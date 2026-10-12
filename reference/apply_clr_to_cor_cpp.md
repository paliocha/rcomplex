# Apply CLR transformation to a correlation matrix

Normalizes a correlation matrix using the Context Likelihood Ratio (CLR)
transformation (Faith et al. 2007). Computes row z-scores, clamps
negatives to zero, and combines row and column z-scores.

## Usage

``` r
apply_clr_to_cor_cpp(cor_matrix, n_cores = 1L)
```

## Arguments

- cor_matrix:

  Symmetric correlation matrix (n x n)

- n_cores:

  Number of OpenMP threads (default: 1)

## Value

CLR-transformed matrix (n x n), symmetric, non-negative

## Details

The CLR transformation:

    z_row(i,j) = (cor(i,j) - mean(row_i)) / sd(row_i)
    z_col(i,j) = (cor(i,j) - mean(col_j)) / sd(col_j)
    CLR(i,j) = sqrt(max(0, z_row)^2 + max(0, z_col)^2)

This matches the original RComPlEx Rmd CLR implementation which uses R's
scale() function (column-wise z-scoring on the symmetric matrix,
equivalent to row-wise z-scoring due to symmetry).
