# Build co-classification matrix with per-pair excess

Build co-classification matrix with per-pair excess

## Usage

``` r
build_coclassification_cpp(memberships, n_genes, return_excess = TRUE)
```

## Arguments

- memberships:

  List of K IntegerVectors (1-based module IDs)

- n_genes:

  Number of genes (N)

- return_excess:

  If TRUE, return excess (adaptive path). If FALSE, return raw
  co-classification (fixed threshold path). Only one N x N matrix is
  returned, saving ~4.6 GB for N = 24k.

## Value

List with one of coclassification or excess_coclassification, plus
expected (length-K per-resolution expected scalars)
