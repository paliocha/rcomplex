# Sparse co-classification restricted to original network edges

Computes co-classification and per-pair excess only for gene pairs
connected in the thresholded co-expression network, reducing memory from
O(N^2) to O(\|E\|).

## Usage

``` r
build_sparse_coclassification_cpp(memberships, n_genes, edges, n_cores = 1L)
```

## Arguments

- memberships:

  List of K IntegerVectors (1-based module IDs)

- n_genes:

  Number of genes (N)

- edges:

  \|E\| x 2 IntegerMatrix of 0-based vertex indices

## Value

List with from, to (0-based), coclassification, excess, and expected
(length-K per-resolution expected scalars)
