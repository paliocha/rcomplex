# Per-gene intramodular statistics (sparse)

Per-gene intramodular statistics (sparse)

## Usage

``` r
module_gene_stats_sparse_cpp(p, i, x, thr, keep, module_members, binary)
```

## Arguments

- p, i, x:

  dgCMatrix slots.

- thr:

  Analysis threshold.

- keep:

  Ascending 0-based indices of the genes to induce on.

- module_members:

  List of integer vectors: local gene indices per module.

- binary:

  Treat every surviving edge as weight 1.

## Value

List with kIM, CC, MAR (one value per induced gene) and scale.
