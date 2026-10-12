# Import a network built outside rcomplex

Wraps given co-expression weights in the network object that
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
returns. Every consumer then reads it unchanged.

## Usage

``` r
as_network(x, density = 0.03, genes = NULL)
```

## Arguments

- x:

  A symmetric numeric matrix, dense or sparse, with gene names as
  dimnames. Or an edge list: a data frame with columns `gene1`, `gene2`
  and `weight`.

- density:

  Fraction of gene pairs kept as edges, in (0, 1).

- genes:

  Character vector of genes, or `NULL` for all genes in `x`. It sets the
  gene order. Genes absent from an edge list get no edges.

## Value

A list as
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
returns with `sparse = TRUE`: `network` (a `dgCMatrix` with both
triangles), `threshold` (the weight at the `density` quantile),
`n_genes`, `params`, `store_density` and `store_threshold`.

## Details

Weights are taken as given, for example a TEA-GCN
`zScore(Co-exp_Str_MR)` or a WGCNA adjacency. `as_network()` does not
recompute mutual ranks; that is the job of
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
on expression. Larger weights mean stronger co-expression. A dense
matrix is stored as
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
stores it, at the `max(density, 0.05)` quantile. A sparse matrix or an
edge list is stored in full, and a pair it omits ranks below every given
weight. The threshold is the weight of the k-th strongest of all
`n (n - 1) / 2` pairs, with k as in
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).
It must be a given weight, so `x` must hold at least k pairs. Functions
that need expression, such as
[`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md),
do not work on the result.

## Examples

``` r
x <- matrix(rnorm(400), 40, 10, dimnames = list(paste0("g", 1:40)))
net <- as_network(compute_network(x)$network)
net$threshold
#> [1] 38.49675
```
