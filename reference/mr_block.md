# Reconstruct mutual-rank values for a gene subset

Computes the exact mutual rank values for a subset of genes. Uses the
stored parameters of the network (`cor_method`, `sign`,
`mr_log_transform`). Recovers values that a sparse network discards.

## Usage

``` r
mr_block(x, genes, net)
```

## Arguments

- x:

  Expression matrix (genes x samples) the network was built from. Must
  contain every network gene as a row; extra rows (e.g. genes removed as
  constant) are ignored.

- genes:

  Character vector of gene identifiers (no duplicates), all present in
  `rownames(net$network)`.

- net:

  Network object from
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  (sparse or dense) with `norm_method = "MR"`.

## Value

A k x k numeric matrix of mutual-rank values with
`dimnames = list(genes, genes)` and zero diagonal, on the same scale as
the network (raw MR, or log-normalized when the network was built with
`mr_log_transform = TRUE`).

## Details

A sparse network (`compute_network(sparse = TRUE)`) discards values
below its `store_threshold`. This function reconstructs them, for
example for a module heatmap, without rebuilding the dense n x n matrix.
The value \\MR\_{ij} = \sqrt{R\_{ij} R\_{ji}}\\ only needs the
correlations of genes i and j ranked over all n network genes, so a k x
n correlation slice suffices.

The gene universe is `rownames(net$network)`: `x` is subset to those
rows first (the constant-gene filter was already applied when the
network was built), and ranks span exactly the network genes. Pass the
same expression matrix that built the network; the reconstructed block
then matches the dense `compute_network(sparse = FALSE)` matrix.

Reconstruction is exact up to ties in the correlation values: exactly
tied correlations get average ranks, and a mathematical tie can be
broken differently by this function
([`stats::cor()`](https://rdrr.io/r/stats/cor.html)) and the network's
correlation backend
([`Rfast::cora()`](https://rdrr.io/pkg/Rfast/man/cova.html)), shifting a
rank by 0.5. Ties are measure-zero for Pearson on continuous data, but
Spearman correlations lie on a finite grid where occasional ties are
possible.

## Examples

``` r
if (FALSE) { # \dontrun{
net <- compute_network(x, density = 0.03) # sparse
# exact values for one module, including sub-threshold entries:
blk <- mr_block(x, module_genes, net)
heatmap(blk)
} # }
```
