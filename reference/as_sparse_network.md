# Convert a dense network object to the sparse representation

Thresholds the dense co-expression matrix at the `store_density`
quantile. Repacks the surviving off-diagonal entries (both triangles) as
a `Matrix::dgCMatrix`, as `compute_network(sparse = TRUE)` does.

## Usage

``` r
as_sparse_network(net, store_density = 0.05)
```

## Arguments

- net:

  Dense network object: output of `compute_network(sparse = FALSE)`, or
  a hand-built `list(network = <matrix>, threshold = )` with gene
  dimnames.

- store_density:

  Fraction of top edges to keep in the store (default 0.05). Must be in
  (0, 1) and, when the network records its density in `params$density`,
  at least that density (the store must contain every analysis edge).

## Value

The network object with `network` as a `dgCMatrix` plus `store_density`
and `store_threshold` fields.

## Details

All other fields of `net` are carried over. `store_density`,
`store_threshold` and `params$store_density` are added. Values below
`store_threshold` are discarded, so downstream analyses at a lower
threshold are refused.

## Examples

``` r
if (FALSE) { # \dontrun{
net <- compute_network(x, density = 0.03, sparse = FALSE)
net_sparse <- as_sparse_network(net, store_density = 0.05)
} # }
```
