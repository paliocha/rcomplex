# Split expression into wiring and deployment layers

The function fits each gene on a per-sample block factor, such as time
point, tree or zone. The block means are the deployment layer, and the
residuals are the wiring layer. Build the network from the wiring layer
with
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).

## Usage

``` r
split_layers(x, block)
```

## Arguments

- x:

  Expression matrix (genes x samples) or a `SummarizedExperiment`, whose
  first assay is used (the default of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)).

- block:

  Per-sample factor (or vector) of length `ncol(x)`, no `NA`, with at
  least two levels and at least two samples per level.

## Value

A list of class `split_layers`:

- wiring:

  Residual matrix, dimnames of `x`. Input for
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md);
  pair it with `null_network(..., block = block)`.

- deployment:

  Genes x levels matrix of block means.

- r2:

  Named numeric, per-gene share of the sum of squares explained by
  `block` (not adjusted).

- df_residual:

  `ncol(x) - nlevels(block)`.

- block:

  List with `factor` (the block factor) and `sizes` (samples per level).

- params:

  List with `n_genes` and `n_samples`.

## Details

**What each layer is for.** The *wiring* layer (the residuals) holds
co-expression beyond the shared course: two genes that only follow the
same time course have no residual correlation. Build the network from it
with
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).
The *deployment* layer (the block means) says where along the course a
gene is expressed. `r2` is the share of each gene's sum of squares the
block explains (the unadjusted variance fraction; see the methods
article).

**Two facts not to miss.** The residuals of one gene sum to zero within
each block, so within a block of size \\n_b\\ they are correlated
\\-1/(n_b - 1)\\ (-1/3 for four replicates). A full shuffle would break
that structure, so the shuffled-expression null for a wiring network
must be `null_network(..., block = block)`. And the wiring layer keeps
only \\n_b - 1\\ degrees of freedom per block (`df_residual` in all, 15
on a 5 x 4 design of 20 samples), so each half of a split-half design
carries about 5 df: split-half replication is not informative there, and
the gate for a wiring network has to be cross-species.

**Which layer to cluster.** On the Pooideae leaf and root time series,
about 75 % of cross-species gene-level conservation (neighbour AUROC
above 0.5) survived in the wiring layer (leaf 0.78, root 0.71), and the
deployment layer conserved weakly. Cluster the wiring layer and carry
`r2` as a per-gene covariate.

## Examples

``` r
set.seed(1)
block <- factor(rep(c("T1", "T2", "T3"), each = 4))
x <- matrix(rnorm(10 * 12), 10, 12,
  dimnames = list(paste0("g", 1:10), paste0("s", 1:12))
)
sl <- split_layers(x, block)
sl
#> Split layers: 10 genes x 12 samples
#>   blocks: T1 (4), T2 (4), T3 (4) 
#>   df_residual: 9 
#>   median r2: 0.124 
net <- compute_network(sl$wiring, density = 0.2)
```
