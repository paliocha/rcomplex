# Shuffled-partner null network

Rebuilds a network after permuting each gene's expression values across
samples. Keeps each gene's marginal distribution and the gene universe.
Destroys all co-expression.

## Usage

``` r
null_network(
  x,
  net,
  seed = NULL,
  n_cores = 1L,
  block_size = NULL,
  block = NULL
)
```

## Arguments

- x:

  Expression matrix (genes x samples, row names are gene identifiers) or
  a `SummarizedExperiment`, whose first assay is used. The same object
  `net` was built from.

- net:

  Network object from
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  built from `x`; its `params` (with `partition`) and `store_density`
  are reused so the null is built the same way.

- seed:

  Integer seed for the permutation, or `NULL` (default) to draw from the
  global RNG. A seeded call draws from a private stream and restores the
  caller's on exit, so it does not displace the caller by however many
  draws the permutation consumed; with `seed = NULL` the draws come from
  the ambient stream and leave it advanced. Same contract as
  `summarize_comparison`.

- n_cores, block_size:

  Passed to
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md);
  pass the `block_size` `net` was built with to keep the null build
  blockwise (it is not recorded in `params`).

- block:

  `NULL` (default) for the full per-gene shuffle, or a per-sample
  grouping (vector or factor of length `ncol(x)`, no `NA`), e.g. time
  point or tree, within which each gene is permuted.

## Value

A sparse network object, as
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
returns, on the same genes as `net`. `params$log10_perms_per_gene` is
the log10 size of one gene's permutation space
(`lgamma(n + 1) / log(10)` for the full shuffle, the sum over blocks
otherwise); with `block`, `params$block` holds the block sizes, named by
level.

## Details

It is the partner network for the empirical calibration of
`method = "rank"`. P-values from a comparison against it are the null
draws that
[`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md)
turns into empirical p-values. The rebuild uses the parameters recorded
in `net`.

With `block`, each gene's samples are permuted independently *within*
each block (one permutation per gene per block). The permutations must
differ between genes: one permutation applied to every gene only
reorders the samples and leaves the correlation matrix unchanged.

**Which null.** The full shuffle (`block = NULL`) is the null of no
co-expression at all. The within-block shuffle keeps each gene's block
means, so the shared time course (or tree effect) survives, and is the
null of no co-expression beyond the shared axis. That is the null a
cross-species statistic on raw networks should be read against: on the
Pooideae time series, cross-species neighbour AUROC is 0.54 for networks
shuffled within time point against 0.50 for the full shuffle, and that
0.04 is the shared time axis, not conserved wiring. The block null
leaves the random-geometric-graph geometry of the shared axis in place,
so within-species module statistics stay uninterpretable under it. After
regressing time out, one gene's residuals correlate -1/3 within a time
point of four samples, which is a second reason to permute within time
point rather than across it; the wiring layer from
[`split_layers()`](https://paliocha.github.io/rcomplex/reference/split_layers.md)
needs this null.
