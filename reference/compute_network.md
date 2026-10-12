# Compute co-expression network

Calculates correlation, applies normalization (Mutual Rank or CLR), and
determines a density-based co-expression threshold. Accepts either a
numeric matrix or a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Usage

``` r
compute_network(x, ...)

# S4 method for class 'matrix'
compute_network(
  x,
  cor_method = c("pearson", "spearman"),
  norm_method = c("MR", "CLR"),
  density = 0.03,
  sign = c("positive", "negative"),
  mr_log_transform = FALSE,
  sparse = TRUE,
  store_density = NULL,
  n_cores = 1L,
  block_size = NULL,
  partition = NULL
)

# S4 method for class 'SummarizedExperiment'
compute_network(x, assay = 1L, ...)
```

## Arguments

- x:

  Expression data: a numeric matrix (genes x samples) with row names as
  gene identifiers, or a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

- ...:

  Arguments passed to methods (see below).

- cor_method:

  Correlation method: `"pearson"` (default) or `"spearman"`.

- norm_method:

  Normalization method: `"MR"` (Mutual Rank, default) or `"CLR"`
  (Context Likelihood Ratio).

- density:

  Fraction of top edges to keep (default 0.03 = 3%).

- sign:

  `"positive"` (default) ranks the strongest positive correlations
  first. `"negative"` negates every correlation between two genes before
  normalization, so the strongest anticorrelations rank first. With
  `partition`, the function negates the correlation of each level before
  it sets negative values to zero.

- mr_log_transform:

  If `FALSE` (default), use the raw MR formula matching the original
  RComPlEx R Markdown. If `TRUE`, use Obayashi & Kinoshita (2009)
  log-normalized formula (values in \[0,1\]).

- sparse:

  If `TRUE` (default), return the network as a sparse
  `Matrix::dgCMatrix` holding both triangles of the thresholded
  co-expression matrix (entries below `store_threshold` are discarded;
  the dense n x n matrix never leaves this function). If `FALSE`, return
  the dense matrix as in rcomplex \< 0.2.0.

- store_density:

  Fraction of top edges to keep in the sparse store (default `NULL` =
  `max(density, 0.05)`). Must satisfy `density <= store_density < 1`;
  the margin above `density` lets downstream sweeps loosen the threshold
  (e.g.
  [`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md))
  without recomputing the network. Only valid with `sparse = TRUE`.

- n_cores:

  Number of threads for parallel computation (default 1).

- block_size:

  `NULL` (default) builds the dense n x n matrix first. A positive whole
  number builds the network a block of genes at a time and never forms
  the n x n matrix (each block holds 8 \* n \* `block_size` bytes, so a
  `block_size` near n saves nothing; values above n are capped at n with
  a message; threads rank the columns of one block, so keep it well
  above `n_cores`, as the 256-1024 benchmarked below): each gene keeps
  only its partners in the top fraction of its correlation ranks, which
  at the default `store_density` is 15% of all pairs held at 8 bytes
  each. On BDIS leaf data (20 samples, n = 20,000, 8 threads) peak
  memory was 1.2-1.4 GB against 5.1 GB for the dense build, in 7.6-7.9 s
  against 8.5 s. The result equals the dense build up to floating-point
  near-ties between correlations: the two builds compute correlations
  with different BLAS calls, so two correlations of one gene that differ
  only in the last bits can be ranked in the other order, shifting an MR
  value. With `cor_method = "spearman"` and few samples, correlations
  that are equal as rationals are common and may round apart, so the
  result can then differ from the dense build in many entries and
  between block sizes. Requires `sparse = TRUE` and
  `norm_method = "MR"`. With `mr_log_transform = TRUE` a pair is kept
  when either gene ranks the other in its top 10% (twice
  `store_density`) and a second correlation pass reads the other rank;
  on the same data peak memory was 1.7 GB against 5.1 GB dense in 9.2 s
  against 8.2 s. If the store threshold cannot prove the kept pairs
  complete the fraction widens and a message says so; peak memory grows
  with it (for log MR it can exceed the dense build), and at all pairs
  the build saves no memory.

- partition:

  `NULL` (default), or a per-sample factor of length `ncol(x)` without
  `NA`, such as tissue or study. The function computes the correlation
  within each level. It sets negative correlations to zero and takes the
  mean over the levels. Then it normalises the result. The function
  drops levels with fewer than 5 samples and tells you which. At least
  two levels must remain. Use the dense build: `block_size` must be
  `NULL`.

- assay:

  Assay name or index to extract from the SummarizedExperiment (default
  1). Requires the SummarizedExperiment package.

## Value

A list with components:

- network:

  Named symmetric matrix (genes x genes) with normalized co-expression
  values. A `dgCMatrix` storing both triangles of the entries at or
  above `store_threshold` (diagonal absent) when `sparse = TRUE`; a
  dense matrix with zero diagonal when `sparse = FALSE`.

- threshold:

  The co-expression threshold at the given density, equal to that of the
  full dense matrix in every mode.

- n_genes:

  Number of genes in the network.

- n_removed:

  Number of constant genes removed before computing correlations.

- params:

  List of parameters used. It holds `partition` when you give one.
  `r_threshold` is the weakest correlation of an edge: the smallest for
  `sign = "positive"`, the largest for `"negative"`. With `partition` it
  is the rectified mean over the levels, negated for `"negative"`. It is
  `NA` when no pair passes.

- store_density, store_threshold:

  Sparse networks only: the stored edge fraction and the value cutoff of
  the store. Analyses at thresholds below `store_threshold` are refused
  (see
  [`as_sparse_network()`](https://paliocha.github.io/rcomplex/reference/as_sparse_network.md)).

## Details

With `sign = "negative"` a gene keeps its own correlation of 1, so it
ranks itself first in both signs. Negative correlations are rarer and
weaker than positive ones in RNA-seq data. The density threshold still
keeps the top fraction of pairs, so a negative network exists at any
sample size even when it holds only noise: at n = 20 samples, r \>= -0.3
is noise. Check the weakest correlation that passed the threshold before
you read a negative network.

## Gene universe

Networks are built on all supplied genes and downstream tests use the
whole network as the hypergeometric population (Netotea et al. 2014).
The canonical ComPlEx implementations restrict the expression tables to
genes with an ortholog before building the networks, which changes
neighbourhoods, thresholds and calls; rcomplex reproduces them only when
the expression matrices are restricted the same way first (see
[`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)).

## Reconstructing sub-threshold values

A sparse network (`sparse = TRUE`) discards values below
`store_threshold`, so visualisations of a gene subset (e.g. module
heatmaps) cannot read them back from `net$network`.
[`mr_block()`](https://paliocha.github.io/rcomplex/reference/mr_block.md)
reconstructs the exact mutual-rank values for a gene subset from the
expression matrix and the network's parameters, without rebuilding the
dense n x n matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a matrix:
net <- compute_network(x,
  cor_method = "spearman",
  norm_method = "mr", density = 0.03
)

# From a SummarizedExperiment:
net <- compute_network(se, assay = "vst", cor_method = "spearman")
} # }
```
