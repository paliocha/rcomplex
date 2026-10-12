# Reduce orthogroups by merging correlated paralogs

Merges the paralogs of each hog whose Pearson correlation exceeds
`cor_threshold`. Clusters them with Ward.D2 agglomeration. Replaces each
merged set by its averaged expression profile.

## Usage

``` r
reduce_orthogroups(
  expr_matrix,
  orthologs,
  gene_col = "gene1",
  cor_threshold = 0.7
)
```

## Arguments

- expr_matrix:

  Numeric matrix (genes x samples) with gene identifiers as row names.

- orthologs:

  Data frame with columns `gene1` (or the column matching gene row
  names), `gene2`, and `hog`, as returned by
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md).

- gene_col:

  Character: which column of `orthologs` contains gene IDs matching row
  names of `expr_matrix` (default `"gene1"`).

- cor_threshold:

  Pearson correlation threshold for merging paralogs within a HOG
  (default 0.7). Higher values are more conservative (fewer merges).

## Value

A list with components:

- expr_matrix:

  Reduced expression matrix (genes x samples) with row names. Merged
  genes have averaged expression.

- gene_map:

  Data frame with columns `original` and `representative`, mapping every
  original gene to its representative in the reduced matrix.

- n_original:

  Number of genes before reduction.

- n_reduced:

  Number of genes after reduction.

- n_merged:

  Number of genes absorbed into representatives.

## Details

This is an optional preprocessing step before
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).
It reduces redundancy from recent duplications where paralogs retain
nearly identical expression patterns, shrinking the expression matrix
and avoiding combinatorial blowup in downstream clique detection.

Subfunctionalized paralogs (distinct expression programs) are preserved
as separate clusters. Zero-variance genes and genes not assigned to any
HOG are kept as-is.

## Examples

``` r
if (FALSE) { # \dontrun{
reduced <- reduce_orthogroups(expr_matrix, orthologs)
reduced$expr_matrix # reduced expression matrix
reduced$gene_map # original -> representative mapping
} # }
```
