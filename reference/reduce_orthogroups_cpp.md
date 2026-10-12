# Reduce orthogroups by merging correlated paralogs

For each HOG with multiple paralogs, computes pairwise Pearson
correlation on the expression matrix, clusters with Ward.D2, and merges
paralogs whose correlation exceeds `cor_threshold` into a single
averaged representative.

## Usage

``` r
reduce_orthogroups_cpp(expr, hog_members, non_hog_idx, cor_threshold)
```

## Arguments

- expr:

  Numeric matrix (genes x samples). Row order matches gene_ids.

- hog_members:

  List of integer vectors. Each element contains 1-indexed row positions
  in `expr` for genes belonging to one HOG.

- non_hog_idx:

  Integer vector of 1-indexed row positions for genes not assigned to
  any HOG (kept as-is).

- cor_threshold:

  Pearson correlation threshold for merging (default 0.7).

## Value

List with: expr_matrix (reduced matrix), out_row_source (1-indexed
original row for each output row), map_from / map_to (1-indexed gene
mapping), n_original, n_reduced, n_merged.
