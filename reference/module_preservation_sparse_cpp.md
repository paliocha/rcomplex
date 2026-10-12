# Module preservation permutation engine (sparse)

Module preservation permutation engine (sparse)

## Usage

``` r
module_preservation_sparse_cpp(
  p,
  i,
  x,
  thr,
  keep,
  module_members,
  ref_kIM,
  ref_CC,
  ref_MAR,
  n_perm,
  n_cores,
  binary,
  store_perm
)
```

## Arguments

- p, i, x:

  dgCMatrix slots.

- thr:

  Analysis threshold.

- keep:

  Ascending 0-based indices of the ortholog-mappable genes.

- module_members:

  List of integer vectors: local gene indices per module, ordered to
  match the reference vectors.

- ref_kIM, ref_CC, ref_MAR:

  Lists of numeric vectors: the reference network's per-gene statistics,
  one vector per module.

- n_perm:

  Number of permutations.

- n_cores:

  Number of OpenMP threads.

- binary:

  Treat every surviving edge as weight 1.

- store_perm:

  Return the raw per-permutation statistic pairs.

## Value

List with observed, perm_mean, perm_sd, n_perm_used, n_exceed, p_value,
p_npc, p_joint, n_joint, scale.
