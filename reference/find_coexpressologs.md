# Find co-expressologs across all species pairs

Runs the comparison pipeline for each species pair. Combines the results
into one edge table for
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md).

## Usage

``` r
find_coexpressologs(networks, ...)

# Default S3 method
find_coexpressologs(
  networks,
  orthologs,
  method = c("hypergeometric", "rank", "permutation"),
  n_cores = 1L,
  use_torch = FALSE,
  pval_combine = c("max", "min"),
  filter_zero = FALSE,
  seed = NULL,
  rho0 = NULL,
  null_networks = NULL,
  ...
)

run_pairwise_comparisons(...)
```

## Arguments

- networks:

  Named list of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  outputs, keyed by species abbreviation.

- ...:

  Additional arguments passed to the default method.

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog` (from
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)).

- method:

  Testing method: `"hypergeometric"` (default, fast), `"rank"` or
  `"permutation"` (rigorous). `"rank"` scores each ortholog pair with
  [`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md)
  and calibrates it against `null_networks` via
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md);
  `rho0` is refused; `power` is the rank-test power (see
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md)).

- n_cores:

  Number of threads (default 1).

- use_torch:

  Logical. Use GPU-accelerated fold-enrichment precomputation for
  permutation method (default `FALSE`).

- pval_combine:

  Hypergeometric and rank methods: how the two directional q-values are
  combined, passed to
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md):
  `"max"` (default; both directions significant – the reciprocal
  criterion of Netotea et al. (2014), the `Max.p.val` filter of the
  original ComPlEx) or `"min"` (permissive; either direction, denser
  edge supply for `find_cliques`).

- filter_zero:

  Analytical method only: passed to
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md).
  `FALSE` (default) keeps every tested ortholog pair, including those
  whose neighbourhood overlap is zero in either direction, so a pair
  that could not have been called still reaches the edge table carrying
  its `power`. `TRUE` removes them before the q-values are computed,
  which is what canonical ComPlEx does – it shrinks the multiple-testing
  set, so every q-value moves. Set by
  [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md).

- seed:

  Integer seed for the call's random draws, or `NULL` (default) to draw
  from the global RNG. Seeding here makes the call reproducible, and
  with it every downstream count thresholded on `q_value`. The seed is
  applied once, for the whole species-pair loop, rather than being
  handed to each
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
  call: one seeded stream flows through the pairs in order, so each pair
  draws its own uniforms instead of every pair reusing the same ones.

  What a seed buys differs by method. Under `"hypergeometric"` it pins
  the randomized-p draws behind pi0, and the result is reproducible at
  any `n_cores`. Under `"permutation"` it pins the per-thread seeds that
  the C++ engines draw from R's RNG, which reproduces the run at
  `n_cores = 1` only: above that the OpenMP guided schedule decides
  which thread takes which HOG, so the pairing of thread streams to HOGs
  still moves between runs.

  A seeded call draws from a private stream and restores the caller's on
  exit; with `seed = NULL` the draws come from the ambient stream and
  leave it advanced. Same contract as `detect_modules` and
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md).

- rho0:

  Analytical method only: reference fold enrichment for the `power`
  column, passed to
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)
  (default `NULL`, the median effect size of the called pairs). The
  permutation path always reports `power = NA`: its HOG-level q-value
  has no per-pair call threshold to measure power against.

- null_networks:

  Specificity method only (required there): named list with one
  [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  object, or a list of them, per species.

## Value

Data frame with columns `gene1`, `gene2`, `hog`, `score`, `evalue`,
`q_value`, `effect_size`, `power`, `species1`, `species2`, `p_value`,
`n_tests`, `jaccard`, `type` (`power` is `NA` under
`method = "permutation"`). `score` is `-log2(p_value)` in bits. `evalue`
is `n_tests * p_value`, the expected number of chance hits this strong.
`p_value` combines the two directions like `pval_combine`. `n_tests`
counts the ortholog pairs tested in that species pair. Ready for
`find_cliques` or `classify_cliques`.

## Details

`score` and `evalue` are BLAST's bit score and E-value, related as
`evalue = n_tests * 2^-score`. Under `method = "permutation"` the
p-value is the HOG's, shared by its pairs, and it cannot fall below
`1 / (n_perm + 1)`, so `score` has a ceiling there (see
[`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)).

## Methods

- hypergeometric:

  Fast path: pair-level Storey q-values via
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md).
  Appropriate when most HOGs are single-copy.

- rank:

  Opt-in:
  [`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md)
  against shuffled-partner networks from
  [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  (argument `null_networks`), q-values via
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md).
  The hypergeometric test counts how many of gene A's translated
  partners are among ortholog B's partners, against an urn of random
  draws. With small neighbourhoods (top-25 lists) that urn let
  technically "sticky" genes, such as genes with many zeros, through as
  false calls; at network density 0.03 it held against a
  shuffled-expression null. The rank test asks where A's translated
  partners sit in B's full co-expression ranking, checks every other
  gene of B's species the same way, and reports B's rank in that
  contest, so a list that thousands of genes recognise as well as B is
  not evidence. The q-value is read against a partner whose expression
  was shuffled within each gene. `effect_size` is
  `sqrt(auroc12 * auroc21)`; `power` is the rank-test power of
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md)
  at its default reference rank.

- permutation:

  Rigorous path: gene-identity permutation via
  [`permutation_hog_test()`](https://paliocha.github.io/rcomplex/reference/permutation_hog_test.md)
  with Besag-Clifford adaptive stopping and Liang discrete q-values.
  Required for multi-copy HOGs where pair-level tests are correlated.

## Gene universe

Each network is built on all genes supplied to
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
and the hypergeometric population is the whole network (Netotea et al.
2014). Canonical ComPlEx implementations restrict expression to ortholog
genes before building the networks; that changes neighbourhoods and
calls, so rcomplex matches them only when the expression matrices are
restricted the same way beforehand (see
[`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
# Fast analytical path; seed pins the randomized pi0 draws
edges <- find_coexpressologs(networks, orthologs, seed = 1)

# Permissive either-direction calls (denser clique input)
edges <- find_coexpressologs(networks, orthologs, pval_combine = "min")

# Rigorous permutation path with GPU
edges <- find_coexpressologs(networks, orthologs,
  method = "permutation", use_torch = TRUE, n_cores = 4L
)
} # }
```
