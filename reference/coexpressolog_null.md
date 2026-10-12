# Degree-preserving edge-swap null for co-expressolog statistics

Tests whether an observed co-expressolog statistic exceeds what network
topology alone produces. Rewires each species network and runs
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
on the result. Compares the statistic with the null distribution.

## Usage

``` r
coexpressolog_null(
  networks,
  orthologs,
  swap_factor = 10L,
  n_cores = 1L,
  seed = NULL,
  ...
)
```

## Arguments

- networks:

  Named list of sparse network objects (`compute_network(sparse = TRUE)`
  outputs), keyed by species abbreviation. Dense networks are rejected;
  convert them with
  [`as_sparse_network()`](https://paliocha.github.io/rcomplex/reference/as_sparse_network.md).

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog`, as for
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md).

- swap_factor:

  Swap trials per permutation, rejected ones included, as a multiple of
  the edge count of each thresholded network (default 10). Must be a
  single finite number \> 0. The trial count is rounded up, so any
  positive factor makes at least one trial on a network with two or more
  edges; a network with fewer than two edges has no swap to make and is
  returned unchanged.

- n_cores:

  Number of parallel workers for the permutation loop (default 1).

- seed:

  Base seed for the run. `NULL` (default) draws one from the ambient RNG
  stream and leaves it advanced by that one draw, so consecutive
  unseeded calls give different nulls and a
  [`set.seed()`](https://rdrr.io/r/base/Random.html) in the caller's
  script reproduces the whole run. With a seed the call draws from a
  private stream and restores the caller's on exit — the package-wide
  contract, see `detect_modules`.

  What a default call guarantees is replayability, not a fixed null. The
  seed actually used — drawn or supplied — is recorded as
  `attr(result, "seed")`, and a drawn one is also announced in a
  message, so any run can be reproduced exactly by passing that value
  back as `seed`. Two default calls return different nulls; the size of
  that difference is what `null_se` and the `p_emp_lo` / `p_emp_hi`
  interval report. The default was `1L` before 0.3.0, which pinned the
  null of every default call while leaving the observed statistic free
  to drift, and hid the Monte Carlo error entirely.

  Two seeds exactly `2^31 - 1` apart (both legal, since a seed may lie
  anywhere in `+/- .Machine$integer.max`) give every rewiring
  permutation the same seed, so the null is identical for the two runs
  even though the observed statistic is not: see the aliasing note under
  `.task_seed()`.

- ...:

  Passed unchanged to
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  for both the observed and every null run (`pval_combine`, ...).
  Zero-overlap ortholog pairs are dropped before the q-values
  (`filter_zero = TRUE`): the statistic reads called edges only, so
  retaining them adds tens of thousands of never-called rows per pair to
  every permutation for nothing.

## Value

Data frame with one row per statistic and columns `statistic`,
`observed`, `null_mean`, `null_sd`, `null_max`, `fold` (observed /
null_mean; `NA` when the null mean is 0), `p_emp`
(`(1 + n_ge) / (n_perm + 1)`), `n_ge` (null draws at or above the
observed value), `null_se` (`null_sd / sqrt(n_perm)`, the Monte Carlo
error of `null_mean`; the relative error of `fold` is
`null_se / null_mean`, which is large whenever `null_mean` is near 0)
and `p_emp_lo` / `p_emp_hi`, an exact Clopper-Pearson 95\\ probability
that `p_emp` estimates — an interval on that probability, not on `p_emp`
itself. The `n_perm x k` matrix of null statistics is attached as
`attr(, "null")`, and the base seed the run used as `attr(, "seed")`.

## Details

Each network is binarised at its analysis threshold and rewired by
degree-preserving edge swaps. The swaps keep every gene's degree and
destroy the match between neighbourhoods and the ortholog mapping. The
null run takes exactly the same arguments (`...`) as the observed run.

Rewiring operates on the network thresholded at `net$threshold` (the
analysis density), not on the stored superset: sparse-store entries with
values at or above the threshold become edges, the rest are dropped. The
rewired networks are unweighted (every stored value is 1, with
`threshold = 1` and `store_threshold = 1`), so only membership-based
consumers are valid downstream; edge weights carry no information after
rewiring. A network with `m` edges gets `ceiling(swap_factor * m)` swap
trials, each the trial of igraph's `keeping_degseq()` rewiring: two
distinct edges are drawn uniformly, a swap that would create a loop or a
multi-edge is rejected, and the rejected trial still counts. Counting it
is what makes the chain sample the realizations of a degree sequence
uniformly. Networks stay simple (no loops, no multi-edges). Adjacency is
held as an `n x n` bit matrix during rewiring, `n^2 / 8` bytes per
network per worker (32 MB at 16 000 genes).

Permutation `b` runs in a worker that seeds itself from the base seed
and its own index (the package-wide per-task seed derivation), so
results are reproducible, independent of `n_cores`, and two neighbouring
base seeds share no rewiring. On Unix the permutations run under
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html); on
Windows they run serially. The seed covers the observed run as well as
the null runs, so settings that draw — the default randomized pi0, for
one — are pinned too. The statistic is the count of rows with
`type == "conserved"` per species pair (named
`paste(species1, species2, sep = "~")`) plus `"total"`, over 100 rewired
permutations. A null run whose rewiring leaves a species pair with no
overlap \> 0 rows records 0 conserved calls for that pair.
`method = "permutation"` in `...` is allowed but slow (a full
permutation test per rewired network).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- coexpressolog_null(networks, orthologs, n_cores = 4L)
res[res$statistic == "total", ]

# replay that exact run, whether or not it was seeded
same <- coexpressolog_null(networks, orthologs,
  n_cores = 4L, seed = attr(res, "seed")
)
} # }
```
