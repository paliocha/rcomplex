# Permutation-based HOG-level conservation/divergence test

Tests each hog for co-expression conservation or divergence. Permutes
gene identities and stops adaptively (Besag & Clifford, 1991).

## Usage

``` r
permutation_hog_test(
  net1,
  net2,
  comparison,
  alternative = c("greater", "less"),
  min_exceedances = 50L,
  max_permutations = 10000L,
  n_cores = 1L,
  use_torch = FALSE,
  seed = NULL
)
```

## Arguments

- net1:

  Network object for species 1 (output of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)).

- net2:

  Network object for species 2 (output of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)).

- comparison:

  Data frame from
  [`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md)
  — raw p-values, **not** FDR-corrected. Used to define HOG structure
  and extract effect sizes.

- alternative:

  Which tail to test: `"greater"` (default) for conservation or `"less"`
  for divergence.

- min_exceedances:

  Besag & Clifford stopping parameter: stop permuting once this many
  permutation statistics exceed (or fall below) T_obs. Higher values
  give more precise p-values but take longer (default 50).

- max_permutations:

  Maximum permutations per HOG (default 10000).

- n_cores:

  Number of threads for parallel computation (default 1).

- use_torch:

  If `TRUE`, precompute the fold-enrichment matrix on GPU via torch,
  then run permutations as fast table lookups. Requires the
  [torch](https://torch.mlverse.org/) package. Default `FALSE`.

- seed:

  Integer seed for the permutation draws, or `NULL` (default) to draw
  from the ambient RNG stream and leave it advanced. The C++ engine
  seeds one RNG per thread from R's stream before the parallel region,
  so a seed reproduces the run at `n_cores = 1` only: above that the
  OpenMP guided schedule decides which thread takes which HOG, and the
  pairing of thread streams to HOGs still moves between runs. A seeded
  call restores the caller's stream on exit — the package-wide contract,
  see
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

## Value

A data frame with one row per HOG, ordered by p-value, with columns:

- hog:

  HOG identifier

- n_pairs:

  Number of ortholog pairs (M x N) in this HOG

- n_sp1:

  Number of unique species-1 genes (M)

- n_sp2:

  Number of unique species-2 genes (N)

- T_obs:

  Observed sum-of-fold-enrichments statistic

- n_perm:

  Number of permutations performed (may be \< max_permutations due to
  adaptive stopping)

- n_exceed:

  Number of permutation statistics exceeding T_obs

- mean_eff:

  Mean geometric-mean effect size across pairs

- p_value:

  Permutation p-value: (n_exceed + 1) / (n_perm + 1)

- q_value:

  Discrete q-value (Liang, 2016) accounting for Besag-Clifford support

## Statistical method

The null hypothesis is that the specific gene identities in a HOG carry
no information about co-expression conservation: replacing the HOG's
genes with randomly chosen genes from the same networks would yield
equally large neighborhood overlap.

For each permutation, M random species-1 genes and N random species-2
genes are drawn (matching the HOG's gene counts), and the
sum-of-fold-enrichments statistic T is computed over all M x N pair x
direction combinations. Each fold enrichment is \\x / E\\ with \\E = m k
/ (n - 1)\\: the anchor gene leaves the ortholog-reachable set (\\k\\)
and the population (\\n - 1\\), as in
[`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md).
The permutation p-value is the fraction of permuted T values that exceed
(or fall below, for divergence) the observed T.

The Besag & Clifford (1991) sequential stopping rule terminates
permutations early once `min_exceedances` permutation statistics exceed
T_obs, providing efficient computation without sacrificing accuracy for
clearly significant or non-significant HOGs.

The internal option `options(rcomplex.force_flag_vector = TRUE)` forces
the C++ backends into the flag-vector intersection mode normally
reserved for networks with more than 100,000 genes (testing hook;
results are identical to the default bit-vector mode). It has no effect
when `use_torch = TRUE`.

## Why not Fisher's method

Fisher's method for combining p-values assumes independent tests. Within
a HOG, pair-level hypergeometric tests share network neighborhoods
(genes co-expressed with one HOG member are often co-expressed with
another) and share the ortholog mapping. This non-independence inflates
Fisher's statistic and produces anti-conservative p-values. The
permutation approach is exact regardless of the dependency structure.

## GPU acceleration

When `use_torch = TRUE`, the fold-enrichment matrix is precomputed on
GPU via GEMM (CUDA, MPS, or CPU fallback). This converts each
permutation from O(M*N*n_words) bit-vector intersections to O(M\*N)
table lookups — typically a 100-300x overall speedup. Requires the
[torch](https://torch.mlverse.org/) package.

On Apple Silicon (MPS), torch uses float32 because Metal does not
support float64. The fold-enrichment matrix is numerically exact in
float32 (all intermediate values are small integers), so
`permutation_hog_test(use_torch = TRUE)` is safe on MPS. However,
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
correlation in float32 can cause rank-swap artifacts in Spearman + MR
normalization. If you need exact Spearman + MR correlation on Apple
Silicon, use `compute_network(use_torch = FALSE)` with
`permutation_hog_test(use_torch = TRUE)` — the 100-300x permutation
speedup comes from the FE precomputation and table-lookup architecture,
not from GPU-accelerated correlation. On CUDA, both functions use
float64 with no precision tradeoff.

## Multiple testing correction

Q-values are computed using the Liang (2016) discrete q-value method via
[`DiscreteQvalue::DQ()`](https://rdrr.io/pkg/DiscreteQvalue/man/DQ.html),
which accounts for the discrete support of Besag & Clifford p-values.
The standard Storey & Tibshirani (2003) pi0-adaptive method cannot be
used here because the adaptive stopping rule produces p-values that are
valid (super-uniform under the null) but **not** uniformly distributed:
the negative binomial stopping concentrates null p-values in the 0.1–0.5
range and depletes the right tail (p \> 0.9), causing Storey's pi0
estimator to severely underestimate pi0. The Liang method estimates pi0
using the discrete support structure and is correctly calibrated (pi0 =
1 on null simulations).

## References

Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
*Biometrika*, 78(2), 301–304.
[doi:10.1093/biomet/78.2.301](https://doi.org/10.1093/biomet/78.2.301)

Liang, K. (2016). False discovery rate estimation for large-scale
homogeneous discrete p-values. *Biometrics*, 72(2), 639–648.
[doi:10.1111/biom.12429](https://doi.org/10.1111/biom.12429)

## Examples

``` r
if (FALSE) { # \dontrun{
hog_results <- permutation_hog_test(net1, net2, comparison,
  n_cores = 8L, use_torch = TRUE
)
significant <- hog_results[hog_results$q_value < 0.1, ]
} # }
```
