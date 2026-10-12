# Summarize a specificity comparison

Computes q-values for the per-direction specificity p-values from
[`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md).
Optionally calibrates them against a shuffled partner null. Returns the
same gene-pair, gene and hog counts as
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md).

## Usage

``` r
summarize_specificity(
  comparison,
  null_p = NULL,
  alpha = 0.1,
  pi0_method = c("storey", "none"),
  species1 = NULL,
  species2 = NULL,
  pval_combine = c("max", "min"),
  p0 = NULL
)
```

## Arguments

- comparison:

  Data frame from
  [`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md).

- null_p:

  `NULL`, or `list(species1 = , species2 = )` of null p-values for the 1
  -\> 2 and 2 -\> 1 directions from the shuffled partner. `NA` entries
  are dropped.

- alpha:

  Significance threshold applied to q-values (default 0.1).

- pi0_method:

  `"storey"` (default) estimates pi0 from the p-values; `"none"` fixes
  pi0 = 1 (Benjamini-Hochberg).

- species1, species2:

  Optional species abbreviations. When both are provided,
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)
  is called internally and the result is returned as a third list
  element `$edges`. This avoids a separate
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)
  call when preparing input for `find_cliques`.

- pval_combine:

  Passed to
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md)
  when `species1` and `species2` are given: `"max"` (default; both
  directions significant – the reciprocal criterion of Netotea et al.
  (2014)) or `"min"` (either direction).

- p0:

  Reference raw p for the `power` column of `$edges`: the rank among its
  candidates at which a conserved ortholog is assumed to sit, as a
  fraction in (0, 1); values below 1e-5 are read as 1e-5, the edge of
  the AUROC grid. `NULL` (default) uses, per direction, the median raw p
  of the pairs called and significant in that direction – the rank-test
  counterpart of `rho0`.

## Value

A list with components:

- results:

  The complete rows of `comparison` with `species1.p.emp` /
  `species2.p.emp` (when `null_p` is given) and `species1.q_value_con` /
  `species2.q_value_con` added.

- summary:

  The counts of
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md),
  `pi0` per direction, `n_null` (null draws used per direction, 0
  without `null_p`) and `n_dropped` (rows with an `NA` p-value).

- edges:

  (Only when `species1` and `species2` are provided.) Edge frame from
  [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md).
  `power` is the probability that the pair would have been called had
  the ortholog ranked at raw p `p0` among its candidates: per direction,
  the anchor's `auroc.grid` (log-linear interpolation) gives the AUROC
  at `p0` and at the largest called raw p, and the ortholog's AUROC over
  its `mapped` genes is taken as normal around the former with the
  Hanley-McNeil (1982) standard error. It is at least 0.5 when `p0` lies
  below the call threshold (always so for the default), and like the
  hypergeometric power it orders pairs by how detectable conservation is
  rather than giving calibrated rates. Directions combine like
  `pval_combine`. A direction with no call gives power 0 for that
  direction (nothing was detectable there, so a miss is not evidence):
  the edge power under `"max"`, while `"min"` takes the other direction.
  The clique classifiers read it through `min_power`.

## Details

Rows where either direction is `NA` (an anchor gene with no mapped
neighbours) are dropped before testing.

## Multiple testing correction

With `null_p`, each direction's p-values are first made empirical
against the p-values of the same comparison run against a shuffled
partner network
([`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)),
pooled over pairs: `p.emp = (1 + #{null p <= p}) / (1 + n_null)`.
Q-values are then Storey q-values on the empirical (or, without
`null_p`, the raw) p-values, with Benjamini-Hochberg as the fallback
when the pi0 fit fails.

No discrete FDR correction is applied. All tests of one direction share
one support, the `1/n2` grid of the rank p-value, so the Liang discrete
q-values of `DiscreteQvalue` gain nothing over Storey, and the
Doehler-Durand-Roquain procedures, which exploit supports that differ
between tests (as the hypergeometric p-values of
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
do), have nothing to exploit. The raw rank p-values are not uniform
under the null; the empirical step against the shuffled partner
calibrates them before Storey is applied.

## References

Hanley, J. A. & McNeil, B. J. (1982). The meaning and use of the area
under a receiver operating characteristic (ROC) curve. *Radiology*,
143(1), 29–36.
[doi:10.1148/radiology.143.1.7063747](https://doi.org/10.1148/radiology.143.1.7063747)

Suresh, H., Crow, M., Jorstad, N., Hodge, R., Lein, E., Dobin, A.,
Bakken, T. & Gillis, J. (2023). Comparative single-cell transcriptomic
analysis of primate brains highlights human-specific regulatory
evolution. *Nature Ecology & Evolution*, 7(11), 1930–1943.
[doi:10.1038/s41559-023-02186-7](https://doi.org/10.1038/s41559-023-02186-7)

Storey, J. D. & Tibshirani, R. (2003). Statistical significance for
genomewide studies. *Proceedings of the National Academy of Sciences*,
100(16), 9440–9445.
[doi:10.1073/pnas.1530509100](https://doi.org/10.1073/pnas.1530509100)
