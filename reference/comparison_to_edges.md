# Build clique edges from pairwise comparison results

Converts output from
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
into the edge format expected by `find_cliques` and `classify_cliques`.
Renames columns, injects species identity, and computes per-pair effect
sizes and classification.

## Usage

``` r
comparison_to_edges(
  comparison,
  species1,
  species2,
  alternative = c("greater", "less"),
  alpha = 0.1,
  pval_combine = c("max", "min"),
  rho0 = NULL,
  p0 = NULL
)
```

## Arguments

- comparison:

  Data frame from
  [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
  (the `$results` element). Must contain columns `gene1`, `gene2`,
  `hog`, plus q-value and effect-size columns from both directions.

- species1:

  Species abbreviation for `gene1` genes (e.g., `"BDIS"`).

- species2:

  Species abbreviation for `gene2` genes.

- alternative:

  Which test direction to use for q-values and classification:
  `"greater"` (conservation, default) or `"less"` (divergence).

- alpha:

  Significance threshold for the `type` column (default 0.1).

- pval_combine:

  How the two directional q-values are combined into `q_value`: `"max"`
  (default; both directions must be significant – the reciprocal
  criterion of Netotea et al. (2014), the `Max.p.val` filter of the
  original ComPlEx) or `"min"` (permissive; a pair is called when either
  direction is significant, yielding a denser edge supply for
  `find_cliques`). A missing directional q-value is ignored in either
  case.

- rho0:

  Reference fold enrichment for the `power` column: how many times more
  often than chance a conserved pair's partners are shared. `NULL`
  (default) takes, per direction, the median `effect_size` of the called
  pairs (about 2.5 on the Pooideae data); a single positive number fixes
  it for both directions.

- p0:

  Rank-test frames only (from
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md),
  carrying the AUROC grid): reference raw p for the `power` column, see
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md).
  `NULL` (default) takes, per direction, the median raw p of the pairs
  called and significant in that direction.

## Value

Data frame with columns:

- gene1:

  Gene identifier (gene1)

- gene2:

  Gene identifier (gene2)

- hog:

  Ortholog group identifier

- score:

  Evidence against chance in bits: `-log2(p_value)`.

- evalue:

  Expected number of pairs this significant by chance:
  `n_tests * p_value`.

- species1:

  Species abbreviation for gene1 (`species1`)

- species2:

  Species abbreviation for gene2 (`species2`)

- p_value:

  The two directional p-values combined like `pval_combine`. Rank-test
  frames use the calibrated p (`*.p.emp`) when present. `NA` when the
  comparison has no p-value columns.

- n_tests:

  Number of ortholog pairs tested, in the direction with more tests.

- q_value:

  Maximum (or minimum, see `pval_combine`) of the two directional
  q-values

- effect_size:

  Geometric mean of directional effect sizes

- jaccard:

  Geometric mean of directional Jaccard indices
  (`sqrt(species1.jaccard * species2.jaccard)`). Range \\\[0, 1\]\\: 1 =
  identical neighborhoods, 0 = disjoint.

- power:

  Detection power: the probability that the pair would have been called
  had its partners been shared at the reference fold enrichment `rho0`,
  given its neighbourhood sizes and the largest p-value called in each
  direction. Directions combine like `pval_combine` (`"max"`: the weaker
  power). Power rises with degree: at typical enrichment a low-degree
  gene's shared partners amount to less than one, which no count test
  can call, so a non-significant low-power pair is absent evidence, not
  evidence of divergence. `classify_gene_cliques` and `classify_cliques`
  read it through `min_power`. `NA` for `alternative = "less"`, when no
  pair is called in a direction, or when the comparison lacks the
  neighbourhood-size columns. For a rank-test frame (one carrying
  `*.auroc.grid`) it is the rank-test power instead, at reference rank
  `p0` (see
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md));
  there a direction with no call gives 0 rather than `NA`, so the
  classifiers read it as uninformative, not as a rejection.

- type:

  `"conserved"` or `"diverged"` if `q_value < alpha`; `"ns"` otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
summary_AB <- summarize_comparison(comparison_AB)
edges_AB <- comparison_to_edges(summary_AB$results, "SP_A", "SP_B")

# Combine multiple species pairs for find_cliques()
edges <- rbind(edges_AB, edges_AC, edges_BC)
} # }
```
