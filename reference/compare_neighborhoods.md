# Compare co-expression neighborhoods across species

For each ortholog pair, tests whether co-expression neighborhoods are
conserved between species using hypergeometric tests in both directions.

## Usage

``` r
compare_neighborhoods(net1, net2, orthologs, n_cores = 1L)
```

## Arguments

- net1:

  Network object for species 1 (output of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)).

- net2:

  Network object for species 2 (output of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)).

- orthologs:

  Data frame with columns `gene1`, `gene2`, and `hog` (output of
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)).

- n_cores:

  Number of threads for parallel computation (default 1).

## Value

A data frame with columns:

- gene1:

  Gene identifier for species 1

- gene2:

  Gene identifier for species 2

- hog:

  Ortholog group identifier

- species1.neigh:

  Number of neighbors of gene1 gene in net1

- species1.ortho.neigh:

  Number of ortholog-mapped neighbors from net2, excluding the gene1
  gene itself

- species1.neigh.overlap:

  Intersection size

- species1.p_value_con:

  Upper-tail hypergeometric p-value for conservation (direction 1); 1
  when the overlap is 0 or 1

- species1.p_value_div:

  Lower-tail hypergeometric p-value for divergence (direction 1)

- species1.p_value_gt:

  \\P(X \> x)\\, ungated (direction 1)

- species1.p_value_eq:

  \\P(X = x)\\ (direction 1)

- species1.effect_size:

  Fold enrichment (direction 1). Values \> 1 indicate conservation, \< 1
  indicate divergence.

- species2.neigh:

  Number of neighbors of gene2 gene in net2

- species2.ortho.neigh:

  Number of ortholog-mapped neighbors from net1, excluding the gene2
  gene itself

- species2.neigh.overlap:

  Intersection size

- species2.p_value_con:

  Upper-tail hypergeometric p-value for conservation (direction 2); 1
  when the overlap is 0 or 1

- species2.p_value_div:

  Lower-tail hypergeometric p-value for divergence (direction 2)

- species2.p_value_gt:

  \\P(X \> x)\\, ungated (direction 2)

- species2.p_value_eq:

  \\P(X = x)\\ (direction 2)

- species2.effect_size:

  Fold enrichment (direction 2). Values \> 1 indicate conservation, \< 1
  indicate divergence.

- species1.jaccard:

  Jaccard index of neighborhood overlap (direction 1): intersection /
  union of anchor and ortholog-mapped neighbors. Range \\\[0, 1\]\\.

- species2.jaccard:

  Jaccard index of neighborhood overlap (direction 2).

- species1.urn:

  Hypergeometric population for direction 1: the number of genes in net1
  less the anchor gene itself. Constant within a call. Downstream
  `power` uses it directly rather than reconstructing it from the effect
  sizes.

- species2.urn:

  Hypergeometric population for direction 2 (genes in net2 less the
  anchor). Differs from `species1.urn` whenever the two networks differ
  in size.

## Details

For an ortholog pair (A in species 1, B in species 2) and direction 1
(anchor A), the neighbourhood of A in net1 (size \\m\\) is tested
against the set of species-1 orthologs of B's neighbours in net2 (size
\\k\\); \\x\\ is their intersection. The urn excludes the anchor gene: A
is never its own neighbour, it is dropped from the mapped set, and the
hypergeometric population is the other \\N - 1\\ network genes.
`species1.p_value_con` is \\P(X \ge x)\\ and keeps the gate of the
original RComPlEx R Markdown (not the paper): it is reported as 1 when
\\x \le 1\\. `species1.p_value_gt` (\\P(X \> x)\\) and
`species1.p_value_eq` (\\P(X = x)\\) are ungated so that
`p_value_gt + U * p_value_eq` is a randomized p-value, exactly uniform
under the null (used by
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md)
to estimate pi0). The effect size \\(x / k) / (m / (N - 1))\\ is
computed for every \\x\\ (1 when \\m = 0\\ or \\k = 0\\). Direction 2 is
symmetric.

## Ortholog preparation for de novo transcriptomes

When working with de novo transcriptome assemblies, it is **critical**
to filter ortholog tables to expressed genes **before** applying any
maximum-paralogs filter. De novo assemblies can produce 10–50x more gene
models than are actually expressed (e.g. 205K models but only 14K
expressed). Unexpressed fragments inflate paralog counts per ortholog
group, causing the paralogs filter to discard groups that would pass if
only expressed members were counted. Failing to do this can reduce the
usable gene set by 40–70% and eliminate all statistical power from the
hypergeometric test.

Recommended workflow:

1.  Build expression matrices for both species

2.  Filter ortholog table to genes present in expression data

3.  Apply maximum-paralogs filter on the reduced ortholog table

4.  Compute networks and run comparison

## Gene universe

Networks are built on all supplied genes and \\N\\ is the number of
genes in the anchor network, as in Netotea et al. (2014). The canonical
ComPlEx implementations (`RComPlEx.Rmd`, `ComPlEx_python`) instead
restrict the expression tables to genes with an ortholog before building
the networks, which changes neighbourhoods, thresholds and therefore the
calls. rcomplex reproduces those implementations only when the
expression matrices are restricted the same way before
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
(`tests/testthat/test-equivalence.R`).

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- compare_neighborhoods(net_A, net_B, orthologs)
head(comparison[, c(
  "gene1", "gene2", "hog",
  "species1.effect_size"
)])
} # }
```
