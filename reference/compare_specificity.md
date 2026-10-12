# Neighbourhood specificity of ortholog pairs

Ranks each ortholog pair against every other gene of the partner
species. This is the engine of `method = "rank"` in
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
and
[`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md).
[`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md)
turns its p-values into q-values.

## Usage

``` r
compare_specificity(
  net1,
  net2,
  orthologs,
  n_cores = 1L,
  directions = c("both", "1to2", "2to1")
)
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

- directions:

  `"both"` (default), `"1to2"` (anchors in species 1) or `"2to1"`.

## Value

A data frame with `gene1`, `gene2`, `hog` and, per requested direction
(`species1.` for 1 to 2, `species2.` for 2 to 1):

- neigh:

  Anchor degree.

- mapped:

  Size of the translated set \\T\\.

- auroc:

  AUROC of \\T\\ in the paired ortholog's ranking.

- p_value:

  Rank p-value of the paired ortholog, as above.

- jaccard:

  As in
  [`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md).

- n.cand:

  Number of candidate genes in the partner network.

- effect_size:

  Equal to `auroc`.

- auroc.grid:

  Matrix, one row per pair and one column per raw-p fraction f (1e-5 to
  1, the column names): the `ceiling(f * n.cand)`-th largest candidate
  AUROC for the anchor, i.e. the AUROC the ortholog needs to reach raw p
  of about f.
  [`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md)
  reads it for the edge `power`. It is computed on every call; a
  comparison against a
  [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  partner, which only needs `p_value`, pays for it too (the package's
  own null runs skip it).

`auroc` and `p_value` are `NA` when \\T\\ is empty or spans every other
partner gene.

## Details

The anchor's co-expression partners are translated to the partner
species. The paired ortholog is scored by how well its own co-expression
ranking recognises that translated list, compared with how well every
other partner-species gene recognises it. The score follows the
co-expression conservation score of Suresh et al. (2023). For anchor
gene \\i\\ in species 1 and direction 1 to 2:

1.  \\T\\ is the set of species-2 genes orthologous, through any copy,
    to a network neighbour of \\i\\ (an entry at or above the analysis
    threshold), minus every species-2 gene in \\i\\'s own ortholog
    group. Removing the own group keeps a pair from scoring itself
    through its paralogs.

2.  For every species-2 gene \\j\\, the AUROC of \\T\\ (without \\j\\)
    among the other \\n_2 - 1\\ genes, ranked by their co-expression
    with \\j\\. Stored entries of column \\j\\ are ranked exactly;
    unstored entries tie at the bottom with their mid-rank. A sparse
    network therefore ranks only its stored top entries, so
    `store_density` sets the resolution.

3.  The p-value of the ortholog pair \\(i, o)\\ is \\(1 + g) / n_2\\,
    where \\g\\ counts the genes \\j \ne o\\ whose AUROC is at least
    that of \\o\\: the rank of the ortholog among all species-2 genes. A
    translated list that thousands of genes recognise as well as \\o\\
    gives a large p-value however high its AUROC.

Direction 2 to 1 swaps the roles. All p-values of one direction lie on
the same \\1 / n_2\\ grid. They are not calibrated on their own: pass
them with the p-values of a comparison against
[`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
to
[`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md).

Runs on the sparse representation; dense networks are converted with
every nonzero entry stored.

## References

Suresh, H., Crow, M., Jorstad, N., Hodge, R., Lein, E., Dobin, A.,
Bakken, T. & Gillis, J. (2023). Comparative single-cell transcriptomic
analysis of primate brains highlights human-specific regulatory
evolution. *Nature Ecology & Evolution*, 7(11), 1930–1943.
[doi:10.1038/s41559-023-02186-7](https://doi.org/10.1038/s41559-023-02186-7)
