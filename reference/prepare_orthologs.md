# Pair orthologs for every species pair

Pairs the genes of each HOG for every pair of species in the long table.
Genes can first be renamed through
[`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.md)
output.

## Usage

``` r
prepare_orthologs(orthologs, reductions = NULL)
```

## Arguments

- orthologs:

  The long table from
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md):
  columns `species`, `gene` and `hog`.

- reductions:

  Named list of
  [`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.md)
  outputs, one per species in `orthologs`, or `NULL` to keep every gene.
  Each element needs a `gene_map` data frame with columns `original` and
  `representative`.

## Value

A data frame with columns `gene1`, `gene2` and `hog`, one row per gene
pair in a shared HOG. With `reductions`, gene names are the
representatives and duplicate rows are removed.

## Details

Species pairs follow the order in which species first appear in
`orthologs`. Rows with `NA` hog are skipped. To build the long table
from `SummarizedExperiment` objects, bind one
`data.frame(species = , gene = rownames(se), hog = rowData(se)$hog)` per
species.

Paralog reduction is a lossy step: correlated paralogs are averaged into
one representative gene, so per-paralog identity is gone by the time the
pairs reach downstream consumers. That trade-off pays for itself for
consumers that need one counterpart per gene (module preservation, the
species-graph clique backend), but the gene-graph clique backend
(`gene_clique_graph` / `classify_gene_cliques`) is built to resolve
which paralog copy is conserved, and needs every original paralog as its
own node to do that. Pass `reductions = NULL` (the default) to skip
paralog reduction entirely and keep every gene at its original identity.

## Examples

``` r
f <- system.file("extdata", "N0.tsv", package = "rcomplex")
ortho <- prepare_orthologs(read_orthologs(f))
head(ortho)
#>       gene1     gene2           hog
#> 1 GeneA_001 GeneB_001 N0.HOG0000001
#> 2 GeneA_001 GeneB_011 N0.HOG0000001
#> 3 GeneA_002 GeneB_002 N0.HOG0000002
#> 4 GeneA_012 GeneB_002 N0.HOG0000002
#> 5 GeneA_003 GeneB_003 N0.HOG0000003
#> 6 GeneA_004 GeneB_004 N0.HOG0000004
```
