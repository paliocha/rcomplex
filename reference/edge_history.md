# Gains and losses of each clique on the species tree

Reads clique presence over species as a binary character. Fitch
parsimony on a rooted species tree places each gain and loss on a
branch. A second table sums the changes per branch.

## Usage

``` r
edge_history(cliques, tree, species)
```

## Arguments

- cliques:

  A
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  table. It needs `clique_id`, `hog`, `missing_species` and
  `missing_reason`. `run` is kept.

- tree:

  A rooted `phylo` from the ape package, with tips named by species.

- species:

  Every species of the analysis. Each must be a tip.

## Value

A list of two data frames. `cliques` has one row per input row, in input
order: `n_present`, `n_absent` (state 0), `n_unknown`, `state_root`,
`root_tied`, `n_gain`, `n_loss`, `n_mpr`, and the branches in `gains`
and `losses`. `branches` has one row per branch, named by the node below
it as in
[`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md):
`n_species`, and `n_gain` and `n_loss` summed over the cliques.

## Details

A member species has state 1. A species with reason `tested_ns` has
state 0: the pair was tested and failed. Every other reason is unknown
(`?`), and so is a tip outside `species`. An unknown tip never counts as
a change. An `absent` species is unknown because gene presence is its
own character, and a gene loss is not a co-expression loss.

The pass is unit-cost Fitch, run once per distinct pattern of missing
species. The tree may hold polytomies. A clique restricted to one clade
costs the same as one gain on its stem or one loss outside it. Such a
root tie resolves to absent, so the clique reads as one gain, as
`lineage_specific` does. `root_tied` flags the tie. Below the root, a
tie takes the change as early as possible (ACCTRAN). `n_mpr` counts the
equally parsimonious reconstructions.

## Examples

``` r
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(text = "((A,B),(C,D));")
  cl <- data.frame(
    clique_id = c("1", "2"), hog = "H1",
    missing_species = c("C,D", "D"),
    missing_reason = c("tested_ns,tested_ns", "tested_ns")
  )
  edge_history(cl, tree, c("A", "B", "C", "D"))
}
#> $cliques
#>   clique_id hog n_present n_absent n_unknown state_root root_tied n_gain n_loss
#> 1         1  H1         2        2         0          0      TRUE      1      0
#> 2         2  H1         3        1         0          1     FALSE      0      1
#>   n_mpr gains losses
#> 1     2 node6       
#> 2     1            D
#> 
#> $branches
#>   branch n_species n_gain n_loss
#> 1      A         1      0      0
#> 2      B         1      0      0
#> 3      C         1      0      0
#> 4      D         1      0      1
#> 5  node6         2      1      0
#> 6  node7         2      0      0
#> 
```
