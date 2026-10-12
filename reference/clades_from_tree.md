# Clades from a species tree

Turns a species tree into the clade list that the classifiers take. Each
internal node gives one clade: the tip labels below it. The root holds
every species, so it carries no contrast and is dropped.

## Usage

``` r
clades_from_tree(phy, min_size = 2L)
```

## Arguments

- phy:

  A tree of class `phylo` from the ape package.

- min_size:

  Smallest clade to keep. Use `1L` to keep the tips.

## Value

A named list of character vectors, one per clade. A clade takes its node
label. An empty or repeated label becomes `"node<k>"`, with `k` the node
number. A tip takes its own label.

## Examples

``` r
if (requireNamespace("ape", quietly = TRUE)) {
  phy <- ape::read.tree(text = "((A,B)AB,(C,D)CD);")
  clades_from_tree(phy)
}
#> $AB
#> [1] "A" "B"
#> 
#> $CD
#> [1] "C" "D"
#> 
```
