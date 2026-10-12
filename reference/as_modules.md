# Build a module assignment from any gene partition

Turns any gene sets into a module assignment. The sets can be regulons,
curated pathways or clusters from another tool.
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md),
[`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
and
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
accept the result without
[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

## Usage

``` r
as_modules(x)
```

## Arguments

- x:

  One of: a named vector mapping gene to module label (`NA` for an
  unassigned gene); a named list mapping module label to a character
  vector of genes; or an existing module assignment (a
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  or `as_modules()` result), which is returned unchanged.

## Value

A list with `modules` (named character vector, gene -\> module label,
unassigned genes omitted), `module_genes` (named list, module label -\>
genes), `n_modules`, `method = "external"` and `params`. `module_genes`
follows the list order for list input, numeric order for numeric labels,
and order of first appearance for other labels; the order changes only
listings such as `coverage`, not statistics.

## Details

The preservation tests give every gene one module label, so the sets
must not overlap. Split overlapping sets (for example regulons that
share genes) into disjoint batches and test each batch separately.

Module labels matter for reproducibility:
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
draws its permutation null module by module in the C-locale (radix)
order of the labels as strings (so "10" before "2"), so under a fixed
seed the permutation p- and q-values are reproducible for a given
labelling but change when the same partition is relabelled. The observed
statistics do not depend on the labels.

## Examples

``` r
mods <- as_modules(list(
  photosynthesis = c("g1", "g2", "g3"),
  ribosome = c("g4", "g5", "g6")
))
mods$modules
#>               g1               g2               g3               g4 
#> "photosynthesis" "photosynthesis" "photosynthesis"       "ribosome" 
#>               g5               g6 
#>       "ribosome"       "ribosome" 
```
