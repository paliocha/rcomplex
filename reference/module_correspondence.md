# Match modules across species by ortholog overlap

Cross-tabulates the modules of two species over a paralog-resolved
ortholog map. Tests each module pair for excess overlap.

## Usage

``` r
module_correspondence(
  modules_ref,
  modules_test,
  map,
  species_ref = NULL,
  species_test = NULL,
  seed = NULL
)
```

## Arguments

- modules_ref, modules_test:

  Module assignments for the two species, from
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  or
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md).

- map:

  Ortholog map from
  [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md),
  with `gene1` in the reference species and `gene2` in the test species.

- species_ref, species_test:

  Optional species labels recorded on the result. `module1` belongs to
  `modules_ref`, and that orientation cannot be recovered from the
  table, so supplying these lets
  [`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
  catch a transposed call.

- seed:

  Integer seed for the randomized-p draws behind pi0, or `NULL`
  (default) to draw from the ambient stream and leave it advanced. A
  seed draws from a private stream and restores the caller's on exit,
  the package-wide contract described under
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

## Value

A list with `pairs` – a data frame of `module1`, `module2`, `size1`,
`size2`, `overlap`, `jaccard`, `p_value` and `q_value`, one row per
module pair – plus the `species_ref` / `species_test` labels. The column
names match what
[`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
expects.

## Details

This answers "which module corresponds to which". That differs from
whether a module's topology is preserved
([`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)).
[`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
needs the correspondence, not the preservation call.

The map assigns one reference gene to each test-species gene, so the
hypergeometric's independence assumption holds. The multi-copy expansion
that made the same test anti-conservative under the retired gene-overlap
engine is gone here: a HOG with three paralogs no longer contributes
three correlated draws to the same urn.

## Examples

``` r
if (FALSE) { # \dontrun{
map <- resolve_ortholog_map(
  ortho, rownames(net_a$network),
  rownames(net_b$network)
)
corr <- module_correspondence(mods_a, mods_b, map)
subset(corr$pairs, q_value < 0.1)
} # }
```
