# Classify hub gene conservation across species and clades

Maps the hub genes from
[`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md)
to hogs. Classifies each hog by its hub pattern across clades (for
example annual and perennial).

## Usage

``` r
classify_hub_conservation(hub_results, ...)

# Default S3 method
classify_hub_conservation(hub_results, clades, module_comparisons = NULL, ...)
```

## Arguments

- hub_results:

  Named list keyed by species name. Each element is the data frame
  output of
  [`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md)
  (with `orthologs` provided so the `hog` column is populated).

- ...:

  Additional arguments passed to the default method.

- clades:

  Named list of species vectors, one per clade, e.g.
  `list(annual = c("SP_A", "SP_B"), perennial = c("SP_C", "SP_D"))`.
  Clades may nest but must not cross. A species in no clade forms its
  own clade. Hubs are counted per top-level clade.

- module_comparisons:

  Optional named list of
  [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
  outputs keyed by alphabetically sorted species pair (e.g.
  `"SP_A.SP_C"`). Required for the conserved_hub vs rewired_hub
  distinction. Two things the caller must now satisfy themselves: each
  element must be built with the alphabetically first species as
  `modules_ref`, and
  [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
  needs a map from
  [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md),
  so the networks are required for the gene universes. Pass
  `species_ref` / `species_test` to
  [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
  and that orientation is checked here instead of taken on trust. A
  module pair absent from the table counts as not corresponding. Modules
  correspond when `q_value < 0.1` and `jaccard >= 0.1`, and a HOG is
  `conserved_hub` when at least half of its cross-clade hub pairs
  correspond.

## Value

A data frame with one row per HOG:

- hog:

  HOG identifier

- classification:

  Conservation category (see Classification waterfall)

- n_species_hub:

  Number of species where the HOG is a hub

- n_species_present:

  Number of species where the HOG has genes. Reported for context; it is
  not the clade denominator.

- hub_trait_groups:

  Comma-separated clades where it qualifies as hub (`NA` for non_hub)

- n_corresponding:

  Cross-clade hub pairs with corresponding modules (`NA` without
  `module_comparisons`)

- n_cross_pairs:

  Total cross-clade hub pairs checked (`NA` without
  `module_comparisons`)

- max_centrality:

  Highest centrality score across species

- best_hub_species:

  Species with highest centrality

## Classification waterfall

For each HOG that appears in at least one species:

- conserved_hub:

  Hub in multiple clades **and** the hub modules correspond across
  clades (checked via `module_comparisons`).

- rewired_hub:

  Hub in multiple clades but in **non-corresponding** modules – the gene
  kept its centrality but changed regulatory context.

- multi_trait_hub:

  Hub in multiple clades; module correspondence unknown
  (`module_comparisons` not provided).

- *trait*\_specific_hub:

  Hub in exactly one clade (e.g. `"annual_specific_hub"`).

- sporadic_hub:

  Hub in some species but does not reach half of a clade.

- non_hub:

  Present in modules but not a hub in any species.

## Examples

``` r
if (FALSE) { # \dontrun{
hub_list <- list(
  SP_A = identify_module_hubs(mods_A, net_A, ortho_A),
  SP_B = identify_module_hubs(mods_B, net_B, ortho_B)
)
clades <- list(annual = "SP_A", perennial = "SP_B")
classify_hub_conservation(hub_list, clades)

# With module correspondence, for the conserved_hub / rewired_hub split.
# The list key must be the alphabetically sorted species pair.
map <- resolve_ortholog_map(
  ortho_AB, rownames(net_A$network), rownames(net_B$network)
)
corr <- list(SP_A.SP_B = module_correspondence(
  mods_A, mods_B, map,
  species_ref = "SP_A", species_test = "SP_B"
))
classify_hub_conservation(hub_list, clades, module_comparisons = corr)
} # }
```
