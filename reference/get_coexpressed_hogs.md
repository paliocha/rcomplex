# Query co-expression partners of a candidate HOG across species

For a given HOG, finds which other HOGs co-express with it in each
species network, then aggregates across species. Useful after
[`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md)
and
[`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
to explore the co-expression neighborhood of a hub gene.

## Usage

``` r
get_coexpressed_hogs(
  candidate_hog,
  networks,
  orthologs,
  clades = NULL,
  min_species = 2L,
  edges = NULL
)
```

## Arguments

- candidate_hog:

  Character string: the HOG ID to query (e.g., `"HOG42"`).

- networks:

  Named list of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  outputs, keyed by species abbreviation.

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog` (from
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)).

- clades:

  Optional named list of species vectors, one per clade. Clades may nest
  but must not cross. A species in no clade forms its own clade. It
  enables the `coexpressed_traits` output column.

- min_species:

  Minimum number of species where co-expression must occur for a partner
  to be reported (default 2).

- edges:

  Optional stacked edge data frame from
  [`run_pairwise_comparisons()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md).
  When provided, attaches conservation statistics for each partner HOG.

## Value

A data frame with columns:

- partner_hog:

  HOG ID of the co-expression partner

- n_species:

  Number of species where co-expression occurs

- coexpressed_species:

  Comma-separated species names

- coexpressed_traits:

  Comma-separated top-level clades of those species (sorted). A species
  in no clade gives its own name. `NA` if `clades` is `NULL`.

- mean_weight:

  Mean co-expression weight across species. Per species, the maximum
  over candidate-HOG copies of the mean edge weight from a copy to the
  partner-HOG genes in that copy's own neighborhood (edge values at or
  above the network threshold), so the value is identical for dense and
  sparse networks

When `edges` is provided, three additional columns are appended:

- partner_conserved:

  Logical: does the partner HOG have any `"conserved"` edge?

- partner_mean_effect:

  Mean effect size from edges for the partner HOG

- partner_min_q:

  Minimum q-value from edges for the partner HOG

Rows are sorted by `n_species` (descending), then `mean_weight`
(descending).

## Examples

``` r
if (FALSE) { # \dontrun{
clades <- list(
  annual = c("SP_A", "SP_B"), perennial = c("SP_C", "SP_D")
)
partners <- get_coexpressed_hogs("HOG42", networks, orthologs,
  clades = clades,
  min_species = 2L, edges = edges
)
# Annual-only partners
partners[partners$coexpressed_traits == "annual", ]
} # }
```
