# Resolve paralog copies in an ortholog map

Reduces the candidate gene pairs of multi-copy hogs to one counterpart
per gene. Uses evidence the pipeline has already produced. Use it before
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md),
which projects module labels between species.

## Usage

``` r
resolve_ortholog_map(
  orthologs,
  genes1,
  genes2,
  species1 = NULL,
  species2 = NULL,
  edges = NULL
)
```

## Arguments

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog` (output of
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)).
  `gene1` / `gene2` hold gene identifiers; species membership is
  resolved against `genes1` / `genes2`.

- genes1, genes2:

  Character vectors giving the gene universes of the two species, e.g.
  `rownames(net$network)`.

- species1, species2:

  Species labels. Required only when `edges` is supplied, to filter
  edges.

- edges:

  Optional coexpressolog edge table from
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
  with columns `gene1`, `gene2`, `species1`, `species2`, `hog` and
  `effect_size`.

## Value

A data frame with columns:

- gene1:

  Species-1 gene identifier

- gene2:

  Species-2 gene identifier

- hog:

  Ortholog group identifier

- source:

  `"coexpressolog"` or `"unresolved"`

## Details

A raw ortholog table implies a paralog expansion, and the projection is
sensitive to it.

## Resolution waterfall

- coexpressologs:

  The mutual-best significant pair: `gene1`'s highest-ranked partner by
  `effect_size` must also rank `gene1` highest.

- unresolved:

  Everything else keeps all its candidate pairs, left for the consumer
  to resolve (see the `source` column).

## Preserved gene set

Resolution never removes a species-2 gene from the map. A candidate pair
survives when it was resolved, or when its species-2 gene was not
claimed by any resolved pair, so the set of mappable species-2 genes is
identical to the set implied by `orthologs` alone. Only which species-1
copy points at each gene changes.

## Examples

``` r
if (FALSE) { # \dontrun{
map <- resolve_ortholog_map(
  ortho, rownames(net_a$network), rownames(net_b$network),
  species1 = "SP_A", species2 = "SP_B", edges = rcx$edges
)
table(map$source)
} # }
```
