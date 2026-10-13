# Find co-expression cliques using C++ two-level decomposition

Finds the maximal species cliques of each hog. Then assigns the best
gene to each species of the clique.

## Usage

``` r
find_cliques(edges, ...)

# Default S3 method
find_cliques(
  edges,
  target_species,
  min_species = length(target_species),
  cost_weights = c(q = 1, effect = 0),
  ...
)
```

## Arguments

- edges:

  Data frame with columns:

  gene1

  :   Gene identifier (first gene in pair)

  gene2

  :   Gene identifier (second gene in pair)

  species1

  :   Species for gene1

  species2

  :   Species for gene2

  hog

  :   Ortholog group identifier

  q_value

  :   q-value for the edge (from pair-level testing)

  effect_size

  :   Numeric effect size

  Optionally includes a `type` column for filtering. Pass the unfiltered
  table (every tested pair, e.g. the output of
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)):
  cliques are built from `"conserved"` rows only, but `intensity` fits
  each edge's weight scale on every row of its species pair, so a
  pre-filtered table changes what it measures. A table whose `type`
  column holds only `"conserved"` rows triggers a warning (class
  `rcomplex_prefiltered_edges`, shown once per session).

- ...:

  Additional arguments passed to the default method.

- target_species:

  Character vector of species abbreviations.

- min_species:

  Minimum number of species per clique (default:
  `length(target_species)`).

- cost_weights:

  Named numeric vector with elements `"q"` and `"effect"` controlling
  the composite cost used for gene-assignment ranking. The cost is
  `q * mean_q - effect * mean_effect` (lower is better). Default
  `c(q = 1, effect = 0)` reproduces the original mean-q-only ranking.

## Value

A data frame with one row per clique:

- hog:

  Ortholog group identifier

- :

  One column per target species (gene ID or NA)

- n_species:

  Number of species in the clique

- mean_q:

  Mean q-value across present clique edges

- max_q:

  Maximum q-value across present clique edges

- mean_effect_size:

  Mean effect size across present edges

- n_edges:

  Number of present edges

- n_missing:

  Number of missing edges (always 0)

- exhaustive:

  `TRUE` when the best-gene search finished within its step budget

- intensity:

  Onnela intensity: geometric mean, across present edges, of each edge's
  ensemble connection probability (in (0, 1); higher = stronger
  conservation). Each edge's association strength (`effect_size`,
  observed overlap over its expectation) is mapped to a probability
  `p = z w / (1 + z w)`, with `z` fitted per species pair by maximum
  entropy, so intensity reads as the per-edge probability that the whole
  clique exists. Not `1 - q_value`: every clique edge already passed
  alpha, so `1 - q` left no range (#11). Not the Jaccard index, whose
  null expectation grows with neighbourhood size, which ranked hub genes
  above equally conserved low-degree ones (#15). `NA` when `edges` has
  no usable `effect_size` column, or when any present clique edge lacks
  a finite value.

- min_effect_size:

  Minimum effect size across present edges (bottleneck enrichment)

- score:

  Sum of the present edges' `score`, in bits. `NA` when `edges` has no
  `score` column.

## Details

The search uses Bron-Kerbosch with pivoting on the species-level
adjacency graph. A backtracker picks the gene per species that minimises
the mean q-value across all present edges.

The clique `score` adds the edge scores. Log-odds add, so the sum is the
log-likelihood ratio of the conserved subnetwork (NetworkBLAST), and it
grows with clique size on purpose.

The best-gene search stops after 2 million steps per species clique. A
cut search keeps its best assignment so far and sets
`exhaustive = FALSE`. Ties at the q floor go to the larger effect size,
then to the gene ids first in radix order. Row order does not matter.

A group of species wired pairwise need not share one set of copies. Then
it reports the largest species subsets that share one set of copies. A
hog gets at most 4,096 subset searches; past that its rows have
`exhaustive = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
cliques <- find_cliques(edges, target_species = c("SP_A", "SP_B", "SP_C"))

# Allow partial species cliques (2 of 3 species)
partial <- find_cliques(edges, target_species, min_species = 2L)
} # }
```
