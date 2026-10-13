# Maximal cliques of the per-orthogroup gene graph

Builds one graph per hog from the co-expressolog calls in `edges`. Nodes
are (species, gene) pairs. Lists every maximal clique of at least
`min_size` nodes, one per paralog combination.

## Usage

``` r
gene_clique_graph(edges, ...)

# Default S3 method
gene_clique_graph(
  edges,
  min_size = 3L,
  alpha_graph = 0.1,
  max_cliques = 1e+05,
  ...
)
```

## Arguments

- edges:

  Data frame of co-expressolog calls, as returned by
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  columns `gene1`, `gene2`, `species1`, `species2`, `hog` and `q_value`.
  The `effect_size` and `score` columns are used when present.

- ...:

  Additional arguments passed to the default method.

- min_size:

  Minimum number of nodes in a reported clique (default 3, matching the
  published workflow).

- alpha_graph:

  Edges with `q_value < alpha_graph` build the graph. Use the calling
  threshold (e.g. 0.1) for complete cliques, a permissive value (0.9)
  for partially-significant cliques, and a value **above** 1 (`Inf`) for
  an unfiltered graph. The comparison is strict, so `alpha_graph = 1`
  drops every edge whose `q_value` is exactly 1 – not a corner case,
  since BH q-values cap at 1 and 20 of 511 module q-values sit there on
  the package's own worked example. Passing several thresholds and
  combining the results is the intended way to feed
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md),
  since a clique that is maximal at one threshold need not be maximal at
  another. Pass the runs to it as a list.

- max_cliques:

  Largest number of maximal cliques one hog may produce. A hog above it
  keeps its wiring classes with the strongest edges until it fits.
  Default 1e5; `Inf` keeps every clique.

## Value

A list of class `gene_cliques` with four tables. `cliques` has one row
per clique:

- clique_id:

  Integer clique id, 1 to the number of cliques

- hog:

  Ortholog group

- n_members:

  Nodes in the clique

- n_species:

  Distinct species in the clique (equal to `n_members` for cross-species
  edge tables)

- n_edges:

  Graph edges inside the clique

- mean_q, max_q:

  q-value summaries over those edges. Twin classes are joined by all
  their gene pairs. A clique edge carries the mean q of those pairs, so
  `mean_q` is the mean over the clique's expansions. `max_q` is the
  largest q of any expansion

- mean_effect_size:

  Present only when `edges` carries `effect_size`. Prefer it over
  `mean_q` for ranking: q-values saturate at the permutation floor and
  cannot separate cliques once they get there

- score:

  Present only when `edges` carries `score`. The sum of the clique
  edges' `score`, in bits. With twins, the mean over the clique's
  expansions, as for `mean_effect_size`

- multiplicity:

  Cliques this row stands for: the product of the members' `n_copies`. 1
  without twins

- mean_q_floor, n_cliques_at_q_floor:

  Smallest clique `mean_q` in the run and how many cliques are tied
  there. A large tie count means `mean_q` cannot rank those cliques at
  all

`members` has one row per clique member: `clique_id`, `species` (a
factor; its levels are the species of `edges` in order of first
appearance), `gene` (the representative of its twin class) and
`n_copies` (the genes with the same called partners, itself included).
`twins` has one row per collapsed gene: `hog`, `species`, `gene` and
`representative`. It has no rows when no two genes are twins. `dropped`
has one row per gene that the clique budget dropped: `hog`, `species`,
`gene`, `n_cliques_full` (the hog's maximal cliques before the drop),
`n_classes` (its twin classes) and `n_classes_kept`. It has no rows when
no hog is over the budget. The two floor diagnostics, plus
`alpha_graph`, `min_size` and `q_floor` (smallest edge q-value in the
graph), are also attributes of `cliques`. Both are always set, empty
results included.

## Details

[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
cliques the species graph instead and returns a single best gene
assignment per species clique. This function is the computation used by
Netotea-style clique workflows.

Because co-expressolog edges are always cross-species, no two genes of
the same species can be adjacent, so a clique carries at most one gene
per species. The `n_species` column records that explicitly rather than
assuming it, so a caller that supplies within-species edges can still
detect the violation instead of getting silently wrong species-pair
arithmetic downstream. Only a row whose two endpoints are the *same*
(species, gene) node is dropped as a self-loop; a within-species edge
between two distinct paralogs is kept, and shows up as
`n_species < n_members`.

The clique `score` adds the edge scores. Log-odds add, so the sum is the
log-likelihood ratio of the conserved subnetwork (NetworkBLAST), and it
grows with clique size on purpose.

Two genes of one species are twins when they have the same partners in
the graph. Twins are never adjacent, and a clique holds at most one of
them. Swapping one twin for another gives another maximal clique, so the
twins of a class collapse to one gene before the cliques are listed.
Each clique stands for every choice of one gene per class. The collapse
is exact: expanding the cliques gives the full set. The representative
is the first gene of its class in radix (C-locale) order.

Paralogs multiply cliques, so a hog may produce millions of them. A cap
on genes per species does not bound that count. One hog with 76 twin
classes, fewer than 10 per species, gave 2.4 million. The budget counts
cliques per hog, after the twin collapse. A hog over `max_cliques` ranks
its twin classes by the sum of their edge scores. It keeps the top class
of each species. It then keeps the most classes, in rank order, that fit
the budget. The genes of the other classes go to `dropped`. Ties go to
the gene id, so row order does not matter.

A hog of 40 or more genes is enumerated one start gene at a time. The
pieces are disjoint, so the cliques are the same. Only the peak memory
drops.

## Duplicate rows

Two rows describing the same undirected pair within one HOG collapse to
the more significant of the two, broken on `effect_size` descending when
the q-values tie – q-values saturate at the permutation floor, so a tie
there must not be settled by input row order. This is a duplicate-row
collapse, not a direction combine:
[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
already merges the two comparison directions with `pval_combine`
(default `"max"`, the reciprocal criterion) and emits one row per pair.
A caller who row-binds the two directional tables instead gets `"min"`
semantics here, so combine upstream if the reciprocal criterion is
wanted.

## References

Rodriguez E, Birkeland S, Chapple ED, et al. (2026). Comparative
regulomics of wood formation across dicot and conifer trees. *Nature
Communications* 17(1).
[doi:10.1038/s41467-026-75624-2](https://doi.org/10.1038/s41467-026-75624-2)

## See also

[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md),
[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)

## Examples

``` r
edges <- data.frame(
  gene1 = c("a1", "a1", "b1"), gene2 = c("b1", "c1", "c1"),
  species1 = c("SP_A", "SP_A", "SP_B"),
  species2 = c("SP_B", "SP_C", "SP_C"),
  hog = "HOG1", q_value = c(0.01, 0.02, 0.03)
)
g <- gene_clique_graph(edges)
g$cliques
#>   clique_id  hog n_members n_species n_edges mean_q max_q multiplicity
#> 1         1 HOG1         3         3       3   0.02  0.03            1
#>   mean_q_floor n_cliques_at_q_floor
#> 1         0.02                    1
g$members
#>   clique_id species gene n_copies
#> 1         1    SP_A   a1        1
#> 2         1    SP_B   b1        1
#> 3         1    SP_C   c1        1
```
