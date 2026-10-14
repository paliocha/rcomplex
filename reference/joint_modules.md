# Find modules jointly across species

Clusters the networks of all species and their hogs as one graph. Labels
are shared: the same label in two species is the same module. The hog
table shows which hogs keep their copies in one module.

## Usage

``` r
joint_modules(
  nets,
  orthologs,
  weight = NULL,
  objective = c("modularity", "CPM"),
  resolution = NULL,
  seed = NULL,
  engine = "star"
)
```

## Arguments

- nets:

  Named list of at least two
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  results, dense or sparse. The names are the species.

- orthologs:

  Long table with columns `species`, `gene`, `hog`, or a file for
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).

- weight:

  Weight of the edge from a copy to its hog. `NULL` uses the median of
  all scaled co-expression edge weights.

- objective:

  Leiden objective: `"modularity"` or `"CPM"`.

- resolution:

  Leiden resolution. `NULL` uses 1 for modularity. For CPM it uses the
  edge density of the joint graph.

- seed:

  Integer seed, or `NULL` to use the current stream.

- engine:

  Clustering engine. Only `"star"` exists.

## Value

A list with:

- modules:

  Named list with one module assignment per species, in the order of
  `nets`. Each has the
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
  shape: `modules` (named integer vector, gene -\> label),
  `module_genes`, `n_modules`, `method = "joint"` and `params`.

- hogs:

  Data frame with one row per hog node, sorted by `hog`: `hog`,
  `n_copies` (genes linked), `n_modules` (distinct labels among the
  copies), `module_main` (the label with most copies; a tie goes to the
  smallest label) and `copies_main` (copies in it).

- params:

  List of `weight` (the value used), `objective`, `resolution`, `engine`
  and `seed`.

- modularity:

  Modularity of the joint graph, hog nodes included.

## Details

The graph is a star expansion. Each species adds its genes and the edges
at or above its `threshold`. Edge weights are divided by that species'
largest weight, so every species lies on the same scale. No edge crosses
species. Each hog with two or more copies among the genes adds one hog
node. An edge of weight `weight` links the hog node to each copy. Hog
nodes carry no other edge. Leiden clusters this graph; the modules then
leave out the hog nodes. `weight` is the coupling between the
co-expression edges and the hogs. A small `weight` gives per-species
modules; a large one pulls the copies of a hog together. The default,
the median co-expression weight, is a guess until a sweep over `weight`
tests it. For CPM, `.leiden_weights()` scales every edge by the graph
maximum, hog edges included; `params$weight` keeps the unscaled value.

## See also

[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
for one species,
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
to test a module in another species.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- joint_modules(list(A = net_a, B = net_b), "orthologs.tsv",
  seed = 1
)
res$hogs[res$hogs$n_modules > 1, ]
pres <- module_preservation(res$modules$A, net_a, net_b, ortho)
} # }
```
