# Line up per-species modules across species

Lines up the modules of each species into shared labels. The same label
in two species is the same module. Modules come from
[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md),
one species at a time.

## Usage

``` r
joint_modules(nets, orthologs, modules = NULL, lambda = 1, seed = NULL)
```

## Arguments

- nets:

  Named list of at least two
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  results, dense or sparse. The names are the species.

- orthologs:

  Long table with columns `species`, `gene`, `hog`, or a file for
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).

- modules:

  Named list with one
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  or
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
  result per species of `nets`, in any order. Each must label every gene
  of its network. `NULL` runs
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  on each species.

- lambda:

  Weight of the hog term against the co-expression term. One number of 0
  or more. At 0 no label spans two species.

- seed:

  Integer seed, or `NULL` to use the current stream.

## Value

A list with:

- modules:

  Named list with one module assignment per species, in the order of
  `nets`. Each has the
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
  shape: `modules` (named integer vector, gene -\> shared label),
  `module_genes`, `n_modules`, `method = "joint"` and `params`.

- members:

  Data frame with one row per species module: `species`, `module` (the
  label in `modules`), `label` (the shared label) and `n_genes`. Rows
  are sorted by `label`, species in `nets` order, then `module`.

- hogs:

  Data frame with one row per hog with two or more copies among the
  network genes, sorted by `hog`: `hog`, `n_copies`, `n_modules`
  (distinct shared labels among the copies), `module_main` (the label
  with most copies; a tie goes to the smallest label) and `copies_main`
  (copies in it).

- params:

  List of `lambda`, `gamma` (1), `seed`, `modules_given` and `quality`
  (the objective of the final labels).

## Details

Each module of each species is one node. Louvain joins these nodes into
shared labels. It never moves a single gene, so labels are species
modules lined up, never genes clustered again. A shared label may hold
several modules of one species. It may hold modules of one species only,
when nothing lines up.

The objective has two terms. Within a species, it is modularity with
that species' own null. Across species, it is the hog overlap of two
modules against the overlap expected for that species pair. A hog with
`c_s` copies in species s and `c_t` in species t gives each copy pair
the weight `1 / (c_s c_t)`. So hogs with many copies count less. The
edges are those at or above each species' `threshold`. Each species'
edge weights are divided by its largest weight.

On wood and leaf, every `lambda` from 0.25 to 16 gave the same labels.
For why family-steered clustering of genes was withdrawn, see
`dev/design-notes/hypergraph-community-literature.md`, sections 10-12.

## See also

[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
for one species,
[`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
for pairs of species with p-values,
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
to test a module in another species.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- joint_modules(list(A = net_a, B = net_b), "orthologs.tsv",
  seed = 1
)
res$members
res$hogs[res$hogs$n_modules > 1, ]
} # }
```
