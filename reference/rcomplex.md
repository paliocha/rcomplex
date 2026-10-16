# Compare co-expression networks across species in one call

`rcomplex()` runs the whole analysis. It builds one network per species,
tests every ortholog pair in every species pair, and classifies the gene
cliques. Give it expression or networks, plus orthologs.

## Usage

``` r
rcomplex(
  expr = NULL,
  orthologs,
  networks = NULL,
  block = NULL,
  clades = NULL,
  density = 0.03,
  sign = c("positive", "negative"),
  method = c("hypergeometric", "rank"),
  alpha = 0.1,
  modules = FALSE,
  null = FALSE,
  n_cores = 1L,
  block_size = NULL,
  seed = NULL
)
```

## Arguments

- expr:

  Named list of expression data, one entry per species: a genes x
  samples matrix or a `SummarizedExperiment`. Give `expr` or `networks`,
  not both.

- orthologs:

  The long ortholog table (`species`, `gene`, `hog`), or a path to an
  ortholog file for
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).

- networks:

  Named list of network objects from
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  or
  [`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md),
  one entry per species.

- block:

  Named list of per-sample factors, one per species in `expr`. Each
  species is split with
  [`split_layers()`](https://paliocha.github.io/rcomplex/reference/split_layers.md).
  The network is built from the wiring layer.

- clades:

  Named list of species vectors, one per clade, as for
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md).
  `NULL` skips the clade tiers.

- density:

  Fraction of gene pairs kept as edges in each network.

- sign:

  `"positive"` or `"negative"`: the correlation sign that makes an edge.
  See
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).

- method:

  `"hypergeometric"` or `"rank"`: the test of each ortholog pair. See
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md).

- alpha:

  The q-value threshold for a call.

- modules:

  If `TRUE`, detect modules in each species and test their preservation
  in every species pair.

- null:

  If `TRUE`, repeat the run with hog labels shuffled among genes of the
  same expression level. The null calls and tiers estimate what
  abundance and network structure alone produce.

- n_cores:

  Number of threads.

- block_size:

  `NULL` (default) builds each species network from the dense gene x
  gene matrix. A positive whole number builds it in blocks of that many
  genes. The full matrix then never forms.

- seed:

  Integer seed, or `NULL` to draw from the global stream. A seed makes
  the run reproducible and leaves your stream unchanged.

## Value

A list of class `rcomplex`:

- networks:

  The network objects, by species.

- edges:

  The edge table from
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
  plus a `sign` column: the sign of every network.

- cliques, members, twins:

  The three tables of
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md):
  one row per clique, one row per clique member, and one row per
  collapsed twin gene. The tables of the three graphs are stacked. `run`
  numbers the graph: 1 strict, 2 loose, 3 all tested pairs.
  `alpha_graph` holds its cut: `alpha`, 0.9 or `Inf`. `clique_id` counts
  within a run, so `run` and `clique_id` together name a clique.

- capped:

  One row per hog and run over the clique budget of
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md),
  with `run` and `alpha_graph`. Its cliques are species blocks. Usually
  empty.

- classification:

  One row per clique, from
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md).
  Its `run` and `clique_id` join `cliques` and `members`.

- modules, preservation:

  With `modules = TRUE`: the
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  results by species, and the
  [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
  result.

- edges_null, classification_null:

  With `null = TRUE`: the edge table and the classification of the
  shuffled ortholog table.

- call:

  The call.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the classification.

## Details

The driver lists cliques on three gene graphs. The strict graph holds
the edges with `q_value < alpha`. It feeds `complete_conserved`,
`lineage_specific`, `partial_present` and `trait_specific`. The loose
graph holds the edges with `q_value < 0.9`. It feeds
`partial_significant`. The third graph holds every tested pair. It
covers only the hogs with no complete or partial tier on the first two
graphs. It feeds `differentiated`, as does the loose graph where every
cross pair has `q_value < 0.9`.
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
classifies the three graphs together, so `hog_class` reads all of them.
A clique maximal on two graphs appears once for each. The third graph is
built only with `clades`, since only the clade tiers read it.

Every clique has at least `min(3, number of species)` genes. Rank edges
are classified at `min_power = 0.9`, hypergeometric edges at 0.8.

`density` and `sign` apply only to networks built from `expr`.
`method = "rank"` and `block` need `expr`. The rank test calibrates
against
[`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
partners, built within each block.

`block_size` passes to
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
and needs `expr`. The blockwise network equals the dense one except at
floating-point near-ties at the density cut. 512 is a good value. A
species with more than 40,000 genes gets a message that suggests it.

[`summary()`](https://rdrr.io/r/base/summary.html) counts each tier in
`cliques`, in `hogs` and in `patterns`. A pattern is one HOG with one
clique species set. The paralog copies of a large HOG can form millions
of cliques, but they share a few patterns. Paralogs with the same
partners are twins: one clique stands for all their combinations. The
`multiplicity` column counts those combinations, and `twins` lists the
collapsed genes.

With `null = TRUE`, the driver shuffles the `gene` column of `orthologs`
within each species. A gene moves only among genes of the same
expression decile, by `gene_mean` from
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).
A network without `gene_mean`, such as one from
[`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md),
forms one decile. The driver then reruns edges, cliques and tiers on the
real networks. Both runs share the rank calibration networks.

With `null = TRUE`, [`summary()`](https://rdrr.io/r/base/summary.html)
reports `calls` and `calls_null` per species pair at `q_value < alpha`.
Their ratio is the empirical false-call rate. The rate is `NA` for a
species pair without calls. The tier table gains `cliques_null`: the
null cliques per tier.

[`print()`](https://rdrr.io/r/base/print.html) counts the cliques of
each graph. It shows, per species, the genes, the samples, the density
and `r_threshold` from
[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md):
the weakest correlation that passed the density threshold. The same
density needs a stronger correlation at few samples. A network from
[`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md)
has no correlation and shows `r ?`, as does a network without edges.

## Examples

``` r
f <- function(x) system.file("extdata", x, package = "rcomplex")
read <- function(x) as.matrix(read.delim(f(x), row.names = 1L))
expr <- lapply(c(SpA = "expr_sp1_small.txt", SpB = "expr_sp2_small.txt"),
  read)
res <- rcomplex(expr, f("orthologs_small.txt"), density = 0.1, seed = 1)
res
#> rcomplex: 2 species, 10 hogs, 0 tiers
#>   SpA  50 genes  10 samples  density 0.1  r >= 0.30
#>   SpB  40 genes  10 samples  density 0.1  r >= 0.23
#>   edges 12 tested, 0 called at q < 0.1
#>   graphs: 0 cliques at q < 0.1, 0 at q < 0.9
#>   cliques 0
summary(res)
#> [1] tier         cliques      cliques_null multiplicity hogs        
#> [6] patterns    
#> <0 rows> (or 0-length row.names)
#> null: not run
write_rcomplex(res, tempfile())
```
