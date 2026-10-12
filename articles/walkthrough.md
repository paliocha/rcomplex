# Walkthrough: the building blocks

The
[quickstart](https://paliocha.github.io/rcomplex/articles/quickstart.md)
runs
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
in one call. This article runs the same analysis step by step on eight
Pooideae grasses. Use it when you need to change one step, to inspect an
intermediate object, or to run the module and trait tests.

| Species                        | Code | Life cycle | Genus pair     |
|--------------------------------|------|------------|----------------|
| *Brachypodium distachyon*      | BDIS | annual     | Brachypodium   |
| *Brachypodium sylvaticum*      | BSYL | perennial  | Brachypodium   |
| *Hordeum vulgare*              | HVUL | annual     | Hordeum        |
| *Hordeum jubatum*              | HJUB | perennial  | Hordeum        |
| *Briza maxima*                 | BMAX | annual     | Briza          |
| *Briza media* (mapped to BMAX) | BMED | perennial  | Briza          |
| *Vulpia bromoides*             | VBRO | annual     | Festuca/Vulpia |
| *Festuca pratensis*            | FPRA | perennial  | Festuca/Vulpia |

The data are VST-normalised RNA-seq for the 1,500 most expressed hogs
from the “Live or let die?” project, stored as `SummarizedExperiment`
objects. Each species has 20 leaf and 20 root samples.

## 1. Inputs

We keep the 20 leaf samples. A network built from pooled leaf and root
samples mostly separates tissues, so use one tissue per network or
remove the tissue with `block` (see the methods article).

``` r
library(rcomplex)
library(SummarizedExperiment)

se_list <- readRDS(system.file("extdata", "pooideae_vignette.rds",
  package = "rcomplex"
))
se_list <- lapply(se_list, function(s) s[, s$tissue == "leaf"])
species <- names(se_list)

clades <- list(
  annual = c("BDIS", "HVUL", "BMAX", "VBRO"),
  perennial = c("BSYL", "HJUB", "BMED", "FPRA")
)
```

`clades` is a named list of species vectors. Every function that takes a
trait or a lineage takes this one shape.
[`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md)
builds it from a species tree.

The ortholog table has one row per gene: `species`, `gene` and `hog`.
Here it comes from
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
For an OrthoFinder or PLAZA file, call
[`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).

``` r
long <- do.call(rbind, lapply(species, function(s) {
  data.frame(
    species = s, gene = rownames(se_list[[s]]),
    hog = as.character(rowData(se_list[[s]])$hog)
  )
}))
long <- long[!is.na(long$hog), ]
tab(head(long, 4))
```

| species | gene                 | hog         |
|:--------|:---------------------|:------------|
| BDIS    | Bradi0012s00100.v3.2 | HOG:0033476 |
| BDIS    | Bradi1g00480.v3.2    | HOG:0031154 |
| BDIS    | Bradi1g00552.v3.2    | HOG:0027601 |
| BDIS    | Bradi1g00560.v3.2    | HOG:0022846 |

### Paralog reduction

Recent duplicates often have near-identical expression. They inflate the
number of cliques without adding information.
[`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.md)
merges the paralogs of a hog that correlate above `cor_threshold` and
keeps divergent copies apart.

``` r
reduced <- lapply(species, function(s) {
  reduce_orthogroups(
    assay(se_list[[s]], 1), long[long$species == s, c("gene", "hog")],
    gene_col = "gene", cor_threshold = 0.7
  )
})
names(reduced) <- species
vapply(reduced, function(r) r$n_merged, integer(1))
#> BDIS BSYL HVUL HJUB BMAX BMED VBRO FPRA 
#>   95   97  121 1235  131  108  104  126
```

[`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)
turns the long table into ortholog pairs and maps the genes through the
reduction.

``` r
orthologs <- prepare_orthologs(long, reduced)
nrow(orthologs)
#> [1] 93823
```

## 2. Networks

[`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
correlates all gene pairs and converts the correlations to mutual rank
(MR). MR asks whether two genes rank each other among their best
partners, so a gene that correlates with everything does not dominate.
The network keeps the top 3% of gene pairs as edges.

The sparse store holds more pairs than the analysis threshold:
`store_density = 0.06` leaves room for the density sweep in section 4.

``` r
networks <- lapply(reduced, function(r) {
  compute_network(r$expr_matrix, density = 0.03, store_density = 0.06)
})
```

Each network records the correlation it needed to pass the threshold. At
20 samples this value is low, and many edges are chance correlations.

``` r
vapply(networks, function(n) n$params$r_threshold, numeric(1))
#>      BDIS      BSYL      HVUL      HJUB      BMAX      BMED      VBRO      FPRA 
#> 0.3506675 0.3186672 0.2894165 0.3056724 0.3620338 0.3456037 0.3967495 0.2942797
```

Networks built elsewhere enter through
[`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md).
Pass a gene by gene matrix of edge weights or a sparse matrix.

## 3. Co-expressologs

[`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
tests every ortholog pair in every species pair. It maps the neighbours
of gene A in species 1 to species 2 and asks how many are neighbours of
the ortholog of A. The test runs in both directions, and both must pass.

``` r
edges <- find_coexpressologs(networks, orthologs, seed = 1)
table(edges$type)
#> 
#> conserved        ns 
#>     23765     69083
tab(head(edges[order(-edges$effect_size), c(
  "gene1", "gene2", "hog", "score", "evalue", "q_value", "effect_size",
  "power"
)], 4))
```

| gene1                               | gene2                        | hog         |  score | evalue | q_value | effect_size | power |
|:------------------------------------|:-----------------------------|:------------|-------:|-------:|--------:|------------:|------:|
| Bradi2g56070.v3.2                   | HORVU.MOREX.r3.3HG0299740.V3 | HOG:0022963 | 10.621 |  1.763 |   0.003 |      16.258 | 0.168 |
| Brasyl.3G236800.v2.1                | gene-Fpra_egapxtmp_024636    | HOG:0011023 |  9.979 |  3.217 |   0.005 |      15.458 | 0.136 |
| Brasyl.8G086300.v2.1                | HORVU.MOREX.r3.5HG0494330.V3 | HOG:0000079 | 38.715 |  0.000 |   0.000 |      14.957 | 0.406 |
| VULPIA.BROMOIDES.r1.ctg234G00302340 | gene-Fpra_egapxtmp_020959    | HOG:0011881 | 62.910 |  0.000 |   0.000 |      14.944 | 0.687 |

- `score` is `-log2(p)` and `evalue` is `n_tests * p`.
- `effect_size` is the fold enrichment of shared neighbours.
- `power` is the chance to call the pair had it been conserved. A
  non-significant pair with low power is weak evidence of divergence.

`type == "conserved"` marks a co-expressolog.

### Hypergeometric or rank

The default test is the hypergeometric test of ComPlEx. It is fast and
matches the published method. `method = "rank"` ranks every gene in
species 2 by its co-expression with the translated neighbours, and
reports the ortholog’s place in that ranking. The rank test needs the
argument `null_networks`: one network per species, built by
[`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
from the expression with the samples shuffled.

### Nulls

[`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md)
rewires each network with the degrees held fixed and repeats the test.
It answers how many calls to expect from the degree sequence alone.

``` r
nul <- coexpressolog_null(networks, orthologs, n_cores = 2, seed = 1)
```

`rcomplex(null = TRUE)` runs the cheaper shuffled-expression null and
reports `calls_null` per species pair.

## 4. Density sweep

A 3% density is a convention. A conserved pair should survive a range of
densities.
[`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md)
repeats the test at a set of multipliers of the threshold. A multiplier
above 1 tightens the network and a multiplier below 1 loosens it. The
loosest multiplier must stay inside the stored edges.

``` r
multipliers <- seq(0.975, 1.025, by = 0.0125)
store_floor <- max(vapply(
  networks, function(n) n$store_threshold / n$threshold, numeric(1)
))
stopifnot(min(multipliers) >= store_floor)

sweep <- density_sweep(networks, orthologs,
  multipliers = multipliers, method = "hypergeometric", seed = 1
)
tab(sweep[, c("multiplier", "eff_density", "n_significant")])
```

| multiplier | eff_density | n_significant |
|-----------:|------------:|--------------:|
|      0.975 |       0.052 |         31923 |
|      0.987 |       0.041 |         28097 |
|      1.000 |       0.030 |         23751 |
|      1.012 |       0.020 |         17918 |
|      1.025 |       0.010 |          9190 |

Raw MR values are compressed near the top. A small multiplier moves the
effective density a lot. Count, for each hog, how many levels keep a
conserved edge:

``` r
survive <- table(unlist(lapply(sweep$edges, function(e) {
  unique(e$hog[e$type == "conserved"])
})))
table(factor(survive, levels = seq_len(nrow(sweep))))
#> 
#>    1    2    3    4    5 
#>    3    7   50  278 1162
```

## 5. Modules

A module is a group of genes that are co-expressed with each other.
[`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
runs Leiden clustering. Check that your modules replicate before you
test them. At 20 samples they often do not.

``` r
modules <- lapply(networks, detect_modules, seed = 42)
vapply(modules, function(m) m$n_modules, integer(1))
#> BDIS BSYL HVUL HJUB BMAX BMED VBRO FPRA 
#>    7    8    9    6    7    7    7    7
```

[`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
builds the same object from any gene partition, such as the output of
another tool.

### Preservation

[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
asks, for each ordered species pair, whether the modules of the
reference species keep their wiring in the partner network. It tests two
statistics against a gene-identity permutation null: the density of the
module (`avg.weight`) and the correlation of intramodular connectivity
(`cor.degree`). A module is preserved only when both are significant.
Each module is `conserved`, `moderate`, `diverged` or `untested`.

The trait test in the next section needs all pairs. Pairs between genera
supply the same-trait contrasts. `n_perm = 500` keeps the example short.
Use the default of 10,000 for results you report.

``` r
pairs <- as.data.frame(t(utils::combn(species, 2)))
names(pairs) <- c("species1", "species2")

trait <- stack(clades)
trait <- setNames(as.character(trait$ind), trait$values)

pres <- preservation_paired(modules, networks, orthologs,
  pairs = pairs, group = trait, edges = edges, n_perm = 500, seed = 42
)
table(pres$classification$classification)
#> 
#> conserved  diverged  moderate 
#>       254        35       110
```

### Is divergence associated with the trait?

[`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
asks one question of the whole matrix: are trait-discordant pairs the
less preserved ones? It compares the mean `Zsummary_std` of concordant
pairs (two annuals or two perennials) with that of discordant pairs. It
refers the difference to a null that relabels the traits of the species.

`block` names the genera. It enables a second null that relabels only
within a genus, which holds the phylogeny fixed. It also drops the
within-genus pairs, where trait and relatedness are confounded.

``` r
genus <- list(
  Brachypodium = c("BDIS", "BSYL"), Hordeum = c("HVUL", "HJUB"),
  Briza = c("BMAX", "BMED"), Festuca = c("VBRO", "FPRA")
)
pmt <- suppressWarnings(preservation_matrix_test(
  pres$classification, trait,
  block = genus, seed = 1
))
tab(pmt$class_means)
```

| class               |   n | mean_z |
|:--------------------|----:|-------:|
| concordant          | 171 | 17.805 |
| annual vs perennial | 171 | 17.903 |

``` r
c(
  observed = pmt$observed, p_free = pmt$p_free, p_blocked = pmt$p_blocked,
  floor_free = pmt$free$p_attainable,
  floor_blocked = pmt$blocked$p_attainable
)
#>      observed        p_free     p_blocked    floor_free floor_blocked 
#>   -0.09869348    0.54285714    0.50000000    0.02857143    0.12500000
```

A positive statistic means discordant pairs are less preserved. On this
data the statistic is near zero, so the leaf modules carry no trait
signal. The two nulls agree.

The floor of a relabelling p-value comes from the species set. Eight
species split 4/4 give `choose(8, 4) = 70` labellings, and two of them
share the maximum, so the floor is `2/70 = 0.029`. Within four genera
the floor is `2/16 = 0.125`. No amount of sampling lowers either floor.
[`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
warns when the floor exceeds 0.05, and the chunk above silences that
warning.

### P-value resolution

A permutation p-value cannot fall below `1 / (n_perm + 1)`. Tests whose
true p-value is smaller all return the floor and tie. Rank on
`Zsummary_std` and keep p and q for the significance call.
[`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)
reports the ties:

``` r
pvalue_resolution(pres$classification$q_value, n_perm = 500)
#> p-value resolution
#>   values          399
#>   distinct        71 (18% of the values)
#>   minimum         0.002 (77 values tied there)
#>   at 1            0
#>   n_perm          500, floor 1/(n_perm + 1) = 0.002
#>   status          permutation-limited: the minimum IS the floor
#>   off grid        259 values are not multiples of the floor
#> 
#>   77 values are tied at 0.002 and carry NO ordering among themselves.
#>   Ranking, weighting or top-k selection over that block reads tie-break
#>   noise, not evidence.
#>   Rank on a continuous effect size instead -- Zsummary_std from
#>   module_preservation() -- and keep the p-value for the significance call
#>   alone.
#>   Separating the tie needs n_perm >= 50000. That is a lower bound: the
#>   true p-values may lie far further below the floor.
```

## 6. Hubs

[`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md)
finds the hub genes of each module. The primary centrality is the
weighted degree. It also reports betweenness and eigenvector centrality.

``` r
hubs <- lapply(species, function(s) {
  identify_module_hubs(modules[[s]], networks[[s]], orthologs)
})
names(hubs) <- species
vapply(hubs, function(h) sum(h$is_hub), integer(1))
#> BDIS BSYL HVUL HJUB BMAX BMED VBRO FPRA 
#>  187  191  184  277  182  184  195  204
```

[`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
labels each hog by where it is a hub. It needs
[`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md),
which tests which module of species 1 matches which module of species 2.
Reuse the ortholog map that
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
resolved. The key is the sorted species pair.

``` r
pairs_genus <- lapply(genus, sort)
correspondence <- lapply(pairs_genus, function(sp) {
  key <- paste(sp, collapse = ".")
  module_correspondence(modules[[sp[1]]], modules[[sp[2]]],
    pres$raw[[key]]$map,
    species_ref = sp[1], species_test = sp[2], seed = 1
  )
})
names(correspondence) <- vapply(
  pairs_genus, paste, character(1), collapse = "."
)

hub_class <- classify_hub_conservation(hubs, clades,
  module_comparisons = correspondence
)
table(hub_class$classification)
#> 
#>    annual_specific_hub          conserved_hub        multi_trait_hub 
#>                     99                     36                      1 
#>                non_hub perennial_specific_hub            rewired_hub 
#>                    630                    124                      7 
#>           sporadic_hub 
#>                    603
```

[`get_coexpressed_hogs()`](https://paliocha.github.io/rcomplex/reference/get_coexpressed_hogs.md)
lists the hogs that co-express with a candidate hub in several species.

``` r
candidate <- hub_class$hog[hub_class$classification == "conserved_hub"][1]
partners <- get_coexpressed_hogs(candidate, networks, orthologs,
  clades = clades, min_species = 2L, edges = edges
)
tab(head(partners[, c(
  "partner_hog", "n_species", "coexpressed_traits", "partner_min_q"
)], 5))
```

| partner_hog | n_species | coexpressed_traits | partner_min_q |
|:------------|----------:|:-------------------|--------------:|
| HOG:0016566 |         8 | annual,perennial   |             0 |
| HOG:0020610 |         8 | annual,perennial   |             0 |
| HOG:0018221 |         8 | annual,perennial   |             0 |
| HOG:0017319 |         8 | annual,perennial   |             0 |
| HOG:0011649 |         8 | annual,perennial   |             0 |

## 7. Cliques

A clique is a set of genes, one per species, that are all co-expressed
with each other. The edges of one hog form a graph, and rcomplex finds
cliques on it in two ways. The two answer different questions.

### Species graph

[`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
cliques the graph of species. It returns one best gene assignment per
clique. Give it the full edge table, so the intensity scale uses every
tested pair.

``` r
cliques <- find_cliques(edges, species, min_species = 2L)
table(cliques$n_species)
#> 
#>    2    3    4    5    6    7    8 
#> 2542 1213  523  303  214  127   90
```

A clique of two species is a single conserved edge. Cliques that span
three or more species need conservation across independent lineages.

[`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md)
removes k species at a time and checks whether each clique persists. A
clique that survives the removal of any one species does not depend on a
single lineage.

``` r
stab <- clique_stability(edges, species,
  clades = clades, all_species = species,
  full_cliques = cliques, max_k = 2L
)
table(stab$stability_class)
#> 
#>    0    1    2 
#> 1558  322 3132
```

[`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md)
assigns each hog to one class: `complete` (all species), `partial` (both
clades, not all species), `differentiated` (a clique in each clade and
none across), `trait_specific` (one clade only) or `unclassified`. The
`underpowered` column flags a `differentiated` or `trait_specific` call
that rests on edges with low power.

``` r
result <- classify_cliques(edges, species, clades,
  min_species = 2L, stability = stab
)
table(result$classification)
#> 
#>       complete differentiated        partial trait_specific   unclassified 
#>             90              6           1315             41             48
```

### Gene graph

[`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
cliques the graph of genes. A hog with paralogs gives one clique per
gene combination. It returns three tables: `cliques` has one row per
clique, `members` one row per clique member, and `twins` one row per
collapsed twin. Twins are paralogs with the same partners; one clique
stands for all their combinations (`multiplicity`).
[`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
applies the six-tier taxonomy of Rodriguez *et al.* (2026):
`complete_conserved`, `lineage_specific`, `partial_significant`,
`partial_present`, `differentiated` and `trait_specific`. A clique that
fits none is `unclassified`. A tier that rests on low power becomes
`underpowered`.

``` r
gene_cliques <- gene_clique_graph(edges)
#> gene_clique_graph: 3 ortholog groups exceeded 10 genes in some species; 16 genes dropped
#> gene_clique_graph: 145 genes collapsed into 71 twin classes; cliques carry multiplicity
gene_class <- classify_gene_cliques(gene_cliques, edges, species,
  clades = clades
)
table(gene_class$classification)
#> 
#> complete_conserved    partial_present     trait_specific       unclassified 
#>                271                 34                 16               7461 
#>       underpowered 
#>                 13
```

Pass the unfiltered edge table. The gap tiers need the pairs that were
tested and failed, so they can refuse them.

Use the species graph for stability and for one clique per hog. Use the
gene graph to see every paralog combination and to follow the published
taxonomy. The species graph finds more `trait_specific` hogs because it
asks for a clique among the species it can see, and the gene graph asks
for a complete clique of one clade with an outside species tested at
adequate power.

## 8. Reading the results

Combine the robustness signals. A hog with a `complete` class, a
stability class of 2 or more, and survival at most density levels is a
strong candidate.

``` r
result$n_densities <- as.integer(survive[result$hog])
result$n_densities[is.na(result$n_densities)] <- 0L
keep <- result$classification == "complete" &
  result$stability_class >= 2L &
  result$n_densities >= nrow(sweep) - 1L
strong <- result[keep, ]
nrow(strong)
#> [1] 87
tab(head(strong[order(strong$best_mean_q), c(
  "hog", "n_species", "best_mean_q", "stability_class", "n_densities"
)], 5))
```

| hog         | n_species | best_mean_q | stability_class | n_densities |
|:------------|----------:|------------:|----------------:|------------:|
| HOG:0028978 |         8 |           0 |               2 |           5 |
| HOG:0020610 |         8 |           0 |               2 |           5 |
| HOG:0018534 |         8 |           0 |               2 |           5 |
| HOG:0011862 |         8 |           0 |               2 |           5 |
| HOG:0021123 |         8 |           0 |               2 |           5 |

Check three things before you believe any result.

1.  Read the correlation each network needed (`r_threshold`).
2.  Read the power of the edges behind a negative call.
3.  Read
    [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)
    for every permutation p-value.

The [methods
article](https://paliocha.github.io/rcomplex/articles/methods.md) gives
the statistics behind each step.
