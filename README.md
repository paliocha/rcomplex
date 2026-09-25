# rcomplex

<!-- badges: start -->
[![R-CMD-check](https://github.com/paliocha/rcomplex/actions/workflows/r.yml/badge.svg)](https://github.com/paliocha/rcomplex/actions/workflows/r.yml)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/paliocha/rcomplex)
<!-- badges: end -->

rcomplex compares gene co-expression networks between species. It asks
which genes, modules and gene groups keep their co-expression partners
across species, and whether that differs between species with different
traits, such as annual and perennial life cycles. You build one network
per species from RNA-seq data, link the networks through ortholog
groups, and test conservation for single genes, for modules (groups of
genes that are co-expressed with each other), and for ortholog groups
across many species at once. The method extends ComPlEx by [Netotea
*et al.* (2014)](https://doi.org/10.1186/1471-2164-15-106).

For the statistics behind each test, see [Methods in
detail](https://paliocha.github.io/rcomplex/articles/methods.html).

## What you need

- One expression matrix per species: genes in rows, samples in
  columns, normalised and on a log-like scale (for example DESeq2 VST).
  A plain matrix or a `SummarizedExperiment` both work.
- Enough samples: aim for at least 15 to 20 samples per species,
  taken under matched conditions (same tissues, time points or
  treatments in every species). A co-expression network built from few
  samples is noisy: many of its edges are chance correlations, and
  conservation tests then have little to find. See
  [pitfalls](#reading-the-results-and-common-pitfalls).
- Ortholog groups: a table that places genes of all species into
  ortholog groups (HOGs), from OrthoFinder, FastOMA, PLAZA or similar.
  See [the file format](#ortholog-file-format).
- For trait tests, several species per trait. A comparison of one
  annual with one perennial species cannot separate the trait from
  everything else that differs between the two. Use several species per
  trait, ideally as phylogenetic pairs (one annual and one perennial in
  each genus).

## Installation

```r
devtools::install_github("paliocha/rcomplex")
```

**System requirements:** C++23 compiler, GNU make. OpenMP is optional
but recommended for parallel permutation and stability tests.

## A first analysis in five steps

This example uses data that ships with the package: VST-normalised
RNA-seq for 2,000 to 4,000 genes in each of four Pooideae grasses, two
annual and two perennial. It keeps only the 20 leaf samples per species
(see [pitfalls](#reading-the-results-and-common-pitfalls) on pooled
tissues).
It runs in under a minute.

### 1. Load expression data and ortholog pairs

```r
library(rcomplex)
library(SummarizedExperiment)

se_list <- readRDS(system.file("extdata", "pooideae_vignette.rds",
                               package = "rcomplex"))
species <- c("BDIS", "BSYL", "HVUL", "HJUB")
se_list <- lapply(se_list[species], function(se) se[, se$tissue == "leaf"])

# Ortholog pairs from the HOG column of rowData(). With an ortholog
# file instead, use parse_orthologs("orthogroups.txt", "BDIS", "BSYL").
orthologs <- prepare_orthologs(se_list)
head(orthologs)
```

`orthologs` has one row per pair of genes that share an ortholog group.
A gene with paralogs appears in several rows.

### 2. Build one network per species

```r
networks <- lapply(se_list, compute_network, density = 0.03)
```

`compute_network()` correlates every pair of genes and turns the
correlations into Mutual Rank (MR) scores. MR asks whether two genes
rank each other among their best partners, so a gene correlated with
everything does not dominate. The network keeps the top 3% of gene
pairs (`density = 0.03`) as edges. Each gene's neighbourhood is the set
of genes it shares an edge with.

### 3. Find co-expressologs

```r
edges <- find_coexpressologs(networks, orthologs, seed = 1)
table(edges$type)
head(edges[order(-edges$effect_size), ])
```

Each row is one ortholog pair in one pair of species. The test asks
whether the gene's neighbours in species A have orthologs among its
partner's neighbours in species B more often than chance (a
hypergeometric test). `type == "conserved"` marks a co-expressolog: an
ortholog pair that kept its co-expression partners in both species.
The test runs in both directions and both must pass.

- `q.value`: the false discovery rate at which this pair would be
  called conserved.
- `effect_size`: fold enrichment, the observed number of shared
  neighbours divided by the number expected by chance. An effect of 5
  means five times more shared partners than chance.
- `power`: the chance that the test would have called this pair had it
  been conserved at a typical effect size. Genes with few neighbours
  have low power, so a non-significant result for them is weak
  evidence of divergence.

### 4. Detect modules and test whether they are preserved

```r
mods <- detect_modules(networks$BDIS, resolution = c(0.5, 1, 1.5, 2),
                       objective_function = "modularity",
                       n_perm_k1 = 20L, seed = 1)
pres <- module_preservation(mods, networks$BDIS, networks$BSYL, orthologs,
                            edges = edges, sp_ref = "BDIS",
                            sp_test = "BSYL", n_perm = 1000L, seed = 1)
classify_preservation(pres)[, c("module", "size", "classification",
                                "Zsummary_std", "q.value")]
```

`detect_modules()` first tests whether the network has module structure
at all, then groups genes into modules with the Leiden algorithm over
several resolutions and keeps the consensus. `module_preservation()`
then asks whether each *Brachypodium distachyon* module keeps its wiring
in the *B. sylvaticum* network: are its genes' orthologs still densely
connected to each other, and are the same genes still the hubs? A module
whose genes all survive but lost their connections counts as diverged.

`Zsummary_std` is the effect size: how many standard deviations the
module's preservation lies above what random gene sets of the same size
show. Values above 10 are strong preservation. `q.value` is the false
discovery rate for calling the module preserved. `classification` is
`conserved` (significant, `Zsummary_std` of at least 10), `moderate`
(significant, weaker), `diverged` (not significant) or `untested`
(a statistic could not be computed). Preservation is directional: run
the reverse direction to ask about *B. sylvaticum* modules in
*B. distachyon*.

With 1,000 permutations no p-value can go below 1/1001. In this example
every module ties at the same smallest q-value, while `Zsummary_std`
ranges from about 8 to 41. Rank modules on `Zsummary_std`, not on q.

### 5. Find conserved and lineage-specific cliques

```r
cliques <- rbind(
  gene_clique_graph(edges, alpha_graph = 0.1, id_prefix = "strict_"),
  gene_clique_graph(edges, alpha_graph = 0.9, id_prefix = "loose_")
)
life_cycle <- c(BDIS = "annual", HVUL = "annual",
                BSYL = "perennial", HJUB = "perennial")
classes <- classify_gene_cliques(cliques, edges, species,
                                 lineage = life_cycle)
table(classes$classification)
```

A clique is a set of genes, one per species, from one ortholog group,
where every pair is a co-expressolog. A clique over all species is an
ortholog group whose co-expression is conserved across the whole
sample. `classify_gene_cliques()` sorts cliques into tiers:

- `complete_conserved`: every species present, every pair conserved.
- `partial_significant`: every species present, a few pairs just miss
  the strict cut-off.
- `partial_present`: conserved in all but one species, and that species
  lacks the gene, was not tested, or had too little power.
- `lineage_specific` or `trait_specific`: conserved within one group
  (here a life cycle), and absent or rejected outside it.
- `differentiated`: conserved within each group but not between them.
- `underpowered`: would be specific or differentiated, but the tests
  that separate the groups had too little power to count as evidence.
- `unclassified`: fits none of the tiers, for example a clique that
  lacks a species whose test was run and failed.

Pass `edges` unfiltered. A tier that claims divergence needs to see
the pairs that were tested and failed. The group-specific tiers need a
clique of at least three species inside one group, so they stay empty
in this four-species example; the tutorial runs eight species.

## Which test should I use?

| Biological question | Level | Functions |
|---|---|---|
| Has this gene kept its co-expression partners in the other species? | Gene pair | `find_coexpressologs()` |
| Is this gene family conserved as a whole, counting all its paralogs? | Ortholog group | `find_coexpressologs(method = "permutation")`, `permutation_hog_test()` |
| Are there more co-expressologs than the network structure alone would produce? | Whole network | `coexpressolog_null()` |
| Does a call depend on the chosen network density? | Gene pair | `density_sweep()`, `coexpressolog_strength()` |
| Does this module keep its wiring in the other species? | Module | `module_preservation()`, `classify_preservation()` |
| Which module in species B corresponds to this module in species A? | Module | `module_correspondence()` |
| Are the hub genes of a module the same across traits? | Gene within module | `identify_module_hubs()`, `classify_hub_conservation()` |
| Do species pairs that differ in the trait preserve fewer modules? | Trait | `preservation_paired()`, `preservation_matrix_test()` |
| Do the same ortholog groups sit in diverged modules in every trait contrast? | Trait | `tag_permutation()` |
| Which ortholog groups are conserved in all species, or in one lineage or trait group only? | Ortholog group | `gene_clique_graph()`, `classify_gene_cliques()` |
| Does a trait-exclusive clique survive when species are left out? | Ortholog group | `find_cliques()`, `clique_stability()`, `classify_cliques()` |
| What does this ortholog group co-express with in each species? | Ortholog group | `get_coexpressed_hogs()` |

## Reading the results and common pitfalls

### Few samples make noisy networks

With 10 samples, a correlation of
0.6 or stronger turns up by chance in about 7% of unrelated gene pairs;
with 20 samples, in 0.5%. Chance edges differ between species, so
conservation looks weaker than it is and modules are not reproducible.
Use 15 to 20 or more samples per species under comparable conditions.

### Network density changes the answer

A denser network gives each gene
more neighbours, which raises power but adds weaker edges. The default
of 3% is a convention. Check that your conclusions hold
over a range with `density_sweep()`, or pick a density from the data
with `suggest_reference_density()`.

### Pooled tissues dominate the network

If leaf and root samples are in
one matrix, the strongest correlations separate leaf genes from root
genes. Most modules then reflect tissue, and species comparisons mostly
compare tissue programs. Build and compare networks one tissue at a
time.

### Permutation tests have a floor

A permutation p-value compares the
real labelling of species with every other possible labelling. With few
species there are few labellings, and the p-value cannot go below a
floor set by the design. Eight species split four annual and four
perennial can be labelled in `choose(8, 4) = 70` ways. Swapping the
labels "annual" and "perennial" gives the same statistic, so at least
two labellings tie at the top and the smallest p-value is 2/70 = 0.029.
If you only allow swaps within each of four genera (to respect
phylogeny), there are 2^4 = 16 labellings and the floor is 2/16 =
0.125, which can never reach 0.05. This is a property of the species
set, not of the number of permutations. `preservation_matrix_test()`
reports the floor, and `pvalue_resolution()` reports how many p- or
q-values are tied at it.

### Paralogs complicate the ortholog map

A gene family with several
copies in one species has several candidate partners in the other.
`resolve_ortholog_map()` picks one copy per gene, first from cliques,
then from mutual best co-expressologs, and leaves the rest unresolved.
`module_preservation()` calls it for you when you pass `edges`. For
recent duplicates with near-identical expression, `reduce_orthogroups()`
can merge them before you build networks.

### Rank on effect size, call on q

Many strong results tie at the
smallest q-value a permutation test can give. Use `q.value` to decide
what is significant, and `effect_size`, `Zsummary_std` or
`mean_effect_size` to rank what is.

## Ortholog file format

`parse_orthologs()` reads a tab-delimited file with at least these
columns:

| Column | Description |
|--------|-------------|
| `species` | Species code of the anchor species |
| `gene_id` | Gene identifier in the anchor species |
| `gene_content` | Members in other species, as `code:gene1,gene2;code:gene3` |

Each row defines one ortholog group from the anchor species' point of
view. `gene_content` lists the members per species, separated by
semicolons. PLAZA writes this format; OrthoFinder and FastOMA output can
be reshaped into it. If your expression data are `SummarizedExperiment`
objects with a HOG column in `rowData()`, use `extract_orthologs()` or
`prepare_orthologs()` instead and skip the file.

## Function index

Full documentation is in the [reference
pages](https://paliocha.github.io/rcomplex/reference/).

**Input and networks**

- [`parse_orthologs()`](https://paliocha.github.io/rcomplex/reference/parse_orthologs.html): read an ortholog group file into ortholog pairs.
- [`extract_orthologs()`](https://paliocha.github.io/rcomplex/reference/extract_orthologs.html): ortholog pairs from the HOG column of two `SummarizedExperiment` objects.
- [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.html): ortholog pairs for many species, optionally after paralog reduction.
- [`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.html): merge paralogs with near-identical expression within each HOG.
- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.html): correlation, MR or CLR normalisation and density threshold for one species.
- [`as_sparse_network()`](https://paliocha.github.io/rcomplex/reference/as_sparse_network.html): convert a dense network object to sparse storage.
- [`mr_block()`](https://paliocha.github.io/rcomplex/reference/mr_block.html): exact MR values for a gene subset, including pairs below the stored edges.
- [`suggest_reference_density()`](https://paliocha.github.io/rcomplex/reference/suggest_reference_density.html): choose a network density from a scale-free fit diagnostic.
- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.html): container that carries species, networks and results through the pipeline.

**Gene and ortholog-group conservation**

- [`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.html): hypergeometric neighbourhood tests for one species pair.
- [`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.html): q-values and summaries for `compare_neighborhoods()` output.
- [`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.html): convert comparison results to the edge table used downstream.
- [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.html): co-expressolog calls for all species pairs (alias `run_pairwise_comparisons()`).
- [`permutation_hog_test()`](https://paliocha.github.io/rcomplex/reference/permutation_hog_test.html): permutation test of conservation for whole ortholog groups.
- [`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.html): rerun the co-expressolog calls at several network densities.
- [`coexpressolog_strength()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_strength.html): edge strength integrated over several densities.
- [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.html): degree-preserving rewiring null for co-expressolog counts.

**Modules and hubs**

- [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.html): Leiden, Infomap or SBM modules, with multi-resolution consensus.
- [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.html): pick one paralog copy per gene where evidence allows.
- [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.html): test whether modules keep their density and hubs in another species.
- [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.html): label modules conserved, moderate, diverged or untested.
- [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.html): match modules across species by ortholog overlap.
- [`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.html): hub genes within each module.
- [`characterize_hubs()`](https://paliocha.github.io/rcomplex/reference/characterize_hubs.html): bridge and betweenness metrics for hub genes.
- [`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.html): hub conservation across trait groups.
- [`get_coexpressed_hogs()`](https://paliocha.github.io/rcomplex/reference/get_coexpressed_hogs.html): co-expression partners of one ortholog group across species.

**Trait tests**

- [`all_species_pairs()`](https://paliocha.github.io/rcomplex/reference/all_species_pairs.html): table of every species pair for `preservation_paired()`.
- [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.html): module preservation over many species pairs, both directions.
- [`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.html): do trait-discordant species pairs preserve less?
- [`tag_permutation()`](https://paliocha.github.io/rcomplex/reference/tag_permutation.html): do the same ortholog groups recur in diverged modules across contrasts?
- [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.html): how many p- or q-values are tied at the permutation floor.

**Cliques**

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.html): all maximal cliques of each ortholog group's gene graph.
- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.html): conservation tiers for gene-graph cliques.
- [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.html): species-graph cliques with one best gene per species.
- [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.html): conservation classes for ortholog groups from species-graph cliques.
- [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.html): leave-k-species-out stability of trait-exclusive cliques.
- [`clique_persistence()`](https://paliocha.github.io/rcomplex/reference/clique_persistence.html): how far each clique's supporting edges sit above the threshold.
- [`clique_threshold_sweep()`](https://paliocha.github.io/rcomplex/reference/clique_threshold_sweep.html): clique survival at stricter densities.
- [`clique_perturbation_test()`](https://paliocha.github.io/rcomplex/reference/clique_perturbation_test.html): clique survival under added noise.
- [`clique_intensity_test()`](https://paliocha.github.io/rcomplex/reference/clique_intensity_test.html): are a clique's edges stronger than matched random edges?

## Citation

If you use rcomplex, please cite the method paper and the package:

- Netotea, S., Sundell, D., Street, N. R. & Hvidsten, T. R. (2014).
  ComPlEx: conservation and divergence of co-expression networks in
  *A. thaliana*, *Populus* and *O. sativa*. *BMC Genomics*, 15, 106.
  [doi:10.1186/1471-2164-15-106](https://doi.org/10.1186/1471-2164-15-106)
- Paliocha, M. rcomplex: Comparative Co-Expression Network Analysis
  Across Species. R package version 0.3.1.
  <https://github.com/paliocha/rcomplex>

Key references for the methods:

- Obayashi, T. & Kinoshita, K. (2009). Rank of correlation coefficient
  as a comparable measure for biological significance of gene
  coexpression. *DNA Research*, 16(5), 249-260.
  [doi:10.1093/dnares/dsp016](https://doi.org/10.1093/dnares/dsp016)
- Jeub, L. G. S., Sporns, O. & Fortunato, S. (2018). Multiresolution
  consensus clustering in networks. *Scientific Reports*, 8, 3259.
  [doi:10.1038/s41598-018-21352-7](https://doi.org/10.1038/s41598-018-21352-7)
- Ritchie, S. C. *et al.* (2016). A scalable permutation approach
  reveals replication and preservation patterns of network modules in
  large datasets. *Cell Systems*, 3(1), 71-82.
  [doi:10.1016/j.cels.2016.06.012](https://doi.org/10.1016/j.cels.2016.06.012)
- Rodriguez, E. *et al.* (2026). Comparative regulomics of wood
  formation across dicot and conifer trees. *Nature Communications*,
  17, 8916.
  [doi:10.1038/s41467-026-75624-2](https://doi.org/10.1038/s41467-026-75624-2)

The full reference list is in [Methods in
detail](https://paliocha.github.io/rcomplex/articles/methods.html#references).

## License

MIT
