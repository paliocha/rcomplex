# Changelog

## rcomplex (development version)

- [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md)
  reads OrthoFinder files whose hog column is `OrthoGroup`, `N1.HOG` or
  any other case of `HOG` or `Orthogroup`. The EVOTREE
  `Orthogroups_20240823_clean_N1.tsv.gz` needs no reshaping.

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  budgets copy combinations as well as maximal cliques. A hog whose
  cliques stand for more than `max_cliques` combinations (summed
  `multiplicity`) is reported as species blocks. `capped` gains
  `n_combinations`.
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  reads a block as one combination: it never splits into variants. The
  clique count stops at ten times `max_cliques`; `n_cliques_full` is
  then `Inf`. Before, igraph counted every clique, and one dense hog of
  the 8-species leaf data did not finish in 16 minutes. On that data,
  the cliques and tiers of
  [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  without `clades` now take 2.4 min and 7 GB, with 17 and 37 hogs as
  blocks. The wood data have no block.

- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  now lists cliques on three gene graphs, as the published workflow did.
  The strict graph holds edges at `q_value < alpha`, the loose graph
  edges at `q_value < 0.9`. The third graph holds every tested pair of
  the hogs with no complete or partial tier on the first two; it is
  built only with `clades`, since only the clade tiers read it. One
  classification covers all graphs, so `partial_significant` and
  `differentiated` now appear. `cliques`, `members`, `twins` and
  `capped` stack the three runs and gain `run` and `alpha_graph`;
  `(run, clique_id)` joins them to `classification`.
  [`print()`](https://rdrr.io/r/base/print.html) counts the cliques of
  each graph. A clique maximal on two graphs appears once for each, so
  clique counts in [`summary()`](https://rdrr.io/r/base/summary.html)
  grow. On the wood data (clades angiosperm, conifer, pinus), hogs per
  tier, before and after:

  - complete_conserved: 938 and 938.
  - partial_significant: 0 and 517.
  - partial_present: 247 and 227.
  - differentiated: 0 and 1,282.
  - lineage_specific: 2,973 and 2,973.
  - trait_specific: 1,433 and 1,099.
  - underpowered: 511 and 551.
  - unclassified: 2,701 and 6,357.

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  replaces the 10-copies-per-species cap by a clique budget,
  `max_cliques = 1e5` per hog. Twins collapse first. A hog over the
  budget is reported as species blocks: one row per group of species
  wired together, with the best copy per species. Gene ids break ties,
  never row order
  ([\#79](https://github.com/paliocha/rcomplex/issues/79)). The result
  gains `capped`, one row per such hog with its full clique count.
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  flags those hogs’ cliques with `capped`;
  [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  stores `capped`. On the wood data no hog is over the budget (2,057
  genes were capped). On the 8-species leaf data fourteen families are
  blocks instead of 214 capped families.

- A group of species wired pairwise need not share one set of copies.
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md),
  [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md)
  and the species blocks of
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  now report its largest species subsets that do. Before, such a group
  gave no row. A hog gets at most 4,096 subset searches. Ties in mean q
  go to the larger effect size, then to the gene id. On the wood data,
  [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md)
  moves 347 hogs: unclassified 2,335 to 2,035, partial 5,423 to 5,689,
  trait_specific 7,832 to 7,862, differentiated 413 to 417.

- [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  and
  [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md)
  drop the 10-copies-per-species cap. The best-gene search order is
  deterministic: gene ids break ties, never row order.
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  gains `exhaustive`, `FALSE` when the search hit its step budget.

- `rcomplex(null = TRUE)` now shuffles hog labels among genes of the
  same expression level and reruns edges, cliques and tiers on the
  shuffled table. `calls_null` counts the calls that abundance and
  network structure alone produce; the tier table gains `cliques_null`.
  The shuffled-expression null gave 0 calls on every data set and is
  gone from the driver;
  [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  stays. `null = TRUE` no longer needs `expr`.

- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  stores `gene_mean`.

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  collapses twin paralogs before it lists the cliques. Twins are genes
  of one species with the same partners. A clique stands for every
  choice of one twin per class, so the collapse is exact. `cliques`
  gains `multiplicity`, `members` gains `n_copies`, and a new table
  `twins` lists the collapsed genes. On the full 8-species leaf data,
  28,343,526 cliques become 3,555,837. The build takes 24 s instead of
  100 s, and 2.9 GB instead of 8.9 GB (max RSS). At 12,000 genes per
  species, 1,452,683 cliques become 42,709.

- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  splits each twin class into copies with the same test status. A clique
  gives one row per choice, with `variant` and `multiplicity`. Every
  combination gets its exact tier. On the full leaf data it classifies
  8,464,243 rows in 169 s at 10.2 GB, instead of 28.3 M rows in 6 min at
  15.0 GB. At 12,000 genes it gives 186,032 rows in 4 s instead of 11 s.

- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  splits a twin class per clique, and only where its copies differ in
  what that clique reads: the species it misses and its other members.
  On the full leaf data, 8,464,243 rows become 3,582,532, in 61 s at 7.8
  GB instead of 169 s at 10.2 GB. At 12,000 genes, 186,032 rows become
  43,039.

- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  stores `twins`, and
  [`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
  writes `twins.tsv`. [`print()`](https://rdrr.io/r/base/print.html)
  shows the distinct cliques, with the count including paralog copies in
  brackets. The [`summary()`](https://rdrr.io/r/base/summary.html) tier
  table gains `multiplicity`.

## rcomplex 0.4.0

rcomplex 0.4.0 is a rewrite of the user surface. The package has 31
exports instead of 57. It is not backward compatible: nothing is
aliased.

### New functions

- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  runs the whole analysis in one call: networks, co-expressologs, gene
  cliques and their tiers.
  [`print()`](https://rdrr.io/r/base/print.html) shows the genes,
  samples, density and `r_threshold` of each species.
  [`summary()`](https://rdrr.io/r/base/summary.html) returns the tier
  table and, with `null = TRUE`, the calls on shuffled networks.
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
  the classification.
- [`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
  writes the tables of an `rcomplex` result as TSV files.
- [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md)
  reads OrthoFinder (`N0.tsv`, `Orthogroups.tsv`), PLAZA and long-table
  ortholog files into the table `species`, `gene`, `hog`. It replaces
  `parse_orthologs()`.
- [`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md)
  wraps a dense matrix, a sparse matrix or an edge list in the network
  object, so a network built elsewhere enters the analysis.
- [`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md)
  returns the clades of a species tree.

### New behaviour

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  returns two tables: `cliques`, one row per clique with an integer
  `clique_id`, and `members`, one row per clique member. The per-clique
  statistics are no longer repeated on every member row. `id_prefix` is
  gone: pass several runs to
  [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  as a list. On the full 8-species leaf data the result takes about 5
  GiB (1.8 + 3.2 GiB) instead of about 20 GB. It is built in 100 s
  instead of 335 s. The peak process memory (max RSS) is 8.9 GB instead
  of 26.5 GB.
- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  is about 40 times faster, with the same tiers. It works on the sparse
  gene x clique incidence. At 12,000 genes per species it classifies
  1.45 M cliques in 11 s instead of 445 s. On the full 8-species leaf
  data it classifies 28.3 M cliques in 6 min. With
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  the peak process memory is 15.0 GB.
- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  stores `members` beside `cliques`, and
  [`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
  writes `members.tsv`.
  [`summary()`](https://rdrr.io/r/base/summary.html) counts each tier in
  cliques, HOGs and species patterns.
- `clades` is a named list of species vectors. Clades may nest. A
  species in no clade forms a clade of its own. It replaces
  `species_trait` and `lineage`. The classifier tables gain a `clade`
  column.
- `compute_network(partition = )` correlates within each level of a
  sample partition, sets negative correlations to zero and averages. It
  aggregates a compendium of mixed designs (after TEA-GCN).
  [`split_layers()`](https://paliocha.github.io/rcomplex/reference/split_layers.md)
  handles a designed experiment with replicates.
- `compute_network(sign = "negative")` builds the network of
  anticorrelation. It replaces `abs_cor`.
- Every edge table has the columns
  `gene1 gene2 hog score evalue q_value effect_size power`, then
  `species1 species2 p_value n_tests jaccard type`. `score` is
  `-log2(p_value)` and `evalue` is `n_tests * p_value`.
  [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  and
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  add a `score`.
- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  records `r_threshold`: the weakest correlation that passed the density
  threshold. [`print()`](https://rdrr.io/r/base/print.html) shows it, so
  a user sees what a sample size bought.
- `rcomplex(null = TRUE)` repeats the comparison on shuffled networks
  and reports `calls_null` and `false_call_rate` per species pair.
- [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md)
  no longer returns a `persistence` column.
- `find_coexpressologs(method = "permutation")` joins HOG q-values to
  pairs with [`match()`](https://rdrr.io/r/base/match.html). It used a
  positional lookup, which gave wrong p and q values for integer hog
  ids.
- `find_coexpressologs(method = "hypergeometric")` is the name of the
  default test. The alias `"analytical"` is gone.

### Renamed

- Columns: `Species1`, `Species2` to `species1`, `species2`; `Gene1`,
  `Gene2` to `gene1`, `gene2`; `HOG` to `hog`; `q.value` to `q_value`;
  `p.value` to `p_value`; `p.calibrated` to `p_calibrated`; the
  per-direction columns to lower case (`species1.p_value_con`).
- Arguments: `sp1`, `sp2` to `species1`, `species2`; `sp_ref`, `sp_test`
  to `species_ref`, `species_test`; `species_trait` and `lineage` to
  `clades`; `abs_cor` to `sign`.
- Functions: `parse_orthologs()` to
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).

### Removed functions

The 17 research probes from the 2026 Orion benchmarks are gone:
`module_auroc()`, `module_auroc_reciprocal()`, `module_replication()`,
`subspace_preservation()`, `as_preservation_matrix()`,
`recurrence_graph()`, `recurrence_modules()`,
`conservation_pattern_table()`, `conservation_lattice()`,
`bicm_species_z()`, `coexpressolog_strength()`,
`suggest_reference_density()`, `clique_persistence()`,
`clique_perturbation_test()`, `clique_intensity_test()`,
`tag_permutation()` and `characterize_hubs()`.

### Demoted to internal

[`compare_neighborhoods()`](https://paliocha.github.io/rcomplex/reference/compare_neighborhoods.md),
[`compare_specificity()`](https://paliocha.github.io/rcomplex/reference/compare_specificity.md),
[`comparison_to_edges()`](https://paliocha.github.io/rcomplex/reference/comparison_to_edges.md),
[`summarize_comparison()`](https://paliocha.github.io/rcomplex/reference/summarize_comparison.md),
[`summarize_specificity()`](https://paliocha.github.io/rcomplex/reference/summarize_specificity.md),
[`permutation_hog_test()`](https://paliocha.github.io/rcomplex/reference/permutation_hog_test.md),
[`run_pairwise_comparisons()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md),
[`mr_block()`](https://paliocha.github.io/rcomplex/reference/mr_block.md),
[`as_sparse_network()`](https://paliocha.github.io/rcomplex/reference/as_sparse_network.md),
`extract_orthologs()` and
[`all_species_pairs()`](https://paliocha.github.io/rcomplex/reference/all_species_pairs.md)
stay in the package but are no longer exported. Call them with
`rcomplex:::`.

### Removed arguments

Every argument that no documented workflow set is removed, and its
default is fixed.

- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md):
  `abs_cor` (now `sign`), `min_var`, `use_torch`.
- [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  `species_pairs`, `alternative`, `alpha`, `min_exceedances`,
  `max_permutations`, `pi0_method`, `out_file`, `p0`.
- [`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md):
  the same, plus `filter_zero`, minus `out_file`.
- [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md):
  `method` (Infomap, SBM), `n_iterations`, `nb_trials`,
  `consensus_threshold`, and the K = 1 test (`test_k1`, `n_perm_k1`,
  `alpha_k1`).
- [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md):
  `map`, `cliques`, `min_module_size`, `binary`, `alpha`,
  `qvalue_method`, `sensitivity`, `copy_draws`.
- [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md):
  `cliques`, `alpha`, `z_conserved`.
- [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.md):
  `alpha`, `z_conserved`, `z_scale`.
- [`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md):
  `statistic`, `exclude_within_block`, `n_perm`, `enum_max`,
  `n_perm_pres`.
- [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md):
  `qvalue_method`.
- [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md):
  `cliques`, `rank_by`, `alpha`.
- [`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md):
  `comparison`, `centrality`, `top_n`, `top_fraction`,
  `min_module_size`.
- [`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md):
  `alpha`, `jaccard_threshold`, `min_trait_fraction`,
  `correspondence_threshold`.
- [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md):
  `max_genes_per_sp`, `max_missing_edges`, `edge_type`, `sweep`,
  `min_stability_class`, `min_persistence`.
- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md):
  `max_gap`, `cross_max`.
- [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md):
  `min_species`, `max_genes_per_sp`, `jaccard_threshold`, `edge_type`,
  `cost_weights`.
- [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md),
  [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md):
  `max_genes_per_sp`, `max_missing_edges`, `edge_type`.
- [`clique_threshold_sweep()`](https://paliocha.github.io/rcomplex/reference/clique_threshold_sweep.md):
  `species_pairs`, `alternative`, `alpha`, `max_genes_per_sp`,
  `max_missing_edges`, `edge_type`, `jaccard_threshold`.
- [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md):
  `statistic`, `n_perm`, `filter_zero`.
- [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md):
  `use_torch`.
- [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md):
  `se_list`, `hog_col`. It reads the long table of
  [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md).
- [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md):
  `min_size`.
- [`get_coexpressed_hogs()`](https://paliocha.github.io/rcomplex/reference/get_coexpressed_hogs.md):
  `species`.

### Other

- The tutorial is now the article `walkthrough`, without the research
  probes. A new `quickstart` vignette shows the driver. `methods` gains
  sections on harmonisation, signed networks, `score` and `evalue`, and
  clades.
