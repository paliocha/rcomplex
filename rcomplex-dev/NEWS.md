# rcomplex 0.4.0

rcomplex 0.4.0 is a rewrite of the user surface. The package has 31
exports instead of 57. It is not backward compatible: nothing is
aliased.

## New functions

- `rcomplex()` runs the whole analysis in one call: networks,
  co-expressologs, gene cliques and their tiers. `print()` shows the
  genes, samples, density and `r_threshold` of each species. `summary()`
  returns the tier table and, with `null = TRUE`, the calls on shuffled
  networks. `as.data.frame()` returns the classification.
- `write_rcomplex()` writes the tables of an `rcomplex` result as TSV
  files.
- `read_orthologs()` reads OrthoFinder (`N0.tsv`, `Orthogroups.tsv`),
  PLAZA and long-table ortholog files into the table `species`, `gene`,
  `hog`. It replaces `parse_orthologs()`.
- `as_network()` wraps a dense matrix, a sparse matrix or an edge list in
  the network object, so a network built elsewhere enters the analysis.
- `clades_from_tree()` returns the clades of a species tree.

## New behaviour

- `clades` is a named list of species vectors. Clades may nest. A
  species in no clade forms a clade of its own. It replaces
  `species_trait` and `lineage`. The classifier tables gain a `clade`
  column.
- `compute_network(partition = )` correlates within each level of a
  sample partition, sets negative correlations to zero and averages. It
  aggregates a compendium of mixed designs (after TEA-GCN).
  `split_layers()` handles a designed experiment with replicates.
- `compute_network(sign = "negative")` builds the network of
  anticorrelation. It replaces `abs_cor`.
- Every edge table has the columns `gene1 gene2 hog score evalue
  q_value effect_size power`, then `species1 species2 p_value n_tests
  jaccard type`. `score` is `-log2(p_value)` and `evalue` is
  `n_tests * p_value`. `find_cliques()` and `gene_clique_graph()` add a
  `score`.
- `compute_network()` records `r_threshold`: the weakest correlation that
  passed the density threshold. `print()` shows it, so a user sees what
  a sample size bought.
- `rcomplex(null = TRUE)` repeats the comparison on shuffled networks
  and reports `calls_null` and `false_call_rate` per species pair.
- `classify_cliques()` no longer returns a `persistence` column.
- `find_coexpressologs(method = "hypergeometric")` is the name of the
  default test. The alias `"analytical"` is gone.

## Renamed

- Columns: `Species1`, `Species2` to `species1`, `species2`; `Gene1`,
  `Gene2` to `gene1`, `gene2`; `HOG` to `hog`; `q.value` to `q_value`;
  `p.value` to `p_value`; `p.calibrated` to `p_calibrated`; the
  per-direction columns to lower case (`species1.p_value_con`).
- Arguments: `sp1`, `sp2` to `species1`, `species2`; `sp_ref`, `sp_test`
  to `species_ref`, `species_test`; `species_trait` and `lineage` to
  `clades`; `abs_cor` to `sign`.
- Functions: `parse_orthologs()` to `read_orthologs()`.

## Removed functions

The 17 research probes from the 2026 Orion benchmarks are gone:
`module_auroc()`, `module_auroc_reciprocal()`, `module_replication()`,
`subspace_preservation()`, `as_preservation_matrix()`,
`recurrence_graph()`, `recurrence_modules()`,
`conservation_pattern_table()`, `conservation_lattice()`,
`bicm_species_z()`, `coexpressolog_strength()`,
`suggest_reference_density()`, `clique_persistence()`,
`clique_perturbation_test()`, `clique_intensity_test()`,
`tag_permutation()` and `characterize_hubs()`.

## Demoted to internal

`compare_neighborhoods()`, `compare_specificity()`,
`comparison_to_edges()`, `summarize_comparison()`,
`summarize_specificity()`, `permutation_hog_test()`,
`run_pairwise_comparisons()`, `mr_block()`, `as_sparse_network()`,
`extract_orthologs()` and `all_species_pairs()` stay in the package but
are no longer exported. Call them with `rcomplex:::`.

## Removed arguments

Every argument that no documented workflow set is removed, and its
default is fixed.

- `compute_network()`: `abs_cor` (now `sign`), `min_var`, `use_torch`.
- `find_coexpressologs()`: `species_pairs`, `alternative`, `alpha`,
  `min_exceedances`, `max_permutations`, `pi0_method`, `out_file`, `p0`.
- `density_sweep()`: the same, plus `filter_zero`, minus `out_file`.
- `detect_modules()`: `method` (Infomap, SBM), `n_iterations`,
  `nb_trials`, `consensus_threshold`, and the K = 1 test (`test_k1`,
  `n_perm_k1`, `alpha_k1`).
- `module_preservation()`: `map`, `cliques`, `min_module_size`, `binary`,
  `alpha`, `qvalue_method`, `sensitivity`, `copy_draws`.
- `preservation_paired()`: `cliques`, `alpha`, `z_conserved`.
- `classify_preservation()`: `alpha`, `z_conserved`, `z_scale`.
- `preservation_matrix_test()`: `statistic`, `exclude_within_block`,
  `n_perm`, `enum_max`, `n_perm_pres`.
- `module_correspondence()`: `qvalue_method`.
- `resolve_ortholog_map()`: `cliques`, `rank_by`, `alpha`.
- `identify_module_hubs()`: `comparison`, `centrality`, `top_n`,
  `top_fraction`, `min_module_size`.
- `classify_hub_conservation()`: `alpha`, `jaccard_threshold`,
  `min_trait_fraction`, `correspondence_threshold`.
- `classify_cliques()`: `max_genes_per_sp`, `max_missing_edges`,
  `edge_type`, `sweep`, `min_stability_class`, `min_persistence`.
- `classify_gene_cliques()`: `max_gap`, `cross_max`.
- `clique_stability()`: `min_species`, `max_genes_per_sp`,
  `jaccard_threshold`, `edge_type`, `cost_weights`.
- `find_cliques()`, `gene_clique_graph()`: `max_genes_per_sp`,
  `max_missing_edges`, `edge_type`.
- `clique_threshold_sweep()`: `species_pairs`, `alternative`, `alpha`,
  `max_genes_per_sp`, `max_missing_edges`, `edge_type`,
  `jaccard_threshold`.
- `coexpressolog_null()`: `statistic`, `n_perm`, `filter_zero`.
- `null_network()`: `use_torch`.
- `prepare_orthologs()`: `se_list`, `hog_col`. It reads the long table
  of `read_orthologs()`.
- `as_modules()`: `min_size`.
- `get_coexpressed_hogs()`: `species`.

## Other

- The tutorial is now the article `walkthrough`, without the research
  probes. A new `quickstart` vignette shows the driver. `methods` gains
  sections on harmonisation, signed networks, `score` and `evalue`, and
  clades.
