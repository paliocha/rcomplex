# rcomplex: Comparative Co-Expression Network Analysis Across Species

Compares gene co-expression networks across species. Maps orthologous
genes, builds one network per species, and tests conservation. The tests
work at the gene, hog, module and clique levels.

## Details

Start with
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md),
which runs the whole analysis. The building blocks below run the same
steps one at a time.

**Run**

- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md):
  the driver, from expression and orthologs to cliques.

- [`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md):
  write the result as tab-separated files.

**Inputs**

- [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md):
  read an ortholog group file.

- [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md):
  pair orthologs for every species pair.

- [`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.md):
  merge correlated paralogs.

- [`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md):
  import a network built outside rcomplex.

- [`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md):
  derive clades from a species tree.

**Harmonise**

- [`split_layers()`](https://paliocha.github.io/rcomplex/reference/split_layers.md):
  split expression into wiring and deployment layers.

**Networks**

- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md):
  build a mutual-rank co-expression network.

**Co-expressologs**

- [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md):
  test every ortholog pair across species pairs.

- [`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md):
  repeat the test at several density thresholds.

**Cliques**

- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md):
  every maximal clique of the gene graph.

- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md):
  conservation tiers for gene cliques.

- [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md):
  cliques of the species graph.

- [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md),
  [`clique_threshold_sweep()`](https://paliocha.github.io/rcomplex/reference/clique_threshold_sweep.md):
  robustness checks.

- [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md):
  conservation pattern per hog.

**Modules**

- [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md):
  find modules with Leiden.

- [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md):
  import modules from any gene partition.

- [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md):
  pick one paralog copy per gene.

- [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md),
  [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.md):
  test and call whether modules are preserved.

- [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md):
  match modules across species.

- [`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md),
  [`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md):
  hub genes.

- [`get_coexpressed_hogs()`](https://paliocha.github.io/rcomplex/reference/get_coexpressed_hogs.md):
  query the partners of a candidate hog.

**Traits**

- [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md):
  preservation over many species pairs.

- [`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md):
  relabelling test on all species pairs.

**Nulls and diagnostics**

- [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md):
  shuffled-partner null network.

- [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md):
  degree-preserving edge-swap null.

- [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md):
  how many distinct p-values a test can give.

## See also

Useful links:

- <https://paliocha.github.io/rcomplex/>

- <https://github.com/paliocha/rcomplex>

- Report bugs at <https://github.com/paliocha/rcomplex/issues>

## Author

**Maintainer**: Martin Paliocha <martin.paliocha@nmbu.no>

Authors:

- Martin Paliocha <martin.paliocha@nmbu.no>
