# Package index

## Run

One call from expression to cliques.

- [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
  : Compare co-expression networks across species in one call
- [`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
  : Write an rcomplex result as tab-separated files

## Inputs

Ortholog groups, clades and networks built elsewhere.

- [`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md)
  : Read an ortholog group file
- [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)
  : Pair orthologs for every species pair
- [`reduce_orthogroups()`](https://paliocha.github.io/rcomplex/reference/reduce_orthogroups.md)
  : Reduce orthogroups by merging correlated paralogs
- [`as_network()`](https://paliocha.github.io/rcomplex/reference/as_network.md)
  : Import a network built outside rcomplex
- [`clades_from_tree()`](https://paliocha.github.io/rcomplex/reference/clades_from_tree.md)
  : Clades from a species tree

## Harmonise

Remove a designed factor before correlation.

- [`split_layers()`](https://paliocha.github.io/rcomplex/reference/split_layers.md)
  : Split expression into wiring and deployment layers

## Networks

- [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  : Compute co-expression network

## Co-expressologs

Pairwise alignment of ortholog pairs.

- [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  [`run_pairwise_comparisons()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  : Find co-expressologs across all species pairs
- [`density_sweep()`](https://paliocha.github.io/rcomplex/reference/density_sweep.md)
  : Sweep density thresholds for robustness assessment

## Cliques

Multiple alignment across species.

- [`find_cliques()`](https://paliocha.github.io/rcomplex/reference/find_cliques.md)
  : Find co-expression cliques using C++ two-level decomposition
- [`clique_stability()`](https://paliocha.github.io/rcomplex/reference/clique_stability.md)
  : Leave-k-out jackknife structural stability for cliques
- [`clique_threshold_sweep()`](https://paliocha.github.io/rcomplex/reference/clique_threshold_sweep.md)
  : Structural survival of cliques across stricter density thresholds
- [`classify_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_cliques.md)
  : Classify HOGs by clique conservation pattern
- [`gene_clique_graph()`](https://paliocha.github.io/rcomplex/reference/gene_clique_graph.md)
  : Maximal cliques of the per-orthogroup gene graph
- [`classify_gene_cliques()`](https://paliocha.github.io/rcomplex/reference/classify_gene_cliques.md)
  : Classify gene-graph cliques into conservation tiers

## Modules

Aligned blocks and their hubs.

- [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  : Detect co-expression modules in a network
- [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
  : Build a module assignment from any gene partition
- [`joint_modules()`](https://paliocha.github.io/rcomplex/reference/joint_modules.md)
  : Line up per-species modules across species
- [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md)
  : Resolve paralog copies in an ortholog map
- [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
  : Test whether co-expression modules are preserved across species
- [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.md)
  : Classify modules as conserved, moderately preserved, or diverged
- [`module_correspondence()`](https://paliocha.github.io/rcomplex/reference/module_correspondence.md)
  : Match modules across species by ortholog overlap
- [`identify_module_hubs()`](https://paliocha.github.io/rcomplex/reference/identify_module_hubs.md)
  : Identify hub genes within co-expression modules
- [`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md)
  : Classify hub gene conservation across species and clades
- [`get_coexpressed_hogs()`](https://paliocha.github.io/rcomplex/reference/get_coexpressed_hogs.md)
  : Query co-expression partners of a candidate HOG across species

## Traits

Is divergence associated with a trait?

- [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
  : Run module preservation across many species pairs
- [`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
  : Trait relabelling test on the all-pairs preservation matrix

## Nulls and diagnostics

- [`null_network()`](https://paliocha.github.io/rcomplex/reference/null_network.md)
  : Shuffled-partner null network
- [`coexpressolog_null()`](https://paliocha.github.io/rcomplex/reference/coexpressolog_null.md)
  : Degree-preserving edge-swap null for co-expressolog statistics
- [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)
  : How much resolution is left in a set of p-values
