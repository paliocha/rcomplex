# Leave-k-out jackknife structural stability for cliques

Leave-k-out jackknife structural stability for cliques

## Usage

``` r
find_cliques_stability_cpp(
  edge_hog,
  edge_g1,
  edge_g2,
  edge_species1,
  edge_species2,
  edge_qval,
  edge_eff,
  n_all_species,
  n_hogs,
  n_genes,
  is_target,
  full_cliques,
  max_k = 3L,
  max_genes_per_sp = 10L,
  jaccard_threshold = 0.8,
  n_cores = 1L,
  w_q = 1,
  w_eff = 0
)
```

## Arguments

- edge_hog, edge_g1, edge_g2, edge_species1, edge_species2:

  0-based edge vectors

- edge_qval, edge_eff:

  per-edge q-value and effect size

- n_all_species:

  total species in the analysis universe

- n_hogs, n_genes:

  total HOGs and genes

- is_target:

  1/0 per species (length n_all_species): target membership

- full_cliques:

  output of find_cliques_cpp on the full dataset

- max_k, max_genes_per_sp, jaccard_threshold, n_cores:

  tuning parameters

- w_q:

  Weight for mean q-value in composite cost (default 1.0)

- w_eff:

  Weight for mean effect size in composite cost (default 0.0)

## Value

List with stability, clique_disruption, stability_class, novel_cliques
