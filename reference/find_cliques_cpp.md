# Find co-expression cliques using two-level decomposition

Find co-expression cliques using two-level decomposition

## Usage

``` r
find_cliques_cpp(
  edge_hog,
  edge_g1,
  edge_g2,
  edge_species1,
  edge_species2,
  edge_qval,
  edge_eff,
  n_target_species,
  min_species,
  n_hogs,
  n_genes,
  max_missing_edges = 0L,
  w_q = 1,
  w_eff = 0,
  max_iterations = 2000000L
)
```

## Arguments

- edge_hog:

  0-based HOG index per edge

- edge_g1:

  0-based gene index for gene 1

- edge_g2:

  0-based gene index for gene 2

- edge_species1:

  0-based species index for gene 1

- edge_species2:

  0-based species index for gene 2

- edge_qval:

  q-value per edge

- edge_eff:

  Effect size per edge

- n_target_species:

  Number of target species

- min_species:

  Minimum species for a valid clique

- n_hogs:

  Total number of distinct HOGs

- n_genes:

  Total number of distinct genes

- max_missing_edges:

  Maximum missing edges tolerated (default 0)

- w_q:

  Weight for mean q-value in composite cost (default 1.0)

- w_eff:

  Weight for mean effect size in composite cost (default 0.0)

- max_iterations:

  Step budget of each best-gene search

## Value

List with: hog_idx (0-based), genes (matrix, 0-based or NA), n_species,
mean_q, max_q, mean_effect_size, n_edges, n_missing, complete (the
search finished within its step budget)
