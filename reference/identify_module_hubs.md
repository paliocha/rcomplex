# Identify hub genes within co-expression modules

Computes within-module centrality for each gene. Flags the top-ranked
genes as hubs. Optionally maps genes to hogs for
[`classify_hub_conservation()`](https://paliocha.github.io/rcomplex/reference/classify_hub_conservation.md).

## Usage

``` r
identify_module_hubs(modules, ...)

# Default S3 method
identify_module_hubs(modules, net, orthologs = NULL, ...)
```

## Arguments

- modules:

  Output of
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

- ...:

  Additional arguments passed to the default method.

- net:

  Output of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).

- orthologs:

  Optional data frame from
  [`prepare_orthologs()`](https://paliocha.github.io/rcomplex/reference/prepare_orthologs.md)
  with columns `gene1`, `gene2`, `hog`. The function auto-detects which
  column matches the gene names in `modules`. If `NULL`, the `hog`
  column in the result is all `NA`.

## Value

A data frame with one row per gene, ordered by module then rank:

- gene:

  Gene identifier

- module:

  Module ID (integer)

- degree:

  Within-module weighted degree
  ([`igraph::strength`](https://r.igraph.org/reference/strength.html))

- betweenness:

  Within-module betweenness centrality

- eigenvector:

  Within-module eigenvector centrality

- mean_edge_weight:

  Mean weight of edges to other module members

- global_degree:

  Weighted degree in the full (thresholded) network

- rank:

  Rank within module by primary centrality (1 = highest; ties use
  `"min"`)

- is_hub:

  `TRUE` if the gene is in the top slice after the tie-breaking cascade

- hog:

  HOG identifier (`NA` if `orthologs` not provided or gene not in the
  ortholog table)

## Details

A hub is in the top 10% by weighted degree, with at least one hub per
module of three or more genes.

## Centrality measures

Centrality is computed on the **within-module subgraph** (edges between
genes in the same module only):

- degree:

  Weighted degree
  ([`igraph::strength`](https://r.igraph.org/reference/strength.html)):
  sum of edge weights to other genes in the same module.

- betweenness:

  Shortest-path betweenness using inverse edge weights as distances.
  Identifies genes that bridge sub-clusters within a module.

- eigenvector:

  Eigenvector centrality
  ([`igraph::eigen_centrality`](https://r.igraph.org/reference/eigen_centrality.html)):
  high for genes connected to other high-centrality genes.

## Tie-breaking cascade

When genes share the same primary centrality score, hub selection uses a
cascade:

1.  Within-module weighted degree

2.  Global weighted degree across the full network

3.  Within-module betweenness

4.  Mean within-module edge weight (strength / degree)

## Examples

``` r
if (FALSE) { # \dontrun{
hubs <- identify_module_hubs(modules, net, orthologs)
hubs[hubs$is_hub, ]
} # }
```
