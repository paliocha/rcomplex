# Detect co-expression modules in a network

Applies the Leiden algorithm (Traag *et al.*, 2019) to a thresholded
co-expression network. Optimises modularity or CPM and guarantees
well-connected communities. The `resolution` argument sets module size.

## Usage

``` r
detect_modules(net, ...)

# Default S3 method
detect_modules(
  net,
  resolution = NULL,
  objective_function = c("modularity", "CPM"),
  seed = NULL,
  n_cores = 1L,
  max_consensus_iter = 10L,
  ...
)
```

## Arguments

- net:

  Network object from
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md).

- ...:

  Additional arguments passed to the default method.

- resolution:

  Resolution parameter for Leiden. `NULL` (default) runs a single
  resolution: for CPM the edge density of the thresholded graph,
  `ecount / choose(vcount, 2)`, so modules must be denser than the
  network average; for modularity 1. The resolved value is returned in
  `params$resolution`. Pass a numeric vector (e.g.,
  `seq(0.5, 2.0, by = 0.5)` for modularity) to run at multiple
  resolutions and produce a consensus partition via co-classification
  (Lancichinetti & Fortunato, 2012); consensus mode has no default.
  Under CPM it is read on the 0–1 weight scale (edge weights divided by
  their maximum, since 0.3.2): a density each module must exceed, so
  values below 1 are the useful range (at 1 no edge is attractive and
  genes stay singletons). Earlier versions read it on raw mutual-rank
  weights with a default of 1, where a resolution of 1 or more returned
  one module containing every gene. On 20,000-gene Pooideae and wood
  networks the density default gave 13–63 modules per species: a few
  large ones (over 1,000 genes) that replicated between sample halves
  and were conserved across species as well as modularity's, plus a tail
  of small communities (under 300 genes) that were not conserved. Drop
  the tail by size before reading modules as units.

- objective_function:

  Leiden objective: `"modularity"` (default since 0.3.2; CPM before) or
  `"CPM"`. The two replicated and conserved large modules equally on
  20,000-gene data, and modularity has no tail of small non-conserved
  communities; see the section on reproducibility. CPM sees edge weights
  rescaled to 0–1; modularity is scale-invariant and sees them
  unchanged. At small sample sizes `"modularity"` replicated better
  across independent samples but also finds modules in noise; see the
  section on reproducibility.

- seed:

  Random seed for reproducibility (default `NULL`).

  With a seed, the result is reproducible and identical at any
  `n_cores`: every parallel task derives its own RNG stream from the
  seed and its task index, so the answer does not depend on how the work
  was distributed. Two limits are worth knowing. Reproducibility holds
  for one machine and one igraph build – the partition is a function of
  what
  [`igraph::cluster_leiden()`](https://r.igraph.org/reference/cluster_leiden.html)
  draws, so an igraph upgrade may move it. And the answer depends on
  `RNGkind` as well as on `seed`: a session that has set
  `RNGkind("L'Ecuyer-CMRG")` gets a different, still reproducible,
  partition. `RNGkind` is deliberately not pinned inside the function.

  With an explicit seed the caller's stream is restored on exit, so the
  call is invisible to anything drawn afterwards: the clustering backend
  draws from a private stream started at `seed`. With `seed = NULL` the
  stream advances by exactly what was taken from it – the one draw used
  to pick a root in consensus mode, whatever the backend consumed in
  single-resolution mode – so consecutive unseeded calls still differ.

- n_cores:

  Number of parallel cores (default 1). Used for `mclapply` Leiden
  sweeps on Unix and OpenMP edge scans in C++. Uses fork-based
  parallelism; avoid combining with active CUDA contexts in the same
  session.

- max_consensus_iter:

  Maximum number of consensus iterations for consensus mode. Default 10.
  Iteration stops when the sweep reproduces its own input (a fixed
  point) or when all resolutions agree. On the eight Pooideae networks
  that stopping rule was measured on it fired at 5–12 iterations, so the
  default cap of 10 cuts two of those eight off before it fires: at the
  default the loop may return a partition it had not finished settling.
  Raise it (20 is above every measured stop) when the returned
  `params$n_consensus_iterations` equals `max_consensus_iter` – a
  truncated run always does, though a run that settles on the last
  allowed iteration does too.

## Value

A list with components:

- modules:

  Named integer vector of module assignments (gene -\> module ID)

- module_genes:

  Named list: module ID -\> character vector of gene names

- n_modules:

  Number of modules detected

- modularity:

  Modularity score of the partition

- graph:

  The igraph graph object used for community detection

- method:

  Method used

- params:

  List of parameters used

When `resolution` is a vector, the output also includes
`resolution_scan`, a data frame with columns `resolution`, `n_modules`,
`modularity`, `ari_next` (Adjusted Rand Index with the next resolution;
NA for the last), and `expected_coclassification` (per-resolution
expected scalar). The `params` list includes `n_consensus_iterations`
(number of iterations until convergence).

## Details

The adaptive consensus path uses sparse co-classification restricted to
the original network's edge set, reducing memory from O(N^2) to
O(\|E\|).

## Reproducibility at small sample sizes

Modules can be stable across seeds and still not replicate across
independent samples, and a correlation network from few samples has
cluster structure even on noise. Measured on Pooideae leaf networks
(2,000 genes, 20 samples) with `resolution = c(0.5, 1, 2)`; the scripts
and tables are in the package's source repository under `dev/bench/`:

- *Replication* (samples split into two halves of 10, two replicates per
  time point each): with `objective_function = "CPM"`, whether modules
  exist at all flipped between the halves – in 4 of 8 species one half
  gave 6–10 modules and the other 1, and two species gave 1 in both
  (median split-half adjusted Rand index 0). With `"modularity"` the
  modules replicated partially (median ARI 0.36, range 0.05–0.44; about
  0 on shuffled expression).

- *Noise*: on random and per-gene shuffled expression (18 data sets, all
  20 samples), CPM returned one module, but on 10-sample halves it gave
  10 modules on shuffled expression in 3 of 8 species. Modularity
  returned 8–14 modules on noise.

These CPM measurements predate 0.3.2 and ran on raw mutual-rank edge
weights, of the order of the gene count, so a resolution of 0.5–2 was
effectively 0 and CPM collapsed to one module whenever the graph was
connected. Since 0.3.2 CPM weights are divided by the maximum edge
weight, so `resolution` is a density on the 0–1 weight scale, and the
default resolution is the network's edge density; the CPM figures above
describe the old scale, not the current one.

The default objective is modularity since 0.3.2 (CPM before). On
20,000-gene Pooideae leaf and wood networks, CPM at the density default
and modularity replicated equally between sample halves and conserved
their large modules equally across species; CPM only added a tail of
small communities that were not conserved. The 2,000-gene `dev/bench/`
numbers above predate the CPM scale fix.

Treat modules as units only after they replicate on independent samples,
well above the same split on shuffled expression; `"modularity"`
replicated better here but needs that check most. To test gene sets from
other sources (pathways, regulons, another tool) use
[`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)
with
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md).

## References

Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to
Leiden: guaranteeing well-connected communities. *Scientific Reports*,
9, 5233.
[doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

Lancichinetti, A. & Fortunato, S. (2012). Consensus clustering in
complex networks. *Scientific Reports*, 2, 336.
[doi:10.1038/srep00336](https://doi.org/10.1038/srep00336)

Jeub, L. G. S., Sporns, O. & Fortunato, S. (2018). Multiresolution
consensus clustering in networks. *Scientific Reports*, 8, 3259.
[doi:10.1038/s41598-018-21352-7](https://doi.org/10.1038/s41598-018-21352-7)

Senbabaoglu, Y. et al. (2014). Critical limitations of consensus
clustering in class discovery. *Scientific Reports*, 4, 6207.
[doi:10.1038/srep06207](https://doi.org/10.1038/srep06207)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single resolution (modularity at resolution 1)
mods <- detect_modules(net)
table(mods$modules) # module sizes

# Multi-resolution consensus (Jeub et al. 2018)
mods_consensus <- detect_modules(net,
  resolution = c(0.5, 1.0, 2.0), objective_function = "modularity",
  n_cores = 4L
)
} # }
```
