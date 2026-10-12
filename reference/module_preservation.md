# Test whether co-expression modules are preserved across species

Projects the modules of a reference species onto a test species through
an ortholog map. Tests, per module, whether its topology survives in the
test species network.

## Usage

``` r
module_preservation(
  modules_ref,
  net_ref,
  net_test,
  orthologs,
  edges = NULL,
  species_ref = NULL,
  species_test = NULL,
  n_perm = 10000L,
  calibrate = c("mixture", "none"),
  n_cores = 1L,
  seed = NULL
)
```

## Arguments

- modules_ref:

  Module assignment for the reference species, from
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  or
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md).

- net_ref, net_test:

  Network objects from
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  for the reference and test species.

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog`.

- edges:

  Optional
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  result, passed to
  [`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md)
  for paralog resolution.

- species_ref, species_test:

  Species labels, required only when `edges` is supplied.

- n_perm:

  Number of permutations (default 10000). The smallest attainable
  p-value is `1 / (n_perm + 1)`, so this sets the floor on how
  significant any module can be: at 1000 permutations every strongly
  preserved module ties at `p = 0.000999` and cannot be ranked by
  q-value. The permutation loop costs O(edges) per iteration – roughly
  0.7 s per 1000 permutations on a 6000-gene network at density 0.03 –
  so the default buys a floor of 1e-4 cheaply. Raise it further when
  many modules sit at the floor.

- calibrate:

  How the combined p-value is calibrated before FDR correction.
  `"mixture"` (default) recalibrates `pmax` toward the empirical joint
  null of the two statistics; `"none"` uses `pmax` unchanged.

  `pmax` is a valid p-value for the intersection-union null, but it is
  calibrated against a bound rather than the actual joint null and runs
  roughly `1/t` conservative when the two statistics are
  near-independent, as they are here – measured on this engine, a
  realised false discovery rate of 1.2e-4 against a nominal 0.05, with
  the smallest attainable q-value 0.10. The mixture blends `pmax` with
  the permutation joint null in proportion to the estimated fraction of
  modules null on *both* statistics, which is a super-uniform bound for
  any dependence structure. The rejection region stays
  `max(p1, p2) <= c`, so a module still cannot be called preserved on
  density alone.

  One caveat, measured rather than argued: when neither statistic shows
  signal anywhere the estimated both-null fraction approaches 1 and the
  procedure calibrates against the intersection null for that contrast.
  Realised FDR there was 0.031 against a nominal 0.05.

- n_cores:

  Number of OpenMP threads (default 1).

- seed:

  Optional RNG seed. Results are independent of `n_cores`. `NULL`
  (default) draws from the ambient stream and leaves it advanced; a seed
  draws from a private stream and restores the caller's on exit, the
  package-wide contract described under
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

## Value

A list with components:

- preservation:

  One row per tested module: `module`, `size`, `size_mapped`, the two
  headline statistics, their permutation p-values, the combined
  `p_value` (raw `pmax`), `p_calibrated` (the mixture recalibration
  described under `calibrate`), `q_value` (Benjamini-Hochberg on
  `p_calibrated`, NOT on `p_value`), `evalue` (`p_calibrated` times the
  number of modules tested: the expected count of modules this preserved
  by chance), `Z.avg.weight`, `Z.cor.degree`, `Zsummary`,
  `Zsummary_null_sd`, `Zsummary_std`, and `medianRank` – the mean of the
  `avg.weight` and `cor.degree` ranks across the tested modules, 1 =
  strongest. No permutation moments enter it, so it avoids the null-sd
  module-size dependence of `Zsummary` – it is a rank *within one run*
  and so is not comparable between runs that tested different numbers of
  modules.
  [`classify_preservation()`](https://paliocha.github.io/rcomplex/reference/classify_preservation.md)
  carries it through.

- observed:

  All six statistics per module, with permutation means and standard
  deviations, plus `n_perm.<stat>` – the number of permutations each
  statistic was actually scored over, which can be fewer than `n_perm`
  when a statistic was undefined – and `n_joint`, the number scorable
  for both headline statistics at once.

- coverage:

  One row per reference module – `size`, `size_mapped`, `tested`, and
  the `reason` it was not – so the tested set reconciles against the
  partition. Modules below 10 mapped genes, and modules whose genes
  never reach the test species, leave the analysis entirely; without
  this the preservation table looks like a complete accounting when it
  is not.

- projection:

  One row per test-species gene that received a module label: the
  reference gene it came from, the module, and which resolution layer
  chose the pair.

- map:

  The ortholog map used.

- params:

  Call parameters, including `calibrate`, the estimated both-null
  fraction `w00`, and `n_joint` per module.

## Details

Unlike gene-overlap tests, this test reports a module as diverged when
it keeps its gene membership but loses its internal wiring.

## Statistics

Two statistics carry the call, the pair NetRep computes when only an
adjacency matrix is available:

- avg.weight:

  `sum(kIM) / (m^2 - m)`, the module's mean within-module adjacency
  ("module density"). On a binary network this is the proportion of
  realised to possible within-module edges.

- cor.degree:

  Pearson correlation between each gene's intramodular connectivity in
  the reference network and in the test network. Measures whether hub
  identity is conserved.

`meanClusterCoeff` and `meanMAR` are reported in `observed` as
diagnostics but take no part in the call; on a hard-thresholded network
the surviving edge weights are nearly constant, which leaves both
statistics with almost no dynamic range.

## Null model

Gene identities are shuffled while edges are held constant, and each
module is handed a contiguous block of the shuffled genes of its own
size. Only ortholog-mappable test-species genes enter the shuffle,
matching NetRep's `"overlap"` null model. P-values are one-sided
(`(exceedances + 1) / (permutations + 1)`) and combined across the two
statistics with `pmax`, so a module is called preserved only when both
are significant – the same reciprocal criterion as
`pval_combine = "max"` elsewhere in the package.

## Paralog resolution

Multi-copy HOGs are reduced toward one counterpart per gene by
[`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md).
Resolution never changes which genes are *mappable* – that is the
invariant
[`resolve_ortholog_map()`](https://paliocha.github.io/rcomplex/reference/resolve_ortholog_map.md)
enforces. It can still change which genes end up *tested*: an unresolved
gene whose candidate module labels tie is dropped by the majority vote,
and resolving its copy rescues it. The tested set therefore does depend
on the resolution.

## References

Ritchie, S. C. et al. (2016). A scalable permutation approach reveals
replication and preservation patterns of network modules in large
datasets. *Cell Systems*, 3(1), 71-82.
[doi:10.1016/j.cels.2016.06.012](https://doi.org/10.1016/j.cels.2016.06.012)

Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. (2011). Is my
network module preserved and reproducible? *PLoS Computational Biology*,
7(1), e1001057.
[doi:10.1371/journal.pcbi.1001057](https://doi.org/10.1371/journal.pcbi.1001057)

## Examples

``` r
if (FALSE) { # \dontrun{
pres <- module_preservation(mods_A, net_A, net_B, ortho,
  edges = rcx$edges,
  species_ref = "SP_A", species_test = "SP_B"
)
classify_preservation(pres)
} # }
```
