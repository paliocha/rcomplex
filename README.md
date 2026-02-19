# rcomplex

Comparative co-expression network analysis across species.

Compares gene co-expression networks across plant species by mapping
orthologous genes (via PLAZA ortholog groups), building co-expression networks
independently per species, then testing whether network neighborhoods are
significantly preserved using hypergeometric tests with q-value correction
and effect sizes.

Based on [Netotea *et al.*, 2014](https://doi.org/10.1186/1471-2164-15-106).

## Installation

```r
# Install from GitHub
devtools::install_github("paliocha/rcomplex")
```

## Workflow

```r
library(rcomplex)

# 1. Parse ortholog groups
orthologs <- parse_orthologs("orthogroups.txt", "Species1", "Species2")

# 2. Build co-expression networks
net1 <- compute_network(expr1, norm_method = "MR", density = 0.03)
net2 <- compute_network(expr2, norm_method = "MR", density = 0.03)

# 3. Compare neighborhoods (pair-level hypergeometric tests)
comparison <- compare_neighborhoods(net1, net2, orthologs)

# 4a. Pair-level summary with q-values (Storey & Tibshirani, 2003)
summary <- summarize_comparison(comparison)

# 4b. HOG-level permutation test (recommended for multi-copy gene families)
#     Uses discrete q-values (Liang, 2016) for multiple testing correction
hog_results <- permutation_hog_test(net1, net2, comparison, n_cores = 4L)
```

## HOG-level permutation test

### The problem with Fisher's method

A naive approach to HOG-level testing would combine pair-level p-values using
Fisher's method. However, Fisher's method assumes independent tests. Within a
HOG, pair-level
hypergeometric tests are correlated because:

- **Shared neighborhoods**: Genes in the same HOG often share co-expression
  neighbors (especially after whole-genome duplication), so their neighborhood
  overlap tests draw from overlapping gene sets.
- **Shared ortholog mapping**: The same ortholog table mediates the
  neighborhood comparison for all pairs in a HOG. Two pairs sharing an sp2 gene
  use identical ortholog-mapped neighbor sets for one direction.

This non-independence inflates Fisher's statistic and produces anti-conservative
(overly optimistic) p-values. The severity depends on HOG size: 1:1 HOGs have
no issue, but multi-copy HOGs (e.g., 5x5 = 25 pairs) can be strongly affected.

### The permutation approach

`permutation_hog_test()` replaces Fisher's method with an exact permutation
test that is valid regardless of the dependency structure.

**Null hypothesis (gene-identity permutation)**: The specific gene identities in
a HOG carry no information about co-expression conservation. Replacing the HOG's
M species-1 genes and N species-2 genes with randomly chosen genes from the
same networks would yield equally large neighborhood overlap.

**Test statistic**: Sum of fold-enrichments across all M x N pair x direction
combinations:

```
T = sum_{i,j} [ x1_ij / E1_ij + x2_ij / E2_ij ]
```

where `x` is the observed overlap and `E = m * k / N_total` is the
hypergeometric expectation. This avoids expensive `phyper()` calls during
permutations while remaining monotonically related to the conservation signal.

**Precomputed reachable sets**: For each gene in each network, the set of
cross-species genes reachable through network neighbors and the ortholog mapping
is precomputed once. This reduces each permutation to set intersections rather
than full neighborhood recomputations.

**Intersection modes**:
- **Bit-vector** (when max(n1, n2) <= 100,000): Each neighbor/reachable set is
  stored as a bit-vector. Intersection is computed via bitwise AND + popcount,
  achieving ~0.15 us per intersection on modern CPUs.
- **Flag-vector** (when max(n1, n2) > 100,000): Sparse set representation with
  a flag array for intersection counting. Avoids the O(n^2/8) memory cost of
  bit-vectors for very large networks.

### Adaptive stopping (Besag & Clifford, 1991)

Rather than running a fixed number of permutations for every HOG, the
implementation uses the sequential Monte Carlo stopping rule of Besag &
Clifford (1991). Permutations continue until `min_exceedances` (default: 50)
permutation statistics exceed (or fall below, for divergence) the observed
statistic, then:

```
p = (n_exceed + 1) / (n_perm + 1)
```

This provides:
- **Efficiency**: Clearly non-significant HOGs stop after ~50 permutations
  (exceedances accumulate immediately). Clearly significant HOGs also stop
  relatively quickly. Only borderline HOGs use the full permutation budget.
- **Precision where it matters**: The p-value estimate has relative precision
  ~1/sqrt(min_exceedances), so with 50 exceedances the coefficient of variation
  is ~14%, sufficient for multiple testing correction.

### Multiple testing correction: discrete q-values (Liang, 2016)

Besag-Clifford adaptive stopping produces p-values that are valid
(super-uniform under the null) but **not** uniformly distributed. The negative
binomial stopping rule concentrates null p-values in the 0.1--0.5 range and
depletes the right tail (p > 0.9). This creates problems for standard
correction methods:

- **BH (Benjamini & Hochberg, 1995)**: Valid but conservative. BH assumes
  nothing about the null distribution beyond super-uniformity, so it correctly
  controls FDR. However, it implicitly assumes pi0 = 1 (all nulls), forfeiting
  the power gain from estimating the true null proportion.
- **Storey q-values (Storey & Tibshirani, 2003)**: Invalid. The pi0 estimator
  relies on null p-values being uniformly distributed in (0, 1) so that the
  right tail (e.g. p > 0.9) is flat. With BC p-values, the depleted right tail
  causes pi0 to be severely underestimated (e.g. median pi0 ~ 0.18 instead of
  ~ 0.4), inflating the number of discoveries by 5--6x.

`permutation_hog_test()` uses the discrete q-value method of Liang (2016) via
`DiscreteQvalue::DQ()`. This method estimates pi0 using the known discrete
support structure of the p-values rather than relying on the right-tail
density. On simulated pure-null BC p-values, Liang correctly returns
pi0 = 1.000 across all seeds, while Storey's estimator fails entirely
(spline returns NA or pi0 ~ 0.2).

The BC p-value support is constructed from the stopping parameters:

- **Early-stopped HOGs** (n_exceed = k): p = (k+1)/(n+1) for
  n = k, k+1, ..., M
- **Max-permutation HOGs** (n_perm = M): p = (j+1)/(M+1) for
  j = 0, 1, ..., k-1

where k = `min_exceedances` and M = `max_permutations`.

### Performance: column-major memory access

The C++ backend (Armadillo, column-major storage) is optimized for sequential
memory access. Since correlation and network matrices are symmetric, column
access (`sim.col(i)`, `net.colptr(i)`) returns the same values as row access
but with sequential rather than strided reads. This applies to:

- **MR normalization** (`mutual_rank_transform_cached_cpp`): ranks are computed
  from `sim.col(i)` and stored by column, avoiding O(n) cache misses per gene.
- **Neighbor list extraction** (`hog_permutation_test_cpp`): `net.colptr(i)`
  gives a sequential scan of gene i's co-expression values.

On a 22K-gene network pair, this gives a ~3x speedup over naive row access
(380s to 126s per species comparison with 8 OpenMP threads).

## References

- Netotea, S. *et al.* (2014). ComPlEx: conservation and divergence of
  co-expression networks in *A. thaliana*, *Populus* and *O. sativa*.
  *BMC Genomics*, 15, 106. [doi:10.1186/1471-2164-15-106](https://doi.org/10.1186/1471-2164-15-106)
- Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
  *Biometrika*, 78(2), 301--304. [doi:10.1093/biomet/78.2.301](https://doi.org/10.1093/biomet/78.2.301)
- Storey, J. D. & Tibshirani, R. (2003). Statistical significance for
  genomewide studies. *PNAS*, 100(16), 9440--9445.
  [doi:10.1073/pnas.1530509100](https://doi.org/10.1073/pnas.1530509100)
- Liang, K. (2016). False discovery rate estimation for large-scale
  homogeneous discrete p-values. *Biometrics*, 72(2), 639--648.
  [doi:10.1111/biom.12429](https://doi.org/10.1111/biom.12429)

## License

MIT
