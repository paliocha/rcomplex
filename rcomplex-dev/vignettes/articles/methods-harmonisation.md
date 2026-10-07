## Harmonising heterogeneous designs

Species rarely share an expression design. A network pooled over leaf and
root samples is a tissue network, and a ten-tissue compendium aligned
against a three-tissue one aligns tissue coverage, not regulation.
rcomplex has one tool per regime.

**A designed experiment with replicates: `split_layers()`.** One OLS
projection of every gene on a per-sample block factor (time point, tree,
zone) gives the block means (the deployment layer) and the residuals
(the wiring layer). Build the network from the wiring layer and pair it
with `null_network(..., block = block)`. Each block level needs at least
two samples. `r2` is the share of each gene's sum of squares the block
explains, the unadjusted variance fraction of Breschi et al. (2016).

No latent factors are removed, on purpose. Only the designed block is
projected out. Removing principal components or other estimated factors
did not improve co-expression network accuracy over unadjusted data (Cote
et al. 2022); principal-component removal lowers false positives without
improving false negatives and is ill-advised when a designed axis exists
(Parsana et al. 2019).

**A compendium: `compute_network(partition = )`.** This is the partition
aggregation of TEA-GCN (Lim et al. 2026) on a partition the user brings
(tissue, study, or a k-means of samples in PCA space); rcomplex does not
cluster samples. The correlation is computed within each level of the
partition, negative correlations are set to zero, and the mean over
levels is normalised by mutual rank (or CLR) and thresholded by density
as usual. Levels with fewer than 5 samples are dropped with a message,
and at least two levels must remain. TEA-GCN's global z-score is not
needed: mutual rank and a density threshold already put networks of
different species on one scale. `null_network()` rebuilds a partitioned
network with the same partition. The partition needs the dense build
(`block_size = NULL`).

### References

Breschi A, Djebali S, Gillis J, et al. (2016). Gene-specific patterns of
expression variation across organs and species. *Genome Biology* 17:151.
doi:10.1186/s13059-016-1008-y

Cote AC, Young HE, Huckins LM (2022). Comparison of confound adjustment
methods in the construction of gene co-expression networks. *Genome
Biology* 23:44. doi:10.1186/s13059-022-02606-0

Lim et al. (2026). TEA-GCN. *Nature Communications* 17:5906.
doi:10.1038/s41467-026-72380-1

Parsana P, Ruberman C, Jaffe AE, et al. (2019). Addressing confounding
artifacts in reconstruction of gene co-expression networks. *Genome
Biology* 20:94. doi:10.1186/s13059-019-1700-9
