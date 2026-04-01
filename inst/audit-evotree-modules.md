# Module Detection Audit: EVOTREE vs rcomplex

*Date: 2026-04-01*

## 1. Executive Summary

EVOTREE and rcomplex take fundamentally different approaches to co-expression module detection:

- **EVOTREE** (`ModuleNetwork.Rmd`): A **hub-centric greedy coverage** algorithm that builds non-overlapping modules around highly connected genes with high expression variance. No optimization objective. Conservation is assessed by the proportion of module genes with significant orthologs.

- **rcomplex** (`R/modules.R`): **Optimization-based community detection** (Leiden, Infomap, or SBM) with optional multi-resolution consensus clustering (Jeub *et al.*, 2018). Conservation is assessed by hypergeometric or Jaccard permutation tests with formal multiple-testing correction.

**Conclusion**: rcomplex's approach is more principled at every level — module detection, resolution handling, and conservation assessment. No features from EVOTREE's module engine need to be incorporated. One optional enhancement is identified: a lightweight hub-gene identification utility for biological interpretability.

---

## 2. EVOTREE ModuleNetwork.Rmd — Plain-Language Walkthrough

Source: <https://github.com/ellendim/EVOTREE/blob/main/GRN/modules/ModuleNetwork.Rmd>

### 2.1 What ModuleNetwork.Rmd does, step by step

**Step 1 — Normalize expression.** Each gene's expression is Z-score standardized across samples. This puts all genes on a comparable scale so that expression variance becomes a meaningful selector for "interesting" genes.

**Step 2 — Build the co-expression network.** Pearson correlation is computed on the transposed expression matrix (samples x genes), then Mutual Rank (MR) normalization is applied: `MR = sqrt(R * t(R))` where `R` is the matrix of column-wise ranks. The diagonal is set to zero. This is the same MR normalization rcomplex uses.

**Step 3 — Threshold the network.** Edges are retained at a 3% density cutoff — only the top 3% of MR-ranked edges survive. Again, this is the same density-threshold approach rcomplex uses.

**Step 4 — Identify hub genes.** Two filters are applied simultaneously:
- **Expression variance > 3** (on the log-scale): selects genes whose expression changes substantially across developmental stages — likely regulators or effectors.
- **Node degree > 20** neighbors (in the thresholded network): selects genes that are highly connected.

Genes passing both filters are ranked by degree (highest first). These are the "hub centroids."

**Step 5 — Build modules by greedy coverage.** Starting from the top-ranked hub:
1. The hub and all its direct neighbors form a module.
2. These genes are removed from the candidate pool.
3. The next highest-ranked hub seeds the next module.
4. Repeat until all hubs are consumed.

This is strictly greedy: the first hub gets the "best" genes; later hubs get whatever remains. There is no reassignment, no objective function, and no iteration.

**Step 6 — Assign leftover genes.** Genes not captured by any hub's neighborhood are assigned to the module whose centroid (mean expression profile) they correlate most strongly with. This uses a relaxed secondary threshold.

**Step 7 — Assess conservation.** For each module, count the fraction of its genes that have at least one significant ortholog (p < 0.05 from the RComPlEx neighborhood comparison). A high fraction means the module is "conserved." There is no formal statistical test of module-level conservation — just a proportion.

### 2.2 Biological rationale

The algorithm assumes that co-expression modules are organized around transcription factor hubs. In the EVOTREE context (wood-forming tissues), this is biologically motivated: SCW (secondary cell wall) cellulose synthase genes cluster around NAC and MYB master regulators. The hub-centric design makes the biological narrative immediate — each module is named by its seed TF.

### 2.3 Limitations

| Limitation | Consequence |
|-----------|-------------|
| No objective function | Partition quality is never measured or optimized |
| Order-dependent | Changing hub ranking changes the entire partition |
| No resolution control | Module granularity is fixed by the variance and degree thresholds |
| No reassignment | A gene assigned to hub #3 might fit better in hub #7's module |
| No K=1 test | No way to know if the network has community structure at all |
| Proportion-based conservation | No formal test of whether module overlap exceeds chance |

---

## 3. EVOTREE Published Method — Manuscript Context

Source: Dimén *et al.*, "An evo-devo resource for wood comparative regulomics across dicot and conifer trees" (Research Square, rs-7656402/v1).

### 3.1 Study design

Six tree species spanning 250 Myr of evolution:
- **Dicots**: aspen (*Populus tremula*), birch (*Betula pendula*), cherry (*Prunus avium*)
- **Conifers**: Norway spruce (*Picea abies*), Scots pine (*Pinus sylvestris*), lodgepole pine (*Pinus contorta*)

522 total samples from wood-forming tissues (secondary phloem, vascular cambium, xylem). Three biological replicates per species. Expression normalized by VST (DESeq2).

### 3.2 Ortholog mapping

OrthoFinder v2.5.2 on longest protein-coding sequences from 27 plant species with a rooted TimeTree species tree. Both standard orthogroups (OGs) and phylogenetic hierarchical orthogroups (HOGs) are used. 5,564 expressed orthogroups across all 6 species.

### 3.3 Co-expression comparison (ComPLEx / RComPlEx)

This is what rcomplex already implements as `compare_neighborhoods()`:
- For each ortholog pair A–B, map co-expression neighbors of A to species B via orthologs.
- Test overlap with actual neighbors of B using a hypergeometric test.
- Do the same in reverse (B's neighbors mapped to species A).
- Both FDR-corrected p-values must be < 0.1 to classify the pair as a "co-expressolog."

### 3.4 Module detection in the manuscript

The manuscript describes an "iterative coverage algorithm applied to the most highly connected genes" — this is the `ModuleNetwork.Rmd` approach described in Section 2.

### 3.5 Clique detection

`igraph::max_cliques()` on cross-species ortholog networks where edges are co-expressolog relationships. Clique categories:
- **Complete conserved**: all 15 ortholog-pair combinations are co-expressologs (FDR < 0.1)
- **Partial conserved**: at least 11 of 15 pairs significant
- **Differentiated**: lineage-specific patterns (all 3 within-lineage pairs significant, at most 4 cross-lineage)
- **Lineage-specific**: 3-membered cliques from dicot-only or conifer-only orthogroups

Result: 2,145 unique cliques from 2,098 orthogroups (38% of expressed orthogroups).

### 3.6 GRN inference (orthogonal to modules)

The manuscript adds a regulatory layer not present in rcomplex:
- ATAC-seq for open chromatin in aspen and spruce
- 808 TF motifs from JASPAR and PlantTFDB mapped to open regions within 10 kb of target genes
- TF–module associations: motif enrichment (Fisher exact, FDR < 0.05) + expression correlation (Pearson > 0.7, FDR < 0.05)
- Three-layer regulatory hierarchy: NAC master regulators → R2R3-MYB TFs → secondary wall NACs/MYBs

This is an entirely separate analysis layer that consumes module/clique output as input. It is not a gap in rcomplex's module detection engine — it is a downstream application.

---

## 4. rcomplex Module Engine — Summary

### 4.1 `detect_modules()` — R/modules.R

Three community detection methods:

| Method | Algorithm | Resolution | Model selection |
|--------|-----------|-----------|-----------------|
| **Leiden** | Modularity or CPM optimization (Traag *et al.*, 2019) | Explicit parameter; scalar or vector | User-specified |
| **Infomap** | Random-walk compression (Rosvall & Bergstrom, 2008) | None (implicit) | Automatic |
| **SBM** | Gaussian stochastic block model, variational EM | None | ICL (automatic) |

**Multi-resolution consensus** (triggered when `resolution` is a vector):
1. Leiden sweep at K resolutions on the original network.
2. Sparse co-classification matrix (C++ `build_sparse_coclassification_cpp`, O(|E|) memory): C[i,j] = fraction of resolutions where genes i and j share a module.
3. Per-pair expected co-classification: E[i,j] = (1/K) * sum_k (s_m(i)/N * s_m(j)/N). This is NOT a scalar mean — it accounts for module-size variation per pair per resolution.
4. Excess = max(C - E, 0). Build consensus graph from excess weights.
5. K=1 null test: permutation of the spectral norm of the excess matrix (degree-preserving rewiring, `sparse_excess_spectral_norm_cpp`).
6. Iterate: Leiden sweep on consensus graph → recompute co-classification → check convergence (ARI > 0.999 across all resolution pairs). Typically 2–5 iterations.

### 4.2 `compare_modules()` — R/modules.R

Two methods for cross-species module comparison:

**Hypergeometric** (R-only, BLAS-accelerated):
- Binary membership matrices M1 (sp1 modules mapped to sp2 space via orthologs) and M2 (sp2 modules).
- Overlap matrix = `crossprod(M1, M2)` — one BLAS call for all pairs.
- Hypergeometric p-value per pair. Q-values via `qvalue::qvalue()` (Storey & Tibshirani, 2003).

**Jaccard permutation** (C++ `module_jaccard_permutation_cpp`):
- Null: Fisher-Yates shuffle of the ortholog mapping (sp2 side), preserving module and network structure.
- Observed Jaccard per module pair. Permuted Jaccard from each shuffle.
- Besag-Clifford adaptive stopping: stop testing a pair once it accumulates `min_exceedances` (default 50).
- P-value = (n_exceed + 1) / (n_perm + 1). Q-values via `DiscreteQvalue::DQ(method = "Liang")` — required because BC p-values have a depleted right tail that causes Storey's method to underestimate pi0.

### 4.3 `classify_modules()` — R/modules.R

Three-way classification based on best-match q-value and Jaccard:

| Classification | Criterion |
|---------------|-----------|
| **Conserved** | q < alpha AND Jaccard >= threshold |
| **Partially conserved** | q < alpha AND Jaccard < threshold |
| **Species-specific** | q >= alpha |

Default: alpha = 0.05, jaccard_threshold = 0.1.

---

## 5. Fundamental Differences

### Comparison table

| Aspect | EVOTREE ModuleNetwork | rcomplex |
|--------|----------------------|----------|
| **Detection philosophy** | Greedy hub-centric coverage | Optimization-based community detection |
| **Objective function** | None (deterministic greedy) | Modularity, CPM, compression, or ICL |
| **Resolution control** | None (fixed by variance/degree thresholds) | Explicit parameter + multi-resolution consensus |
| **Robustness** | Sensitive to hub ranking order | Consensus across resolutions with ARI convergence |
| **K=1 null test** | No | Spectral norm permutation |
| **Conservation test** | Proportion of significant orthologs | Hypergeometric or Jaccard permutation + q-values |
| **Multiple testing** | None at module level | Storey FDR or Liang discrete FDR |
| **Hub identification** | Central to algorithm | Emergent (not explicit) |
| **Leftover genes** | Assigned to nearest centroid | All connected genes assigned by the algorithm |
| **GRN inference** | ATAC + JASPAR + TF correlation | Out of scope |
| **Implementation** | Pure R notebook | R + C++/Armadillo (OpenMP, sparse) |

### 5.1 Hub-centric vs community detection

EVOTREE assumes modules are organized around transcription factor hubs — genes with high expression variance and high connectivity. This is a strong biological prior. When it holds (e.g., SCW regulons around NAC/MYB TFs), the modules are immediately interpretable.

Community detection makes no such assumption. Leiden, Infomap, and SBM discover modules from network topology alone. Modules without clear hubs (metabolic pathways, housekeeping gene clusters) are detected just as well. The biological narrative requires post-hoc enrichment analysis, but the partition itself is unbiased.

**Verdict**: Community detection is more general. Hub-centric is a useful biological lens but should not drive the partition.

### 5.2 Greedy coverage vs global optimization

EVOTREE's greedy assignment means the first hub captures the densest neighborhood. Subsequent hubs receive whatever genes remain. If two hubs share many neighbors, the first one "wins" and the second gets a depleted module. There is no mechanism to redistribute genes.

Leiden optimizes a global objective (modularity or CPM) and can move genes between modules to improve the overall partition. The consensus procedure further stabilizes this by averaging over multiple resolutions.

**Verdict**: Global optimization produces higher-quality partitions. Greedy coverage is fast but fragile.

### 5.3 Resolution sensitivity

EVOTREE has no resolution parameter. Module granularity is entirely determined by the variance > 3 and degree > 20 thresholds. Changing these thresholds changes the number and identity of hubs, which changes the entire partition. There is no way to systematically explore different scales.

rcomplex's multi-resolution consensus explicitly addresses this: sweep Leiden across a range of resolutions, compute co-classification across all of them, extract structure that is robust across scales. The iterative Jeub *et al.* (2018) procedure with per-pair null subtraction ensures that the consensus is not biased by module size.

**Verdict**: Multi-resolution consensus is the principled solution to the resolution limit problem.

### 5.4 Conservation assessment rigor

EVOTREE assesses module conservation by counting "what fraction of module genes have significant orthologs (p < 0.05)." This has several problems:
- It uses gene-level p-values, not a module-level test.
- The threshold (p < 0.05) is arbitrary and not adjusted for the number of modules tested.
- A module where 60% of genes are significant is called "conserved," but there is no test of whether 60% exceeds what would be expected by chance given module sizes and ortholog density.

rcomplex tests each module pair directly:
- Hypergeometric: does the overlap between sp1-module (mapped to sp2) and sp2-module exceed chance?
- Jaccard permutation: does the observed Jaccard exceed what you get by shuffling the ortholog mapping?
- Both methods apply q-value correction across all module pairs.

**Verdict**: rcomplex's module-level testing with FDR control is strictly more rigorous.

---

## 6. What Hub-Centric Captures That Community Detection Does Not

An honest assessment of what rcomplex currently lacks:

### 6.1 Explicit hub identification

EVOTREE's modules are defined by their hubs. Each module has a named "seed gene" — often a known TF — that anchors the biological interpretation. In rcomplex, hubs exist within detected modules but are not explicitly identified. A user would need to compute within-module degree or betweenness centrality post-hoc to find them.

### 6.2 Variance-driven gene pre-selection

EVOTREE's variance > 3 filter selects biologically variable genes before module detection. This focuses the analysis on genes likely to be developmentally regulated. rcomplex's `compute_network()` has a `min_var` parameter that removes near-constant genes, but it does not identify "highly variable" genes as a distinct class.

### 6.3 Immediate biological narrative

When a module is built around a known TF hub (e.g., "the VND7 module"), the biological story is immediate. Community detection modules require downstream enrichment analysis to interpret. This is a presentation difference, not an algorithmic one — but it matters for biological communication.

---

## 7. Cross-Species Clique Architecture

| Feature | EVOTREE | rcomplex |
|---------|---------|----------|
| Algorithm | `igraph::max_cliques()` on co-expressolog graph | C++ Bron-Kerbosch with pivoting + backtracking gene assignment |
| Graph type | Cross-species ortholog network | Within-species networks + HOG mapping |
| Gene assignment | All possible (combinatorial) | One best per species clique (minimize mean q-value) |
| Partial cliques | Manual filtering (>= 11/15 edges) | `max_missing_edges` parameter |
| Stability testing | None | Leave-k-out jackknife (`clique_stability()`) |
| Threshold robustness | None | Density sweep persistence (`clique_persistence()`) |
| Classification | 4 categories (complete, partial, differentiated, lineage-specific) | 5-class HOG waterfall (complete, partial, differentiated, trait_specific, unclassified) |

rcomplex's clique engine is strictly more capable. The jackknife stability analysis and threshold sweep persistence have no analog in EVOTREE.

---

## 8. Recommendations

### 8.1 No EVOTREE module detection features need incorporation

The hub-centric iterative coverage algorithm is a less principled approach to module detection than Leiden consensus. It has no optimization objective, is sensitive to hub ranking order, lacks resolution control, has no K=1 test, and uses an arbitrary proportion threshold for conservation. rcomplex already implements the superior approach at every level.

### 8.2 Optional: `identify_module_hubs()` utility

A lightweight post-processing function that identifies hub genes per module from `detect_modules()` output. This would provide the biological interpretability that EVOTREE gets from its hub-centric design, without compromising module detection quality.

Possible signature:
```r
identify_module_hubs(modules, net,
                     method = c("degree", "betweenness", "eigenvector"),
                     top_n = 5L)
```

Returns the top-N hub genes per module by within-module centrality. Low priority — nice-to-have for downstream reporting.

### 8.3 GRN inference is out of scope

EVOTREE's ATAC-seq + JASPAR motif + TF correlation pipeline is an entirely separate analysis layer. It consumes module and clique output as input. If rcomplex users need GRN inference, they should use dedicated tools (SCENIC, GRNBoost2, or custom ATAC-seq pipelines). This is not a gap in rcomplex's module engine.

### 8.4 No changes needed to conservation assessment

rcomplex's `compare_modules()` with either hypergeometric or Jaccard permutation testing, combined with proper q-value correction (Storey or Liang), is strictly more rigorous than EVOTREE's proportion-of-significant-orthologs approach.

---

## 9. Summary

| Feature | EVOTREE | rcomplex | Action needed |
|---------|---------|----------|---------------|
| Module detection | Hub-centric greedy | Leiden/Infomap/SBM + consensus | None |
| Multi-resolution | No | Jeub *et al.* (2018) iterative | None |
| K=1 null test | No | Spectral norm permutation | None |
| Module conservation | Proportion threshold | Hypergeometric or Jaccard perm + q-values | None |
| Hub identification | Built into algorithm | Not explicit | Optional utility |
| GRN inference | ATAC + JASPAR + TF | Not in scope | Out of scope |
| Cross-species cliques | `igraph::max_cliques` | C++ BK + backtracking + stability | None |
| Clique stability | No | Leave-k-out jackknife | None |
| Threshold robustness | No | Density sweep persistence | None |
| HOG classification | 4 categories | 5-class waterfall | None |

---

## References

- Dimén *et al.* (2025). An evo-devo resource for wood comparative regulomics across dicot and conifer trees. Research Square, rs-7656402/v1.
- Jeub *et al.* (2018). Multiresolution consensus clustering in networks. *Scientific Reports*, 8, 3259.
- Lancichinetti & Fortunato (2012). Consensus clustering in complex networks. *Scientific Reports*, 2, 336.
- Netotea *et al.* (2014). ComPlEx: conservation and divergence of co-expression networks in *A. thaliana*, *Populus* and *O. sativa*. *BMC Genomics*, 15, 106.
- Rosvall & Bergstrom (2008). Maps of random walks on complex networks reveal community structure. *PNAS*, 105(4), 1118--1123.
- Storey & Tibshirani (2003). Statistical significance for genomewide studies. *PNAS*, 100(16), 9440--9445.
- Traag *et al.* (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
