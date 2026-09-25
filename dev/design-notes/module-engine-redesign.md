# Module engine: assessment and redesign options

Date: 2026-09-23. Status: design note, nothing implemented. Audience:
Martin and an implementing agent. Martin: "the module engine is
fundamentally flawed", asked whether PR #4 is a starting point and what
the literature offers, "from the computational/mathematical side rather
than bioinformatics" too, and "I do not need to shoehorn Leiden in if
there are better algorithms at hand, I just used it because it is super
fast."

Companion notes: [network-sparsification-plan.md](network-sparsification-plan.md)
(the gated work-package pattern this note reuses),
[mdl-engine.md](mdl-engine.md) (MDL for *edge* selection; this note is
MDL for *module* selection, a different question),
[clique-classifier-consolidation.md](clique-classifier-consolidation.md).

## 1. Recommendation

The engine is flawed at the level of the **unit**, not the algorithm.
Independently detected per-species partitions of n = 20 networks are
not reproducible objects, so no pairwise preservation statistic
downstream can rescue them. Three things follow.

0. **Fix the null and the unit before the algorithm.** At n = 20 every
   gene is a unit vector in R^19, so the correlation graph is a random
   geometric graph under H0 and every configuration-model or MDL null
   sees structure in it (Section 5.1). Shuffled expression is the only
   valid null, and the only statistic with a clean null is a
   cross-species one, because the eight noise geometries are
   independent under permutation. Within-species module p-values go,
   whichever engine is chosen.
1. **Detect modules jointly across species**, with orthology as the
   coupling, so that a module's evidence is pooled over eight species
   instead of being matched post hoc between two unstable partitions.
   The published objective is OrthoClust's (Section 4); the same
   objective is Mucha's multislice modularity with ortholog inter-layer
   edges, and Peixoto's layered SBM is the MDL version with model
   comparison built in (Section 5.3, design B).
2. **Gate it** exactly as the sparsification plan gates T1 to T3: one
   probe script, one go/no-go statistic (split-half replication of the
   joint partition on real minus shuffled expression, versus the
   per-species baseline), no package code before the gate passes.

PR #4 is not the starting point. Its instinct (ortholog-anchored local
structure instead of global per-species partitions) is right and is
kept here as one of the candidate designs (Section 7, design C), but the
branch is an analysis script with a ZDS interpretation on top, not an
engine. Salvage the idea, not the branch.

Added 2026-09-23 after Martin supplied Tanay 2002 and Ben-Dor 1999:
design D (Section 7), significance-scored biclustering on the gene x
sample graph with orthology as a shared row set, is the one candidate
that removes the geometric-null problem instead of working around it,
because its null is a random bipartite graph with fixed gene degrees,
which is the shuffled-expression null in closed form. It should be
probed beside A and B, not after them.

## 2. What is broken, with the measurements

All numbers from the 2026-09-15 leaf-only diagnostics (n = 20 per
species: 5 time points x 4 replicates; MR density 0.03; scripts were in
a session scratchpad, not the repo).

| observation | value | consequence |
|---|---|---|
| Leiden seed stability (HVUL, modularity, res 1) | ARI 0.97 | the algorithm is deterministic enough |
| Leiden split-half stability (10 v 10 replicates) | ARI 0.15 HVUL, 0.06 FPRA | the partition is not a property of the biology |
| edge Jaccard between halves | 0.107 / 0.056 (chance 0.015) | the network itself is mostly noise at n = 20 |
| edge FDR of top-3 % Pearson edges (gene-wise shuffle null) | HVUL 0.04, FPRA 0.37 | species differ tenfold in data quality |
| Leiden on shuffled expression | 8 to 11 modules, Q 0.25 | modularity finds structure in noise |
| transitivity of shuffled-data graph | 0.18 (rewired 0.03) | n = 20 correlation graphs are geometric |
| nested SBM (graph-tool 2.98) on shuffled data | 205 blocks (1 on rewired) | MDL does not fix it: a geometric graph is compressible by blocks, so these are real blocks, not modules (5.1) |
| `detect_modules()` K = 1 test | rejects on pure noise | its rewiring null destroys the geometry every correlation graph has |
| MR degree, 5 to 95 % (HVUL) | 278 to 1057 | flat: hub-based and local-search methods have nothing to grip |
| gene-level AUROC conservation (Crow 2022, HVUL v BDIS) | 0.63 real, 0.50 shuffled, ceiling 0.78, per-gene reliability r 0.71 | the cross-species signal is real and replicates at the gene level |

Two conclusions. The signal survives at the level of *gene
neighbourhoods compared across species*; it does not survive the step
"partition one species, then compare". And the null model question is
prior to the algorithm question: any method whose null is the
configuration model (modularity, degree-corrected SBM, the K = 1 test)
will call geometric noise "structure".

Two further known confounds are unchanged by any engine: pooled
leaf + root networks are dominated by the tissue contrast (PC1 85 to
96 % of variance; run per tissue), and the trait test's label space on
four genera of two is 2^4 under the blocked null (`pvalue_resolution()`).

## 3. PR #4 assessed

Branch `experiment/clique-module-deployment`, 11 commits, 50 files,
46 k added lines (44 of the files are CSV/HTML outputs), CI failing,
diverged from `main`. Package code touched: none. The method lives in
`analysis/probe-clique-module-deployment.R` (2445 lines).

What it does:

1. `find_cliques()` at `min_species = 4`; rank cliques; keep the top 12
   HOGs as anchors (6 multi-copy, 6 "strongest remainder").
2. Per anchor and species, take the anchor gene's 1-hop MR neighbourhood
   (`threshold` and above), map neighbours to HOGs.
3. "Core" = HOGs present in the neighbourhood in at least
   `max(3, ceiling(0.5 * n_represented))` species. This is the local
   module family.
4. Compare edge profiles among core HOGs across species (topology);
   score the module eigengene against time and tissue (deployment);
   screen annual v perennial rewiring; the last eight commits are
   heterochrony / ZDS interpretation.

What is right:

- modules anchored on cross-species evidence rather than on one
  species' partition;
- paralogs resolved by the clique's gene assignment, so a multi-copy
  HOG contributes one gene per species per anchor;
- conservation expressed as recurrence across eight species, not
  pairwise matching.

This is COMODO reinvented (Zarrineh et al. 2011, Section 4): seed,
expand per species, keep what orthologs share.

What disqualifies it as a starting point:

- **Twelve hand-picked anchors.** Coverage is whatever the clique
  ranking puts first; it is a case study, not detection.
- **No null on the core.** With MR degree flat at 278 to 1057 the
  neighbourhoods are large and overlapping, so "in 50 % of species" is
  uncalibrated. The pairwise version of this recurrence is exactly
  `compare_neighborhoods()`, which does have a null.
- **Circular.** Anchors are chosen for conservation, the core is the
  conserved part of their neighbourhood, and divergence is read off
  the remainder.
- **Not package code**, and the branch would have to be rebuilt on
  0.3.1 anyway.

## 4. Published cross-species module methods (bioinformatics side)

Verified against paper or repository unless marked "(abstract only)".

| method | citation | coupling | paralogs | software | fit |
|---|---|---|---|---|---|
| OrthoClust | Yan et al. 2014, Genome Biol 15:R100, doi:10.1186/gb-2014-15-8-r100 | sum of per-species modularities + kappa x ortholog-edge agreement | explicit: each ortholog edge weighted down by copy number at both ends; unweighted version collapsed paralog families into one giant module (their Fig S4) | Julia, github.com/gersteinlab/OrthoClust, 31 commits, no release, ~3 h for 35 k genes | **the objective to use**; software not |
| fastOC | github.com/mzinkgraf/fastOC; used in Zinkgraf et al. 2020, New Phytol 228:1811, doi:10.1111/nph.16819 (abstract only) | same objective, approximated by `igraph::cluster_louvain` on one merged graph, ortholog weight `(1/n_A + 1/n_B)/2` | as OrthoClust | R, 44 commits, dormant, WGCNA dependency | closest plant precedent (13 tree species, conserved v lineage-specific wood modules); global null model is an approximation |
| SCSC | Cai et al. 2010, PLoS Comput Biol 6:e1000707, doi:10.1371/journal.pcbi.1000707 | probabilistic: orthologs *encouraged*, not forced, into one cluster | soft | none maintained | likelihood analogue of OrthoClust; idea only |
| COMODO | Zarrineh et al. 2011, NAR 39:e41, doi:10.1093/nar/gkq1275 | seed module expanded per species until ortholog sharing is statistically optimal | via orthology tables | KU Leuven web page, dormant, <= 3 species | what PR #4 rediscovered; its stopping rule is the null #4 lacks |
| BiTSC | Sun, Zhou, Li 2021, Bioinformatics 37:1225, doi:10.1093/bioinformatics/btaa741 | orthology as bipartite graph, expression as node covariates, bipartite spectral clustering | designed for many-to-many and orphans | Python, github.com/edensunyidan/BiTSC, two species; n-partite sketched, not built | **borrow the subsample-consensus tight clustering** (Section 6) |
| ManiNetCluster | Nguyen, Blaby, Wang 2019, BMC Genomics 20(S12):1003, doi:10.1186/s12864-019-6329-2 | manifold alignment into a shared latent space, k-medoids | 0/1 correspondence matrix, many-to-many never discussed | R over Python, last commit 2019 | two networks only; no |
| IsoRankN / Ficklin & Feltus | Liao et al. 2009, Bioinformatics 25:i253; Ficklin & Feltus 2011, Plant Physiol 156:1244, doi:10.1104/pp.111.173047 | network alignment, homology weight alpha | forces a node mapping | C++, old | grass precedent (maize v rice), but alignment tops out at 5 to 6 networks and does not yield modules |
| Stuart metagenes | Stuart et al. 2003, Science 302:249 | 1:1 RBH orthologs, joint P over species | none | none | historical |
| Bergmann / ISA; Piasecka | Bergmann et al. 2004, PLoS Biol 2:E9; Piasecka et al. 2012, BMC Genomics 13:124 | per-species biclusters compared through orthologs | none | none | historical |
| WGCNA lineage | Oldham 2006 PNAS; Miller 2010 PNAS; Langfelder et al. 2011, PLoS Comput Biol 7:e1001057 | per-species modules, Zsummary preservation | 1:1 | WGCNA / NetRep | the lineage the current engine descends from and whose flaw it inherits |
| ComPlEx | Netotea et al. 2014, BMC Genomics 15:106 | pairwise hypergeometric neighbourhood tests | many-to-many, gene level | rcomplex | gene level only; no modules |
| CoCoCoNet / Crow 2022 | Lee et al. 2020, NAR 48:W566; Crow et al. 2022, NAR 50:4302, doi:10.1093/nar/gkac276 | aggregate rank networks, neighbour-voting AUROC | 1:1 OrthoDB | R/web | the gene-level score that *did* replicate on our data; not a module method |
| MVBC | Sun et al. 2016, Bioinformatics 32:i137, doi:10.1093/bioinformatics/btw278 | joint sparse rank-one factorisation | 1:1 common gene set | R/C++ | excludes the paralog-rich HOGs; no |
| MEGENA, DiffCoEx, DINGO, multiWGCNA | Song & Zhang 2015; Tesson 2010; Ha 2015; Tommasini & Fogel 2023 (Bioconductor) | shared gene set across conditions | none | R | only after collapsing to HOGs; no |
| Juxtapose, GenePlexusZoo | Ovens et al. 2021, BMC Bioinf 22:125; Mancuso et al. 2024, PLoS Comput Biol 20:e1011773 | random-walk / node2vec on a union graph with ortholog cross-edges | GenePlexusZoo: degree-weighted many-to-many eggNOG | Python | embeddings, not modules; same union-graph object as OrthoClust |
| Russell et al. 2023 | PLoS Comput Biol 19:e1011616, doi:10.1371/journal.pcbi.1011616; github.com/russell-madison/corr_comm_detection | multilayer modularity on correlation matrices, correlation-matrix null, data-driven omega, GenLouvain + consensus | identity coupling (tissue layers) | MATLAB | **borrow the correlation null and the generalist/specialist readout** |
| EVOTREE | Rodriguez et al. 2026, Nat Commun 17:8916, doi:10.1038/s41467-026-75624-2 | per-orthogroup gene graph, cliques | one gene per species per clique | R | clique-based; rcomplex reproduces it already |
| TEA-GCN, CoNekT, Plant Atlas Viewer | Lim et al. 2026, Nat Commun, doi:10.1038/s41467-026-72380-1; Proost & Mutwil 2018 | per-species networks / clusters | none | Python / web | network construction and lookup, not joint detection |
| Ruprecht 2017; Julca 2021 | Plant J 90:447 (abstract only); Nat Plants 7:1143 | per-species clusters mapped through families; orthogroup-level organ programs | family level | none | not joint detection |

Review: Ovens, Eames & McQuillan 2021, Front Genet 12:695399,
doi:10.3389/fgene.2021.695399, tabulates IsoRankN, OrthoClust, WGCNA,
BiNA, SCHype, COMODO, ManiNetCluster, Juxtapose and concludes alignment
methods top out at 2 to 6 species and that multilayer approaches are
the open direction. No 2024 to 2026 paper was found that jointly detects
modules across Poaceae with paralog-explicit coupling.

Not found or unverifiable: "phylo-WGCNA", "MultiSpecies WGCNA",
"Sherlock", "SimBic".

Statistical framing for the trait step (unchanged by the engine):
Dunn et al. 2018, PNAS 115:E409, doi:10.1073/pnas.1707515115 (pairwise
cross-species comparisons are non-independent); EVE (Rohlfs & Nielsen
2015, Syst Biol 64:695), CAGEE (Bertram et al. 2023, MBE 40:msad106) and
`phytools::phylANOVA` take one continuous statistic per species per
module, so a per-layer module statistic can be fed in, with residual df
about 6 on eight tips; phylogenetically structured permutations
(Lapointe & Garland 2001, J Classif 18:109; Adams & Collyer 2015,
Evolution 69:823, `RRPP`) make the free null honest but do not lift the
2^4 floor.

## 5. Computational side: MDL, SBM, hierarchy, stability, correlation-native nulls

Verified against Crossref, arXiv, publisher or repository unless marked
"unverified". Read 5.1 first: it is the finding that orders everything
else.

### 5.1 The null is a random geometric graph, and that is the whole problem

With n = 20 samples each standardised gene is a unit vector in R^19,
Pearson r is a cosine, and a thresholded correlation graph is a
spherical random geometric graph (RGG) in 19 dimensions *even under
H0*. Mutual rank is a per-node monotone transform of the same cosines,
so it stays geometric. This is the mechanism behind every number in
Section 2: transitivity 0.18 on shuffled data, 8 to 11 Leiden modules,
205 nested-SBM blocks.

- Devroye, György, Lugosi & Udina 2011, Electron J Probab 16,
  doi:10.1214/EJP.v16-967: clique number of high-dimensional RGGs is far
  above Erdős–Rényi when d is fixed and p grows.
- Bubeck, Ding, Eldan & Rácz 2016, Random Struct Alg 49:503,
  doi:10.1002/rsa.20633: geometry is detectable from signed triangle
  counts whenever d is small relative to p^3; at d = 19, p = 20 000 it
  is overwhelming.
- Krioukov 2016, Phys Rev Lett 116:208302 ("clustering implies
  geometry"); Boguñá et al. 2021, Nat Rev Phys 3:114,
  doi:10.1038/s42254-020-00264-4.
- Guimerà, Sales-Pardo & Amaral 2004, Phys Rev E 70:025101 (modularity
  of random graphs is large); Good, de Montjoye & Clauset 2010, Phys Rev
  E 81:046106 (the modularity landscape is degenerate: the split-half
  ARI 0.06 to 0.15).
- Schaub, Delvenne, Yaliraki & Barahona 2012, PLoS ONE 7:e32210,
  doi:10.1371/journal.pone.0032210: the field-of-view limit. Modularity
  *and* Infomap see only clique-like communities; geometric structure at
  scale looks like many small modules to both. So Infomap is not an
  escape.
- Peixoto 2023, *Descriptive vs. inferential community detection*,
  Cambridge Elements, doi:10.1017/9781009118897: inferential methods
  report only compressible structure. That does not rescue us: an RGG
  *is* compressible by blocks (spatial neighbourhoods), so the 205 SBM
  blocks on noise are real blocks, they are just not modules.
- MacMahon & Garlaschelli 2015, Phys Rev X 5:021006,
  doi:10.1103/PhysRevX.5.021006: modularity on correlation matrices
  with a random-matrix (Marchenko–Pastur) null. At p/n = 1000 the bulk
  edge is lambda_+ = (1 + sqrt(1000))^2 ~ 1064 while a rank-19
  correlation matrix of trace 20 000 has 19 non-zero eigenvalues
  averaging ~1053: essentially the whole spectrum sits at the bulk edge
  and the RMT null declares nearly everything noise. That is the honest
  verdict on within-species module significance at n = 20, not a reason
  to avoid the method.

Consequences, which Sections 6 to 8 build on:

1. Shuffled expression is the RGG null. Degree-preserving rewiring is
   not, and any statistic with a configuration-model null (modularity,
   DC-SBM, OSLOM, conductance, the K = 1 test) rejects "no structure" on
   noise. This is a property of the data, not of Leiden.
2. **The one thing shuffling cannot fake is cross-species agreement.**
   Each species' noise geometry is independent under permutation. So a
   statistic of the form "genes that co-cluster jointly are
   co-expressed across species beyond what permuted expression gives"
   has a clean null; a within-species statistic at n = 20 does not.
   This is exactly why the gene-level AUROC replicated (Section 2) and
   the per-species partitions did not.
3. Spatial null models for modularity (Expert, Evans, Blondel &
   Lambiotte 2011, PNAS 108:7663) need an external geometry; for a
   correlation graph the "space" is the data itself, so they are
   circular here. Geometric block models (Galhotra et al. 2018, AAAI;
   Abbe, Baccelli & Sankararaman 2018, SODA; arXiv 2403.02802) are
   theory without maintained software (unverified).

### 5.2 MDL and inferential community detection

- Karrer & Newman 2011, Phys Rev E 83:016107 (degree-corrected SBM);
  Newman 2016, Phys Rev E 94:052315 (modularity at fixed resolution is
  DC-SBM maximum likelihood, so Leiden already is an SBM fit with K and
  gamma chosen for it).
- Peixoto 2014, Phys Rev X 4:011047 (nested SBM, MDL); Peixoto 2017,
  Phys Rev E 95:012317 (the microcanonical MDL actually implemented).
  Agglomerative heuristic O(N ln^2 N), "reliable results for networks in
  excess of 10^7 edges"; graph-tool's `multilevel_mcmc_sweep` is
  O(E ln^2 N). One species at 20 k x 6 M edges is inside the tested
  range (tens of minutes to hours per fit); the eight-layer union at
  ~50 M edges is beyond the paper's largest example but near-linear in
  E, so day-scale per fit (estimate). Top-k at 8 x 200 k edges: minutes.
  Posterior sampling multiplies by the number of sweeps. Python only;
  graph-tool 2.98 is installed on Martin's Mac; R via reticulate.
- Zhang & Peixoto 2020, Phys Rev Research 2:043271 (assortative
  "planted partition" SBM with MDL; graph-tool `PPBlockState`). The
  inferential replacement for modularity when modules, not arbitrary
  blocks, are wanted. The 205 nested-SBM blocks on shuffled data were
  general blocks; the assortative-only MDL fit is the cleaner question
  and has not been run.
- Peixoto 2018, Phys Rev E 97:012306 (weighted SBM, `rec_types`
  real-normal / real-exponential): fits the top-k *weighted* graph
  without a hard threshold.
- Peixoto 2021, Phys Rev X 11:021003 (partition modes: consensus *and
  dissensus* of a partition ensemble; graph-tool `PartitionModeState`,
  `ModeClusterState`). Usable today on the existing Leiden ensembles:
  split-half ARI 0.1 is dissensus, and this says how many alternative
  module systems there are and which genes are stable across modes.
- Rosvall & Bergstrom 2007, PNAS 104:7327 (compression view);
  Rosvall & Bergstrom 2011, PLoS ONE 6:e18209 (hierarchical Infomap);
  De Domenico, Lancichinetti, Arenas & Rosvall 2015, Phys Rev X
  5:011027 (multilayer Infomap); Edler, Bohlin & Rosvall 2017,
  Algorithms 10:112, doi:10.3390/a10040112 (state-node formulation).
  Infomap v2.15.1 (2026), C++ with Python bindings, R via r-universe and
  `infomapecology` (Farage et al. 2021, Methods Ecol Evol 12:778).
  Multilayer input allows links between *different* physical nodes in
  different layers, and several state nodes per physical node: set
  physical node = HOG, state nodes = (species, copy), intra-layer links
  = co-expression, and paralogs are first-class with no coupling
  constant (a relax rate instead). Different objective (flow), so not
  comparable with modularity results, and field-of-view applies.
- Newman & Reinert 2016, Phys Rev Lett 117:078301 (Bayesian K for the
  DC-SBM): single layer, superseded by nested-SBM MDL in practice.

### 5.3 Multilayer SBMs where layers are different node sets

Only one published model handles "different nodes per layer, coupled by
an explicit inter-layer edge set", and it does so by construction:

- **Peixoto 2015, Phys Rev E 92:042807, doi:10.1103/PhysRevE.92.042807**
  (layered SBM; graph-tool `LayeredBlockState`, `LayeredNestedBlockState`,
  edge covariate `ec`, `independent = TRUE/FALSE`). Layers are edge
  covariates on one node set, so take the union of all (species, gene)
  nodes, give each species' co-expression edges its own layer, and put
  the ortholog edges in a ninth layer; a node absent from a layer has
  degree zero there. `clabel` / `pclabel` = species forces
  within-species blocks, after which the block-block edge counts of the
  ortholog layer *are* the module-correspondence matrix, inferred
  jointly with the modules. Leaving the constraint off allows
  cross-species blocks. The description-length difference between fits
  with and without the ortholog layer is a direct test of whether
  orthology explains co-expression block structure; against shuffled
  expression it is a clean one (5.1, consequence 2).

The rest assume a shared node set with identity coupling, so they
apply only after collapsing to HOG nodes (losing paralogs): Stanley,
Shai, Taylor & Mucha 2016, IEEE Trans Netw Sci Eng 3:95 (strata
multilayer SBM: layers grouped into strata sharing one SBM; the idea
"strata = annual v perennial" is attractive but identity-coupled);
Vallès-Català et al. 2016, Phys Rev X 6:011036; Paul & Chen 2016,
Electron J Stat 10:3807 and 2020, Ann Stat 48:230; De Bacco, Power,
Larremore & Moore 2017, Phys Rev E 95:042317 (MULTITENSOR, mixed
membership); Bazzi et al. 2020, Phys Rev Research 2:023100 (generative
benchmark, useful for synthetic tests; node-set flexibility
unverified).

### 5.4 Correlation-native methods (no thresholding)

- MacMahon & Garlaschelli 2015 (5.1): RMT null; MATLAB reference, any
  Louvain with a custom B. Dense 20 k^2 per species is feasible; the
  verdict at n = 20 is "almost nothing".
- Masuda, Kojaku & Sano 2018, Phys Rev E 98:012312 (maximum-entropy
  configuration model for correlation matrices); Kojaku & Masuda 2019,
  Proc R Soc A 475:20190578 (Scola, github.com/skojaku/scola, Python):
  an edge only where the correlation is unexpected under a null chosen
  by a BIC-type criterion that knows n. The correlation-native
  replacement for "MR + hard threshold". O(p^2) memory (3.2 GB at
  20 k), iterative lasso; scalability to 20 k nodes unverified; expect
  very sparse output at n = 20, and that sparseness is information.
- Russell et al. 2023 (Section 4): Masuda 2018 null per layer,
  GenLouvain + consensus, 203 genes, identity coupling.
- Bazzi et al. 2016, Multiscale Model Simul 14:1 (temporal multilayer
  modularity on correlation networks with RMT nulls; the null plumbing
  for GenLouvain).
- Hoffmann, Peel, Lambiotte & Jones 2020, Sci Adv 6:eaav1478,
  doi:10.1126/sciadv.aav1478: communities inferred end to end from the
  node time series, no edges ever formed, full posterior. Hundreds of
  nodes; not scalable to 20 k without re-engineering. Related: Peixoto
  2019, Phys Rev Lett 123:128301 (dynamics, not iid samples).
- Bongiorno, Miccichè & Mantegna 2022, Physica A 593:126933: bootstrap
  replicas of the dissimilarity matrix give a p-value per clade of a
  hierarchical clustering. Correlation-native, but bootstrapping 20
  samples is thin.

### 5.5 Markov stability and scale selection

Delvenne, Yaliraki & Barahona 2010, PNAS 107:12755; Lambiotte,
Delvenne & Barahona 2014, IEEE Trans Netw Sci Eng 1:76; Arnaudon et al.
2024, "Algorithm 1044: PyGenStability", ACM Trans Math Softw 50(2),
doi:10.1145/3651225 (Python, Louvain or Leiden, robust scale selection
by NVI across scales and across runs). Custom constructors are
supported, so a supra-adjacency with eight intra-species blocks and an
ortholog coupling block can be supplied; Mucha's multislice modularity
is itself a Laplacian-dynamics derivation, so this is well defined.
`linearized` constructor at 6 M edges; the continuous ones need a
matrix exponential. Caveat: robust scales on an RGG are still
geometric scales. NVI fixes the split-half instability, not the null.

### 5.6 Significance and stability of communities

- Lancichinetti, Radicchi, Ramasco & Fortunato 2011, PLoS ONE 6:e18961
  (OSLOM: local significance of a cluster via order statistics of its
  worst member; the principled "seed, expand, stop when no longer
  significant"; C++, R wrapper `bioregion::netclu_oslom`). Null is
  configuration, so geometric clusters pass.
- Zhang & Moore 2014, PNAS 111:18144 (belief propagation, retrieval
  modularity as a structure-vs-null test). Same null caveat.
- Rosvall & Bergstrom 2010, PLoS ONE 5:e8694 (bootstrap edge weights,
  recluster, report significance cores). Bootstrapping *samples*
  instead of edges turns this into the split-half test already run.
- Decelle, Krzakala, Moore & Zdeborová 2011, Phys Rev E 84:066106
  (detectability threshold): our problem is the opposite regime, too
  much detectable structure.
- Lancichinetti & Fortunato 2012, Sci Rep 2:336; Jeub, Sporns &
  Fortunato 2018, Sci Rep 8:3259 (both already in `detect_modules()`).
- Monti et al. 2003, Mach Learn 52:91 (consensus clustering); von
  Luxburg 2010, Found Trends Mach Learn 2:235 (stability can be high at
  the wrong K; it is not significance); **Şenbabaoğlu, Michailidis &
  Li 2014, Sci Rep 4:6207, doi:10.1038/srep06207** (consensus
  clustering finds stable structure on null data unless calibrated
  against a null; PAC score); **Tseng & Wong 2005, Biometrics 61:10,
  doi:10.1111/j.0006-341X.2005.031032.x** (tight clustering: subsample,
  keep only tight stable cores, leave the rest unassigned; what BiTSC
  borrowed). Ballouz, Verleyen & Gillis 2015, Bioinformatics 31:2123
  (sample-size guidance for co-expression).
- Ghasemian, Hosseinmardi & Clauset 2020, IEEE TKDE 32:1722 (held-out
  link prediction as model selection across 16 methods).
- Peel, Larremore & Clauset 2017, Sci Adv 3:e1602548,
  doi:10.1126/sciadv.1602548: metadata is not ground truth; their
  BESTest is the right shape for "module divergence v annual/perennial"
  on eight tips.

### 5.7 Local, seed-based and overlapping communities

- Palla, Derényi, Farkas & Vicsek 2005, Nature 435:814 (clique
  percolation): the original clique-anchored method and PR #4's closest
  relative; on a high-dimensional RGG cliques percolate massively, so it
  returns giant components on noise.
- Lancichinetti, Fortunato & Kertész 2009, New J Phys 11:033015 (local
  fitness expansion; OSLOM is its significance-based successor).
- Andersen, Chung & Lang 2006, FOCS, doi:10.1109/FOCS.2006.44
  (PageRank-Nibble, conductance guarantee via local Cheeger); Kloster &
  Gleich 2014, KDD (HK-relax); Whang, Gleich & Dhillon 2016, IEEE TKDE
  28:1272 (NISE seed-set expansion); Yang & Leskovec 2013, Knowl Inf
  Syst 42:181 (conductance among the best scoring functions against
  ground truth); Coscia et al. 2012, KDD (DEMON). Software:
  LocalGraphClustering (github.com/kfoynt/LocalGraphClustering, Python;
  ACL, HK, MQI, FlowImprove, SimpleLocal; 100 M-edge graphs on a
  laptop). Conductance has no null, so calibrate against shuffled
  expression; the cross-species comparison of per-species expansions
  from one ortholog-clique seed is the quantity with a clean null.
- Ahn, Bagrow & Lehmann 2010, Nature 466:761 (link communities): cost
  ~ sum of k^2 over nodes, fine at 200 k edges, heavy at 6 M.

### 5.8 Hierarchy beyond nested SBM

Ravasz & Barabási 2003, Phys Rev E 67:026112 (descriptive signature);
Clauset, Moore & Newman 2008, Nature 453:98 (HRG, thousands of nodes at
most); Sales-Pardo et al. 2007, PNAS 104:15224 (co-classification
against the modularity-fluctuation null; small networks); Bonald et al.
2018, arXiv 1806.01664 (Paris: agglomerative, near-linear,
parameter-free; scikit-network, cdlib); Lyzinski et al. 2017, IEEE
Trans Netw Sci Eng 4:13 (HSBM via spectral embedding; no official
software); Schaub, Li & Peel 2023, Phys Rev E 107:054305 (what
"hierarchical" should mean, with a spectral test). Multilayer
hierarchical options that exist today: graph-tool nested + layered;
Infomap multilevel on multilayer input; Jeub 2018 hierarchical consensus
over any partition ensemble.

### 5.9 Multilayer modularity optimisers and R coverage

- Mucha, Richardson, Macon, Porter & Onnela 2010, Science 328:876,
  doi:10.1126/science.1184819. GenLouvain (MATLAB + MEX, v2.2 2019)
  accepts an arbitrary modularity matrix B, so ortholog inter-layer
  edges are encodable.
- leidenalg (Traag; Traag, Waltman & van Eck 2019, Sci Rep 9:5233):
  `optimise_partition_multiplex` needs all layers on one vertex set but
  allows a different partition type, resolution and weight per layer.
  Hand-built construction: vertex set = all (species, gene) pairs, one
  `RBConfigurationVertexPartition` per species (its own null), one
  `CPMVertexPartition(resolution_parameter = 0, node_sizes = 0)` layer
  holding the ortholog edges, layer weight = kappa. Undocumented but
  uses only documented API.
- igraph R: single graph, but the CPM `vertex_weights` trick reproduces
  the per-layer null (Section 8, P2, verified).
- multinet 4.3.4 (CRAN 2026-03; Magnani, Rossi & Vega 2021, J Stat
  Softw 98(8)): `glouvain_ml(gamma, omega)`, `infomap_ml`,
  `clique_percolation_ml`; the data model accepts inter-layer edges
  between *different* actors, but whether `glouvain_ml` uses them
  rather than only the omega identity coupling is unverified.
- muxViz (De Domenico, Porter & Arenas 2015, J Complex Netw 3:159; R,
  Infomap-based, GUI-first); MolTi (Didier, Brun & Baudot 2015, PeerJ
  3:e1525; same node set only); leidenAlg, leidenbase, netmem, cdlib:
  single layer or no multiplex section.

### 5.10 Biclustering on the gene x sample graph

Two papers Martin supplied on 2026-09-23 as "highly relevant", both read
in full, plus the cross-species precedent they lead to.

- **Tanay, Sharan & Shamir 2002 (SAMBA)**, Bioinformatics 18:S136,
  doi:10.1093/bioinformatics/18.suppl_1.S136. Expression becomes a
  bipartite graph, conditions U x genes V, with an edge when the gene's
  standardised level in that condition is above 1 or below -1 (a signed
  variant finds *consistent* biclusters by a reduction that doubles the
  graph). A bicluster is a heavy subgraph. Weights are log-likelihood
  ratios: an edge scores `log(p_c / p_uv)`, a non-edge
  `log((1 - p_c) / (1 - p_uv))`, where `p_uv` is the probability of that
  edge in a random bipartite graph with the observed degree sequence
  (estimated by Monte Carlo) and `p_c > max p_uv` is the bicluster's
  constant edge probability (0.9 in practice). So a bicluster's weight
  *is* its significance against a degree-preserving null. Polynomial
  algorithm under bounded gene degree d: hash every subset of each
  gene's neighbourhood of size N1 to N2, keep the k heaviest bicliques
  per gene, local add/remove improvement, greedy overlap filter. A
  p-value by the Liapunov CLT on the weight of the best subgraph for a
  fixed condition set, Bonferroni over condition subsets. Validation: a
  random bipartite graph with the same degree sequence gives
  significance values well separated from the real data (their Fig.
  3c). 15 000 genes x 500 conditions in minutes on a 2002 PC. Software:
  SAMBA shipped inside EXPANDER (Shamir lab, Java GUI); EXPANDER 7.2
  (October 2017) is the last release and its page highlights ISA, not
  SAMBA; no library, no CLI. A reimplementation is a few hundred lines.

  Why it matters here (5.1): the null sits on the gene x sample
  *response* graph, and gene-wise shuffling of expression is exactly a
  random bipartite graph with fixed gene degrees. SAMBA's null is our
  shuffled-expression null in closed form, and a gene x gene
  correlation graph is never formed, so the RGG structure never
  arises. With 20 samples a gene's neighbourhood has at most 20
  conditions, so the bounded-degree algorithm is exact, not heuristic.
  The price is binarisation: a bicluster says "these genes all respond
  in these samples", not "these genes are correlated". On a 5 time
  point x 4 replicate design that is close to what a co-expression
  module means, and the signed variant keeps direction.

- **Ben-Dor, Shamir & Yakhini 1999 (corrupted clique, CAST)**, J Comput
  Biol 6:281, doi:10.1089/106652799318274. Generative model: the true
  clustering is a clique graph (disjoint union of cliques) and the
  observed similarity graph flips every edge and non-edge independently
  with probability alpha < 1/2. PCC (theoretical) recovers the
  clustering with high probability in O(n^2 log^c n) by sampling a small
  core and classifying every other vertex by attraction; CAST
  (practical) grows one cluster at a time by an affinity threshold t
  with add *and remove* steps, then moves vertices to their
  highest-affinity cluster until stable. Matlab, 1999; the model, not
  the code, is the asset.

  Why it matters: (i) it is the generative model behind rcomplex's
  clique layer, a conserved module being a clique corrupted by noise,
  and the per-species alpha is exactly the edge FDR measured on shuffled
  data (HVUL 0.04, FPRA 0.37), so it says which species' networks are
  recoverable at all (alpha near 1/2 is not); (ii) its independent-error
  assumption is what the RGG violates (transitivity 0.18), so the
  guarantee does not transfer to a thresholded correlation graph, which
  is 5.1 stated a third way, while on SAMBA's bipartite graph it does
  hold; (iii) CAST's add/remove clean-up is the step design C's seed
  expansion lacks.

- **Multi-species cMonkey** (Waltman et al. 2010, Genome Biol 11:R96,
  PMC2965388; full text read via Europe PMC). The cross-species
  biclustering precedent. Phase 1, shared space: biclustering over an
  "orthologous core" (InParanoid families, paralogs allowed), where an
  ortholog pair's membership score is a logistic combination of the two
  species' single-species cMonkey scores, `pi_ik ∝ exp(b0 + b1 (g_U +
  g_V))`, "easily extended to more than two organisms"; each species
  keeps its own condition space, so a bicluster is a shared gene core
  with a per-species condition set. Phase 2, elaboration: per-species
  addition of species-specific genes with the core locked, and several
  paralogs per family allowed. Three Firmicutes. R code was at
  meatwad.bio.nyu.edu (not re-verified; cMonkey2, Reiss et al. 2015,
  NAR, is the maintained successor and single-species). No explicit
  null: scores are cMonkey's per-data-type likelihood P-values combined
  under annealing, not SAMBA's statistical model.

- **Hochbaum 1998**, "Approximating clique and biclique problems",
  J Algorithms 29:174, doi:10.1006/jagm.1998.0964 (Martin, 2026-09-23,
  "interesting, too"; read in full). The complexity map behind SAMBA's
  design choices, and the exact tool for the step SAMBA does
  heuristically. Results: (i) the *node*-deletion biclique problem on a
  bipartite graph (delete the minimum weight of nodes so that what is
  left is complete bipartite) is polynomial, by one minimum cut, since
  it is maximum-weight independent set on the bipartite complement
  (Yannakakis 1981); the general-graph version without the independence
  requirement is also polynomial (monotone IP2), and equals node
  connectivity of the complement; (ii) the *edge*-deletion biclique
  problem is NP-hard (Dawande et al. 1997, by reduction from maximum
  clique, which is why Tanay 2002 needs the bounded-degree restriction)
  but 2-approximable by a single minimum cut on a bipartite network
  whose node weight is half the adjacent edge weight (their Fig. 5);
  (iii) edge-deletion to a clique is NP-hard and 2-approximable via a
  vertex cover in which edges cover non-edges; (iv) every NP-hard
  variant here is MAX SNP-hard, so 2 is the floor unless vertex cover
  falls. All via Hochbaum's IP2 half-integrality framework, minimum cut
  in T(n, m) ~ O(nm log(n^2/m)).

  Why it matters for design D: SAMBA's exactness rests on hashing every
  subset of a gene's condition neighbourhood, which is fine at 20
  conditions per species but not for a joint bicluster whose condition
  side is 8 species x 20 samples. Hochbaum supplies the cleaning step
  for the joint object: given a candidate (HOG set, sample set across
  species), form the bipartite graph HOG x (species, sample) with an
  edge when some copy of the HOG responds in that sample, and delete the
  minimum weight of HOGs and samples so that what remains is a biclique.
  That is CAST's add/remove step done *optimally* and in polynomial time
  by one minimum cut, with node weights taken from SAMBA's
  log-likelihood contributions; paralogs enter through the OR over
  copies, so no copy choice is needed. The edge-deletion 2-approximation
  is the fallback when a block is allowed to keep a few non-responding
  cells. For the clique layer the edge-deletion-clique result is the
  formal name of what `find_cliques()`'s backtracker does (fewest missing
  edges first); at <= 8 species exhaustive search is fine and the
  2-approximation is not needed.

- Related, none with SAMBA's degree-corrected significance: Bergmann
  2004 ISA across six organisms (Section 4); R `isa2` (CRAN) is ISA;
  `biclust` (CRAN: Cheng-Church, Bimax, Plaid, Xmotifs, Quest,
  spectral); `QUBIC` (Bioconductor; discretised expression, gene graph
  weighted by shared levels, greedy expansion, C, fast); `fabia`
  (Bioconductor; factor-analysis biclustering); Bi-EB (2022, empirical
  Bayes cross-species multi-omics biclustering, PMC9690013,
  unverified). None installed locally.

### 5.11 Summary table

| method | different-node-set layers via explicit inter-layer edges | resolution-free / MDL | significance or stability built in | correlation-native | software | 8 x 200 k edges | 8 x 6 M edges |
|---|---|---|---|---|---|---|---|
| nested + layered SBM (Peixoto 2014/2015) | yes: union graph + ortholog layer + `clabel` | MDL | posterior modes (Peixoto 2021); null still SBM | no (weighted variant on top-k) | graph-tool, Python | minutes | hours to a day per fit (est.) |
| assortative SBM `PPBlockState` (Zhang & Peixoto 2020) | no (single layer) | MDL | no | no | graph-tool | minutes | hours |
| Infomap multilayer (De Domenico 2015; Edler 2017) | yes: state-node links between different nodes | no gamma; MDL of flow | no (Rosvall 2010 bootstrap separate) | no | C++; Python; R via r-universe | seconds to minutes | minutes to hours |
| leidenalg multiplex | yes: hand-built ortholog CPM layer | no (gamma per layer) | no | no | Python | seconds | minutes |
| igraph `cluster_leiden` CPM trick | yes, but with a cross-layer null term: a different objective (Section 8, P2) | no | no | no | R, installed | seconds | minutes |
| GenLouvain (Mucha 2010) | yes: arbitrary B | no | no | with RMT null (Bazzi 2016) | MATLAB | minutes | memory-bound |
| multinet `glouvain_ml` | data model yes; use in objective unverified | no | no | no | R (CRAN) | minutes | untested |
| OrthoClust (Yan 2014) | yes: kappa-coupled ortholog edges | no | no | no | Julia 0.4, dead | port needed | no |
| MULTITENSOR / sMLSBM / Vallès-Català / Paul–Chen | no (identity coupling) | mixed | no | no | Python / MATLAB / none | n/a | n/a |
| BiTSC (Sun 2021) | bipartite orthology, 2 species | K by tightness | tight clustering | expression as covariates | Python | fine | n/a |
| Hoffmann 2020 no-edge SBM | no | Bayesian K | full posterior | yes | Python | hundreds of nodes | no |
| MacMahon–Garlaschelli RMT modularity | no (multilayer via Bazzi 2016) | no | RMT null | yes | MATLAB reference; any Louvain with custom B | dense 20 k^2 feasible | same |
| Scola / Masuda null | no | null-model BIC | edge test | yes | Python | O(p^2), unverified at 20 k | same |
| PyGenStability (Arnaudon 2024) | yes via custom constructor | scale selection by NVI | stability (NVI) | with custom B | Python | minutes | hours (linearized) |
| OSLOM | no | free | local significance, configuration null | no | C++; R `bioregion` | minutes | hours (est.) |
| Peixoto 2021 partition modes | post hoc, any partitions | n/a | dissensus | n/a | graph-tool | fast | fast |
| tight clustering / consensus + null (Tseng–Wong; Şenbabaoğlu) | post hoc | stability-chosen | stability, needs null calibration | either | R `tightClust` / own code | cheap | cost = detector x runs |
| ACL / HK-relax / NISE seed expansion | per layer, seeds from ortholog cliques | local, no gamma | conductance, no null | no | LocalGraphClustering, Python | instant | fine |
| clique percolation (Palla 2005) | no | k | no | no | igraph-based | fine | percolates on noise |
| Paris (Bonald 2018) | no | dendrogram | no | no | scikit-network, cdlib | fast | fast |
| link communities (Ahn 2010) | no | dendrogram | no | no | R `linkcomm` (not re-verified) | fine | heavy |
| SAMBA (Tanay 2002) on gene x sample graph | orthology as shared HOG rows (design D), not edges | free; LLR weights | yes: degree-preserving bipartite null, CLT p-value | bypasses correlation entirely | EXPANDER 7.2 (2017, Java GUI); reimplement | trivial (degree <= 20) | n/a (no gene graph) |
| multi-species cMonkey (Waltman 2010) | shared ortholog core + per-species conditions; paralogs allowed | annealed scores | no explicit null | expression-native | R 2010, unmaintained | n/a | n/a |

### 5.12 Suresh et al. 2023 (primate MTG, Gillis lab): aggregation, conservation AUROC, expressolog

Martin (2026-09-24 evening): read the paper (Nat Ecol Evol 7:1930) and the
code (github.com/hamsinisuresh/Primate-MTG-coexpression, six R scripts).
What the code actually does:

1. **Consensus states across species** (MetaNeighbor one_vs_best AUROC on
   within-species clusters, reciprocal best hits and AUROC > 0.6 in at
   least one pair) gives 57 cell types shared by five primates. Not in
   this repo (AllenInstitute/Great_Ape_MTG), but the idea is what a
   data-driven alignment of the wood gradient needs: match the
   published per-species section clusters (`DATA/sampleClusters/cct_*`)
   across species by replicability instead of by annotation boundaries
   or by my two-landmark warp.
2. **Aggregate network** (`get_coexpression_network.R`): per state,
   pseudobulk samples of 20 cells; Spearman correlation; the upper
   triangle rank-standardised to [0, 1]; the 57 state networks summed
   and rank-standardised again. Genes restricted to 4,500 by expression
   breadth. This is CoCoCoNet-style meta-analysis: composition effects
   cancel because each state contributes an equal-weight rank matrix.
   The analogue here is per-tree networks (20 to 28 sections each) on
   wood and per-tissue networks on Pooideae, rank-averaged; it is the
   answer to the tree-split Jaccard of 0.04 in 11.7 that does not
   throw away samples.
3. **Co-expression conservation** (`get_coexpression_conservation.R`):
   top-10 neighbours of gene i in species 1 as a 0/1 row, multiplied by
   species 2's rank matrix, gives for every gene j of species 2 the
   rank-sum of i's neighbours in j's row; the analytic AUROC of that
   sum; the score is the ortholog's AUROC, and the **specificity** is
   the rank of the ortholog's AUROC among all j (`roc_bg_predict`),
   bidirectional and averaged. That is the co-expressolog test with a
   continuous, self-calibrated null: no hypergeometric, no discrete
   FDR, and paralogs surface as the best-ranked j rather than being
   excluded. It also matches the AUROC conservation that already
   worked on the leaf diagnostics.
4. **Expressolog score** (`get_expressolog_score.R`): Pearson
   correlation of the two orthologs' profiles across the 57 matched
   states, ranked against all other genes of the other species, as an
   AUROC; computed at class, subclass and cluster level and averaged.
   Needs matched states; on our data those are the five time points or
   the eight gradient bins, so it is coarse but defined.
5. **Divergence** = expressolog < 0.55 in a class and co-expression
   conservation human-vs-mammal below mammal-vs-mammal (one-sided
   Wilcoxon over species pairs, BH) on 22 CoCoCoNet species. The test
   is one species against many; with six or eight species and a
   balanced contrast the pairs are as dependent as our labellings, so
   the form transfers and the power does not.

Limits for us: 1:1 orthologs only (14,131), bulk meta-networks from
thousands of samples, and a single-species divergence question.

## 6. What to borrow

Ordered by how much of the problem each removes.

1. **The null and the unit (5.1).** Shuffled expression is the only
   valid null; the reported statistic must be cross-species. Every
   candidate engine below is evaluated that way, and within-species
   module p-values are dropped from the package whichever engine wins.
   This is a reframing, not software, and it applies to the current
   engine today.
2. **Ortholog-coupled joint detection with paralog down-weighting**
   (Yan et al. 2014; fastOC). Inter-layer edge weight
   `kappa * (1 / n_A + 1 / n_B) / 2`; the unweighted version collapses
   paralog families, which with HJUB at 86 % multi-copy is decisive.
3. **Per-layer null** (Mucha 2010), so that the densest species does
   not set the null for the sparsest. Verified in R via the CPM
   `vertex_weights` trick (Section 8, P2).
4. **MDL model comparison with and without the ortholog layer**
   (Peixoto 2015), real against shuffled. This is what an SBM backend
   buys that a modularity backend cannot: a description-length gain
   that is a test, not a score, and a module-correspondence matrix
   inferred jointly with the modules.
5. **Tight-core consensus calibrated against null co-clustering**
   (Tseng & Wong 2005; Şenbabaoğlu 2014; BiTSC's use of it). Subsample
   the *samples* (not the edges), re-detect, keep only gene pairs whose
   co-clustering frequency exceeds the frequency on shuffled data, and
   leave everything else unassigned. Attacks the split-half instability
   directly, reuses the existing consensus machinery in `R/modules.R`,
   and costs only detector runs.
6. **Dissensus, not only consensus** (Peixoto 2021; Rosvall 2010).
   Report how many alternative module systems the ensemble contains and
   which genes are stable across them. Runs on the Leiden ensembles the
   package already produces.
7. **Per-layer contribution and the generalist/specialist readout**
   (Russell et al. 2023) as the per-module per-species statistic that
   replaces `Zsummary_std` and feeds `preservation_matrix_test()`.
8. **Seed expansion with a conductance or significance stop** (Andersen
   2006; Kloster & Gleich 2014; OSLOM's criterion; COMODO's stopping
   rule) as the principled version of PR #4's neighbourhood step.
9. **Paralogs as state nodes** (Edler 2017; Infomap multilayer) if an
   Infomap backend is ever tried: the only formulation where copies
   need no weighting rule.
10. **BESTest** (Peel et al. 2017) for the trait step: is the
    annual/perennial labelling of species informative about the module
    profile, against the correct label space.
11. **SAMBA's null placement** (Tanay 2002): score modules on the
    bipartite gene x sample response graph against the degree-preserving
    random bipartite null, which is the shuffled-expression null in
    closed form. The only item in this list that removes 5.1 instead of
    working around it. Design D.
12. **Corrupted-clique alpha** (Ben-Dor 1999) as the per-species
    recoverability parameter: report the shuffled-data edge FDR per
    species next to every module statistic, and use CAST's add/remove
    clean-up in any seed expansion.
13. **cMonkey's two phases** (Waltman 2010): a shared core across
    species first, then per-species elaboration with the core locked.
    That is "conserved core plus species-specific departure" as an
    algorithm, and it is the shape Section 9 asks for.


Added 2026-09-24 evening from Suresh et al. 2023 (5.12):

14. **Rank-aggregated networks across replicate units** (per tree, per
    tissue): equal-weight rank matrices summed and re-ranked. Removes
    tree and composition effects without splitting samples.
15. **Conservation AUROC with a specificity rank** as the co-expressolog
    score: continuous null, paralog-aware, replaces the hypergeometric
    and the discrete-FDR item in Section 9.
16. **Replicability-matched states** (MetaNeighbor one_vs_best) to align
    the wood gradient across species from the published section
    clusters, in place of annotation boundaries or landmark warps.
17. **Expressolog profile score** over matched states as the per-gene
    complement of the program-level shape coherence.

## 7. Candidate designs

All three consume the same cached union graph (Section 8, P1) and are
judged by the same gate (P3). They are not exclusive: C is a local
complement to A or B.

### A. Multilayer modularity with ortholog coupling (OrthoClust objective, Leiden optimiser)

Objective: sum over species of per-layer modularity, plus kappa times
ortholog-edge agreement with paralog down-weighting. R-native today via
`igraph::cluster_leiden` with the CPM trick; exact via leidenalg
multiplex through reticulate if the cross-layer null contamination
turns out to matter; fastOC is the published R precedent (13 tree
species).

- Pro: an afternoon to run; seconds per fit, so subsample consensus
  (borrow 5) and kappa sweeps are cheap; one partition gives
  `module x species` membership; per-layer Q_s(c) is the divergence
  statistic (borrow 7).
- Con: two free parameters (kappa, gamma) with no internal criterion;
  no significance of its own, so everything rests on borrow 1 and 5;
  field-of-view limit (5.1).
- kappa selection: the smallest kappa at which split-half replication
  on real data separates from replication on shuffled data (Section 8,
  P3), not an external gold standard as in OrthoClust.

### B. Layered nested SBM on the union graph (graph-tool)

Eight co-expression layers plus one ortholog layer as edge covariates
on the union node set, `clabel = species`, nested, MDL. One fit returns
joint modules, the module-correspondence matrix (ortholog-layer block
counts), per-species block matrices, a hierarchy, and a description
length to compare with the no-ortholog-layer fit and with the shuffled
fit.

- Pro: no kappa, no gamma, no K; model comparison is built in; posterior
  modes (borrow 6) come from the same machinery; `PPBlockState` gives
  the assortative-only variant when modules rather than blocks are
  wanted.
- Con: Python dependency (reticulate, graph-tool install is heavy;
  Orion has Apptainer, the Mac has 2.98); minutes at top-k, hours to a
  day per fit at full 3 % density, and subsample consensus multiplies
  that; the SBM null is still not geometric, so the 205-blocks-on-noise
  behaviour will recur within species and only the ortholog-layer
  description-length gain against shuffled expression is a test.
- If B is adopted the package boundary is a graph export plus a result
  import, not an embedded solver; the note in the sparsification plan
  about not porting PANINIpy applies here too.

### C. Seed expansion from co-expressolog cliques (the principled PR #4)

Seeds = cliques from `find_cliques()` or `gene_clique_graph()`; per
species, expand the seed gene with a PageRank-Nibble or heat-kernel
push on that species' graph, stop by conductance (or OSLOM's
significance), map the expansion to HOGs, and compare expansions across
species against shuffled expression. Overlapping by construction;
covers only seeded HOGs.

- Pro: reuses the clique machinery, which is the part of the package
  that is validated (EVOTREE reproduction); milliseconds per seed; the
  cross-species comparison of per-species expansions is the statistic
  with the clean null; answers the question #4 actually asked
  ("what travels with a conserved clique") without pretending to be a
  partition.
- Con: not a module engine; no global partition, so
  `preservation_matrix_test()` gets a seed x species matrix rather than
  a module x species one; needs a small Rcpp port of ACL push or Python
  LocalGraphClustering.

### D. Joint significance-scored biclustering on the gene x sample graph (SAMBA + orthology)

Rows = (species, gene), columns = (species, sample). A joint bicluster
is a set of HOGs H and, for every species s that carries it, a copy set
from H and a sample set S_s. Its weight is the sum over species of
SAMBA's log-likelihood-ratio weight on species s's (copies, S_s)
subgraph; orthology enters as the shared H (cMonkey's phase 1), not as
coupling edges, so there is no kappa. Seeds: per-species SAMBA
bicliques (exact at degree <= 20), mapped to HOGs, joined across
species where the HOG overlap is significant; CAST-style add/remove on
the joint weight; then per-species elaboration with the core locked
(cMonkey's phase 2).

- Pro: the null is right at n = 20 and within-species significance is
  honest without a shuffled control (5.1 removed, not worked around);
  paralogs are rows and need no weighting rule; conservation is which
  species contribute a copy set, divergence is a species whose copies
  respond in no sample subset; the per-species sample sets say *when* a
  module is deployed, which is what PR #4's deployment analysis wanted
  and what no partition method gives, and they absorb tissue and time
  point without a covariate model.
- Con: no maintained software (SAMBA in EXPANDER 7.2, 2017, GUI only;
  MSCM R code from 2010); binarisation at |z| > 1 discards correlation
  magnitude; output is overlapping and non-exhaustive, so consumers
  built for partitions get a bicluster x species matrix instead of a
  module x species one; the joint weight's p-value needs the CLT plus
  Bonferroni argument redone for a sum over species (same form).
- Cost: SAMBA's core is a few hundred lines (hash the subsets of each
  gene's <= 20-condition neighbourhood, keep the k heaviest, local
  improvement); the joint step reuses `gene_clique_graph()`'s HOG
  bookkeeping, and its clean-up is Hochbaum's node-deletion biclique by
  minimum cut (5.10), for which igraph's `max_flow()` suffices. Gate: P3 unchanged (split-half replication of bicluster
  membership per species, real minus shuffled), plus SAMBA's own
  calibration against a degree-preserving random bipartite graph.

### Not carried forward

Infomap multilayer (field-of-view, incomparable objective; keep as a
fallback if A and B both fail the gate for reasons that look like
resolution); PyGenStability (scale selection is a refinement of A once
A passes); correlation-native nulls at n = 20 (Scola / RMT: the answer
is known to be "almost nothing" within species, which is borrow 1 said
differently; revisit if a compendium with n in the hundreds appears);
network alignment, ManiNetCluster, MVBC, identity-coupled multilayer
SBMs (Section 4 and 5.3).

## 8. Probe plan (the gate)

One script under `prepare_data/` (gitignored, like the other validation
scripts), leaf and root separately, all eight species. No package code
until P3 passes.

Script: `prepare_data/probe-module-engine/p1_probe_multilayer.R`
(gitignored with the rest of `prepare_data/`; args `<tissue> <outdir>
[k_top] [gene_cap] [n_cores] [engine]`). First run 2026-09-23, leaf,
k = 25, 20 000 most variable HOG-mapped genes per species, both
engines, output under `prepare_data/probe-module-engine/out-2026-09-23/`.

**P1. Build the union graph.** Node = (species, gene) for every gene in
the species' MR network (top-k sparsified, k in {25, 50}, so about
8 x 0.5 M edges; the 3 % stores are 8 x 6 M and only needed to check
that top-k did not change the answer). The script caps each species at
the 20 000 most variable HOG-mapped genes (HJUB has 44 000) so eight
dense MR builds fit a 64 GB laptop; the gene universe is fixed on the
full data so the ortholog edges are identical across halves and
shuffles. Intra-layer edges = MR edges,
weight 1 or the MR-derived weight. Inter-layer edges = every
cross-species HOG pair, weight `kappa * (1 / n_A + 1 / n_B) / 2` where
`n_A`, `n_B` are the copy numbers in the two species (OrthoClust's
down-weighting; the unweighted version is known to collapse). Write the
graph once, cache it.

**P2. Run the candidate optimisers on the same graph.** At minimum:

- multilayer modularity with per-layer null (design A), exact form:
  Python leidenalg `optimise_partition_multiplex` through reticulate
  (`py_require()` pulls leidenalg 0.12 into reticulate's uv-managed
  Python; verified 2026-09-23 on a toy two-layer graph), with one
  `RBConfigurationVertexPartition` per species holding that species'
  edges and a `CPMVertexPartition(resolution_parameter = 0)` layer
  holding the ortholog edges at layer weight kappa. Two Leiden
  iterations per run: at full size (158 049 nodes, 2.36 M intra-layer
  and 1.92 M ortholog edges) graph construction takes under 20 s and
  one optimisation takes 194 s at `n_iterations = 2` (11 joint modules
  at kappa = 1), while run-until-stable had not returned after 35
  minutes and was killed. Seed-to-seed ARI is reported per run so the
  cost of stopping at two iterations is visible. The kappa grid is
  split across two processes (`kappas` argument).
- the same objective's R stand-in, engine `cpm`:
  `igraph::cluster_leiden(objective_function = "CPM", vertex_weights =
  k_is / sqrt(2 * m_s))`, which reproduces Mucha's per-layer
  configuration null exactly for within-species pairs. It also adds a
  null term `v_i v_j` to every *cross-species* pair in a module, which
  Mucha's objective does not have. Per pair it is tiny, but it sums over
  all `n_A x n_B` pairs of two merged modules to a size-product penalty
  of about `f_s f_t sqrt(4 m_s m_t)` (f = the fraction of a layer's
  edges inside the module), so coupling only takes effect above a
  threshold kappa* ~ f x mean degree / p_ortholog. The smoke run
  (1500 genes, top-10) showed exactly that: ortholog-edge agreement
  0.01 at kappa = 1, 0.75 at kappa = 2, on shuffled data too. So the
  trick is a *different* objective (multilayer modularity with a
  cross-layer configuration null), not an approximation of OrthoClust's.
  The exact optimiser is Python leidenalg
  `optimise_partition_multiplex` (Section 5.9), which the probe script
  uses; the R trick is kept as a second engine so the effect of the
  cross-layer null is itself measurable. On the smoke run the exact
  engine's agreement rose smoothly (0.27 at kappa = 0.5 on real data
  against 0.08 shuffled) with no threshold. Checked
  2026-09-23 on igraph
  2.3.3: on one graph the trick returns the identical partition and
  modularity as `objective_function = "modularity"` (0.4707 both,
  NMI 1); on two disjoint layers of different density the global null
  merges the sparse layer into one module while the per-layer weights
  recover its three planted blocks (1/2 v 3/2 modules per layer). This
  is why fastOC's single-graph Louvain is only an approximation: the
  global `m` lets the dense species set the null for the sparse ones.
  Sweep kappa over
  {0, 0.25, 0.5, 1, 2, 4}; kappa = 0 is the current engine on the union
  graph and the baseline (agreement saturates by kappa = 4).
- the MDL candidates from Section 5 that passed the scale check, on
  the top-k graph.

**P3. Go / no-go statistic.** For each species, split the 20 samples
2 v 2 replicates within each of the 5 time points, 10 v 10 (as in the
2026-09-15 diagnostic),
rebuild both halves' networks, rebuild the union graph twice, run the
optimiser on each, and report per species: ARI between halves,
restricted to that species' genes. Baseline is the kappa = 0 column
(measured 0.06 to 0.15). Also report modules found on gene-wise
shuffled expression (baseline 8 to 11) and the module-size
distribution. **Go** if coupling raises split-half ARI materially in
most species at some kappa without the shuffled-expression module count
rising. **No-go** if not; then the joint direction is dead and no
engine choice matters, and the honest fallback is gene-level
conservation only (`compare_neighborhoods()` plus the Crow 2022 AUROC
score), with modules dropped from the trait test.

One trap in P3, found while writing this note: the ortholog edges are
identical in both halves (HOG membership does not depend on the
samples), so at large kappa the joint partition follows orthology alone
and split-half ARI rises *trivially*, on shuffled expression as much as
on real. The gate statistic is therefore split-half ARI on real minus
split-half ARI on shuffled expression, per species and per kappa, and
"go" means that gap widens with kappa before the shuffled ARI itself
climbs. For design B the equivalent is the description-length gain of
the ortholog layer on real minus on shuffled expression. Rewiring nulls
are never used anywhere in the gate (5.1).

**P4. Only after go**: kappa selection rule (Section 6), per-layer
statistics, and the consensus / tight-cluster wrapper. Then a plan for
the package change, as a separate note.

## 9. Downstream consequences in the package

If the gate passes, the change is a contraction, not an addition.

| today | after |
|---|---|
| `detect_modules()` per species, then `preservation_paired()` over `choose(8, 2)` contrasts | one joint call returning `module x species` membership; the per-species runs become the kappa = 0 special case |
| `module_correspondence()` to match modules between species | unnecessary: the joint label *is* the correspondence |
| `resolve_ortholog_map()` to choose which paralog copy carries a label | mostly unnecessary: every copy is a node and gets its own label; keep only for consumers that still need a 1:1 projection |
| `module_preservation()`'s `avg.weight` and `cor.degree` under gene-identity permutation | keep the statistics as **descriptive** per-module per-species scores (they are exactly the per-layer contribution); drop the pairwise permutation p-value, which would be circular after joint detection (the coupling already pulled the genes together) |
| `classify_preservation()` on `Zsummary_std` per contrast | a per-module *profile* over species: which layers carry the module (Russell 2023's generalist v specialist), read against the trait |
| `preservation_matrix_test()` | unchanged in mechanics; the input matrix becomes module x species instead of module-direction x contrast, which removes the all-pairs non-independence Dunn 2018 warns about |
| `tag_permutation()` | unchanged |
| K = 1 test (`test_k1`) | retire: wrong null (Section 2); replaced by the shuffled-expression control and by MDL model selection where an MDL backend is used |

The pair-level discrete-FDR item recorded here (Martin, 2026-09-24:
per-pair hypergeometric supports, with the Döhler, Durand and Roquain
step-up as the fix) is retired by the specificity score
(`compare_specificity()`, Section 11.9). Its p-values share one `1/n2`
support per direction, so there is no support heterogeneity to exploit,
and they are calibrated empirically against a shuffled partner before
Storey q-values are taken.

Circularity, stated once: with kappa > 0 a module is found partly
*because* its orthologs co-cluster, so "is module c preserved in species
s" cannot be tested by a within-species permutation of the same genes.
What can be tested is (a) whether c's within-species density in s is
above what a species-s-only partition would give (the kappa = 0
contrast), and (b) whether the per-species profile of c is associated
with the trait under relabelling. Divergence becomes "departure from a
well-supported core", which is the better-posed question.

## 10. Limits no engine lifts

- **n = 20 per species per tissue.** Coupling pools evidence for
  modules present in several species; a lineage-specific module still
  rests on one species and stays as fragile as today. Expect the engine
  to be good at conserved cores and honest about the rest.
- **Species data quality differs tenfold** (edge FDR 0.04 v 0.37).
  A species with a noisy network will look "diverged" in every module.
  Report per-species edge reliability next to every per-species module
  score; consider weighting layers by it.
- **The trait test's label space** is 2^4 under the blocked null
  whatever the module engine produces; see `pvalue_resolution()`.
- **Tissue confound**: per-tissue networks, or species x tissue layers
  with identity coupling across tissues within a species and ortholog
  coupling across species. The latter is the multilayer generality
  paying off, but it is scope for after the gate.
- **Geometric graphs.** At n = 20 every configuration-model or SBM
  null sees structure in a correlation graph (5.1), and the RMT null
  says the honest within-species answer is "almost nothing". No engine
  changes that; the shuffled-expression control and the cross-species
  statistic are the only defences and must stay in the gate whatever
  the backend.

## 11. First probe results (2026-09-23, leaf, exact engine)

Settings: leaf, top-25 MR neighbours, 20 000 most variable HOG-mapped
genes per species (158 049 nodes, 2.36 M intra-layer edges, 1.92 M
ortholog edges from 17 839 multi-species HOGs), leidenalg multiplex,
two iterations, kappa in {0, 0.25, 0.5, 1, 2, 4}, halves = 2 v 2
replicates per time point. Merged table:
`prepare_data/probe-module-engine/out-2026-09-23/gate_leaf_exact_merged.tsv`.

| species | split-half ARI, kappa 0 | kappa 1 | kappa 2 | kappa 4 | shuffled, kappa 4 | seed ARI, kappa 4 | Q_s, kappa 0 -> 4 |
|---|---|---|---|---|---|---|---|
| BDIS | 0.081 | 0.066 | 0.086 | 0.176 | 0.008 | 0.33 | 0.60 -> 0.23 |
| BMAX | 0.081 | 0.073 | 0.153 | 0.186 | 0.008 | 0.32 | 0.63 -> 0.27 |
| BMED | 0.033 | 0.024 | 0.032 | 0.155 | 0.009 | 0.32 | 0.57 -> 0.15 |
| BSYL | 0.119 | 0.159 | 0.177 | 0.188 | 0.008 | 0.34 | 0.62 -> 0.24 |
| FPRA | 0.079 | 0.088 | 0.108 | 0.193 | 0.006 | 0.33 | 0.55 -> 0.21 |
| HJUB | 0.138 | 0.164 | 0.144 | 0.208 | 0.009 | 0.39 | 0.66 -> 0.24 |
| HVUL | 0.122 | 0.135 | 0.224 | 0.184 | 0.011 | 0.32 | 0.65 -> 0.29 |
| VBRO | 0.103 | 0.116 | 0.169 | 0.191 | 0.009 | 0.35 | 0.64 -> 0.24 |

Ortholog-edge agreement on the full data, real / shuffled: kappa 0.25
0.11 / 0.06, kappa 1 0.15 / 0.09, kappa 2 0.28 / 0.96, kappa 4 0.95 /
1.00. Modules of at least 10 genes per layer: 9 to 13 at kappa <= 1, 8
at kappa 2, 7 at kappa 4 (real and shuffled alike). Split-half ARI on
shuffled expression: <= 0.005 at kappa <= 1, 0.05 to 0.07 at kappa 2,
<= 0.011 at kappa 4.

What it says:

1. **The kappa = 0 baseline reproduces the 2026-09-15 diagnostic**:
   per-species split-half ARI 0.03 to 0.14 (median 0.09) against
   about 0 on shuffled expression.
2. **Weak coupling does nothing.** At kappa <= 1 only 10 to 15 % of
   ortholog edges fall inside a module, per-layer Q and module counts
   are unchanged, and the gap moves by at most 0.01. The per-species
   partitions resist the pull.
3. **The transition is sharp and sits between kappa 2 and 4.** On
   shuffled expression the partition follows orthology already at
   kappa 2 (agreement 0.96); on real data only at kappa 4 (0.95). Real
   within-species structure resists the coupling that noise does not,
   which is itself evidence that the structure is real. Kappa 2 is the
   worst place to be: mixed regime, shuffled split-half artefact at its
   maximum (0.05 to 0.07), gap down in three species.
4. **At kappa 4 replication doubles, uniformly.** Real split-half ARI
   0.16 to 0.21 in every species (median gain +0.086 over kappa 0),
   shuffled 0.01, so the gap is 0.15 to 0.20 and is not the coupling
   artefact of Section 8's trap paragraph: on shuffled data the
   orthology-driven partition does *not* replicate between halves,
   because the ortholog graph alone has many equivalent cuts and noise
   picks different ones. What stabilises the cut on real data is the
   pooled within-species structure of eight species. BMED, the species
   with the weakest own structure (Q_s 0.15 at kappa 4), gains most
   (0.03 -> 0.16): its modules there are the other species' modules
   projected through orthology.
5. **The price**: at kappa 4 the partition is orthology-dominated, 7
   modules per layer of roughly 2 500 genes, and per-layer Q_s falls
   from about 0.6 to 0.15 to 0.29. The joint partition fits each
   species' wiring poorly. There is no kappa that gives both decent
   Q_s and improved replication at resolution 1; the objective flips
   from per-species to orthology with nothing usable in between.
6. **Half of what is left is optimiser noise.** Seed-to-seed ARI on the
   full data is 0.3 to 0.5 at every kappa (two Leiden iterations), so
   the split-half ARI is bounded by the optimiser, not only by the
   data; gap / seed ARI is 0.05 to 0.28 at kappa 0 and 0.45 to 0.56 at
   kappa 4. Consensus over seeds and subsamples (Section 6, item 5)
   is the obvious next lever and should raise both numbers.

Verdict against the P3 rule: **conditional go**. Coupling raises
split-half replication materially in all eight species at kappa 4
without the shuffled artefact rising, so the joint direction is not
dead. But what replicates is a coarse shared partition, not
per-species modules, at ARI 0.18 in absolute terms, and it is bought
by giving up within-species fit. The result supports design A as a
way to find conserved *cores* and argues against reading its
per-species Q_s as divergence. It does not lift the n = 20 ceiling on
within-species structure, which is exactly what 5.1 predicted, and it
does not change the case for design D, whose null is placed where
this ceiling does not apply.

Next steps in order: (a) consensus over seeds and sample subsamples at
kappa 4 with the null calibration of item 5, to see how much of the
0.18 is optimiser noise; (b) a resolution sweep at kappa 4 (gamma
below 1 was not tried; 7 modules is the resolution limit of modularity
on a 4 M-edge graph); (c) done, Section 11.1; (d) done, Section 11.2; (e) the design D probe, which needs
SAMBA's scoring reimplemented and is independent of all of the above.

### 11.1 The cpm engine: what the cross-layer null term does

Same graph, same halves, `igraph::cluster_leiden` CPM with vertex
weights `k_is / sqrt(2 m_s)`, run until stable. Merged table:
`gate_leaf_cpm_merged.tsv`.

| kappa | agreement real / shuffled | modules >= 10 per layer | Q_s (range) | split-half real (range) | split-half shuffled (range) | seed ARI |
|---|---|---|---|---|---|---|
| 0 to 1 | 0.000 / 0.000 | 9 to 16 | 0.56 to 0.66 | 0.03 to 0.16 | <= 0.005 | 0.36 to 0.63 |
| 2 | 0.009 / 0.362 | 10 to 15 | 0.54 to 0.66 | 0.03 to 0.15 | ~0.01 | 0.36 to 0.76 |
| 4 | 0.976 / 0.869 | **548 to 648** | 0.08 to 0.13 | 0.28 to 0.40 | **0.09 to 0.23** | 0.66 to 0.72 |

Three things the exact engine does not show:

1. **No coupling at all below kappa 2.** Agreement is exactly zero up
   to kappa 1: the size-product penalty of the cross-layer null term
   blocks every cross-species merge, as the smoke run predicted.
2. **When coupling arrives it fragments.** At kappa 4 the partition has
   about 600 modules of at least 10 genes per layer instead of 7. The
   penalty grows with the product of the two module sizes, so the
   optimum under this objective is many small cross-species modules,
   which is the ortholog graph's fine structure (HOG families), not
   co-expression modules. Q_s is 0.08 to 0.13.
3. **The coupling artefact in full.** Split-half ARI on shuffled
   expression is 0.09 to 0.23 at kappa 4, against 0.01 for the exact
   engine, because a fine orthology-dictated partition *is*
   reproducible between halves whatever the expression. The real minus
   shuffled gap (median gain +0.088) happens to match the exact
   engine's (+0.086), but it sits on top of an artefact and describes a
   different object.

Conclusion: the R trick is retired as an engine. It stays in the
script only as the measurement above. Any future R-native
implementation must carry the per-layer null without a cross-layer
term, which means leidenalg's multiplex bookkeeping in C++ (Section
5.9), not a vertex-weight trick.

### 11.2 Root (exact engine, same settings)

Root networks: 2.03 M ortholog edges from 18 107 multi-species HOGs,
2 v 2 halves per time point, run 16:46 to 18:52. Merged table:
`gate_root_exact_merged.tsv`.

| species | split-half ARI, kappa 0 | kappa 2 | kappa 4 | shuffled, kappa 4 | seed ARI, kappa 4 | Q_s, kappa 0 -> 4 |
|---|---|---|---|---|---|---|
| BDIS | 0.097 | 0.111 | 0.223 | 0.004 | 0.51 | 0.62 -> 0.35 |
| BMAX | 0.078 | 0.120 | 0.220 | 0.003 | 0.50 | 0.61 -> 0.30 |
| BMED | 0.027 | 0.033 | 0.187 | 0.004 | 0.49 | 0.63 -> 0.15 |
| BSYL | 0.087 | 0.121 | 0.224 | 0.003 | 0.50 | 0.70 -> 0.31 |
| FPRA | 0.061 | 0.078 | 0.210 | 0.002 | 0.50 | 0.57 -> 0.25 |
| HJUB | 0.131 | 0.238 | 0.210 | 0.007 | 0.50 | 0.66 -> 0.35 |
| HVUL | 0.029 | 0.053 | 0.215 | 0.003 | 0.50 | 0.67 -> 0.32 |
| VBRO | 0.161 | 0.136 | 0.219 | 0.003 | 0.50 | 0.67 -> 0.27 |

Agreement real / shuffled: kappa 1 0.13 / 0.06, kappa 2 0.28 / 0.83,
kappa 4 0.93 / 0.98. Modules >= 10 per layer: 10 to 14 at kappa <= 1,
7 to 9 at kappa 2, 6 at kappa 4. Median gap change over kappa 0: 0.00
up to kappa 1, +0.014 at kappa 2, **+0.140 at kappa 4** (leaf +0.086).

Same shape as leaf, three differences worth recording:

1. **The gain is larger and cleaner.** Gap 0.18 to 0.22 in every
   species at kappa 4, shuffled split-half <= 0.007 at every kappa
   (leaf had a 0.05 to 0.07 bump at kappa 2; root's is <= 0.009), and
   the transition on shuffled data is less abrupt (agreement 0.83 at
   kappa 2 against 0.96 in leaf).
2. **The weakest species gain most, again by projection.** HVUL's own
   root network replicates at 0.03 (leaf: 0.12) and reaches 0.22 under
   coupling; BMED 0.03 -> 0.19. Their Q_s at kappa 4 (0.15, 0.32) says
   how little of their own wiring the joint modules use. Species data
   quality differs by tissue as much as by species.
3. **Optimiser noise is lower** (seed ARI 0.50 at kappa 4 against 0.33
   in leaf), so gap / seed ARI is 0.37 to 0.44: closer to the ceiling,
   same conclusion that consensus is the next lever.

Verdict unchanged: conditional go for design A on both tissues, with
the same caveats (coarse orthology-driven partition, 6 modules per
layer, within-species Q_s not to be read as divergence). Two tissues
agreeing on the shape of the curve rules out the leaf result being a
property of one dataset.

### 11.3 Consensus over seeds (Orion array, 2026-09-23 evening)

600 runs on Orion (job 1368015 plus two reruns): both tissues, all six
conditions, {kappa 0, 2, 4} at gamma 1 and kappa 4 at gamma 0.5 and
0.25, ten seeds each, two Leiden iterations, 1.1 min per run at 200
concurrent. Analysis `p4_consensus.R`: co-classification frequency
f_ij over seeds on the union-graph edges; consensus partition = Leiden
on edges with f_ij >= 0.5 weighted by f_ij; per-node stability s_i =
mean f_ij over the node's intra-layer edges; stable core = s_i above
the 95th percentile of s_i on the shuffled full data; split-half ARI
between the consensus partitions of the two halves. Tables:
`consensus_leaf.tsv`, `consensus_root.tsv`, per-node
`consensus_nodes_<tissue>.tsv.gz`.

Split-half ARI per species, gamma 1, real (shuffled in brackets):

| tissue | statistic | kappa 0 | kappa 2 | kappa 4 |
|---|---|---|---|---|
| leaf | single seed | 0.03 to 0.15 (<= 0.005) | 0.04 to 0.20 (0.03 to 0.05) | 0.12 to 0.16 (0.01) |
| leaf | consensus, all nodes | 0.01 to 0.25 (0.00 to 0.48, degenerate) | 0.04 to 0.28 (0.04 to 0.05) | **0.22 to 0.28** (0.03 to 0.04) |
| leaf | consensus, stable core | 0.02 to 0.38, cores 30 to 61 % | 0.04 to 0.29, cores 96 to 99 % | **0.24 to 0.35**, cores 38 to 63 % (shuffled 0.02 to 0.06) |
| root | single seed | 0.02 to 0.15 (<= 0.012) | 0.02 to 0.27 (<= 0.02) | 0.23 to 0.26 (0.01) |
| root | consensus, all nodes | 0.02 to 0.19 (0.00 to 0.11) | 0.02 to 0.29 (~0) | **0.25 to 0.27** (0.03 to 0.05) |
| root | consensus, stable core | **0.69 to 0.90 in BDIS, BMAX, FPRA, HJUB, VBRO; 0.11 to 0.22 in BMED, BSYL, HVUL**; cores 0.6 to 6 % (shuffled cores empty or ~0) | 0.44 to 0.78 / -0.03 to 0.33, cores 5 to 14 % | **0.28 to 0.36**, cores 28 to 66 % (shuffled ~0) |

Seed-to-seed ARI on the full real data: leaf 0.40 to 0.60 at kappa 0,
0.48 to 0.51 at kappa 4; root 0.34 to 0.63 at kappa 0, 0.24 to 0.25 at
kappa 4. On shuffled data at kappa 4 it is 0.03 (leaf) and 0.08
(root). Gamma 0.25 and 0.5 at kappa 4 collapse every layer into one
giant module in both tissues (single-seed ARI undefined, consensus ARI
0.87 on shuffled data from the giant cluster alone): resolution below
1 is out.

What it says:

1. **Optimiser noise explained part of leaf, little of root.** Ten-seed
   consensus at kappa 4 lifts leaf from about 0.14 to 0.25 and root
   from 0.24 to 0.27. Both tissues end at ARI 0.25 to 0.28, which is
   the data ceiling of the joint partition on this design, not the
   optimiser's. Section 11's caveat "half of what is left is optimiser
   noise" was right for leaf and wrong for root.
2. **Stable cores at kappa 4 are large and only modestly better.**
   Co-classification stability separates real from shuffled cleanly
   (real cores 28 to 66 % of genes, shuffled 5 % by construction and
   with ARI ~0), but restricting to the core raises ARI only to 0.24 to
   0.36. The joint modules are reproducible as a whole at that level,
   not as a small hard nucleus plus noise.
3. **At kappa 0 in root there *is* a hard nucleus.** The 95th percentile
   of shuffled stability saturates at 1.0, so the core is the set of
   genes whose every intra-layer edge is co-classified in all ten
   seeds: 0.6 to 6 % of genes (about 340 in BDIS). Between the halves
   those cores replicate at ARI 0.69 to 0.90 in five species, and the
   shuffled cores are empty or at ~0. The three species where the real
   core fails (BMED, BSYL, HVUL, ARI 0.11 to 0.22) are exactly the three
   whose *shuffled* root networks carry large seed-stable cores (5, 39,
   22 % of genes), i.e. whose data behave most like noise; their
   per-species baseline was also the lowest (0.02 to 0.03). This is
   Tseng and Wong's tight clustering doing what it promises: a small
   reproducible set per species, the rest unassigned. Leaf does not
   show it (threshold 0.74, cores 30 to 61 %, core ARI 0.02 to 0.38).
4. **Consensus partitions are finer than single-seed ones.** At kappa 4
   the consensus graph at tau 0.5 splits the 6 to 7 modules per layer
   into 100 to 170 pieces per species (400 to 550 on shuffled data), so
   the ARI above compares fine partitions. At kappa 0 the shuffled
   consensus collapses to one giant module, which makes its ARI (up to
   0.48 in BMAX leaf) a giant-cluster artefact rather than a null.

Verdict, updated: design A's joint partition replicates at about 0.25
after consensus in both tissues and no cheaper lever remains
(resolution below 1 collapses, more seeds do not help). The
per-species tight cores in root are the more interesting object: a
reproducible nucleus of a few hundred genes per species that exists
without any coupling, in the species whose data are clean enough. That
is design C/D territory (a seed set, not a partition), and it argues
for the next probe being D (SAMBA-style scoring on the gene x sample
graph), with sample-subsample consensus (the A/B halves generalised)
as the replication test.

Next steps, revised: (a) done, this section; (b) done, negative;
(c), (d) done; (e) done, Section 11.4; (f) design B on the cached
graph; (g) sample-subsample consensus (Tseng and Wong proper) for the
root kappa 0 cores, to see whether the nucleus survives more than one
split; (h) the unit decision at the end of 11.4.

### 11.4 Design D probe (2026-09-23 evening, both tissues)

Scripts: `p5_bicluster.R` (sample-subset biclusters),
`p6_timepoint_biclusters.R` and `p6b_timepoint_auroc.R` (time-point
variant). Same halves and shuffles as the module probe. Martin: "Begin
with D, then continue with D."

**Sample-subset biclusters (SAMBA as published, applied to 20 samples).**
Gene x sample response graph (|z| > z_thr), Tanay log-likelihood
weights against a Chung-Lu degree null, every sample subset of size 3
to 6 enumerated exactly (60 249 on the full data, 582 on a half),
gene set = positive summed weight (loose) or a complete biclique
(strict), top 50 per sign after an overlap filter, significance
against the shuffled subsets of the same size.

| variant | full data, real | split-half best Jaccard, real (shuffled) | cross-species HOG Jaccard, real (shuffled) |
|---|---|---|---|
| loose, z > 1 (leaf) | 50 per sign, median 1 500 genes, all above the shuffled maximum | 0.05 to 0.11 (0.03); none at >= 0.5 | 0.10 (0.06) |
| complete, z > 1 (leaf) | 50 per sign, median 215 genes, all above the shuffled maximum | **0.02 to 0.03 (0.02 to 0.03)**: chance | 0.03 (0.02) |
| complete, z > 1.5 (leaf) | 14 to 31 per sign, median 18 genes, 90 to 100 % above the shuffled maximum; 346 of 348 are triples of single samples mixing time points | nothing survives on a 10-sample half, real or shuffled | 0.02 (no shuffled biclusters) |

Verdict on the unit: dead on this design. Every real bicluster is
"significant" against the degree null, yet the objects do not
replicate between halves at all, because a subset of three to six
*individual samples* out of twenty is a sample-specific object (the
strict biclusters are mostly co-outlier triples), and a half has
different samples. SAMBA's null was placed correctly (5.1); what fails
is the condition side being twenty samples rather than hundreds of
conditions. The root strict runs were stopped unrun: two leaf variants
at chance and one with no objects settle it.

**Time-point response sets (the design-aware unit).** Condition = time
point; a gene responds at t when |z| > 1 in at least 3 of 4 replicates
(2 of 2 on a half); biclusters = complete bicliques over the 31
time-point subsets. Because time points are shared by design, halves
and species are compared for the *same* condition.

| | leaf | root |
|---|---|---|
| genes responding per single time point, real (shuffled), median over species | T1 899, T2 717, T3 296, T4 272, T5 788 (161 to 169) | T1 636, T2 196, T3 152, T4 338, T5 432 (125 to 130) |
| multi-time-point sets | empty | empty |
| split-half Jaccard, single time points, real (shuffled) | 0.01 to 0.38 (0.01); T1 and T5 best | 0.00 to 0.41 (0.01) |
| **split-half AUROC** (genes called in one half ranked by the other half's mean z), real (shuffled) | **0.72 to 0.91 at T1, 0.67 to 0.89 at T5**, 0.38 to 0.86 at T3; median 0.77 / 0.79 down / up (0.50) | 0.66 to 0.91 at T1; median 0.72 / 0.73 (0.50) |
| cross-species HOG Jaccard, single time points, real (shuffled) | 0.05 at T1 and T5, 0.01 to 0.04 elsewhere (0.01) | 0.05 at T1 (0.01) |
| **cross-species AUROC** (HOGs called in species 1 ranked by species 2's mean z at the same time point), real (shuffled) | 0.57 at T1, 0.56 at T5, 0.51 to 0.54 elsewhere (0.50) | 0.60 at T1, 0.56 at T2, 0.50 to 0.55 elsewhere (0.50) |
| best time point in the other species (heterochrony check) | max over tp2 0.54 to 0.58, same-time 0.51 to 0.58; matrix diagonal-dominant | max 0.55 to 0.60, same-time 0.50 to 0.60 |

What it says:

1. **Within a species the time-point programs are real and reproducible
   in rank terms.** Called sets are 2 to 6 times the shuffled size and
   the other half ranks them at AUROC about 0.8 at the ends of the time
   course. The Jaccard of two hard-thresholded two-replicate sets
   (0.1 to 0.4) understates this badly; the threshold, not the biology,
   is what fails to replicate.
2. **Across species the same programs are only weakly shared.** AUROC
   0.55 to 0.60 at T1 and T2, near 0.50 at T3 and T4, no species
   receiving above 0.57, and allowing the other species to respond at a
   different time point adds at most 0.01. So this is not
   heterochrony masking conservation; the T1 response (and T5 in leaf)
   is partly shared, the middle of the course is not.
3. The sets never combine across time points, so the "bicluster" here
   is a per-time-point differential-expression set, not a
   co-expression module. D on this design reduces to: which orthologs
   respond at the same point of the course in which species.

**Verdict on D**: the null placement is right and the within-species
half of the promise holds (AUROC 0.8 against 0.25 ARI for any
partition), but the cross-species half is weak, and the unit is a
response program, not a wiring module. It does not replace the module
engine; it is a cleaner cross-species *readout* than module
preservation: per HOG and time point, the rank of its orthologs'
response in each species, with shuffled expression as the null and no
partition anywhere. That readout is cheap (seconds), exact, and gives
the trait test a HOG x species x time-point matrix with a clean null,
which is closer to what `tag_permutation()` wants than anything the
module engine produces.

**Overall after the D probe.** Three units have now been gated on the
same halves and shuffles: joint modules (design A) replicate at ARI
0.25 after consensus and are coarse; per-species tight cores (kappa 0,
root) replicate at 0.7 to 0.9 but cover 1 to 6 % of genes and only in
the cleaner species; time-point programs replicate at AUROC 0.8 within
species and 0.55 to 0.60 across. None of the three lifts the n = 20
ceiling on cross-species inference; each measures a different thing
honestly. The next decision is which object the trait test should be
built on, and that is a scientific choice, not an engineering one:
(i) design A cores for "conserved wiring", (ii) time-point programs for
"conserved response", or (iii) both, with the gene-level
neighbourhood conservation (`compare_neighborhoods()`, the Crow AUROC)
kept as the per-gene statistic. Design B (layered SBM) and the
subsample consensus of the root cores remain unrun.

### 11.5 Anchored cliques (Martin's proposal, 2026-09-24)

Martin: co-expressolog cliques as anchors, within-species co-expression
cliques around the anchor gene as the nucleus, compared across species
through HOGs; and per anchor, a species x species matrix of the
similarity of those nuclei. Scripts `p7_anchored_cliques.R` (probe)
and `p7b_trait_readout.R` (matrix readout). Same top-25 graphs, halves
and shuffles as before; every statistic cross-species or cross-half.

Design: anchors selected on **half A only**, by a lightweight
co-expressolog test (hypergeometric overlap of two genes' neighbour-HOG
sets, single-copy HOGs, BH per species pair, q < 0.05, overlap >= 3)
and a clique of >= 4 species in the significant-pair graph. Nucleus per
species = union of the maximal cliques of size >= k containing the
anchor gene; the other top-25 neighbours are the non-clique
neighbours. Tests on **half B**: recurrence of nucleus HOGs across the
anchor's species (nucleus3 = HOGs in >= 3 of them), against the same
anchors on shuffled B and against degree-matched random anchors;
replication of the nucleus between halves; and the pairwise Jaccard of
nuclei across the anchor's species on the full data.

| | leaf, k = 4 | leaf, k = 6 | root, k = 4 |
|---|---|---|---|
| anchors (half A) | 32 (21 x 4 species, 10 x 5, 1 x 7) | 32 | 43 (25 x 4, 13 x 5, 5 x 6) |
| anchor-species pairs with no such clique, real (shuffled) | 0 (0) | 0.035 (**0.326**) | 0 (0) |
| nucleus genes per species, median | 30 | 24 | 31 |
| nucleus3 per anchor, real (random anchors; shuffled) | **4.3 (0.06; 0)**, 81 % of anchors >= 1 | 3.6 (0.00; 0), 75 % | **6.5 (0.02; 0)**, 98 % |
| nucleus3 by anchor scope, 4 / 5 / 6 or 7 species | 2.5 / 7.1 / 14 | 1.9 / 6.2 / 13 | 5.0 / 8.2 / 10 |
| nucleus HOGs in >= 2 species, real (random; shuffled) | 13.4 (0.9; 0.8) | 11.3 (0.3; 0.1) | 16.4 (0.9; 0.4) |
| within-anchor contrast, fraction of HOGs recurring in >= 2 species, clique vs non-clique | 0.128 vs 0.000 | 0.142 vs 0.013 | 0.157 vs 0.000 |
| replication A vs B, median Jaccard, real (shuffled) | 0.118 (0.000) | 0.121 (0.000) | 0.132 (0.000) |
| pairwise nucleus Jaccard on full data, median, real (non-clique; random; shuffled) | **0.121 (0; 0; 0)** | 0.129 (0; 0; 0) | **0.159 (0; 0; 0)** |
| species pairs sharing >= 3 nucleus HOGs, full data | 92 % | 90 % | 98 % |

Trait readout on the pairwise matrices (`p7b`): concordant minus
discordant mean Jaccard over between-genus pairs, aggregated over
anchors, exact relabelling nulls: leaf -0.0005 (p 0.60 free, 0.63
blocked), root +0.002 (p 0.37 free, 0.50 blocked). Nothing, at floors
of 1/70 and 1/16 with 31 to 43 anchors.

What it says:

1. **The idea works, and it is the cleanest cross-species object so
   far.** Around a conserved anchor there is a nucleus of a few HOGs
   (2 to 14, growing with the anchor's scope) that recur in at least
   three of its species on data the anchor was not selected on;
   degree-matched random genes have essentially none, and shuffled
   expression has none. Every species pair shares nucleus HOGs on the
   full data. Specificity is as good as it gets on this design.
2. **Absolute similarity is modest, as everywhere.** Nucleus Jaccard
   between species is 0.12 to 0.16 and between halves 0.12; a quarter
   of one half's nucleus is even among the other half's top-25
   neighbours. The nucleus is a reproducible *core* of a few HOGs
   inside a noisy neighbourhood, which is the same shape as the
   kappa = 0 tight cores of Section 11.3, now with a cross-species
   identity and a null.
3. **Size-4 cliques are not selective within a species; size 6 is.**
   On the top-25 graphs every gene sits in 4-cliques, real or
   shuffled (geometry, 5.1), so the "clique" at k = 4 is the dense
   neighbourhood and the within-anchor contrast is uninformative (the
   non-clique remainder is 2 HOGs per anchor). At k = 6 a third of the
   shuffled anchor-species pairs have no clique while 3.5 % of the real
   ones do, at a cost of about 20 % of the nucleus. The clique
   criterion buys within-species selectivity; the cross-species
   nulls carry the inference either way.
4. **Few anchors.** 32 and 43, because the lightweight test on a
   10-sample half is conservative. The package's own co-expressolog
   calls on the full data would give hundreds, at the price of
   selecting and testing on the same samples; a proper version splits
   the design once for selection and once for testing, or uses
   subsamples.
5. **The trait readout is exactly as powered as the note predicted.**
   Per-anchor matrices feed the relabelling test directly, and with
   two labellings tied at the top of a 16-element space nothing under
   about 40 anchors of strong effect can reach p < 0.05.

Verdict: adopt as the core layer of the engine sketch (Section 10's
successor in the next note): anchors from co-expressologs, nuclei as
k-cliques around the anchor with k chosen where shuffled data lose
theirs, nucleus recurrence across the anchor's species against
matched-random and shuffled nulls, per-anchor species x species
matrices as the trait-test input. It replaces design C and makes PR
#4's question answerable; the multilayer partition (design A) becomes
the frame that groups anchors, and the time-point programs (design D)
the response readout beside it.

### 11.6 What the anchored nuclei are (annotation, 2026-09-24)

Martin: "Does this surface any biology in the example data?" Nuclei
recomputed on the full data (`p7c_nucleus_members.R`, k = 4, HOGs in
>= 3 of the anchor's species), annotated through the Brachypodium
member of each HOG via NCBI Gene (RefSeq descriptions, 747 of 785
genes; barley has almost no UniProt cross-references in Ensembl,
21 of 650). Anchors were grouped by connected components of the
anchor-to-nucleus links among anchors: leaf 32 anchors in 15 groups,
root 43 in 9. Full listings: `anchored_nucleus_members_<tissue>_k4.tsv`
and `bdis_annotation_ncbi.tsv`.

The groups are canonical regulons, named here by their members:

| group (tissue) | anchors | nucleus members (Brachypodium RefSeq names) | reading |
|---|---|---|---|
| leaf 1 / root 2, 7 species | SPX5, SPX6, SQD2, NIGT1, VIP1, inorganic pyrophosphatase 1, UGPase 3, LTI6A, FadD26 | SPX1, MGD2, GDPD1/2, sn1-DAG lipase, PAP15/22/23, NPC-type PI-PLC X, U-box 33, SPX membrane protein | **phosphate-starvation response**: PHR1/SPX signalling with the VIP1 InsP8 sensor, NIGT1, and the phospholipid-to-galacto/sulfolipid remodelling enzymes; the same nucleus in both tissues |
| root 1, 8 species | expansin A17, XTH26, extensin, AGP30, peroxidase 1/7, blue copper protein, IQD14, SFH3 (COW1), PBL23 | endoglucanases, ROP-GEFs, XTH12, peroxidase 5, WRKY25, MIZU-KUSSEI 1 | **root hair / cell expansion**: ROP-GEF and SFH3 tip growth, wall loosening, hydrotropism |
| root 3 | GPAT5, BODYGUARD 3, GPAT6, peroxidase 11, cytochrome b561 DOMON, ADIPOR1 | CASP-like 1C1, GDSL lipases, KCS1, LTP, laccases, MYB93, PELPK1 | **suberin / Casparian strip** endodermal barrier, with its regulator MYB93 |
| root 4 | RAP2-3 and ERF071 (group VII ERFs), prolyl 4-hydroxylase 6, stearoyl-ACP desaturase | plant cysteine oxidases 1/2/3, RBOH, PHOS32/34, LOB42, CYP73A | **low-oxygen response**: the ERF-VII / PCO N-degron oxygen sensing module, anchored in HVUL in all five |
| root 7 | G6PD2, 6PGD2, FNR root isozyme | nitrite reductase, nitrate reductase, GS1, APR1, MDAR5, root ferredoxin, PGI | **root nitrate assimilation** with the OPPP reductant supply |
| root 8 | MTR-1-P isomerase, DEP1, DMAS1 | NAS, NAAT, ARD, MTK, MTN, YSL9, ZIF1, ZTP29, FIT | **Strategy II iron uptake**: phytosiderophore synthesis fed by the methionine salvage cycle, grass-specific |
| root 6 | LTI65 | dehydrins DHN3/4, Rab16B, Rab21, LEA6/14, DC-8, PM19L, aldose reductase | **ABA / dehydration** LEA module |
| root 5 | four uncharacterised | PTM, BIG, PIE1, LSD1-like, BRM, MED12, UPL3, PRP8 | large nuclear regulators; co-expressed by size and expression class rather than function |
| leaf 2 | PSI subunits II, VI, XI; CP26; OEE2 | LHCII, PSI-O/IV/III/V/psaK, PSII 5 kDa, OEE1, beta-amylase | **light harvesting**; anchored only in BMAX, BSYL, FPRA, BMED (and once VBRO), not in BDIS, HVUL, HJUB |
| leaf 3, 9 / root 9 | HSA32, HSP70-8, BAG6, LIFEGUARD 2 | ClpB, sHSPs, DnaJ, HOP | **heat shock**; BMAX, BSYL, FPRA, VBRO |
| leaf 4, 6 | RPL26, RPS17, RPS4 | cytosolic ribosomal proteins, EF1 subunits | cytosolic translation |
| leaf 7, 8 | plastid RPL21, RP3 | plastid ribosomal proteins | plastid translation |
| leaf 10, 12, 14 | plastocyanin, FBPase, an uncharacterised gene | RbcS, transketolase, GAPDH A, SBPase, PGK, glycine cleavage H, CP41 | **Calvin cycle / photorespiration** |
| leaf 11 | chalcone-flavonone isomerase 3 | CHS, PAL, CHI, C-glucosyltransferase, CYP93G2, MYB P | **flavone biosynthesis** with its MYB, in BDIS, BMAX, HJUB, HVUL |
| leaf 13, 15 | uncharacterised; FLA16 | PTAC10/12, tRNase Z; RLK, GAUT-like, GT | plastid transcription; wall glycosylation |

What it says:

1. **The unit finds biology, and the biology is conserved regulons.**
   Every group with more than one anchor is a known co-regulated
   program (phosphate starvation, root hair growth, suberin, hypoxia,
   nitrate assimilation, iron uptake, dehydration, light harvesting,
   heat shock), recovered from 20 samples per species by cross-species
   recurrence alone, with no annotation used in the selection. The
   phosphate module is the same nucleus in leaf and root across seven
   species, down to the InsP8 sensor and the lipid-remodelling enzymes.
2. **Anchors are not independent; they are cores of the same module.**
   The 32 and 43 anchors reduce to about a dozen programs, and within a
   program the anchors anchor each other's nuclei. That is the
   structure the engine's joint frame (design A) should organise: one
   module per program, many anchors per module.
3. **Species coverage of a program is a result in itself.** Light
   harvesting anchors only in BMAX, BSYL, FPRA and BMED; heat shock in
   BMAX, BSYL, FPRA and VBRO; hypoxia in HVUL every time; iron uptake
   rarely in the Hordeum pair. Whether these are biology (sampling
   stage, annual leaves senescing), design (which species carry a
   co-expressolog clique on half A) or data quality (edge FDR) is
   exactly what the per-species reliability and the split-selection
   version of the anchor step have to separate before any of it is
   read against the trait.
4. The nuclear-regulator group (root 5) is the reminder that
   co-expression also groups genes by expression class; a function
   annotation is not evidence that a nucleus is a regulon.

### 11.7 Wood replication (EVOTREE data, 2026-09-24 evening)

Martin: "Redo the same module detection and visualisation with the
wood data." Same anchored-clique pipeline on the Rodriguez et al. 2026
wood-formation series (aspen, birch, cherry; Norway spruce, Scots pine,
lodgepole pine), scripts `prepare_data/evotree/module-probe/w0-w4`.
Samples are tangential cryosections per tree (Asp 106 samples / 4
trees, Birch 66 / 3, Cher 82 / 3, Nor 65 / 3, Scots 84 / 3, Lodge
84 / 3); halves are odd versus even trees, so half B of spruce has 14
samples. The gradient is warped per species into the four published
section zones (phloem-cambium, expansion, secondary wall, maturation)
and binned into eight steps for the page. Parameters as on Pooideae
except `min_species = 3` (six species instead of eight) and the
lineage contrast angiosperm versus conifer with its ten unordered
3-versus-3 relabellings (floor p = 0.1).

| | wood, k = 4 |
|---|---|
| co-expressolog-lite tests on half A (significant) | 59,002 (1,827; 3.1 %) |
| anchors (clique >= 3 species) | 103 (101 of size 3) |
| anchor scope | 93 conifer-only, 6 angiosperm-only, 4 mixed |
| nucleus3 on half B: real / shuffled / random | 0.86 / 0 / 0.01 |
| anchors with >= 1 recurrent nucleus HOG on B | 44 % |
| replication A vs B, median Jaccard: real / shuffled | 0.136 / 0 |
| pairwise nucleus Jaccard (full): real / B / shuffled / random | 0.094 / 0.043 / 0 / 0 |
| programs (merged groups) | 115 (12); 76 with >= 3 HOGs |

Three things differ from Pooideae.

1. **The anchors are almost all conifer cliques.** The two pines and
   spruce share far more co-expressologs than aspen, birch and cherry,
   which belong to three orders; 93 of 103 anchors are Lodge/Nor/Scots
   cliques. So a wood program is wired in conifers by selection, and
   the angiosperm side of every comparison is unselected. The page
   states this and adds a within-lineage **shape coherence** (mean
   Pearson correlation of the binned courses between the species of a
   lineage) so a reader can see whether the angiosperm orthologs still
   form a program at all; the "coherent divergence" ranking multiplies
   the split by the smaller coherence. Half of the 62 scored programs
   sit at the expression floor p = 0.1 against 6 expected, which is
   this selection bias, not lineage biology.
2. **The nuclei are weaker but still clean.** nucleus3 is 0.86 against
   2.4 to 4.1 on Pooideae leaf and root, and only 44 % of anchors carry
   any recurrent nucleus HOG on half B; the shuffled and random
   controls are at zero as before. Half B has three trees per species,
   and in spruce 14 sections, so the replication half is thin.
3. **The one program anchored in all six species is the secondary-wall
   program** (CESA4 / CESA8 / TRM30 group, 18 HOGs), and it is the most
   conserved course on the page (|split| 0.15, p 0.5): the positive
   control the wood series was published on.

What the coherent-divergence ranking surfaces (both lineages coherent,
r >= 0.8 within each, courses differ): the AtMC9 (metacaspase 9,
xylem cell death) nucleus peaks at secondary wall in conifers and one
zone later, at maturation, in all three angiosperms; the POK2 / kinesin
(cytokinesis) nucleus stays high through expansion in angiosperms and
falls earlier in conifers; the pectinesterase / AtDMP7 nucleus and the
GUT1 (IRX10) / FLA15-16 / GPK1 nuclei show the same one-zone offsets
in opposite directions. Whether these are heterochrony or where each
species' published zone boundaries fall is not separable on this
warping; the offsets run both ways, so they are not one systematic
boundary shift. Artifact: "Wood Regulon Gradient",
https://claude.ai/artifact/83Lyauv8DC3YkieCvkRYan.

Diagnosis of the asymmetry (scratchpad `w_diag.R`, same half-A test
per species pair): the significant co-expressolog-lite pairs are
Scots-Lodge 1,016 of 6,714 shared single-copy HOGs (15 %), Nor-Scots
244 and Nor-Lodge 243 (4 %), Birch-Cher 174 (3 %), Asp-Birch 54 and
Asp-Cher 56 (1.5 %), and every cross-lineage pair 0 to 7. The conifer
anchors are the Pinus congeners plus spruce; the angiosperm trio's
weakest links run through aspen, which has 5,946 single-copy HOGs
against 8,000 to 9,200 elsewhere (salicoid WGD). And the tree split
itself is thin: the top-25 edge sets of half A and half B overlap at
Jaccard 0.04 in five species (0.015 in spruce with 14 sections on B),
0.21 in aspen with two trees a side. The half-B nucleus tests were run
against networks that barely reproduce half A's edges, which is why
nucleus3 is 0.86 here against 2.4 to 4.1 on Pooideae.

For the engine: the anchor step's coverage bias is now measured on a
second dataset and is the design's main exposure. An anchor set that
is 90 % one lineage cannot be read against that lineage's trait; the
engine needs either a per-lineage anchor quota or, better, the
within-lineage co-expression clique as an alternative anchor when the
cross-lineage co-expressolog clique is absent.

### 11.8 Version 2 of the anchored engine on both datasets (2026-09-24 evening)

Martin: "DO that for both the Pooideae and wood data", after the
improvement list in 11.7. Shared engine `engine2.R` (Pooideae drivers
`p10_regulons_v2.R` per tissue, `p11_score_v2.R`; wood driver
`w5_regulons_v2.R`; all under `prepare_data/`, gitignored). What
changed, in the order of 11.7:

1. **Species-level cross-fitting instead of sample halves.** Networks
   on all samples. Validation is leave-one-species-out: anchors are
   re-selected from the pairs that do not involve the held-out
   species, and its nucleus around the ortholog is scored by the HOGs
   it shares with the reference nuclei (HOGs in at least two of the
   anchor's scope nuclei), against the held-out species' shuffled
   network and against ten degree-matched random genes.
2. **Blocks.** Congeners count once: an anchor clique must span
   `min_blocks` blocks (Pooideae: four species over three genera; wood:
   two species over two blocks with Scots and Lodge one block, so a
   conifer anchor needs spruce plus a pine). Nulls are enumerated as
   unordered partitions where possible; on Pooideae the 16 genus flips
   keep the complement tie, so the floor is 2/16 (and 2/70 free).
3. **Best paralog pair.** Every copy pair up to five copies a side is
   tested, the minimum p is Sidak-corrected for the pairs tried, BH per
   species pair. Tests rose from 59k to 110k on wood and significant
   pairs from 1.8k to 6.5k; on Pooideae 284k / 272k tests with 12.9k /
   19.4k significant (leaf / root).
4. **Scores on member genes, paralogs averaged per HOG**, so copy
   number no longer weights a HOG; per-replicate bin means (plants on
   Pooideae, trees on wood) with a 200-fold bootstrap band and CI, and
   |split| also in units of the pooled replicate SD.
5. **Wiring at HOG level** (any copy pair is an edge) with a z against
   100 random HOG sets matched for copy number.
6. **Wood landmarks.** Each tree is re-warped on the peak of a cambium
   set (PXY, WOX4, ANT, ATHB8, LBD1, LBD4) and of a secondary-wall set
   (CESA4/7/8, COBL4, IRX9, GUT1), running mean over three sections,
   mapped to 0.31 and 0.60; all 19 trees had the cambium peak before
   the wall peak. Per-tree peaks of every program (loess over sections)
   are drawn as a strip and summarised as the lineage peak offset.
7. **Names.** 1,283 further Pooideae HOGs named by NCBI Gene on the
   Brachypodium locus tag; 754 wood HOGs without an Arabidopsis member
   named from the spruce eggNOG annotation (marked with an asterisk).

Grouping changed too: merged groups are the connected components of
**mutual** anchor-nucleus links only, and a merged program is emitted
only when its union has at most 30 HOGs. Connected components of
one-directional links had swallowed the map (a 397- and a 740-HOG
"group" on Pooideae, 3,485 on wood).

| | Pooideae leaf | Pooideae root | wood |
|---|---|---|---|
| anchors | 431 | 668 | 2,175 |
| anchors spanning both sides | 407 | 599 | 167 |
| LOSO tests (anchor x held-out species) | 2,488 | 3,887 | 6,634 |
| shared HOGs, real / shuffled / random | 5.50 / 0.05 / 0.10 | 5.04 / 0.07 / 0.10 | 0.61 / 0.00 / 0.02 |
| held-out species on the anchors' own side | 5.53 (82 % >= 1) | 5.10 (80 %) | 1.08 (43 %) |
| held-out species on the other side | 2.41, n = 27 (70 %) | 2.72, n = 100 (76 %) | 0.15, n = 3,360 (11 %) |
| weakest held-out species | BMED 3.16 | BMED 1.83 | Asp 0.36, Nor 0.38 |
| programs with >= 4 HOGs (merged) | 1,063 over both tissues (39 merged) | | 2,095 (145 merged) |

Readings.

1. **The cross-fit confirms the unit on Pooideae and qualifies it on
   wood.** A held-out grass recovers five to seven of the reference
   HOGs around the ortholog against a tenth of one for random genes,
   in 80 % of anchors, and the trait side of the anchor does not
   matter (2.4 to 2.7 shared HOGs when the held-out species is from
   the other trait). On wood a held-out tree species recovers 0.6, and
   0.15 across the lineage boundary, ten times the random level but a
   nucleus that only rarely exists on the other side. Spruce is the
   weakest conifer (0.38) and aspen the weakest angiosperm (0.36).
2. **Anchor coverage on wood is repaired but not symmetric.** 167
   anchors now span both lineages (was 4) and the angiosperm pairs
   carry 430 to 600 significant co-expressologs (was 55 to 174), but
   the Pinus pair still carries 2,702 against 900 for spruce-pine and
   the cross-lineage pairs 39 to 71. 1,014 of 1,593 scored wood
   programs sit at the blocked floor (1/4) against 398 expected, which
   is the one-sided anchoring, not lineage biology.
3. **On Pooideae the trait is not the axis that organises the
   programs.** Of 1,061 scored programs, 62 (leaf) and 73 (root) put
   the true annual/perennial labelling at the blocked floor of 2/16,
   against 133 expected under exchangeability. Fewer than expected
   means the within-genus annual-minus-perennial differences carry
   mixed signs across genera more often than not, so the trait
   grouping cancels them where a genus flip does not. Programs that do
   divide the traits exist (galactoside fucosyltransferase, ELF3-like,
   peroxidase 2, CML31, LHCII), each with a coherent course on both
   sides, but there is no excess of them.
4. **Wood divergences with both sides coherent** are of the kind 11.7
   suggested, now with tree-level peaks: AtSS2 (starch synthase 2),
   WVD2/WDL1 and APX3 nuclei are angiosperm-anchored, coherent in both
   lineages (r 0.7 to 0.97), and peak a zone apart; FLY2/FLY1 is
   anchored in aspen, spruce and lodgepole pine, coherent on both sides
   at r 0.95 / 0.91, and the angiosperms peak 0.12 axis units later
   than the conifers at the wall stage. The AtMC9 one-zone offset of
   11.7 is not among the top coherent divergences on the landmark
   axis.

Artifacts: "Pooideae Regulon Course" version 4,
https://claude.ai/artifact/ESPjqBzfs2hXrZQn1h2qxn; "Wood Regulon
Gradient" version 3, https://claude.ai/artifact/83Lyauv8DC3YkieCvkRYan.
Both pages: rank by expression, wiring, both or coherent divergence;
scope filter for anchors on both sides; hash parameters
`#rank=coh&show=merged&scope=both&ctx=root`.

For the engine, three things carry over. The anchor step should
select on all samples and validate by leaving a species out; sample
halves were the wrong holdout for data with three replicates. Anchor
scope must be reported per side and the blocks must enter the clique
rule, or a congener pair dominates. And the trait readout has to
start from the number of programs at the attainable floor against its
expectation, before any program is read.

### 11.9 Two Suresh borrows tried: rank aggregation and the specificity score (2026-09-24 night)

Martin: "Try the first two (network agg. and conservation score)".
Script `p13_agg_cons_probe.R <pooideae|wood> [n_genes] [skip_loto]
[pair_from]`, 8,000 top-variance HOG-mapped genes per species,
Spearman networks rank-standardised as in `get_coexpression_network.R`,
six species pairs per dataset. Outputs under `p13_wood/` and
`p13_pooideae/` (the first run died when macOS revoked the directory
during a write; the rerun writes to the session scratchpad and the
files were copied back).

**Aggregation.** Wood, leave one tree out: the network of the held-out
tree against a network built from the other trees either by pooling
their sections or by rank-averaging their per-tree networks (Suresh's
aggregate). Reproducibility as the mean AUROC of the training top-25
neighbourhood in the test tree's rank rows, and the top-25 Jaccard.

| species (trees) | pooled AUROC / Jaccard | aggregate AUROC / Jaccard | shuffled test |
|---|---|---|---|
| Asp (4) | 0.981 / 0.205 | 0.981 / 0.206 | 0.500 / 0.002 |
| Cher (3) | 0.929 / 0.062 | 0.925 / 0.058 | |
| Lodge (3) | 0.915 / 0.073 | 0.914 / 0.073 | |
| Scots (3) | 0.907 / 0.072 | 0.911 / 0.075 | |
| Birch (3) | 0.888 / 0.043 | 0.880 / 0.034 | |
| Nor (3) | 0.835 / 0.035 | 0.831 / 0.032 | |
| all 19 trees | 0.913 / 0.088 | 0.911 / 0.086 | |

Aggregation over trees gives nothing over pooling, and for the
cross-species step it costs calls (below). Two other things the table
says: the fourth aspen tree lifts the held-out AUROC from 0.93 to 0.98,
so sample count, not tree effects, limits the wood networks; and a
neighbourhood AUROC of 0.9 coexists with a top-25 Jaccard of 0.09, the
same lesson as the leaf diagnostics, set overlap is the wrong
reproducibility measure at this depth. On Pooideae the aggregate is
over tissues (leaf and root rank-averaged), where it did help (below).

**Conservation score.** For gene i of species 1, its top-25 neighbour
HOGs mapped to every copy in species 2 as a 0/1 row, times species 2's
row-rank matrix, gives an analytic AUROC for every species-2 gene j;
the specificity is the rank of the ortholog's AUROC among all j,
averaged over both directions; the best copy pair per HOG is kept.
Calibration: the same score against a partner network built from
per-gene sample-permuted expression, calls at the threshold where
shuffled calls are at most 5 % of real calls. Same universe and same
top-25 sets for the hypergeometric-lite test (best copy pair, Sidak,
BH q < 0.05, overlap >= 3), also run against the shuffled partner.
Means over six pairs:

| network | spec calls at 5 % FDR | hypergeometric calls (on shuffled) | Jaccard of call sets | best copy not the top-variance copy |
|---|---|---|---|---|
| Pooideae leaf | 455 | 453 (112) | 0.33 | 0.65 |
| Pooideae root | 989 | 493 (8) | 0.33 | 0.69 |
| Pooideae leaf+root aggregate | 1,137 | 465 (0) | 0.32 | 0.67 |
| Pooideae leaf+root pooled | 1,589 | 751 (0) | 0.36 | 0.66 |
| wood pooled sections | 697 | 434 (0) | 0.27 | 0.70 |
| wood tree aggregate | 578 | 380 (1) | 0.23 | 0.73 |

Per pair on wood (pooled): Scots-Lodge 1,700 against 1,423, Nor-Scots
420 against 557, Asp-Birch 882 against 255, Birch-Cher 1,117 against
344, Asp-Nor 5 against 1, Cher-Lodge 57 against 21.

Readings.

1. **The top-25 hypergeometric is anti-conservative on the n = 20 leaf
   networks** (correction, 2026-09-25: this is the probe's test on top-25
   neighbour lists; the package's `compare_neighborhoods()` at density 0.03
   is not, see 11.12): 112 of 453 calls per pair recur on a partner whose
   expression was shuffled (25 % empirical FDR), against 8 of 493 on
   root and none on the aggregates. Ties and sparse genes survive a
   permutation of samples, and the hypergeometric has no way to see
   them; the specificity score is calibrated against exactly that.
2. **At matched FDR the specificity score calls two to three times as
   many pairs** on root, on the tissue aggregate and on the angiosperm
   tree pairs, the same number on leaf, and fewer only on Nor-Scots.
   The call sets overlap at Jaccard 0.3: the hypergeometric calls are
   nearly all high-specificity (0.97), the score adds pairs whose
   overlap is small in count but specific in rank.
3. **Paralogs matter**: in two thirds of multi-copy calls the best
   pair is not the two top-variance copies, on both datasets and all
   network types.
4. **Tissue aggregation helps on Pooideae, tree aggregation does not on
   wood.** The leaf+root rank-average calls 1,137 against 455 and 989
   for the single tissues with zero false hypergeometric calls, and
   sits below the pooled 1,589 that carries the tissue axis (the
   confound in the memory note). Aggregating over units that differ in
   state removes the state axis; aggregating over replicate trees only
   throws away sample size.
5. **Cross-lineage wood pairs stay near zero** under either score
   (5 to 57 calls), so the anchor asymmetry of 11.7 is in the data,
   not in the test.

For the engine: replace the hypergeometric co-expressolog test by the
specificity score with shuffled calibration (this also retires the
discrete-FDR item in Section 9), keep best-copy pairs, and offer
rank-aggregation across contexts (tissues) but not across replicates.

### 11.10 Package implementation against the probe (2026-09-25)

`method = "specificity"` on branch `feature/specificity-score`
(`compare_specificity()`, `null_network()`, `summarize_specificity()`),
cross-checked by `prepare_data/probe-module-engine/p14_package_crosscheck.R`
and `p14b_definition_check.R` on the p13 genes (8,000 per species,
Spearman correlation, MR network at density 0.03, one null network per
species). Outputs in `out-2026-09-23/p14/`.

| pair / tissue | Spearman probe vs package | probe calls | package calls | Jaccard | analytical calls | Jaccard analytical vs package |
|---|---|---|---|---|---|---|
| BDIS-HVUL leaf | 0.780 | 788 | 624 | 0.47 | 1,532 | 0.40 |
| BDIS-HVUL root | 0.785 | 1,233 | 2,866 | 0.42 | 2,442 | 0.85 |
| FPRA-VBRO leaf | 0.774 | 475 | 996 | 0.39 | 1,738 | 0.55 |
| FPRA-VBRO root | 0.816 | 1,231 | 1,774 | 0.54 | 2,361 | 0.74 |

The gap is definitional. The probe's own score recomputed on the package's
MR network with density neighbourhoods (mean 228 genes) agrees with the
package at Spearman 0.962 and picks the same best copy pair in 96 % of
HOGs; with top-25 neighbourhoods on the same network it agrees at 0.80.
Runtime: about 1 s per pair for both directions plus both null directions
at 8 cores, after 5 to 6 s of network builds.

Two readings. At density 0.03 the package's specificity and analytical
calls overlap more with each other (Jaccard 0.40 to 0.85) than either does
with the top-25 probe, so the neighbourhood definition matters as much as
the test; the n = 20 leaf case, where the hypergeometric was
anti-conservative, is where they part most. Paralog choice is the least
stable quantity across definitions (53 to 59 % agreement on multi-copy
HOGs against the probe).

Bug found on the way, not fixed on this branch: `compute_network()`'s
`min_var` filter lets a gene that is constant within the samples through
(floating-point variance of about 1e-30), and Spearman then errors with
"sim contains NaN". It bites on tissue subsets.

### 11.11 The anchored regulons rebuilt on the rank engine (2026-09-25)

The anchor pair test of 11.8 (hypergeometric-lite on top-25 neighbour HOG
sets, q < 0.05, overlap >= 3) replaced by the package's
`find_coexpressologs(method = "rank")` on MR networks at density 0.03 (store
0.05), one `null_network()` per species, every copy pair up to five a side,
call at q < 0.10. Everything downstream (blocks, LOSO, nuclei, programs,
scores, wiring) unchanged. Scripts `engine2.R::pair_tests_rank()`,
`p10_regulons_rank.R`, `p11_score_rank.R`, `w5_regulons_rank.R`; outputs
`v3r_*`, payloads `regulon_course_v4_rank.json`,
`regulon_course_wood_v3_rank.json`. Both pages republished (Pooideae
version 5, wood version 4).

| | Pooideae leaf | Pooideae root | wood |
|---|---|---|---|
| copy pairs tested / called at q < 0.10 | 736,338 / 120,715 (16 %) | 737,759 / 226,776 (31 %) | 296,808 / 59,011 (20 %) |
| anchors (11.8 lite test) | 3,625 (431) | 6,972 (668) | 8,419 (2,175) |
| anchors spanning both sides | 3,445 | 6,718 | 1,011 (167) |
| LOSO shared HOGs, real / shuffled / random | 1.12 / 0.02 / 0.04 | 0.94 / 0.03 / 0.04 | 0.33 / 0.00 / 0.01 |
| LOSO, held-out species from the other side | 0.28 | 0.42 | 0.17 |
| programs with >= 4 HOGs | 1,698 over both tissues | | 2,709 |

Wood calls per species pair: Scots-Lodge 59 % (Lodge is mapped on the Scots
genome), Birch-Cher 32 %, spruce-pine 26 %, aspen-birch 23 %, aspen-cherry
20 %, then spruce-angiosperm 4 to 12 % and pine-angiosperm 1 to 3 %.

Readings. The rank test calls five to fifteen times as many copy pairs as
the lite test, so there are many more anchors, and each anchor's nucleus is
smaller and replicates less across held-out species (LOSO 1 against 5 on
Pooideae), while shuffled and random controls stay near zero. The extra
calls are orthologs whose partner lists are recognised better than a
shuffled partner would recognise them, which includes broad conserved
programs (a shared time or tissue axis), not only tight regulons: a call
rate of 31 % in root says the calibration is honest about the null, not that
a third of genes sit in regulons. For anchoring, the useful filter is now the
nucleus (recurrent HOGs across the scope), not the pair test. The lineage
asymmetry on wood is repaired further (1,011 cross-lineage anchors), and
spruce links to the angiosperms more than the pines do. The top wood
program is the same as under the lite test (AtSS2 nucleus, starch synthase,
maturation in angiosperms and cambium in conifers), which is the stability
check that matters.

### 11.12 Correction: the package hypergeometric at density 0.03 is not anti-conservative (2026-09-25)

The regulon ablation (scripts `abl_phase1.R`, `abl_phase2.R`, outputs in
`out-2026-09-23/ablation/`) ran the package's pair tests at q < 0.1 on MR
networks at density 0.03 and repeated each with species 2 replaced by its
`null_network()` (expression shuffled within genes), three species pairs
per dataset:

| arm | leaf calls, real / null | wood calls, real / null |
|---|---|---|
| hypergeometric, raw MR | 10,312 / 0; 10,504 / 0; 7,242 / 0 | 11,339 / 0; 14,837 / 0; 4,398 / 0 |
| hypergeometric, log MR | 10,485 / 0; 10,692 / 1; 7,377 / 0 | 11,837 / 0; 15,353 / 0; 4,459 / 0 |
| rank, raw MR | 5,847 / 0; 6,204 / 0; 963 / 0 | 5,516 / 0; 8,078 / 0; 1,535 / 0 |
| rank, log MR | 5,889 / 0; 6,088 / 0; 1,235 / 0 | 6,036 / 0; 8,719 / 0; 1,640 / 0 |

Across all copy pairs the hypergeometric calls 35 % (leaf) and 43-45 %
(wood), the rank test 16-17 % and 20-22 %; log versus raw MR moves the
counts by 1-9 %. So the anti-conservativeness of 11.9 belonged to the
top-25 neighbour-list test, whose small neighbourhoods let ties and sparse
genes dominate; at density 0.03 (about 600 neighbours) the urn test holds
against the same null. The case for the rank test rests on calibration by
construction and per-copy ranking, not on fewer false calls. Regulon-level
results of the ablation follow when phase 2 completes.

### 11.13 Regulon ablation: pair test x MR form x regulon definition (2026-09-25)

Question (Martin): are the noisy regulons of 11.11 an artefact of how
regulons are defined? Call threshold fixed at q < 0.1. Arms: A0 legacy
top-25 hypergeometric-lite; A1/A2 package hypergeometric on raw/log MR at
density 0.03; A3/A4 package rank test on raw/log MR. Definitions: D0
current (clique anchors, 4-clique nuclei on the top-25 graph, recurrence
>= 2); D1 D0 restricted to the N strongest anchors (N = smallest anchor
count across arms: 433 leaf, 2,224 wood); D2 D0 without the clique step
(recurring top-25 neighbours); D3 nuclei from the pair test's own network
(density 0.03, recurrence >= 60 % of scope). Scores: leave-one-species-out
shared HOGs (real / shuffled / random), known-program recovery
(photosynthesis and cytosolic ribosome in leaf, secondary wall in wood),
within-species wiring z. Scripts `abl_*.R`, table
`out-2026-09-23/ablation/ablation_metrics.tsv`.

Held-out shared HOGs per anchor, real (shuffled):

| | leaf D0 | leaf D1 | leaf D2 | wood D0 | wood D1 | wood D2 |
|---|---|---|---|---|---|---|
| A0 legacy lite | 5.49 (0.05) | 5.49 (0.05) | 5.32 (0.10) | 1.08 (0) | 1.08 (0) | 1.06 (0.03) |
| A1 hyper raw | 3.14 (0.03) | 3.36 (0.03) | 3.12 (0.07) | 0.84 (0) | 0.84 (0) | 0.81 (0.02) |
| A2 hyper log | 3.28 (0.03) | 3.43 (0.03) | 3.24 (0.07) | 0.92 (0) | 0.92 (0) | 0.89 (0.02) |
| A3 rank raw | 2.61 (0.02) | 2.68 (0.03) | 2.58 (0.06) | 0.50 (0) | 0.50 (0) | 0.48 (0.02) |
| A4 rank log | 2.77 (0.03) | 2.87 (0.03) | 2.73 (0.06) | 0.49 (0) | 0.49 (0) | 0.48 (0.02) |

Readings.

1. **The pair test moves replication most, and the rank test is last.**
   Legacy lite > package hypergeometric > rank, on both datasets; the
   rank test's regulons replicate at about half the legacy level in leaf
   and wood. Every arm is 25-60 times above its shuffled and random
   controls, so none is noise; the difference is effect size.
2. **The regulon definition is not the artefact, with one exception.**
   Restricting to the strongest anchors (D1) barely changes replication
   (+0.1 to +0.2) but tightens programs (leaf wiring z 17 -> 25, wood
   10 -> 14) at a small cost in control recall; dropping the clique step
   (D2) changes nothing. D3 (nuclei on the test's 600-neighbour network)
   is the exception: programs of 60-155 HOGs, real/null ratio 4-5, and
   control precision collapses (ribosome 0.015). Large neighbourhoods do
   not make regulons.
3. **MR form is a small effect**: log MR adds 0.1-0.15 shared HOGs in
   leaf, nothing in wood.
4. **Known programs are recovered by every arm under D0-D2**: ribosome
   recall 0.79-0.88 (random 0.0001), photosynthesis 0.30-0.44, secondary
   wall 0.67-1.0. Photosynthesis precision is highest for the legacy test
   (0.28 against 0.15 hypergeometric, 0.20 rank).
5. **Cross-trait replication** (held-out species from the other trait,
   leaf): legacy 2.4, rank 0.45, hypergeometric 0.15-0.2.
6. The legacy lite test is slightly anti-conservative in leaf against the
   shuffled null (0-1 % of calls), clean in wood.

Interpretation. What separates the arms is the neighbourhood scale of
anchor selection. The legacy test selects genes whose top-25 lists are
conserved and then builds nuclei on the same top-25 graph; the package
tests select on conservation of about 600-gene neighbourhoods (density
0.03) and the nuclei are built on the top-25 graph, a scale the selection
never looked at. The rank test is the most selective about broad
neighbourhood identity, which is the least related to tight nuclei. The
next test is therefore the rank test at a tight density (about 25-50
neighbours per gene, density ~0.0015-0.0025) feeding D1 nuclei: calibrated
calls at the scale regulons live at.

### 11.14 Tight neighbourhoods, fastOC borrowings, and the pipeline switch (2026-09-25/26)

**Tight-density arms** (`abl_phase1_tight.R`, same scoring as 11.13; q < 0.1):
A5 rank at density 0.0015 (~30 neighbours per gene), A6 rank at 0.0025
(~50), A7 hypergeometric at 0.0025. Held-out shared HOGs per anchor, D1
(strongest anchors, same N as 11.13):

| arm | leaf | leaf, other trait | wood | wood, other lineage | leaf wiring z |
|---|---|---|---|---|---|
| A0 legacy lite (top-25) | 5.49 | 2.41 | 1.08 | 0.12 | 29.2 |
| A1 hypergeometric, 0.03 | 3.36 | 0.19 | 0.84 | 0.03 | 25.4 |
| A3 rank, 0.03 | 2.68 | 0.45 | 0.50 | 0.06 | 23.2 |
| A5 rank, 0.0015 | 3.80 | 1.06 | 0.85 | 0.28 | 32.9 |
| A6 rank, 0.0025 | 3.67 | 1.05 | 0.81 | 0.20 | 30.0 |
| A7 hypergeometric, 0.0025 | 4.64 | 1.26 | 0.80 | 0.05 | 31.2 |

All tight arms made no calls against shuffled partners (six pairs each).
D3 (nuclei on the test's own network) becomes usable at tight density
(7-10-HOG programs, 20-30x above shuffled) where at 0.03 it collapsed.
Known-program recovery is unchanged (A5 D1: photosynthesis 0.39, ribosome
0.85, secondary wall 0.87). The scale of anchor selection was the artefact
of 11.13: selecting at the scale the nuclei are built at recovers most of
the gap to the legacy test with calibrated calls.

**fastOC borrowings** (Zinkgraf et al. 2018, 2020; github.com/mzinkgraf/fastOC,
OrthoClust on Louvain with top-5 kNN graphs, ortholog weight
(1/cA + 1/cB)/2, co-appearance over 100 runs; `abl_phase2b.R`):
1. Copy-number weighting as an anchor-ranking penalty: multi-copy share of
   the top anchors falls from 65 % to 6 % (wood) and 52 % to 1 % (leaf),
   but replication falls 0.85 -> 0.49 (wood) and 3.80 -> 1.64 (leaf) and
   secondary-wall recall 0.87 -> 0.40. The best anchors are multi-copy
   families. Rejected.
2. Tight neighbourhoods: adopted (above).
3. Continuous membership score (fraction of scope nuclei holding a HOG):
   ranks control members above other members at chance (AUC wood 0.55,
   leaf photosynthesis 0.59, ribosome 0.51). Rejected.

**Pipeline switch.** `engine2.R::pair_tests_rank(density = 0.0015)`,
drivers `p10_regulons_tight.R`, `p11_score_tight.R`, `w5_regulons_tight.R`
(outputs `v4t_*`, payloads `regulon_course_v5_tight.json`,
`regulon_course_wood_v4_tight.json`); definition D0 unchanged. Held-out
replication over all anchors (loso_test): Pooideae leaf 1.12 -> 2.34,
root 0.94 -> 1.69; wood 0.33 -> 0.26 (other lineage 0.17 -> 0.12). Wood
gains on the matched top-500 sample (11.13 protocol) but not over all
anchors: the tight test calls 5,628 anchors there, more weak ones than the
top-500 sample sees. Pages republished (Pooideae version 6, wood version 5).
Next candidates: D1-style anchor restriction in the pipeline for wood, and
the hypergeometric at tight density as an alternative anchor test.

## 12. Sources and provenance

- Two literature surveys run 2026-09-23 by subagents in this session,
  one on the bioinformatics side (Section 4) and one on the
  network-science side (Section 5). Each citation was checked against
  Crossref, arXiv, publisher page, PMC or the software repository
  unless marked "(abstract only)" or "unverified" in place. Items the
  surveys could not find or confirm are listed as such rather than
  omitted.
- Local measurements: the 2026-09-15 leaf-only diagnostics (Section 2;
  scripts were in a session scratchpad, not the repo, numbers recorded
  in the maintainer's memory notes) and the 2026-09-23 igraph 2.3.3
  check of the CPM `vertex_weights` trick (Section 8, P2; planted SBM,
  600 and 1200 nodes).
- PR #4 read from `origin/experiment/clique-module-deployment` at
  `e6f32aa` (2026-09-23).
- Data shape from `prepare_data/data/*_se.rds`: 8 species x 40 samples
  (38 for VBRO), 14 211 to 16 144 HOGs per species, 16 to 30 % of HOGs
  multi-copy (HJUB 86 %), maximum copies per HOG 78 to 329.
- Corrections made while verifying: PyGenStability is ACM TOMS 2024
  (Algorithm 1044), not JOSS; Yang & Leskovec is Knowl Inf Syst 2013
  (42:181); Stanley et al. is IEEE TNSE 3:95 (2016); ManiNetCluster is
  BMC Genomics 2019, not Nat Commun.
