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
