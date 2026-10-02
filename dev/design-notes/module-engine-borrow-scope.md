# Module engine: what to borrow (scope, 2026-10-02)

Companion to `module-engine-redesign.md` (the survey and probe record,
sections 1-12). That note's verdict stands: the engine's problem is the
unit and the null, not the optimiser. This note records a second, wider
crawl -- fastOC and OrthoClust read from source, Breschi et al. 2016, ten
package and concept crawls (Opus agents, source read, software verified by
fetching on 2026-10-02), an Asta/OpenAlex literature pass for 2019-2026 --
and turns it into gated work packages. Everything marked *sim* was run on
simulated data only, *data* on the Pooideae leaf/root expression, and
*untested* on neither. Probe scripts the agents left are copied to
`prepare_data/probe-module-engine/borrow-scope-2026-10-02/` (gitignored).

## 1. Verdict

Three things change relative to `module-engine-redesign.md` section 1.

1. **One shared cross-species test for any engine.** A module is a gene
   set from species A, whatever produced it (per-species Leiden, the
   multilayer objective, a time-program set, a clique nucleus). Its
   conservation in species B is the neighbour-voting AUROC of its
   translated HOG set in B's real network against degree- and
   copy-matched random HOG sets (Section 4, EGAD / CoCoCoNet). The null
   is valid by construction because B's noise is independent of A's
   detection. Per-species z-scores replace `Zsummary_std` in
   `preservation_matrix_test()`. No within-species module p-value
   survives. Work package 1.
2. **The unit has two layers.** Time explains R² 0.34-0.54 of expression
   per gene (*data*), and PC1 is the time axis in HVUL (R² 0.95). But
   about 75 % of cross-species gene-level conservation survives
   residualising time (neighbour AUROC 0.63-0.69 raw, 0.59-0.66 after
   time, 0.50 shuffled; *data*). So modules are not only "the time
   axis": a *deployment* layer (the five time-point means) and a
   *wiring* layer (correlation of the within-time residuals, 15 df) are
   separate objects with separate nulls. Work package 0.
3. **Kappa-free cross-species engines exist, and kappa has a criterion
   where it stays.** The HOG-pair summary-count graph with a
   Poisson-binomial null (CODENSE / MULE / BioNERO's consensus idea) puts
   paralog copy number into the null probability instead of an edge
   weight; subspace agreement between species' spectral subspaces is a
   partition-free statistic. Where the OrthoClust objective is kept
   (design A), star expansion removes the copy-number weight formula and
   the supra-Laplacian lambda_2(kappa) transition gives kappa a
   data-internal choice. Work packages 2, 3, 5.

One reopening. Section 5.1 of the redesign note closed the within-species
question. A crawl of geometric-graph theory found a statistic that
separated a spherical random geometric graph from the same graph with
planted modules at n = 20 (*sim*): whiten the gene vectors by Tyler's
M-estimator, compute per-edge excess common neighbours given edge length
and degrees, cluster the z > 4 residual graph with HCS. It is conditional
(iid genes under the null; modules smaller than about p/19; not run on
mutual rank or real data) and it does not replace the cross-species test.
Work package 4.

Borrowed from fastOC (Section 2): the multilayer object, the
co-appearance consensus, and the lesson that the per-species
`dynamicTreeCut` on co-appearance is where its modules come from. Not
borrowed: top-5 unweighted kNN layers, Louvain, copy-number edge weights
(rejected on data, redesign note 11.14), the absence of a null.

## 2. fastOC and OrthoClust, read from source

`github.com/mzinkgraf/fastOC` (R, v0.99.2, 814 lines in `R/functions.R`,
GPL-3, dormant; chapter: Zinkgraf, Groover & Filkov 2018, CCIS 940:3,
doi:10.1007/978-3-030-00825-3_1; application: Zinkgraf et al. 2020, New
Phytol 228:1811, 13 tree species, 291,375 genes, 4.6 M edges).

| step | fastOC | source |
|---|---|---|
| per-species layer | Pearson on RPKM after a variance filter; each gene's top 5 neighbours; unweighted, symmetrised | `getEdgelist()`, `weighted2rankList()` (`functions.R:255-350`) |
| ortholog layer | every ortholog pair from a pairwise table; weight `(1/c_A + 1/c_B)/2 * couple_const`, where c is the gene's number of partners in that table | `getOrthoWeights()` (`:375-409`) |
| optimiser | `igraph::cluster_louvain` on the merged edge list, 100 runs, edge order shuffled per run | `louvain()` (`:433-475`) |
| consensus | drop communities with < 10 members; per species, co-appearance = `tcrossprod(occurrence) / nRuns`; average-linkage `hclust` on `1 - coappearance`; `cutreeDynamic` per species with a hand-set `minModuleSize` and `cutHeight` | `filterCommunityAssign()`, `multiSppHclust()`, `multiSppModules()` (`:42-76`, `:707-762`) |
| conserved vs lineage-specific | read off the cross-species blocks of the co-appearance heatmap | `plot_MultiSpp()` (`:552-615`) |
| null | none | -- |

Three consequences for rcomplex.

- The Louvain runs are a seed ensemble on one global modularity null
  (igraph's), not OrthoClust's per-layer null. The redesign note (11.1,
  11.9) measured what that cross-layer null term does: agreement exactly
  0 up to kappa 1, then HOG-family fragmentation. fastOC's "100 runs plus
  tree cut" hides this because the modules are re-derived per species
  from co-appearance, so the global partition only has to be *stable*,
  not right. Leiden on the exact multiplex objective (leidenalg
  `optimise_partition_multiplex`, or an R port of its bookkeeping) is the
  replacement; `igraph::cluster_leiden` with the CPM vertex-weight trick
  is not (11.1).
- The co-appearance consensus is the useful part, and rcomplex already
  has the machinery (`build_sparse_coclassification_cpp()`, Jeub 2018
  excess subtraction). What fastOC lacks and the redesign note already
  asked for (borrow 5) is a null-calibrated threshold: keep gene pairs
  whose co-appearance exceeds the co-appearance on shuffled expression.
- The ortholog weight was tested on our data as an anchor penalty and
  rejected (11.14: replication 3.80 -> 1.64 on leaf). Two replacements
  from this crawl: star expansion (one auxiliary HOG node, weight kappa
  per copy; the penalty is the number of defecting copies, the optimum
  placement is the majority vote `resolve_ortholog_map()` already takes;
  Section 5) and expression-learned weights pruned relative to the best
  copy (SAMap / EPIPHITES, Section 6; circular unless recomputed inside
  the shuffled null).

Two facts the redesign note's Section 4 has wrong. An R OrthoClust
exists (`OrthoClust_1.0.tar.gz` inside
`github.com/LiLabAtVT/CompareTranscriptome`, Plant Methods 2019,
doi:10.1186/s13007-019-0440-x; two species, greedy spin flips, per-layer
null); its `get_energy()` (89 lines) is an exact two-species oracle for
testing an rcomplex implementation of the objective. And no Leiden-based
OrthoClust successor was found anywhere (OpenAlex: 54 citers, about 20
relevant, none a multi-species method; Europe PMC full text, "multilayer
modularity" AND ortholog AND coexpression, 2019-2026: 0 hits), so design
A would be a first implementation.

## 3. Breschi et al. 2016 (Genome Biol 17:151)

Per-gene additive linear model `y = mu + organ + species + e` on
log-cRPKM, one sample per organ x species, 1:1 orthologs, fractions of
the sum of squares; SVG/TVG = organ + species >= 0.75 with a twofold
ratio. Their point: whether transcriptomes cluster by organ or by species
"may not be a property of the transcriptomes, but rather a consequence of
the dominant behavior of a subset of genes".

What transfers: within-species, per gene copy, R²(time) and
R²(harvest day | time), closed-form OLS, under a second for 20k x 20.
Carry them into every module and co-expressolog table, and label a
module by its mean R²(time) as "conserved deployment" or "conserved
within-time wiring". What does not transfer: the species term. Our
expression is VST per species on separate genomes with paralog
aggregation, so a species fraction would be technical.

## 4. Package crawl

Source read, not abstracts. Licences matter: rcomplex is MIT, so GPL
packages contribute ideas only.

| package | verdict | what exactly | plugs into | effort |
|---|---|---|---|---|
| BioNERO (Bioc 1.21, GPL-3) | no as code; idea | collapses paralogs to orthogroup median (`exp_genes2orthogroups`); dense TOM (about 50 GB for 8 x 20k); consensus branch for > 3 sets errors (`do.call` on a list of TOMs; tests cover 2 sets); vignette's own rice-maize run preserved no module. Idea: consensus graph as an elementwise low quantile across species, rebuilt on HOG pairs with an analytic null (Section 5, summary-count graph) | WP 2 | -- |
| multiWGCNA (Bioc 1.9, GPL-3) | partial | one species, no orthology, WGCNA. Borrow its null's *design*: split halves balanced within each time point (2 of 4 reps), modules re-detected on one half and scored on the other, size-matched module as reference. For us: a per-species **replication ceiling** that cross-species scores are read against; normalises the tenfold edge-FDR gap between species. Also `bidirectionalBestMatches()` as a reciprocal-best flag for `module_correspondence()` | WP 1, 6 | 1-2 d; flag 1 h |
| CoSIA (Bioc) | no | six model animals hard-wired, no networks or modules | -- | -- |
| EGAD (Bioc, GPL-2) | yes, reimplement | `neighbor_voting()`: score = `W 1_train / degree`, 3-fold CV over the set, analytic rank-sum AUROC; "degree null" = AUROC of ranking by degree alone (hub diagnostic, not a p-value). Calibration (*sim*, 4,000 iid genes, n = 20): random-set AUROC mean 0.50, SD 0.09 / 0.06 / 0.03 at set size 20 / 50 / 200 -- 1.4x the Mann-Whitney SD, so the analytic p is anti-conservative; naive leave-one-out is biased low (0.38 at size 20). About 15 lines | WP 1 | 1 d |
| MetaNeighbor (Bioc, MIT) | partial | transposable: genes as cells, species as studies, modules as cell types, votes from `W_B`. `one_vs_best` + reciprocal hits + connected components = cross-species meta-modules; null = permute A's labels | replaces `module_correspondence()` | 1-2 d |
| CoCoCoNet / Crow 2022 | idea | gene-set conservation = EGAD AUROC on the set in A and its ortholog set in B + degree AUROC; 1:1 OrthoDB, no null. An unpublished Gillis-lab notebook (`Family_level_coexpression`) does it N:M on Arabidopsis WGCNA modules into 9 plant species, "real" = mean AUROC > 0.75, no null | WP 1 | in EGAD |
| lionessR (Bioc, MIT) | algebra only | `e_q = N(r_all - r_-q) + r_-q`; no null; dense n² x samples (about 70 GB at 20k x 20). At n = 20 the per-sample edge is 0.88-0.94 correlated with `z_iq z_jq` and 83 % of its energy is rank one (*sim*, two independent checks), so a single-sample network is the outer product of that sample's z-profile. Useful identity: `r = (1/N) sum_q z_iq z_jq` splits r into per-time-point contributions for free | WP 0 | hours |
| PRANA (CRAN, GPL-3) | idea | jackknife pseudo-values of a connectivity statistic regressed on covariates. Package: ARACNE per group, infeasible at 20k, and an indexing bug (`rnaseqdatB[-j, ]` with global `j`) makes group-B pseudo-values wrong. Idea: pseudo-values of per-module `avg.weight` over 20 leave-one-out rebuilds as a covariate-adjusted divergence score | after WP 1 | 2-3 d |
| DCENt (arXiv 2605.30577; `samozm/DCENt`, MIT) | idea | mixed model with correlated random slopes per node and subject; needs repeated measures on one subject and a dense 2k x 2k covariance (12.8 GB at 20k). Idea: split r into between-time and within-time components | WP 0 | -- |
| CoReg (two unrelated tools) | no | LiLab 2018 TF-target Jaccard; Lee/Pan/Chen 2026 dense-peeling factor covariates | -- | -- |
| sva `sva_network` | no | blind PC removal; `num.sv` asks for 5 PCs on all three species tested (*data*), which deletes the time axis; `dat - colMeans(dat)` recycles down columns (mis-centred) | -- | -- |
| RUVcorr (Bioc) | no | no design matrix; k and nu by eye; unmaintained (2015) | -- | -- |
| variancePartition (Bioc 1.42) | idea | `get_prediction(fit, ~ 0 + (1 | plant))` equals BLUP removal exactly; for a fixed factor the residualisation is one OLS projection, so write it in base R rather than import lme4/pbkrtest. Fractions `R²(time)`, `R²(day | time)` likewise | WP 0 | 10-20 lines |
| Cote 2022; Parsana 2019 | lesson | "no form of data correction substantially improved the accuracy of co-expression networks as compared to unadjusted data" (Cote); known-covariate ~ RUVcorr ~ none; PC/PEER over-correct; Parsana: PC removal lowers FDR, "no improvement on false negative rates", and warns against it when a designed axis exists | protocol | -- |
| permute (CRAN) | vocabulary | `how(within = Within("free"), blocks = tp)`; `numPerms()` for floors. Trap: `shuffleSet(..., check = TRUE)` truncates to the enumeration (1023 rows instead of 20,000 for tp x day blocks); use `check = FALSE`. One permutation applied to all genes leaves cor unchanged; network nulls permute each gene independently. Two lines of base R, no dependency | WP 6 | hours |
| permuco (CRAN) | no | tests a factor on a per-observation response; an edge or module statistic is one number per network. Only use: eigengene ~ time point | -- | -- |
| rmcorr; misty; mlVAR (CRAN) | rmcorr idea | rmcorr = Pearson on time-point-centred data, df = N - k - 1 = 14. The "within-plant" correlation does not exist here (one sample per plant per tissue). misty and mlVAR unidentifiable at 5 clusters / no series | WP 0 | -- |

Design facts established by the crawl (*data*, `prepare_data/data/*_se.rds`):
harvest is destructive, so within a tissue plant = sample and a plant
random effect is not identifiable; each time point is split 2 + 2 over
two harvest days; 34-50 % of genes respond to harvest day beyond time,
and in FPRA the top residual PC after time has R² 0.88 with harvest day
(leaf-root plant residual r 0.18 vs 0.01 permuted; BDIS, HVUL about
0.04). VBRO has 19 plants (T3 = 3). Residuals of one gene are correlated
-1/3 within a time point after OLS on time, so the shuffle null must
permute within time point; a time-residual network has 15 df (null r SD
0.258 vs 0.229 raw) and its split halves 5 df, so the gate for residual
networks has to be cross-species.

## 5. Graph-theory concepts

Four crawls: dense-subgraph and multi-network mining; bipartite,
biclique and concept lattices; hypergraphs, tensors and multilayer
objects; connectivity-defined clusters and geometric graphs. The
redesign note's Section 5 has none of these.

One structural fact runs through three of the four. With 8 species, every
HOG, clique or HOG-pair edge has a species profile with at most 2^8 = 256
patterns (3^8 with a tested-negative state). CODENSE's second-order
graph, closed frequent edge sets (gSpan after MULE's ortholog
contraction), the recurrent-heavy-subgraph tensor, maximal bicliques and
concept lattices all collapse into bookkeeping over those patterns. The
algorithms were built for 39-130 networks; here they are a `table()` over
packed bit patterns, and the whole cost is in defining the cell and its
null.

| concept | verdict | what exactly | plugs into | effort |
|---|---|---|---|---|
| Summary-count graph + Poisson-binomial null (CODENSE, Hu et al. 2005; MULE, Koyuturk 2004) | **yes** | HOG-pair edge present in species s if any copy pair is an MR edge (ortholog contraction). Under shuffled expression the species are independent and `p_s(A,B) = 1 - (1 - d_s)^(c_A c_B)` puts copy number in the null, not in a weight; "present in >= k species" is Poisson-binomial. At 20k HOGs and d = 0.03 the noise count is about 1e4 edges at k >= 4 and 250 at k >= 5 (*untested*). Within-species triangles break independence between edges, not each edge's marginal; check with shuffles. Then Leiden or an s-of-8 core on the significant edges, into `as_modules()`. Global version of the anchored nuclei (redesign 11.5), no kappa | WP 2 | 1 d |
| CODENSE second-order graph | yes | nodes = recurrence-graph edges, joined when their 8-vector species profiles correlate (binary: <= 256 classes exact; weighted: MR percentile per species); dense components within a profile group = modules wired together in the *same* species subset. The profile is the trait readout, straight into relabelling | WP 2 | with above |
| Query-anchored densest subgraph (Goldberg min-cut, anchor forced to the source side) | yes | parameter-free nucleus on the weighted recurrence subgraph; exact; `igraph::max_flow`; statistic = density gap real vs shuffled | WP 2; anchored engine | 0.5 d |
| k-truss / (r,s)-nucleus (Sariyuce) | yes, later | polynomial stand-in for k-clique nuclei: nucleus = anchor's t-truss component, t where shuffled anchors lose it (as k = 6 was chosen). igraph has `coreness` only; Rcpp triangle peeling on the dgCMatrix with sorted-vector intersection | anchored engine | 1-2 d |
| gamma-quasi-cliques, k-plexes; cross-graph quasi-cliques (Pei, Jiang & Zhang 2005) | definition | gap tolerance for one missed MR edge at n = 20; "nucleus in >= m species" is formally a frequent cross-graph quasi-clique. No public code; Rcpp branch-and-bound on ego graphs only | nucleus definition | 1-2 d |
| RHS tensor (Li et al. 2011) | partial | anchor-seeded truncated power iteration on the 8 sparse slices gives a continuous per-species deployment vector y (replaces nucleus Jaccard). Its rewiring null is invalid here | per-anchor species matrix | 1 d |
| Filtration persistence (Giusti et al. 2015) | partial | within one species `clique_persistence()` is the analytic birth time (min edge MR) and needs no sweep; the cross-species call is not monotone in density, so stability theorems do not cover `clique_threshold_sweep()`. Betti curves of the anchor ego order complex (`ripserr`) as an RGG diagnostic | docs; diagnostic | hours |
| HOG/clique x species pattern table (3-valued `s+` / `s-` / `?`) | yes | one row per *gene clique* from `gene_clique_graph()` (OR over paralogs hides the trait-specific case: copy 1 conserved in annuals + copy 2 in perennials reads "all 8"); cells from `classify_gene_cliques()`'s `missing_reason` | `conservation_pattern_table()` | 50 lines |
| Iceberg concept lattice + implications (fcaR 2.1.0 for exploration; `arules::apriori(target = "closed frequent itemsets")` in production) | yes | species-level tiers become named intents in a lookup; `partial_significant` / `differentiated` count pairs and need the 28-column species-pair context; waterfall precedence stays as lookup order. Adds support per intent, unnamed high-support species sets, implications with confidence (phylogenetic nesting check). Null: BiCM as standardiser, relabelling over same-composition intents for inference (70 free / 16 blocked; `pvalue_resolution()`) | beside `classify_gene_cliques()` | 1-2 d |
| BiCM species-pair projection (Saracco et al. 2017) | yes | margin-corrected species x species z (15 lines, Newton over <= 9 degree classes x 8 species); separates FPRA's 0.37 edge FDR from shared conservation. Standardiser only: cells are dependent (pairwise construction, phylogeny), so within-genus pairs all validate | `preservation_matrix_test()` input | 0.5 d |
| Bimax / FCA on HOG x (anchor, species) nuclei | yes, anchors | the one bipartite object where enumeration is non-trivial (about 320 columns); groups anchors into programs against shuffled nuclei. CRAN `biclust` was archived 2025-12-19 (redesign 5.10 is stale); 40 lines of Bimax | anchored engine | 1 d |
| MBEA / iMBEA / LCM-MBC; bipartite modularity; bipartite SBM | no | trivial at 8 columns; species modules = genera; igraph Louvain/Leiden use a unipartite null on bipartite input | diagnostics | -- |
| Star expansion of a HOG hyperedge | yes | one auxiliary HOG node with zero null mass, one edge of weight kappa per copy; the penalty is the number of defecting copies, the inner optimum is the majority vote; drops the copy-number formula and shrinks the coupling layer about 12x (1.9 M ortholog edges -> <= 158k). Kappa stays. Clique expansion is just the weight `1/(d - 1)` plus within-species paralog edges to drop | design A | low |
| Supra-Laplacian lambda_2(kappa) (Radicchi & Arenas 2013) | yes | lambda_2 grows with coupling then kinks at the structural transition; the redesign note's jump at kappa 2 -> 4 (11.1) looks like it. Curve real vs shuffled, RSpectra shift-invert on 158k nodes, minutes for a 20-value grid: a data-internal kappa | design A (P3) | low |
| Subspace agreement (co-regularised multi-view spectral read as a statistic; Kumar & Daume 2011) | **yes** | `S_AB = norm(U_A^T P_AB U_B)_F^2 / K`, mean squared cosine of principal angles between A's top-K Laplacian subspace and B's mapped through the row-normalised ortholog map P. Partition-free, so per-species module instability never enters; null shuffled expression; K and P's normalisation calibrated on the null. Seconds (`arma::eigs_sym` already in `src/`). Output = all-pairs species matrix, the input shape of `preservation_matrix_test()` | WP 5 | 1 d |
| Hypergraph modularity (Kaminski 2019; h-Louvain); Kumar IRMM; hypergraph SBM (Chodrow 2021; Hy-MMSBM) | no | the weight becomes a choice of `w(d,c)` plus per-size weights (the HOG-vs-2-edge weight *is* kappa) and the volume null reintroduces the cross-layer term of 11.1; IRMM pulls diverged paralogs back into family modules; the SBM treats orthology as an observation and collapses families or ignores them by edge-count ratio (a hidden kappa). `wdc = "majority"` is a cheap HOG-coherence readout real vs shuffled | readout only | -- |
| Multilayer k-core / densest (Galimberti 2017, 2020) | yes | "dense in >= s of 8" = union over `choose(8, s)` uniform-k peelings, seconds in C++; community-search variant = principled nucleus. Needs the OR-projection's copy inflation handled (mean projection or per-pair thresholds) | WP 2 | low-medium |
| Sparse symmetric NN-CP on projected HOG networks | later | per-module per-species loadings; dense `rTensor` / `nnTensor` / `tensorly` infeasible (25.6 GB), sparse MTTKRP fine | after WP 2 | medium |
| Tyler / angular-central-Gaussian whitening (latent space with observed positions) | **probe** | the shuffled null is an *isotropic* spherical RGG; real data is *anisotropic* (sample PCs). "Rejects shuffled" can mean "has PCs". Three spikes alone give transitivity 0.53 vs 0.17 (*sim*). Tyler's M-estimator (19 x 19 fixed point) whitens; plain covariance whitening fails (24k false edges). Caveat: modules large enough to form PCs (about >= p/19) become "geometry" by definition | WP 4 | low |
| Excess common neighbours (Ricci-curvature proxy) | **probe** | `CN - E[CN | r_ij, deg_i, deg_j]`, E fitted on an isotropic null draw; exact Ollivier-Ricci is one OT problem per edge (days), and on sparse RGGs curvature is a function of edge length plus triangles anyway. Sparse `A A` restricted to edges, cost sum deg² about 1e9 per species | WP 4; new kernel | medium |
| HCS (Hartuv & Shamir 2000) | **probe** | cluster = edge connectivity > n/2 (implies diameter <= 2, density >= 0.5); a definition with a size floor, can return "no clusters"; `igraph::min_cut`, 30 lines. Only on the residual graph: the raw MR graph is one dense core (max coreness 255, every node in the 60-core) | WP 4 | low |
| CLICK; Gomory-Hu; spectral eigengap; SPC; MCL | partial / no | CLICK: exact `f0(r) ∝ (1 - r²)^8` at n = 20 as a log-odds edge weight (a third `norm_method`); Gomory-Hu as the nested hierarchy of the residual graph (R igraph exposes only `igraph:::gomory_hu_tree_impl`); eigengap says K about 20 on pure noise (spherical-harmonic multiplicities 1, 19, 189: a geometry diagnostic); SPC and MCL no | -- | -- |

Random-geometric-graph theory itself gives only a negative result: no
theorem covers RGG vs RGG + planted; global statistics (transitivity,
max core) do not move (0.173 vs 0.176; 255 vs 247, *sim*); the test has
to be local and length-conditioned. One correction to redesign 5.1: its
RMT figure sets the *mean* of the 19 eigenvalues (1053) against the bulk
edge (1064); the null band is `(sqrt(p/n) ± 1)²` about 935-1064 and real
PCs are spikes far outside it, so within-species RMT detects anisotropy,
not "almost nothing".

## 6. Literature 2019-2026 (Asta + OpenAlex + Europe PMC; code fetched)

| paper | verdict | what exactly | plugs into |
|---|---|---|---|
| CroCoNet (Termeg et al., bioRxiv 2025, doi:10.1101/2025.11.18.689002; R, `Hellmann-Lab/CroCoNet`, active 2026-03) | **yes** | regulator-anchored modules (same shape as redesign 11.5) on a phylogeny-weighted consensus network; preservation between every pair of *replicate* networks within and across species; NJ tree per module; total tree length regressed on within-species diversity; conserved / diverged = outside the prediction interval, studentised-residual BH, jackknife weights, matched random-module filter. Liftoff onto one genome, no paralogs. Borrow: cross-species distance read against each module's own within-species sample-half distance -- the same fix multiWGCNA's null design points at | WP 1, 6 |
| SAMap (Tarashansky 2021 eLife) + EPIPHITES (Passalacqua & Gillis 2024 Nat Plants) | yes | homology edges re-weighted by cross-species expression correlation, pruned below `0.25 x` the gene's best edge; EPIPHITES trims many-to-many plant families to co-expression proxies. Replaces `(1/c_A + 1/c_B)/2`; circular unless recomputed inside the shuffled null | design A ortholog layer |
| Arboretum 2.x / Muscari (Roy lab; NAR 2021 six plants with OrthoFinder gene trees; Genome Biol 2021 five cichlids) | comparator | module labels propagate down the species tree through per-branch transition matrices; at a duplication node the copies draw independently, so paralogs are explicit with no weighting rule. Gaussian mixture on profiles over matched conditions (our 5 time points x 2 tissues qualify). C++/GSL; no null; k fixed. The only phylogenetic multi-species module model with explicit paralogs, and the largest gap in redesign Section 4 | external comparator |
| Melo, Pallares & Ayroles 2024 (PLoS Comput Biol; `ayroles-lab/SBM-tools`) | partial | weighted nested DC-SBM on a BH-thresholded Spearman graph, weight `2 atanh(rho)` real-normal; about 10x more blocks, many non-assortative, equally enriched. Recipe for design B; at n = 20 the Fisher-z variance is 1/17 and the geometry stays | design B |
| OrthoClust R (2014; used by CompareTranscriptome 2019 and the Hydra study Nat Commun 2025 at kappa 2, r >= 0.975, RBH) | yes | `get_energy()` as the exact two-species oracle | design A tests |
| Pembroke et al. 2021 (Genome Biol) | partial | divergence = within-species cross-study preservation minus cross-species, CIs from permuting studies; "asymmetric" (reference species). CroCoNet's idea, cruder | `preservation_paired()` |
| Mutwil kingdom-wide stress atlas (bioRxiv 2026) | diagnostic | per-species TEA-GCN + Louvain, modules matched by orthogroup overlap against a gene-to-orthogroup shuffle (the per-species design that failed here). Paralog co-localisation: whole-genome-duplicate pairs in the same module 57.6 % vs 19.3 % permuted -- an engine sanity statistic (HJUB is 86 % multi-copy) | engine probe |
| Masuda, Boyd, Garlaschelli & Mucha 2025, Phys Rep, "Introduction to correlation networks: beyond thresholding" | read | the review of the component papers Section 5.1 cites | -- |

Already in the redesign note, nothing new: BiTSC, ManiNetCluster,
Juxtapose, Ovens 2021, Russell 2023, fastOC, CoCoCoNet / Crow 2022,
Suresh 2023, GenePlexusZoo, the EVOTREE precursor preprint.

## 7. Corrections to `module-engine-redesign.md`

- 4: OrthoClust has an R implementation (two species); software is not
  "dead", only the Julia port is.
- 5.1: the RMT comparison is mean-vs-edge; the honest reading is
  "anisotropic", and Tyler whitening is the one-component null.
- 5.1 / 5.6: restricted nulls (permute within time point) keep the
  geometry -- transitivity 0.28 vs 0.17 under full shuffle (*sim*) --
  and raise the bar from "noise" to "shared time course"; they do not
  make a within-species module statistic interpretable.
- 5.10: CRAN `biclust` archived 2025-12-19.
- 8 / 11: the sample design is destructive harvest with two harvest days
  per time point; "plant" is not a repeated-measures factor and harvest
  day is a candidate batch axis (FPRA).
- Section 4's catalogue misses Arboretum / Muscari and CroCoNet.

## 8. Work packages

Gated as before: no exported function changes until its gate passes.
Each package names its files, its deliverable, its acceptance check and
its dependencies; packages without a dependency arrow can run in parallel
in worktrees (agents must commit). The ladder applies: one kernel per
package, no new object slots, no refactor off the path.

**WP 0 -- two-layer decomposition (data probe; 0.5 d).**
Script `prepare_data/probe-module-engine/p17_two_layers.R`. Per species
and tissue: OLS-residualise VST on `time_point` (fixed; base R `qr`), keep
the five time-mean profiles; `compute_network()` on residuals (wiring)
and compute per-gene `R²(time)`, `R²(real_day | time)`. Cross-species
neighbour AUROC (k = 50, single-copy HOGs, as the confound agent ran it)
on raw, wiring and deployment layers against shuffled-within-time-point.
*Acceptance:* table `two_layers_<tissue>.tsv` reproduces the 0.63-0.69 /
0.59-0.66 / 0.50 pattern on all 28 pairs, plus the harvest-day falsifier
(AUROC after also removing `real_day`; if it drops to shuffled in a
species, that species' wiring layer is flagged batch). Depends on nothing.
Blocks WP 2-4's choice of input layer.

**WP 1 -- `module_auroc()` (package; 1 d).** `R/module_auroc.R`, tests
`tests/testthat/test-module-auroc.R`, one C++ helper only if the sparse
`W %*% L` product in R is too slow (it should not be: one product per
species, L = n x (M + n_null) indicator columns). Inputs: an
`as_modules()` object from species A, the species-B network, the HOG
map. Steps as in Section 4's EGAD row: translate to all B copies, drop
within-HOG edges, 3-fold CV over the set's HOGs, rank-sum AUROC, null =
`n_null` random HOG sets matched on set size, B copy-number distribution
and B degree decile, same pipeline; return per-species AUROC, z, p, and
the degree-only AUROC as a hub flag. `seed` through `.seed_scope()`; add
the entry point to the RNG-contract table. *Acceptance:* `devtools::test()`
green, and a fixture test that modules detected on shuffled-A expression
give z with mean within ±0.1 and SD within 0.9-1.1 over B (the
calibration gate), and that split-half within A is the same statistic
with the other half as B. Depends on nothing; WP 2-4 are scored by it.

**WP 2 -- summary-count graph with Poisson-binomial null (probe, then
package; 1 d + 1 d).** Probe `p18_recurrence_graph.R`: ortholog
contraction of the eight MR stores to HOG pairs; count per pair; `p_s`
with copy number; Poisson-binomial p per pair (`poibin` or the DFT
convolution in 10 lines); significant-edge graph; CODENSE profile
grouping (256 classes) and Leiden on the significant edges; also the
anchored densest subgraph for the 11.5 anchors. Nulls: shuffled
expression per species (independence of the eight layers is the point),
and the predicted noise count. *Acceptance:* on shuffled data the number
of significant pairs at FDR 0.05 is within the Poisson-binomial
prediction (the dependence check); on real data the Leiden modules score
under WP 1 with median z > 3 in >= 4 species and the shuffled-A gate
holds. If it passes, package as `recurrence_graph()` returning an igraph
plus the edge table, feeding `as_modules()`. Depends on WP 1 (scoring)
and WP 0 (which layer).

**WP 3 -- design A with star expansion and the lambda_2 criterion (probe;
2 d).** Extend `p1_probe_multilayer.R` (exact leidenalg engine): replace
the weighted ortholog layer by star-expanded HOG nodes with zero null
mass; compute lambda_2 of the supra-Laplacian over a kappa grid on real
and shuffled; run the kappa-4 consensus with the null-calibrated
co-appearance threshold (fastOC's consensus, Jeub excess, keep pairs
above the shuffled co-appearance). *Acceptance:* the lambda_2 kink sits at
the agreement transition of 11.1; the star layer's modules score at least
as well as the weighted layer's under WP 1, with the HJUB family
fragmentation of 11.9 gone. Then decide between WP 2 and WP 3 as the
engine by WP 1 z and split-half replication; keep one. Depends on WP 1.

**WP 4 -- within-species clean-up (sim, then data probe; 2 d).**
`p19_tyler_hcs.R`: (a) 20k-gene simulation with an angular-central-
Gaussian background from the real Tyler fit of BDIS leaf, planted modules
of 20-500 genes at r 0.3-0.8, 1,500 outlier/zero-inflated genes, MR at
density 0.03; (b) BDIS leaf real vs shuffled. Pipeline: Tyler whitening ->
MR -> excess common neighbours per edge against an isotropic draw ->
z > 4 residual graph -> HCS with a size floor set where shuffled gives
zero clusters. *Acceptance:* (a) recovery vs false nodes and zero
clusters on the null; (b) split-half ARI over called genes only above the
shuffled split, and the called clusters score under WP 1. If (b) fails,
the within-species question stays closed and this note's Section 5
records why. Depends on WP 1 for scoring only.

**WP 5 -- `subspace_preservation()` (package; 1 d).** `R/subspace.R`:
top-K eigenvectors of each species' normalised Laplacian from the sparse
store (`arma::eigs_sym` path already in `src/`), row-normalised ortholog
map P from `prepare_orthologs()`, `S_AB` for all pairs, K and the P
normalisation calibrated on shuffled expression; returns the all-pairs
species matrix in the shape `preservation_matrix_test()` reads.
*Acceptance:* shuffled expression gives `S_AB` at its null level on all
28 pairs; the real matrix correlates with the WP 1 per-species z of the
chosen engine. Independent of WP 1-4; can run first.

**WP 6 -- nulls and the replication ceiling (package; 1 d).**
`null_network(block = )`: per-gene independent shuffle within block
(within time point, `check = FALSE` semantics, two lines), with
`numPerms()`-style floors reported; a split-half helper balanced within
time point (2 of 4 reps) that returns the within-species sample-half
score of size-matched modules as the reference every cross-species score
is read against (multiWGCNA's null design, CroCoNet's internal
reference). *Acceptance:* RNG-contract test passes for the new argument;
on the 11.5 anchors the ceiling is reported next to the cross-species z.
Depends on WP 1 for the score.

**WP 7 -- clique-side pattern table and lattice (package; 2 d).**
`conservation_pattern_table()` (one row per gene clique, 3-valued species
and species-pair cells from `classify_gene_cliques()`'s bookkeeping),
`arules` closed itemsets as the iceberg lattice, named intents mapped to
the existing tiers, BiCM species-pair z as a new `preservation_matrix_test()`
input. *Acceptance:* every `classify_gene_cliques()` species-level tier
is reproduced as a named intent on the Pooideae edge tables (exact
match), and the BiCM z-matrix on shuffled expression is null. Independent
of WP 0-6.

Order: WP 5 and WP 7 are independent and cheap; WP 0 and WP 1 first;
WP 2 before WP 3 (if the kappa-free engine passes, WP 3 is a comparison,
not a requirement); WP 4 last. Not scheduled: hypergraph engines,
NN-CP, truss kernel, Arboretum port, PRANA pseudo-values.

## 9. Sources

fastOC: `github.com/mzinkgraf/fastOC`; Zinkgraf, Groover & Filkov 2018,
doi:10.1007/978-3-030-00825-3_1; Zinkgraf et al. 2020, New Phytol 228:1811,
doi:10.1111/nph.16819. OrthoClust: Yan et al. 2014, Genome Biol 15:R100;
R package inside `github.com/LiLabAtVT/CompareTranscriptome`
(doi:10.1186/s13007-019-0440-x). Breschi et al. 2016, Genome Biol 17:151,
doi:10.1186/s13059-016-1008-y. Cote et al. 2022, doi:10.1186/s13059-022-02606-0;
Parsana et al. 2019, doi:10.1186/s13059-019-1700-9. Kuijjer et al. 2019,
lionessR, doi:10.1186/s12885-019-6235-7 (method: iScience 14:226).
PRANA: doi:10.1186/s12864-023-09787-3. DCENt: arXiv 2605.30577. EGAD:
Ballouz et al. 2017; CoCoCoNet: Lee et al. 2020, doi:10.1093/nar/gkaa348;
Crow et al. 2022, doi:10.1093/nar/gkac276; Passalacqua & Gillis 2024,
doi:10.1038/s41477-024-01738-4. CODENSE: Hu et al. 2005, Bioinformatics
21:i213, doi:10.1093/bioinformatics/bti1049. MULE: Koyuturk et al. 2004/2006.
Cross-graph quasi-cliques: Pei, Jiang & Zhang 2005 KDD; Jiang & Pei 2009,
doi:10.1145/1460797.1460799. RHS: Li et al. 2011, doi:10.1371/journal.pcbi.1001106.
Nucleus decomposition: Sariyuce et al. 2015/2017. Giusti et al. 2015,
doi:10.1073/pnas.1506407112. Multilayer cores: Galimberti et al. 2017,
doi:10.1145/3132847.3132993; 2020, doi:10.1145/3369872. BiCM: Saracco et al.
2017, New J Phys 19:053022. Hypergraph modularity: Kaminski et al. 2019,
doi:10.1371/journal.pone.0224307; Kumar et al. 2020, doi:10.1007/s41109-020-00300-3;
Chodrow, Veldt & Benson 2021, doi:10.1126/sciadv.abh1303. Supra-Laplacian:
Radicchi & Arenas 2013, doi:10.1038/nphys2761; De Domenico et al. 2013,
doi:10.1103/PhysRevX.3.041022. Co-regularised spectral: Kumar, Rai & Daume
2011, NeurIPS. HCS: Hartuv & Shamir 2000, IPL 76:175. CLICK: Sharan &
Shamir 2000, ISMB. Ricci flow: Ni et al. 2019, Sci Rep 9:9984. Latent
space: Hoff, Raftery & Handcock 2002, JASA 97:1090. CroCoNet:
doi:10.1101/2025.11.18.689002. SAMap: doi:10.7554/eLife.66747. Arboretum
applications: doi:10.1093/nar/gkaa1041; doi:10.1186/s13059-020-02208-8.
Melo et al. 2024, doi:10.1371/journal.pcbi.1012300. Pembroke et al. 2021,
doi:10.1186/s13059-020-02257-z. Masuda et al. 2025, doi:10.1016/j.physrep.2025.06.002.
Asta Scientific Corpus MCP: `asta-tools.allen.ai/mcp/v1`.
