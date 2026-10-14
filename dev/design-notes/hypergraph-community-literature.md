# Hypergraph community detection: literature for WP15 (2026-10-07)

Companion to `sharpen-plan.md` WP15 / WP15a / WP15b and to
`module-engine-borrow-scope.md` (sections 2, 5 and 10.10, which already probed
star expansion and the multiplex coupling on Pooideae data). This note surveys
hypergraph modularity, its algorithms, nulls and validation, and multilayer
coupling, and says what each finding means for the plan. No package code
changed. Method and failures are at the end.

## 0. Verdict

1. The plan's `q(alpha)` loses co-expression at alpha = 1 or gives it a
   species-blind null, and scores the HOG 2-section with a global null term
   (risks 1 and 3).
2. Hypergraph modularity beats the 2-section by small margins, at moderate
   noise, and only for hyperedges of size >= 3 (Kaminski 2024).
3. No C++ hypergraph Leiden was found; `h-louvain` has no licence file. At
   n = 3,000 Bethe Hessian spectral beats h-Louvain (Li 2026).
4. The field validates against planted partitions, not by permutation.
5. Which HOGs "split" depends on tau, not only on the data (Li 2026).
6. In HSBM terms unweighted HOGs add well under 1 % of the signal against
   about 600 co-expression edges per gene. WP15b's port should wait on a
   three-engine probe (section 7). Peixoto, Peel & Gross 2026 put the burden
   of proof on the hypergraph (section 8).
7. Probe outcome (sections 10-11): E1 fails on wood (no coupling window);
   E2, the species-blocked objective, couples the species but ties the
   two-step partition matched by hog overlap on planted data and on wood
   split halves. Recommendation: withdraw WP15b, ship a matching step, do
   not build E3.

## 1. Hypergraph modularity definitions, resolution, degeneracy, nulls

**Kaminski family.** Kaminski, Poulin, Pralat, Szufel & Theberge 2019 (PLOS
ONE 14:e0224307) define strict and general hypergraph modularity against a
Chung-Lu hypergraph null: a hyperedge of size d falls inside part A with
probability `(vol(A)/vol(V))^d`. Kaminski, Misiorek, Pralat & Theberge 2024 (J
Complex Netw 12:cnae041; arXiv:2406.17556) generalise this to weights
`eta_{c,d}` for a hyperedge of size d with c > d/2 members in one part. The
tau-family sets `eta = (c/d)^tau`: majority (tau = 0), linear (1), quadratic
(2), strict (tau -> Inf). Their default is tau = 2 because a d-edge with c
members in one part contributes `c(c-1)/(d(d-1)) ~ (c/d)^2` of its weight to
the 2-section modularity. The difference from the 2-section is that a
hyperedge with c <= d/2 in every part counts for nobody, and a hyperedge
counts for at most one part. They note the usual resolution parameter gamma on
the degree tax carries over unchanged, and that values of different modularity
functions must not be compared across tau.

**All-or-nothing (AON).** Chodrow, Veldt & Benson 2021 (Sci Adv 7:eabh1303)
derive hypergraph modularity as the likelihood of a degree-corrected
hypergraph SBM (DCHSBM). In the AON case a hyperedge counts only if all its
nodes share a cluster, which is Kaminski's strict variant with size-specific
resolution parameters. Their generalised Louvain optimises it; AON beats
clique-expansion Louvain when hyperedges are mostly pure and loses when they
are not.

**Splitting functions.** Veldt, Benson & Kleinberg 2022 (SIAM Rev 64:650)
treat every such choice as a splitting function: the penalty a hyperedge pays
depends on how it is split. Clique expansion, star expansion and AON are
special cases; not every splitting function is graph-reducible. For WP15 this
frames the whole decision: star expansion (WP15a) is a linear penalty in the
number of defecting copies; tau-modularity is a step at c = d/2 followed by
`(c/d)^tau`.

**Reductions.** Clique-expansion weight `w/C(d,2)` keeps total weight,
`w/(d-1)` keeps degree (the plan's choice; Kumar et al. 2020, Appl Netw Sci,
who also reweight by hyperedge balance, IRMM).

**Unified comparison.** Poda & Matias 2024 (Peer Community J 4:e37) put these
variants in one framework and benchmark the codes: the objectives answer
different questions.

**Resolution limit and degeneracy.** No hypergraph-specific resolution limit
theorem was found. The graph results carry over because every tau-modularity
is a sum of a coverage term and a global volume null: Fortunato & Barthelemy
2007 (PNAS 104:36) for the limit, Good, de Montjoye & Clauset 2010 (Phys Rev E
81:046106) for the plateau of near-optimal, mutually dissimilar partitions.
Traag, Van Dooren & Nesterov 2011 (Phys Rev E 84:016114) show CPM is
resolution-limit-free; no hypergraph CPM was found. Vidiella, Duran-Nebreda &
Valverde 2026 (Commun Phys 9:198) is an ecology perspective with no new
method. It repeats that modularity has resolution limits and degeneracy,
points to inference (Peixoto 2023; Contisciani 2022) instead, and shows
nestedness and modularity arising from different projections of one
hypergraph: the projection decides the structure you find.

**Nulls.** Chodrow 2020 (J Complex Netw 8:cnaa018) gives configuration models
for hypergraphs (vertex- and stub-labelled, MCMC sampling) and shows that
answers depend on whether you randomise the hypergraph or its projection. The
Chung-Lu null of the Kaminski family is species-blind: it lets any hyperedge
member land on any gene. A HOG's species composition is fixed, so the right
null for an ortholog hyperedge draws each member from its own species' volume
(a product of per-species probabilities). The redesign note (11.1) measured
what a global null on cross-species terms does to the multiplex objective. For
d >= 3 the strict tax is tiny (`sum_A p_A^d`), so this matters mostly for d =
2 and for majority/linear tau.

## 2. Algorithms and implementations

**Local move.** h-Louvain (Kaminski 2024) is Louvain on the blend `alpha q_H +
(1 - alpha) q_{H[2]}`, alpha rising by the `(p_b, p_c)` schedule, tuned by
Bayesian optimisation over 10 seeds. No Leiden refinement. Chodrow 2021 and
Veldt's HyperModularity.jl have a generalised Louvain for AON /
splitting-function objectives. HyperNetX `last_step()` is a
hypergraph-modularity local-move refinement run after a 2-section Kumar or
Louvain start. Leiden itself (Traag, Waltman & van Eck 2019, Sci Rep 9:5233)
guarantees gamma-connected communities on graphs; no paper found states it for
a hypergraph objective.

**Spectral and inference.** Chodrow, Eikmeier & Haddock 2023 (SIAM J Math Data
Sci 5:251) give nonbacktracking spectral clustering for nonuniform
hypergraphs. Li, Schaub & Peel 2026 (Sci Adv 12:eaef2184; arXiv:2601.10502)
derive a Bethe Hessian for nonuniform HSBMs with model selection by negative
eigenvalues. Ruggeri, Lonardi & De Bacco 2024 (J Stat Mech 043403;
arXiv:2312.00708) give message-passing detectability bounds; Ruggeri et al.
2023 (Sci Adv 9:eadg9159) scale a mixed-membership HSBM (Hy-MMSBM) to large
hypergraphs; Contisciani, Battiston & De Bacco 2022 (Nat Commun 13:7229) fit
overlapping communities and predict hyperedges (Hypergraph-MT). del Genio 2024
(arXiv:2412.06935) defines a spectral hypermodularity matrix. Kovacs, Benedek
& Palla 2025 (Sci Rep 15:36032) adapt clique percolation to hyperedges.

**Software (checked 2026-10-07).**

| code | language | licence | objective / algorithm |
|---|---|---|---|
| `pawelwm/h-louvain` | Python (hypernetx 2.3.5) | **none** | tau-modularity, Louvain + BO; last commit 2024-06 |
| HyperNetX (pnnl) | Python | PNNL "Other" | `modularity()`, `kumar()`, `last_step()` |
| SimpleHypergraphs.jl | Julia | MIT | Kaminski modularity, Julia original of h-Louvain |
| HyperModularity.jl | Julia | MIT | Chodrow/Veldt AON and splitting functions, Louvain |
| hypergraphx | Python | BSD-style | Hy-MMSBM, Hypergraph-MT, spectral |
| Hypergraph-MT | Python | GPL-3 | mixed-membership EM |
| graph-tool | C++/Python | LGPL | SBM on the bipartite incidence graph only |
| KaHyPar / Mt-KaHyPar | C++ | GPL-3 / MIT | balanced partitioning (cut, km1), not modularity |
| libleidenalg | C++ | GPL-3 | Leiden on graphs and multiplex layers; no hypergraph |
| igraph C, NetworKit | C / C++ | GPL-2 / MIT | no hypergraph community detection |
| XGI | Python | BSD-3 | no community detection |
| CRAN HyperG, Bioc `hypergraph` | R | GPL-2 / Artistic | no modularity |
| CRAN multinet | R/C++ | Apache-2 | multilayer clustering, not hypergraph |

`euxhenh/hypergraph-leiden` is a single-cell notebook project running graph
Leiden. **Verdict: the plan's claim stands.** No C++ hypergraph Leiden was
found, and no CRAN or Bioconductor package implements hypergraph modularity.
Both claims rest on absence of evidence (GitHub and CRAN searches by name and
keyword).

**Licence consequence.** `h_louvain.py` carries no licence, which means all
rights reserved. WP15b may call it from a dev-only probe, but its
"bookkeeping" must be re-derived from the paper (section 4.4 and appendix
pseudo-code), not translated from the file. libleidenalg (GPL-3) is the only
code the plan should port, as G9b already says.

## 3. When does the hypergraph beat the 2-section, and lift-off

**Evidence for.** Kaminski 2024: on h-ABCD (n = 300-1,000, sizes 2-8),
h-Louvain at tau 2-3 beats 2-section Louvain and Kumar in AMI, but the gain is
concentrated at moderate noise; at low noise all methods agree and at high
noise all fail. On real data the gains are slight (primary-school contacts:
small tau best; cora co-citation: larger tau best). Strict modularity loses
when many community hyperedges are impure. Their rule of thumb: cluster the
2-section first, tabulate hyperedge purity (c, d), and choose tau from it.
Chodrow 2021: in their experiments AON wins when hyperedges are pure, clique
expansion otherwise. Ruggeri 2024 compare a hypergraph's entropy with its
clique expansion's: detection gets easier when hyperedges overlap heavily on
node pairs.

**Evidence against or neutral.** Pellegrin, Fesser & Weber 2025
(arXiv:2502.09570, a seed) find that graph-level models on hypergraph
expansions often outperform hypergraph-level models, in learning tasks. Li,
Schaub & Peel 2026 argue that every hypergraph method must choose which orders
and split shapes to keep, and show their spectral method keeps higher-order
hyperedges whole and prefers near-homogeneous splits (3-1 over 2-2), even
where lower-order edges support another partition as strongly. Kovacs 2025
find hyperedge and projected clique percolation differ.

**Lift-off.** Kaminski 2024: from singletons no single move completes a
hyperedge of size >= 4, so the edge term stays 0 while the tax grows; with
mixed sizes, small edges dominate early merges. Their fix is the alpha
schedule `alpha_i = 1 - (1 - p_b)^(i-1)`, switched when the community count
first falls to `n p_c^(i-1)`. The best `(p_b, p_c)` depends on the data and on
tau; avoid both near 0 or both near 1; good settings have `p_b + p_c ~ 1`.
They tuned it by Bayesian optimisation because no default was reliable. In
WP15 the lift-off is different in kind: HOG hyperedges are disjoint
(hyperdegree 1), so q_H never supplies any merge signal between HOGs;
co-expression does all merging at every alpha. The schedule only decides when
copies start to be pulled toward their HOG majority.

## 4. Multilayer coupling as the alpha = 0 limit; cross-species work

Mucha, Richardson, Macon, Porter & Onnela 2010 (Science 328:876) add
inter-layer edges of weight omega between copies of a node, with no null term
on them. Bazzi et al. 2016 (Multiscale Model Simul 14:1) analyse how omega
trades intra-layer modularity against persistence of labels across layers,
with decoupled layers at small omega and one shared partition at large omega.
Pamfil, Howison, Lambiotte & Porter 2019 (SIAM J Math Data Sci 1:667) show
multilayer modularity equals maximum likelihood in a multilayer SBM, with
omega set by a layer-persistence parameter, so omega can be estimated rather
than swept. Weir et al. 2017 (CHAMP, Algorithms 10:93) prune the parameter
plane to the partitions optimal somewhere. This matches what the borrow-scope
note measured (10.10): coupling is a phase transition in kappa, star expansion
shifts it about 8x, and in the coupled regime co-localisation is 0.99, so a
star module cannot show subfunctionalisation.

The plan's degree-preserving 2-section of HOG hyperedges (`w/(d-1)`) is the
multilayer coupling for paralog-aware ortholog maps; for 1:1 orthologs between
two species it has the same edges as Mucha's categorical coupling, but the
plan scores it as a modularity, with a null term Mucha does not have.
Tau-modularity differs in how defection is charged: a copy that leaves its
HOG's majority costs `(c/d)^tau - ((c-1)/d)^tau` while c - 1 > d/2, then the
whole HOG term at the majority line. Star expansion charges kappa per
defecting copy. So tau = 2 with large d tolerates one or two defecting
paralogs cheaply; strict does not.

Cross-species, ortholog-anchored detection: OrthoClust (Yan et al. 2014,
Genome Biol 15:R100) is the two-species multilayer objective with a per-layer
null; fastOC (Zinkgraf et al. 2020, New Phytol 228:1811) runs Louvain on the
merged graph without one. The borrow-scope note (section 2) found no
Leiden-based or hypergraph successor among 54 OrthoClust citers. Eriksson et
al. 2021 (Commun Phys 4:133) compare unipartite, bipartite and multilayer
random-walk representations of the same hypergraph for Infomap; the
representation changes the number, size and overlap of communities; the
bipartite form is the analogue of a gene-HOG incidence graph. Russell et al.
2023 (PLoS Comput Biol 19:e1011616) run multilayer community detection on
co-expression networks of several tissues, the closest biological analogue of
species layers; it is still pairwise coupling of identical genes. The Asta
citation walk over 390 papers citing the five seeds found no cross-species,
ortholog-anchored or hypergraph network-alignment community method. Feng et
al. 2023 (Proc ACM Manag Data 1:4, 10.1145/3617335) build a modularity null
(PIC) for hyperedges that are not all-or-nothing; it is the nearest published
alternative to the tau-family for HOGs that split.

## 5. Significance and validation

Most hypergraph papers validate against planted partitions: AMI against h-ABCD
ground truth (Kaminski 2024; h-ABCD, arXiv:2210.15009), or recovery and
detectability thresholds under an HSBM (Chodrow 2021, Ruggeri 2024, Li 2026).
Model selection, where it exists, comes from likelihood (DCHSBM), from
negative Bethe Hessian eigenvalues (Li 2026), or from description length
(Peixoto 2023, *Descriptive vs. inferential community detection in networks*,
Cambridge Elements). Musciotto, Battiston & Mantegna 2021 (Commun Phys 4:218)
filter hyperedges against an analytic null (statistically validated
hypergraphs); Young, Petri & Peixoto 2021 (Commun Phys 4:135) reconstruct
hyperedges only where the evidence supports them. Consensus over seeds
(Lancichinetti & Fortunato 2012; Jeub et al. 2018) is used on graphs, and
Grassetti & Mastrandrea 2026 (arXiv:2602.21838, STAR) pick one representative
partition out of the degenerate high-modularity set; Kaminski 2024 run 10-100
seeds and report the best or the mean, without consensus. No paper found tests
a detected hypergraph community by permutation or split-half replication.

Implication: rcomplex keeps its own validation. Joint detection couples
species by construction, so a joint module "preserved" across the species it
was detected on is circular. It must replicate across sample halves, beat
shuffled expression, and pass the held-out AUROC test (borrow-scope verdict
1).

## 6. Design implications and risks for WP15a / WP15b

1. **Keep co-expression outside the alpha blend** (Kaminski 2024, eq. 3; the
   plan's `q(alpha)`). In Kaminski's setting every 2-edge is also a hyperedge
   of q_H, so alpha = 1 still sees them. In the plan the co-expression edges
   sit in q_2 with a per-species null; if they stay only there, q_H at alpha =
   1 sees disjoint HOGs and is maximised by one module per HOG; if they are
   also put in q_H, they get the global volume null that redesign 11.1
   rejected. Write the objective as `Q = sum_s modularity_s(coexpr, gamma) +
   lambda [alpha q_H^HOG + (1 - alpha) q_2^HOG]`, so alpha changes only the
   form of the HOG penalty and lambda keeps its role as the coupling (kappa).
   The brute-force acceptance test should include a case where this matters.
2. **Lambda is omega, and it has a phase transition** (Mucha 2010; Bazzi 2016;
   Pamfil 2019; borrow-scope 10.10). `lambda = NULL` = median co-expression
   weight is a guess on an 8x-sensitive scale. Pick lambda by the borrow-scope
   criterion (smallest lambda with real agreement near 1 and shuffled
   split-half at its floor), or estimate it as a persistence parameter per
   Pamfil 2019. Sweep it in the Orion probe as planned.
3. **Do not put a global null on the HOG terms** (Mucha 2010; Kaminski 2019
   Chung-Lu; Chodrow 2020; redesign 11.1). The plan's q_2 scores the HOG
   2-section as a modularity at gamma, which is the cross-layer null term that
   11.1 measured: no coupling up to kappa 1, then fragmentation into HOG
   families. Mucha's coupling has no null, and WP15a's star expansion avoids
   one by giving the HOG nodes zero null mass. For q_H, replace `Pr(Bin(d,
   vol(A)/vol(V)) = c)` by the probability that c of the HOG's members fall in
   A when each member is drawn from its own species' volume (a
   Poisson-binomial). The species-blind null taxes placements that cannot
   happen; the effect is largest at d = 2 and for small tau.
4. **Tau decides what "split HOG" means** (Kaminski 2024 section 5.1; Li,
   Schaub & Peel 2026). Choose tau from the purity table of a WP15a run (their
   EDA rule), report subfunctionalisation candidates only when they hold over
   tau in {1, 2, Inf}, and calibrate the split count against shuffled copy
   labels within HOG (the `p_copy` idea of `resolve_ortholog_map()`). Expect
   no tau effect at all for two-species, 1:1 designs.
5. **Expect a small gain; let the gate decide** (Kaminski 2024; Poda & Matias
   2024; Pellegrin 2025; borrow-scope 10.10, where the recurrence engine beat
   every coupled multiplex run, z cross 31.4 against 23.0 on leaf). G9 should
   require WP15b to beat WP15a on split-half and cross-species replication,
   not on modularity. Do not tune `(p_b, p_c)` by Bayesian optimisation of q:
   that maximises the objective, and a higher q is not more replicable
   modules.

## 7. Objective and optimiser: is Leiden shoehorned in?

Li, Schaub & Peel 2026 (Sci Adv 12:eaef2184), read in full, define
`SNR_BH = [sum_k (k-1)(d_in^k - d_out^k)]^2 / sum_k (k-1) d^k` for a symmetric
non-uniform HSBM, where d^k is a node's mean number of size-k hyperedges and
d_in^k / d_out^k count the pure and the boundary-crossing ones. Bethe Hessian
(BH) spectral clustering, with K = number of negative eigenvalues, reaches the
belief propagation (BP) limit for one hyperedge size and falls short of it
for mixed sizes. Their supplement has BH beating h-Louvain in AMI at
n = 3,000; at n = 30,000 h-Louvain would not run. The code
(`eggplantisme/HyperGraphBetheHessian`) is notebooks without a licence. Li &
Peel 2026 (arXiv:2604.18565) show three phases above the detectability
threshold: small communities merged into large ones, then separated as one
group, then resolved. BH needs a stronger signal than BP for the last phase,
and the K rule (negative eigenvalues, free energy, MDL) moves all three
boundaries.

**(a) Maximisation or inference.** Neither family fits rcomplex as stated. The
HSBM has no hyperedge weights and no species. Co-expression edges never cross
species, so an unblocked HSBM or BH finds the species first; an inference
engine needs species-blocked degree correction (Pamfil 2019's multilayer SBM,
extended to HOGs). The layers are not sparse in SBM terms: at density 0.03 on
20k genes a gene has about 600 co-expression edges and at most one HOG. The
null the redesign note measured is a geometric graph, not a block model. For
inference: rcomplex's problem is the minority one, since conserved modules of
300-1,000 genes are 2-5 % of a layer. Inference also brings K from the
data, where modularity brings gamma and a resolution limit (Fortunato 2007).
For Leiden: it runs at 8 x 20k genes, and gamma and lambda are already
calibrated against shuffled data (borrow-scope 10.10). Leiden is not
shoehorned in, but its objective is descriptive and misspecified, so
replication stays the judge (section 5).

**(b) The SNR from rcomplex data.** d^k comes straight from the data: the
co-expression degree per species at the analysis density gives k = 2, and a
gene in a size-k HOG adds 1 to d^k. d_in and d_out need a partition; use the
existing `detect_modules()` or the WP15a partition. That gives two numbers
before any engine exists. The co-expression signal per node is
`d_in^2 - d_out^2`, on a degree of about 600. The most HOGs can add is
`sum_k (k-1) d^k` with every HOG pure, a few units. Unweighted HOGs therefore
change SNR_BH by well under 1 %, and everything rests on lambda, which lies
outside the published model. Two species with 1:1 orthologs give k = 2 only:
the hypergraph is a graph, tau-modularity is modularity, and hypergraph BH is
graph BH.

**(c) Probe (Orion, before WP15b).** Engines:

- E1: star expansion on igraph Leiden (WP15a).
- E2: corrected h-modularity (risk 1 objective, risk 3 tax), a 200-line
  probe-only serial Louvain on 4k-gene subsamples. Written from the paper,
  not from h-Louvain's code.
- E3: species-blocked BH: per-species graph BH plus the lambda-weighted HOG
  projection `A^(k)`, K from negative eigenvalues. BP via hypergraphx (BSD)
  on planted data only.

Data:

- Planted: a species-blocked h-SBM, 8 species, degree 600 and 10, 2-5 %
  minority modules, HOG sizes from the real table. Each copy sits in its
  HOG's module with probability `1 - delta`, delta in {0, 0.1, 0.3}.
  h-ABCD adds degree heterogeneity.
- Real: leaf and wood split halves, scored as in borrow-scope 10.10.

Metrics: AMI, minority best-match Jaccard, and delta-copy precision and recall
(planted); split-half ARI against shuffled, z cross on the held-out AUROC
test, and co-localisation (real).

Decision rule:

- E2 earns WP15b only if both hold: it beats E1 on delta-copy recall at equal
  AMI, and its real-data gains exceed seed noise (0.7-2.2 in z).
- E3 earns its own WP if it beats E1 on minority recovery at degree 600 and
  its K replicates across halves.
- If E1 ties both, WP15 stops at WP15a.

**(d) Recommendation: re-scope WP15b and wait for the probe.** Build WP15a
(it is E1). Hold the libleidenalg port and the GPL-3 switch it brings. The
SNR budget says unweighted HOGs add almost nothing; the only comparison at
this scale favours BH over h-Louvain (Li 2026); and the recurrence engine
already beat every coupled run (borrow-scope 10.10). Run (b) first; it takes
an afternoon. If the HOG share stays negligible at the lambda where coupling
happens (star kappa 16-32), E2 differs from E1 only in how it charges
defecting copies. The probe then tests one splitting function, which a
per-HOG cost on E1 might emulate (Veldt 2022 show this for cuts, not for
modularity).

## 8. Second seed batch

On topic, three seeds:

- Contisciani 2022 (already in section 2) reports Hypergraph-MT faster than
  its clique expansion when hyperedges are large. On gene-disease and
  contact data it predicts hyperedges better than the expansion, and on
  Congress, Walmart and Trivago it does no better.
- Kritschgau et al. 2024 (Sci Rep 14:6933) fit a degree-corrected
  microcanonical hypergraph SBM by simulated annealing. They report detection
  near the conjectured sparse thresholds but stop short of them, with K
  given. Their gains over the multi-edge clique projection are modest. MDL
  picks too few clusters (6 against 9 school classes).
- Kovacs 2025 (already in section 2).

Off topic, two seeds:

- Lotito et al. 2025 (Commun Phys 8:43): directed-hypergraph motifs and
  reciprocity. Microscale only; HOGs are undirected.
- Miyashita et al. 2025 (Sci Rep 15:20729): a clustering coefficient for
  hypergraphs. A local density measure with no communities or nulls.
- Felippe, Kirkley & Battiston 2026 (Sci Adv 12:eaec5619): a normalised
  mutual information for hypergraph similarity, covering cross-order
  overlap and node coarse-graining. It compares two hypergraphs on one
  labelled node set and does not detect communities, so there was no walk.
  One possible later use: an all-pairs species similarity after
  coarse-graining genes to HOGs, as a partition-free input to
  `preservation_matrix_test()`. Untested.

The 2024 Commun Phys collection (19 articles) has nothing on communities,
nulls or alignment.

The Asta walk from the on-topic seeds and Li & Peel 2026 returned 136 citing
papers. Read beyond the title:

- Peixoto, Peel & Gross 2026 (arXiv:2602.16937) argue that graph models
  represent multibody interactions fully, and that hypergraphs are a
  constrained special case of them. They find no evidence for the claimed
  broad advantage of hypergraphs. This backs section 7's default to E1
  unless the probe says otherwise.
- Ni, Deng & Mu 2025 (arXiv:2505.04967) give an SBM over several
  hypergraphs. It infers communities together with the edges between
  hypergraphs (their example is genes and their protein products). It is
  the nearest published model to "species layers joined by orthology", and
  a candidate E3 variant.
- Kirkley, Felippe & Malizia 2026 (arXiv:2606.00893) prune redundant
  hyperedges non-parametrically. This is a possible HOG filter before
  detection.

## 9. Hypergraph encodings and multi-copy hogs (2026-10-12)

All numbers come from the 12,000-gene Pooideae leaf table: 8 species,
1,452,683 gene cliques, 3,742 hogs. Scripts: `dev/probes/encodings-frc.R`
and `dev/probes/encodings-twins.R`.

### 9.1 The paper

Pellegrin, Fesser & Weber 2025, "Enhancing the Utility of Higher-Order
Information in Relational Learning", NeurIPS 2025 (arXiv:2502.09570; code
`github.com/Weber-GeoML/Hypergraph_Encodings`, Python 3.11, torch 2.0, Julia
for Ollivier-Ricci). They compute node-level encodings from a hypergraph and
feed them to GNNs. At hypergraph level the models are UniGCN, UniGAT, UniGIN
and UniSAGE. At graph level they are GCN and GPS on expansions.

| family | definition |
|---|---|
| H-k-RWPE | return probabilities after 1..k steps of a hypergraph walk |
| H-k-LAPE | eigenvectors of `B1^T B1` or `I - Dv^{-1/2} B1 De^{-1} B1^T Dv^{-1/2}` |
| HCP | per node, min, max, mean, median and sd of curvature over incident hyperedges |
| H-LDP | degree profile |

The walk comes in three kinds: equal-node, weighted-edge and equal-edge. The
equal-edge walk picks an incident hyperedge uniformly, then a co-member
uniformly: `P(i -> j) = 1/|{e : i in e}| * sum_{e containing i and j}
1/(|e| - 1)`. HCP uses Forman-Ricci, `F(e) = sum_{k in e} (2 - d_k)`, or
Ollivier-Ricci (Coupette, Dalleiger & Rieck 2023, ORCHID).

The main result: graph-level models on hypergraph expansions beat hypergraph
architectures. The encodings help the graph-level models and not the
hypergraph ones. Section 3 already counts the paper as evidence against or
neutral.

Verdict: do not adopt the code. The definitions are one page of sparse
products on an incidence matrix. The `members` table of `gene_cliques` is
that matrix (shape B, `clique-table-shape.md`).

### 9.2 Forman curvature of the clique-lifted hypergraph

Each hyperedge is one maximal gene clique. `d_k` is the number of cliques
holding gene k. One `tabulate()` and one grouped sum over `members` give
`F(e)`.

| measure | value |
|---|---|
| cliques with F > 0 (every member in at most two cliques) | 0.27 %; 2,898 of 3,742 hogs have at least one |
| HOG:0025026 | 1,400,795 cliques (96 % of all), median F about -1.3 million |
| multi-copy (hog, species) cells where copies differ in mean F | 1,460 of 2,653 (55 %) |

Here F is minus the clique multiplicity. A copy in fewer cliques is either
diverged or more specific. Without edge weights the sign cannot be read.
Keep F as a probe column, not a tier rule. In several hog cells, copies had
exactly the same clique count and mean F. That led to 9.3.

### 9.3 False twins: what no encoding can resolve

Every topological encoding gives identical values to nodes with identical
neighbourhoods. Same-species paralogs with the same cross-species
co-expressolog partners in the gene graph at `alpha_graph` are false twins.
They are never adjacent and share one neighbour set. They multiply the
cliques. Twin class = hash of the sorted neighbour set per gene. Cliques are
compared as multisets of twin classes.

| measure | value |
|---|---|
| multi-copy (hog, species) cells | 5,657; 17,263 copies -> 15,420 twin classes |
| cells where all copies are twins | 1,233 (22 %) |
| cells with at least one twin pair | 1,620 |
| cliques | 1,452,683 -> 145,476 distinct up to twin substitution (10x) |
| HOG:0025026 | 1,400,795 -> 120,383; copies per species 8,9,8,6,10,15,13; twin classes 8,4,7,5,5,14,10 |
| HOG:0024976 | 11,178 -> 604 |
| HOG:0031783 | 6,383 -> 108 |

The reduction is lossless. A quotient clique stands for every choice of one
representative per twin class. Its multiplicity is the product of the class
sizes. Two consequences:

1. The network cannot resolve twins, so do not enumerate them. Collapse each
   class to one node before `igraph::max_cliques()`. Report `n_copies` per
   member.
2. Non-twin copies differ in topology. They are the real resolution targets.
   There, curvature, random-walk encodings or plain edge scores apply.

Any hypergraph community method in WP15 inherits the same degeneracy.

The probe in 9.6 shows twins are twins only in the gene graph. Within their
own species the two copies' co-expression neighbourhoods differ (hog-projected
W1 median 0.52; 6 of 2,082 twin pairs below 0.01). So a twin class is an
enumeration identity, not a biological class. The copies pass the same
cross-species calls, nothing more. The reduction in cliques stays exact. The
reading "twins = redundant copies" does not hold.

### 9.4 Two ideas under probe

Results are in 9.6.

1. **Equal-edge walk on the mixed hypergraph.** Hyperedges are hogs across
   species plus within-species co-expression edges. The 4-step return is
   gene -> co-expressed neighbour -> that neighbour's hog -> ortholog's
   neighbour -> own hog -> gene. It is a per-copy co-expressolog statistic.
   The `1/(|e| - 1)` factor normalises it by copy number. It is a soft,
   multi-step form of the hypergeometric overlap. Test: does it separate
   gene-graph false twins? Their within-species neighbourhoods differ even
   when their cross-species partners do not.
2. **Ollivier-Ricci curvature of a hog hyperedge.** Use a 0/1 ground metric
   (same hog = 0). On hog-projected neighbourhood measures,
   `W1(mu_g, mu_g') = 1 - sum_h min(nu_g(h), nu_g'(h))` in closed form. Then
   `kappa(hog) = 1 - mean pairwise W1` across species. The transport plan
   between two species' copy sets is a soft copy map. It is an alternative
   to the single-layer rule in `resolve_ortholog_map()`.

### 9.5 Recommendation

WP-c: twin reduction inside `gene_clique_graph()`, the next shape step
(issue #72). Curvature and the transport copy map stay in `dev/probes/`; 9.6
shows neither earns a package column. Do not adopt the repository.

### 9.6 Probe results (2026-10-12)

Scripts: `dev/probes/encodings-common.R`, `encodings-rw.R`, `encodings-orc.R`.

Setup: 8 species x 12,000 genes = 96,000 nodes. 15,796 hog hyperedges (max
581). 17.3 M co-expression edges at density 0.03 (degree 360 everywhere).
P has 36.1 M non-zeros. 2,082 twin pairs; 66,049 non-twin paralog pairs.
Tiers at 12k: complete_conserved 4, partial_present 19, trait_specific 65,
unclassified 3,647. With n = 4, AUROCs against complete_conserved are noise.

**Equal-edge walk.** A gene's one hog edge competes with 360 co-expression
edges. The hog step has weight 1/361. The hog-mediated share of the 4-step
return is about 0.1 %.

| Q1: abs(dlog p4) | twin | non-twin paralog | random |
|---|---|---|---|
| all paths, median | 0.130 | 0.317 | 0.307 |
| hog part, median | 0.048 | 0.239 | 0.540 |
| degree gap, median | 0.116 | 0.297 | 0.258 |

| Q2: agreement with the most-called copy | h4 | p4 | highest degree | baseline |
|---|---|---|---|---|
| 9,745 cells | 0.442 | 0.554 | 0.560 | 0.410 |

p4 is degree. h4 is at baseline. Run time 93 s, 5.8 GB.

**Hog-projected Ollivier-Ricci.** kappa median 0.040 (IQR 0.032-0.054).
Neighbourhoods of different species overlap about 4 %.

| Q4: tier | n | 1 - mean W1 | 1 - max W1 | hog size |
|---|---|---|---|---|
| complete_conserved | 4 | 0.056 | 0.001 | 57 |
| partial_present | 19 | 0.129 | 0.033 | 13 |
| trait_specific | 65 | 0.061 | 0.016 | 8 |
| unclassified | 3,647 | 0.063 | 0.022 | 6 |

AUROC, any tier vs unclassified: 0.532 (mean), 0.403 (max). The max version
is a size effect: a big hog always holds one distant pair. Per copy, the
lowest mean cross-species W1 picks the most-called copy in 0.760 of 9,745
cells (baseline 0.410). Run time 341 s, 12.9 GB.

| Q5: within-species W1 | n | hog-projected median | share < 0.01 |
|---|---|---|---|
| twin | 2,082 | 0.519 | 0.003 |
| non-twin paralog | 20,000 | 0.955 | 0 |
| random | 20,000 | 0.968 | 0 |

**Copy map.** 51,562 (hog, A, B) triples have a multi-copy side. Optimal
transport (OT) by `lpSolve`, uniform marginals.

| Best-q called pair vs the plan | n | heaviest cell | argmin C | base | in OT support | base | RBP | base |
|---|---|---|---|---|---|---|---|---|
| both multi | 8,448 | 0.730 | 0.797 | 0.124 | 0.897 | 0.446 | 0.849 | 0.166 |
| one side 1 | 9,370 | 0.919 | 0.919 | 0.461 | 1 | 1 | 0.919 | 0.461 |

| `resolve_ortholog_map()` pairs vs the plan | pairs | in OT support | base | RBP on C | base |
|---|---|---|---|---|---|
| both multi | 12,056 | 0.880 | 0.446 | 0.735 | 0.166 |
| one side 1 | 9,370 | 1 | 1 | 0.898 | 0.461 |

The LP returns a vertex, so plans are sparse by construction. Argmin C
beats the heaviest OT cell (0.797 vs 0.730). OT adds nothing over
reciprocal best pairs (RBP).

**What separates copies.** Hog-projected W1 to the other species does
(argmin C, RBP, lowest mean W1: 0.76-0.85 against baselines of 0.17-0.41).
But W1 re-scores the neighbourhood overlap the hypergeometric test already
uses. It confirms the calls. It is not new evidence.

**What does not.** The equal-edge walk: one hog edge drowns among 360
co-expression edges. Curvature: kappa is about 0.04 in every tier and tracks
hog size. OT over reciprocal best pairs. Twin classes as biology: twins
differ within their species.

**Worth a package probe.** None of the scores. At most, use the lowest mean
cross-species W1 as a tie-breaker for `source == "unresolved"` copies in
`resolve_ortholog_map()`, after a held-out test. A walk first needs a weight
on the hog layer.

## 10. WP15a built; star expansion fails on wood (2026-10-14)

`joint_modules()` landed as PR #88 (engine E1 of section 7c). Two findings,
both from the species-blind null of igraph's modularity (risk 3 of section 6).

**Fixture.** On a two-species fixture with four planted modules of 40 genes
(loading 0.7), the default weight (median scaled co-expression weight, about
1) recovers each species (ARI 0.98 / 0.99) but gives no shared labels: the
null term `k_i k_j / 2m` is charged on every cross-species gene pair, so a
module of species A and its module in B only merge when the hog links pay
that tax. Coupling starts between weight 1 and 2, holds at 2 to 8, breaks at
16. At loading 0.5 the joint run is worse than `detect_modules()` on one
species (ARI 0.29 against 0.61). Tests run at weight 4.

**Wood (six species, `prepare_data/evotree/e9_snr.R`, cache
`prepare_data/cache/snr/`).** The same null from the other side. At
resolution 1 the joint `2m` is six times one species' `2m`, the null
deflates, and Leiden returns species blobs (13 modules; Birch and Cherry one
module each) at every weight. Resolution S = 6 restores the per-species null
up to the spread in layer size (ARI 0.74 against per-species Leiden). The
weight sweep at resolution 6, weight as a multiple of the median:

| x median | 1 | 4 | 16 | 32 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|---|---|---|
| joint modules | 43 | 46 | 71 | 114 | 157 | 311 | 2,242 | 6,319 |
| hogs pure | .02 | .02 | .02 | .03 | .05 | .11 | .41 | 1.00 |
| ARI vs two-step (median) | .74 | .73 | .72 | .72 | .66 | .54 | .27 | .03 |

Up to 32x nothing couples. Beyond, modules split rather than join: the
cross-species tax on merging two modules is `n_A n_B d^2 / 2m_s` (quadratic
in module size), the hog gain is `n_h w` (linear), so Leiden shrinks modules
to the size the weight can pay for. At 512x the hog link equals a gene's
whole co-expression strength (about 540 scaled units) and every gene follows
its hog. There is no weight at which per-species structure and shared labels
coexist. Section 7b's SNR: hogs are 1.5 % of `sum_k (k-1) d^k` (9.4 of 632);
under the per-species partitions they are mostly impure, so unweighted hogs
subtract from the numerator. Leiden on the joint graph takes 25 s; the load
and the edge extraction take the minutes.

**What it means.** E1 is not a tie, it is a failure: star expansion on a
global-null optimiser cannot give joint modules on this data, and section
7c's rule ("if E1 ties both, stop at WP15a") does not apply. The only
candidates left are objectives with a per-species co-expression null and a
hog term with its own null (risk 1's `Q = sum_s modularity_s + lambda q_HOG`,
E2) or the species-blocked spectral route (E3). Neither exists in igraph, so
the probe must build E2 before WP15b can be judged, as section 7d said. WP15a
stays as the alpha = 0 reference and as a working two-step module detector
with shared labels at weight 0. The default `weight` stays the documented
guess; the user docs say to sweep it.

## 11. E2 probe: the species-blocked objective ties the matching step (2026-10-14)

Engine E2 of section 7c, built as a probe (`dev/probes/joint_louvain.{cpp,R}`,
`joint_louvain_check.R`, `joint_planted.R`, `joint_wood_halves.R`; outputs
and `REPORT.md` under `prepare_data/cache/probe_e2/`, gitignored).

**Objective.** Risk 1's form with risk 3's blocked null at tau = 2:

```
Q = sum_s (1 / 2m_s) sum_{i,j in s} (A_ij - gamma k_i k_j / 2m_s) delta(c_i, c_j)
  + (lambda / 2m_bar) sum_{i in s, j in t, s != t} (w_ij - h_i^t h_j^s / H_st) delta(c_i, c_j)
```

Co-expression edges are `.joint_edges()` scaled per species; each species
has its own `2m_s`; `2m_bar` is their mean. A hog contributes every
cross-species pair of copies at weight `w_ij = 1 / (c_s c_t)` (one vote per
hog per species pair); `h_i^t` is a gene's hog weight into species t and
`H_st` the species-pair total, so the null is a configuration model blocked
by species pair. Cross-species pairs carry no co-expression null and
same-species pairs no hog null; `lambda` is in scaled-edge units, as E1's
`weight`. Serial Louvain in Rcpp with one refinement pass on the base graph
(Rotta & Noack 2011; plain Louvain could not release a defecting copy from
its hog's supernode). Checks: C++ quality equals the R formula to 2e-16; no
single move improves the returned partition on a toy; at `lambda = 0` the
per-species terms equal igraph's modularity to 3e-15 and match Leiden within
0.7 %; on `jm_fixture()` at `lambda = 4` ARI 1.00 per species, hog purity
0.99, the planted defector outside `module_main`. Wood (121,385 genes,
75.6 M co-expression entries): 30 s per run.

**Control.** Two-step + match: `detect_modules()` per species, then the same
kernel with genes frozen in their modules, so only whole modules move. It
isolates the one thing E2 adds, letting hogs move single genes.

**Planted.** Six species, 4,000 genes each, 12 modules of 2-5 % of a layer,
hog copy profile from the wood table covering half the genes, each copy in
its hog's module with probability `1 - delta`. Mean degree 60 (per-species
Leiden AMI 0.90) and 10 (AMI 0.49). Each engine at its best lambda by AMI,
three seeds (`planted_best.tsv`):

| d | delta | engine | lambda | AMI | modules found (Jaccard > 0.5) | defector recall |
|---|---|---|---|---|---|---|
| 10 | 0.1 | match | 0.25 | 0.40 | 0.06 | 0.79 |
| 10 | 0.1 | E2 | 0.25 | 0.38 | 0.11 | 0.77 |
| 10 | 0.3 | match | 4 | 0.41 | 0.06 | 0.77 |
| 10 | 0.3 | E2 | 0.25 | 0.37 | 0.06 | 0.76 |
| 60 | 0.1 | match | 0.25 | 0.73 | 0.89 | 0.95 |
| 60 | 0.1 | E2 | 0.25 | 0.72 | 0.86 | 0.93 |
| 60 | 0.3 | match | 1 | 0.73 | 0.89 | 0.91 |
| 60 | 0.3 | E2 | 0.25 | 0.74 | 0.75 | 0.88 |

AMI is capped near 0.73 for every joint engine because the truth keeps
each species' background apart and the hogs merge it. E2 never beats the
match on defector recall at equal AMI. The one place E2 does something the
match cannot is weak layers with near-complete conservation: at d = 10,
delta 0 to 0.1, E2 at lambda 4-16 recovers 0.72-0.97 of the planted modules
against 0.03-0.06 for the match, by pooling six weak signals through the
hogs; the price is defector recall 0.05-0.20 (hogs drag defectors along),
and at delta 0.3 the gain is gone. E1 at weight 4 keeps per-species
modules (shared labels 0-0.08) and from 32 fragments (560-1,730 modules),
as on wood.

**Wood split halves.** Samples split per species (seed 1; 32-53 per half),
networks at density 0.03 on all genes, modules of 20 or more genes, three
seeds, held-out AUROC (`dev/probes/module_auroc.R`, 100 nulls, raw z) of
half-A modules on half-B networks:

| engine | lambda | joint modules | per-species split-half ARI | joint ARI | hog purity | hog pairs in one module | ARI vs two-step | z cross | z within |
|---|---|---|---|---|---|---|---|---|---|
| two-step | 0 | 43 per-species | 0.48 | 0.55 | 0 | 0 | 1.00 | 12.8 | 73.6 |
| match | 0.25-16 | 8 | 0.48 | 0.37 | 0.10 | 0.25 | 1.00 | 12.8 | 73.7 |
| match | 64 | 8 | 0.48 | 0.36 | 0.11 | 0.26 | 1.00 | 13.2 | 72.8 |
| E2 | 0.25 | 9 | 0.49 | 0.38 | 0.10 | 0.24 | 0.88 | 12.9 | 73.4 |
| E2 | 4 | 9 | 0.47 | 0.34 | 0.10 | 0.25 | 0.83 | 13.1 | 75.7 |
| E2 | 16 | 8 | 0.51 | 0.33 | 0.12 | 0.27 | 0.86 | 13.4 | 73.6 |
| E2 | 64 | 9 | 0.49 | 0.40 | 0.23 | 0.34 | 0.66 | 14.1 | 66.9 |
| E1 | 32x | 46 per-species | 0.53 | 0.56 | 0.00 | 0.00 | 0.69 | 11.8 | 66.9 |

Both species-blocked engines couple the species at every lambda (8-10 joint
modules from 0.25 on), which E1 never did: a cross-species merge costs
nothing in co-expression, so any positive hog term joins the modules. E2
reshapes within-species modules (ARI 0.83-0.88 against two-step up to
lambda 16; two-step against itself over seeds is 0.85-0.95) without gain:
every E2 advantage over the match (per-species ARI +0.03, joint ARI +0.04,
z cross +0.9) is inside the three-seed spread (0.07, 0.04, 1.4). On its own
objective E2's Louvain finds no more than the match up to lambda 16 (Q 3.856
against 3.855-3.897); at 64 it does (4.13 against 4.03) and pays with the
within-species ceiling (z 66.9) and ARI 0.66 against two-step. Hog purity is
0.10 for both engines: with 8-10 joint modules, one multi-species hog in ten
has all copies in one module, and a quarter of ortholog pairs share a module
(chance about 0.11). The wood data sit below the planted delta = 0.3 (pair
agreement about 0.5) where even E2's pooling gain was gone.

**Verdict.** E2 fails both halves of 7c's rule against the honest control:
no defector-recall gain at equal AMI on planted data, no wood gain beyond
seed noise. The species-blocked objective's optimum is the two-step
partition with modules matched by hog overlap. Recommendation: withdraw
WP15b (h-Leiden port, GPL-3 switch), since nothing a kernel adds on this
data is worth a licence change, and make the product a matching step:
per-species modules, then joint labels from hog overlap against the
species-pair null, which the frozen kernel computes in one second. E3
(species-blocked Bethe Hessian) would optimise the same blocked structure
by another route; its one open question, minority modules in weak layers
with near-complete conservation, is the d = 10 regime, which wood is not
in at density 0.03. Recommendation: do not build it. Both are Martin's
calls. Consequence for `compact-null.md` section 0: joint modules by
matching are the per-species modules, so the compact-null rerun on them
is the run already done (STOP stands). Open knob: `gamma` was 1 throughout; finer
modules may hold more conserved families, but every engine would need the
same sweep and the control would move with it.

## Reading list

| rank | citation | DOI / arXiv | why |
|---|---|---|---|
| 1 | Kaminski, Misiorek, Pralat & Theberge 2024, J Complex Netw | 10.1093/comnet/cnae041; arXiv:2406.17556 | the WP15b objective, tau-family, alpha schedule, lift-off |
| 2 | Li, Schaub & Peel 2026, Sci Adv | 10.1126/sciadv.aef2184; arXiv:2601.10502 | SNR budget for mixed sizes; BH with K from the spectrum; beats h-Louvain; split bias |
| 3 | Chodrow, Veldt & Benson 2021, Sci Adv | 10.1126/sciadv.abh1303 | modularity = DCHSBM likelihood; AON; when it beats clique expansion |
| 4 | Pamfil, Howison, Lambiotte & Porter 2019, SIAM J Math Data Sci | 10.1137/18M1231304 | omega (lambda) as an estimable persistence parameter; multilayer SBM |
| 5 | Li & Peel 2026, arXiv | arXiv:2604.18565 | minority communities: three phases; BH weaker than BP; K rule moves boundaries |
| 6 | Peixoto, Peel & Gross 2026, arXiv | arXiv:2602.16937 | graph models already cover multibody interactions; the burden of proof sits with the hypergraph |
| 7 | Poda & Matias 2024, Peer Community J | 10.24072/pcjournal.404 | side-by-side benchmark of hypergraph modularities and codes |
| 8 | Veldt, Benson & Kleinberg 2022, SIAM Rev | 10.1137/20M1321048 | splitting functions: star, clique, AON, tau in one frame |
| 9 | Ni, Deng & Mu 2025, arXiv | arXiv:2505.04967 | SBM over several hypergraphs joined by inter-hypergraph edges (gene-protein) |
| 10 | Kaminski, Poulin, Pralat, Szufel & Theberge 2019, PLOS ONE | 10.1371/journal.pone.0224307 | Chung-Lu hypergraph null, strict modularity |
| 11 | Pellegrin, Fesser & Weber 2025, NeurIPS | arXiv:2502.09570 | node encodings from hypergraphs; graph-level models on expansions win; section 9 |
| 12 | Coupette, Dalleiger & Rieck 2023, ICLR | arXiv:2210.12048 | ORCHID: Ollivier-Ricci curvature for hyperedges; basis for the hog-curvature probe |

Dropped from the top 10, still cited in the text: Bazzi 2016, Chodrow 2020
and Mucha 2010.

## Method

Lead (Opus) read both seed PDFs (Vidiella 2026; St-Onge et al. 2022, Commun
Phys 5:25, contagion, nothing on communities), Kaminski 2024 and Li 2026 in
full text, the abstracts of the other arXiv seeds (2609.12175, temporal
hypergraph motifs; 2502.09570; 2304.10031, a deep-learning survey; none on
communities) and the plan and borrow-scope notes; all DOIs checked on
Crossref. Four Sonnet agents: citation walk, software survey (gh, GitHub,
CRAN), Nature collection ("Higher-order interaction networks 2021", Commun
Phys, 20 articles, 4 relevant, none on hypergraph modularity), x.com skim.
Asta worked as an MCP JSON-RPC endpoint (`asta-tools.allen.ai/mcp/v1`, key
as a header only; REST paths 404); it walked papers citing Kaminski 2019 and
2024, Chodrow 2021, Contisciani 2022 and Li 2026 (390) plus keyword searches.
Semantic Scholar returned 429 and was not needed. Not done: a forward walk
from Mucha 2010; section 4 rests on known papers checked by DOI. Title-only
hits are not used. x.com and Bluesky `site:` searches found no relevant
posts, so nothing was fetched. Second batch (section 8): five seed PDFs, the
2024 Commun Phys collection (Sonnet skim), and one more Asta walk; abstracts
read through arXiv and Crossref, and every new DOI checked on Crossref.
