# A compactness-matched null for module preservation (WP-k)

2026-10-14. Branch `docs/compact-null`, base `8a1afb3` (main after the
driver-tiers note). Source: Lee H, Lee H, Baek Y, Choi A, Shin H (2026)
MoBC-Net, Bioinformatics Advances, doi 10.1093/bioadv/vbag279 (accepted
manuscript, 11 pages), and its code, github.com/mogam-ai/MoBC-Net at
`df11d88` (2026-08-12), MIT licence, copyright 2026 MOGAM. Martin: plan
a second null for `module_preservation()` that asks "is this module
more preserved than an equally sized, equally hub-rich, equally
cohesive random gene set from the same species?". The current null
stays.

## 0. Verdict

- Borrow the idea, not the distance. ABin (adaptive degree bins) and
  RandSDM (degree-matched sequential sampling that favours compact
  sets) port directly. Their Gaussian kernel on shortest-path distance
  does not: at density 0.03 a module's neighbourhood covers most genes
  after a few dozen draws, so every candidate sits one step away and
  the weights go flat. Replace the distance by the share of the chosen
  set a candidate is wired to, `w = exp(beta * f)`, and calibrate
  `beta` per module so the random sets are as cohesive as the module.
- Draw the random set in the reference species, map it with the same
  ortholog map, score it in the test species. Matching cohesion in the
  test species would match `avg.weight` by construction and remove the
  signal the test measures.
- It is a second null beside the gene-identity null, not a
  replacement. The gene null asks "more preserved than random genes";
  this one asks "more preserved than any cohesive set of the same
  degrees". New columns `Zsummary_compact`, `p_compact`, `q_compact`.
- Gate first. A probe on the wood data (section 4) must show that the
  compact null removes the correlation of the score with module size
  and hubness, and that it changes the module ranking. If either fails,
  build nothing.
- Three work packages after the probe: kernel (WP-k1), wiring into
  `module_preservation(compact = TRUE)` (WP-k2), and a gated call
  (WP-k3). The surface guard allows a 12th formal on
  `module_preservation()`; that is gate G-k2.

### Probe outcome (2026-10-14): STOP, deferred to after WP15

The probe (`prepare_data/evotree/e8_compact_null.R`; tables in
`prepare_data/cache/compact_probe/`) ran on the six wood directions
at the default module resolution and at resolution 8, 200 draws.
Medians over directions, resolution 8 (39-55 modules per direction;
the default resolution has 6-8 modules of about 19 % of the mappable
genes each, too few and too large to test):

| Condition (section 4)                  | Result                 | Gate       |
|----------------------------------------|------------------------|------------|
| Bias present: rho(size), rho(hubness)  | 0.358, 0.667           | pass       |
| Removed by arm 3                       | 0.186, 0.214           | fail by 0.007 on size |
| Ranking changed: rho, top-quartile     | 0.789, 0.757           | pass       |
| Feasible: share cohesion_ratio >= 0.9  | 0.377                  | fail (0.75) |

- The bias is real and the kernel nearly removes it, but matching is
  infeasible: no degree-matched gene set outside a module is as
  cohesive as the module, at any `beta` (checked to 1024; the 64 cap
  is not the cause). Wood modules are about twelve times denser than
  degree alone predicts. Only dropping degree matching reaches the
  module's cohesion, which defeats the purpose.
- `beta = 0` alone does not remove the bias. The expression stratum
  is droppable (arm 4 within 0.09 of arm 3 on every rho).
- WP-k1, WP-k2, WP-k3 are not built. Martin (2026-10-14): rerun the
  probe after WP15a and WP15b land. Joint modules across species are a
  different unit, built on the mapped genes and finer than one
  species' Leiden modules, so section 4 should be rerun on them before
  the verdict is revisited. The scoring kernel and spec in section 3
  stand as written.

## 1. What the paper did

### 1.1 Setting

MoBC-Net tests whether two disease modules are closer in a STRING PPI
network than expected, then ranks link genes on the shortest paths
between them. Networks: 15,867 proteins and 235,783 interactions
(TCGA-BRCA), 12,543 and 169,320 (mouse cachexia, CAC). Mean degree
about 30. Modules: Walktrap on the DEG-induced subnetwork, 8 modules
over 20 genes (BRCA), 9 over 10 genes (CAC). Distance: unweighted
shortest paths in the giant component (section 2.1.3).

Closest distance between modules U and W (section 2.2):

```
D(U, W) = 1 / (|U| + |W|) * ( sum_{u in U} min_{w in W} d_uw
                            + sum_{w in W} min_{u in U} d_wu )
```

### 1.2 ABin (section 2.3.1; code `abinning_estimate_deg_bag_prob()`,
`R/Degree_random_v2.r`)

Input: degree `k(v)` of every gene, the module sets `M_m`, the draw
count `n_cv`, `kappa` (default 0.07).

1. One bin `B_d` per distinct degree, ascending.
2. Demand per bin: `p_{d,m} = |B_d n M_m| / |B_d|`;
   `n_d = n_cv * max_m p_{d,m}`. Threshold `T = kappa * n_cv`.
3. While some bin has `n_d >= T`: take the highest-degree such bin and
   merge it into its lower neighbour (the lowest bin merges upward).
   Recompute `n_d` for the merged bin.

`n_cv` cancels: the rule is `max_m |B_d n M_m| / |B_d| < kappa`. No
module may take more than 7 % of any bin. In the code a bin holding no
module gene starts at `n_d = Inf`, so it is always merged: every final
bin holds at least one module gene. The code shares one set of bins
across all modules (the max over `m`); the introduction says bins are
built per module. The code's default is `ratio = 0.1`; the paper uses
0.07, chosen to match fixed bins of 100 genes (Guney 2016).

A module larger than `kappa` of the pool cannot satisfy the rule in any
bin and collapses to one bin: no degree matching at all.

### 1.3 RandSDM (section 2.3.2; code `modularity_sampling_multi()`)

Per module `m`, per draw:

1. Take the degrees of the module's genes and shuffle them.
2. For each degree in that order, the candidates are the unused genes
   of its ABin bin. The first gene is drawn uniformly. Every later gene
   is drawn with weight

   ```
   w(v) = exp( - min_{u in U_t} d_uv^2 / 2 )
   ```

   where `U_t` is the set drawn so far (sigma fixed at 1 in the code).
   The chosen gene leaves the pool.
3. The drawn set is the random module. In the code one draw builds all
   modules in turn from one shared pool, so the random modules of one
   draw are disjoint. 1,000 draws per module pair.

The distance matrix is `igraph::distances()` over all gene pairs.

### 1.4 Score (section 2.3.3)

`z = (D_obs - mean(D_null)) / sd(D_null)`; empirical
`p = #{D_null < D_obs} / n` (no +1 in the code), then BH. Monte Carlo:
start at 1,000 draws, add batches of 250 until "prespecified stability
criteria" (supplement, not seen) or 5,000. The GitHub code has no batch
loop: it runs a fixed `random` count.

### 1.5 Results that matter here

- Table 1, MoBC distance against hubness (mean degree of the top 5 %
  genes of the pair) and pair size. CAC: raw -0.46 / -0.67; after
  ABin + RandSDM -0.09 / -0.31. BRCA: raw -0.42 / -0.14; after 0.11 /
  0.52 (the size sign flips: large pairs become conservative).
- Table 2: the null mean and sd themselves track hubness and size (null
  mean -0.82 / -0.62 BRCA, null sd -0.47 / -0.94 CAC). Pair-specific
  standardisation is needed because the null moves with size.
- Figure 4: RandSDM gives random modules with higher modularity and a
  smaller diameter than degree-only sampling (RandSD). Cohesion is
  approached, not matched.
- Section 3.3: in CAC the closest pair changes from M1-M2 (the largest,
  most hub-rich modules) to M4-M8 after standardisation.
- Section 3.3 and the Discussion: on a WGCNA co-expression network the
  modules were larger and more consolidated, and no pair was robustly
  significant. Dense, redundant wiring hides local bridges; the authors
  suggest sparsification and did not test it.

## 2. What rcomplex does today

- `module_preservation()` (`R/module_preservation.R:180-185`), 11
  named formals. Reference modules are projected onto test genes by
  `resolve_ortholog_map()` (`:207`) and `.pres_project()` (`:346`).
- The pool is every ortholog-mappable test gene (`:220-225`, NetRep's
  "overlap" null). Reference per-gene statistics are computed once and
  held fixed (`:302-312`).
- The kernel shuffles the pool once per permutation and hands each
  module a contiguous block of its own size
  (`src/module_preservation.cpp:527-541`). One pass over the induced
  edges scores every module (`score_modules()`, `:179`). Under this
  null `cor.degree` falls to about 0, because the reference vector is
  paired with random test genes.
- Matched: mapped size. Not matched: degree, cohesion, expression
  level, copy number. The reference side never varies, so the null
  never asks how a cohesive reference set fares in the test species.
- Score (`.pres_assemble()`, `:427-543`): one-sided
  `(exceed + 1) / (n_perm + 1)`, `pmax` of the two statistics (`:446`),
  the mixture calibration with `w00` (`:459-466`), BH, and
  `Zsummary_std = ((Z.avg.weight + Z.cor.degree) / 2) /
  (sqrt(2 + 2 rho) / 2)` (`:474-512`).
- `classify_preservation()` (`:661-747`): `q_value < 0.1`, then
  `Zsummary_std >= 10` for "conserved". Its own documentation
  (`:618-634`) records the size problem MoBC-Net's Table 2 shows: over
  110 modules in 12 Pooideae contrasts the null sd of `cor.degree`
  tracks mapped size at Spearman -0.996, that of `avg.weight` at
  -0.873. `medianRank` is the size-free workaround.
- `preservation_paired()` forwards `...` to `module_preservation()`
  (`:936`, `:1026-1035`). `preservation_matrix_test()` reads
  `Zsummary_std` (`R/preservation_matrix.R:335`).
- Surface guard (`tests/testthat/test-surface.R`): exports over 8
  formals must be in `wide_ok` (`:2-4`), with no cap on their count
  (`:34-38`). The formal counts are a snapshot (`:24-32`).
- RNG: `.seed_scope(seed)` at `:200`; the kernel draws one seed per
  permutation from `R::runif` before the parallel region.
- WP-d (NEWS, development version) shuffles hog labels within
  expression deciles for the driver null. The removed `module_auroc()`
  (0.3.2) matched random HOG sets on size, copy class and degree
  decile; its gate failed until the copy class was added
  (`module-engine-borrow-scope.md` 10.1).

## 3. Spec

### 3.1 Where the random set is drawn

Source: the reference species. The module was detected there, so it is
cohesive there by construction. The question is whether that cohesion
travels further than the cohesion of other reference gene sets of the
same size and degrees. So: draw a cohesive random set in the reference
network, map it, score it in the test network, exactly as the module
is scored.

Not the test species. A test-side set matched to the module's test
cohesion has the module's `avg.weight` by construction, so
`Z.avg.weight` goes to 0 whatever the biology.

Per module and per draw:

- Pool: mappable reference genes (at least one test partner in the
  resolved map), minus the module's own genes. This is the package's
  self-excluded urn: a compact draw would otherwise regrow the module
  and test it against itself.
- Draws are independent per module. Modules of one draw are not made
  disjoint. With a full partition, disjoint draws leave the last
  modules whatever remains, and no compactness is possible there. The
  test side scores one module at a time, which costs the same as one
  pass for all (each gene's edges are read once per module it is in).

### 3.2 Strata and bins

Each reference gene gets a stratum: copy class (1 or 2+ test partners
in the map) crossed with `gene_mean` quintile over the pool (when the
network has `gene_mean`; else one level). Within each stratum, ABin
(1.2) runs per module on degree, with `kappa = 0.07`, computed over the
pool. Degree is counted in the reference graph induced on the mappable
genes, the same graph the sampler walks.

- Copy class keeps the mapped size and the multi-copy share of the
  random set equal to the module's. `module_auroc()` failed its
  calibration gate without it.
- Quintiles, not deciles: ten strata already halve the room ABin needs.
- Report per module `n_bins` (sum over strata). One bin per stratum
  means degree is not matched there; the probe reports how often.

`kappa` is an internal constant, not an argument.

### 3.3 Sampler

Choice (b) below. Per module of `s` reference genes, per draw:

1. List the module's (stratum, bin) labels, one per gene; shuffle.
2. The first gene is uniform among the pool genes of its (stratum,
   bin). Each later gene `v` among the unused genes of its (stratum,
   bin) gets weight

   ```
   f_U(v) = |N(v) n U| / |U|        (share of the chosen set wired to v)
   w(v)   = exp(beta * f_U(v))
   ```

   with `N(v)` the neighbours of `v` in the thresholded reference
   graph. Edges count as 0 or 1: at density 0.03 the stored weights are
   near-constant (kernel comment, `src/module_preservation.cpp:390-
   393`), and integer counts keep the sampler exact and cheap.
3. Keep an integer array `c(v) = |N(v) n U|`. Adding gene `u` adds 1
   to `c` of each neighbour of `u`: O(degree). At the end `c(u)` for
   `u` in `U` is the reference intramodular degree of the draw, which
   `cor.degree` needs. No second pass.

The knob is `beta >= 0`. `beta = 0` is RandSD (degree and stratum
matched, uniform within bins: Guney 2016). Large `beta` grows the set
greedily around the first gene.

Cohesion of a gene set `S` in the reference network is its
`avg.weight` there, `sum_{u != v in S} a_uv / (|S|^2 - |S|)`, the
statistic the module is tested on. Calibration per module: bisect
`beta` on [0, 64], eight steps of 25 pilot draws, until the median
cohesion of the draws is within 5 % of the module's. Record `beta` and
`cohesion_ratio` (median draw cohesion / module cohesion over the real
draws). `cohesion_ratio < 0.9` at `beta = 64` means no equally cohesive
set exists outside the module; the row is flagged, not dropped.

Alternatives considered:

- (a) Published: Gaussian kernel on hop distance. Cost: all-pairs
  BFS, O(n E), about 20,000 x 6 M edges, and a 1.6 GB distance matrix
  per species. Fails: with mean degree 450-700 a random set of 50
  genes is adjacent to about 80 % of all genes (1 - 0.97^50), so nearly
  every candidate sits at distance 1 and the weights go flat. The
  paper's own WGCNA run found nothing.
- (b) Share of the chosen set wired to the candidate, `exp(beta f)`.
  Cost: O(s (b + k)) per draw. Chosen: graded where (a) is binary, one
  knob, calibrated to the module.
- (c) Random walk or snowball from a seed, within the bins. Cheap.
  Fails: a walk visits genes in proportion to degree, which breaks the
  degree match, and on a dense graph it mixes in one or two steps, so
  it is barely more cohesive than (b) at `beta = 0`.
- (d) Uniform degree-matched draws, rejected unless cohesion matches.
  Cost: O(s) per try. Fails: acceptance near 0, because modules sit
  far above the 0.03 background density and a uniform set does not.

`s` is module size, `b` the mean pool size of a bin, `k` the mean
degree.

### 3.4 Mapping and statistics

- Precompute the resolved map as a reference-to-test adjacency (sorted
  integer vectors). A random set's test block is the union of its
  genes' test partners. A test gene reached from two drawn reference
  genes takes the smaller reference index, as `.pres_project()` takes
  the smallest `gene1`.
- Per test gene of the block, the paired reference value is `c(u)` of
  its reference gene, so a reference gene mapped to two paralogs counts
  twice, as in the observed run (`:308-310`).
- In the test network: `avg.weight` and `cor.degree` of the block, by
  the same code as the observed module. `meanClusterCoeff` and
  `meanMAR` are not computed for the compact null; they take no part in
  the call.
- A block under 10 test genes, or with constant `kIM`, gives `NA` for
  that draw and leaves the joint count, as now (`n_joint`).

### 3.5 Scores and columns

Per module and statistic: `Z = (obs - mean(null)) / sd(null)`; one-sided
`(exceed + 1) / (n_perm + 1)`; `pmax`; the same mixture calibration and
BH, reusing `npc_combine()` and `.pres_w00()`. New columns in
`preservation`, `NA` when `compact = FALSE`:

- `Z.avg.weight_compact`, `Z.cor.degree_compact`
- `Zsummary_compact`: `((Z1 + Z2) / 2) / (sqrt(2 + 2 rho) / 2)`, rho the
  null correlation of the two statistics over the compact draws, clamped
  to [0, 1] as for `Zsummary_std`.
- `p_compact` (calibrated), `q_compact` (BH).

New table `compact`, one row per tested module: `module`, `beta`,
`cohesion_ref` (module), `cohesion_null` (median draw), `cohesion_ratio`,
`n_bins`, `size_mapped_null` (mean test block size), `n_joint`.

Read the two nulls side by side. `Zsummary_std` and `q_value`: the
module's wiring survives better than random genes'. `Zsummary_compact`
and `q_compact`: it survives better than other cohesive sets of its
degrees. A module high on the first and near 0 on the second is
preserved as any dense region of the reference network would be.
`classify_preservation()` does not change in WP-k2 (gate G-k3).

### 3.6 Expression level

Match it: yes, as the quintile stratum in 3.2. It costs one integer
per gene and no extra pass. WP-d measured that abundance doubles the
false call rate of the pair test (high-high cells 18-24 % against
7-14 %), and highly expressed genes have more reliable correlations,
so a highly expressed module is likely to be preserved for that reason
alone. Degree tracks expression, but not fully. The probe keeps an arm
without the stratum; if the score's correlation with module mean
expression moves by less than 0.1 between the arms, drop the stratum
and say so.

### 3.7 Monte Carlo

Fixed `n_perm`, the same argument as the gene null, no batch stopping.
The package reads attainable p-values from `n_perm`
(`pvalue_resolution()`), and ranks on Z, which needs only the null
mean and sd. The paper's stopping criteria are in its supplement and
absent from its code. If the compact null at the default 10,000 draws
is too slow (3.8), the fix is a faster sampler, not a second count.

### 3.8 Cost

Per direction: `M * n_perm * (s (b + k) + s k_test)` simple operations,
plus calibration (`8 * 25 = 200` pilot draws per module, 2 % at 10,000).
Wood: 15,000-23,000 genes, mean degree 450-700. With `M = 20` modules
of `s = 500`, `b = 1,500`, `k = 600`: about 3e7 per draw, 3e11 per
direction at 10,000 draws, under a minute at 8 threads. The probe
measures the real figure. Upgrade path, if needed: group a bin's
candidates by their count `c(v)`, so a draw costs O(s) per step, not
O(b). Mark it `ponytail:` in the kernel.

C++ (`src/compact_null.cpp`, new): reference neighbour lists from
`induced_weighted_lists_dense()` / `_sparse()` (`src/neighbor_lists.h`)
on the mappable reference genes, binarised; the test side reuses
`WeightedNeighbors` and a one-module variant of `score_modules()`
factored out of `src/module_preservation.cpp`. Integer indices, sorted
vectors, no string maps. One `std::mt19937` per draw, seeded from
`R::runif` before the parallel region, so the result does not depend
on `n_cores`. OpenMP pragmas carry `if(n_cores > 1)`. No R-level
forking: `preservation_paired()` runs its directions serially.

### 3.9 Surface

`module_preservation()` gains `compact = FALSE`, its 12th named formal.
The guard passes: the function is in `wide_ok` by name and the budget
test has no count cap. The surface snapshot changes by one number.
`preservation_paired()` reaches it through `...`, no change there.

Cheaper routes, rejected:

- Run both nulls always (no formal): doubles the cost for every call
  and changes the default output.
- Overload `calibrate` or `n_perm`: two meanings in one argument.
- A separate export: 32 exports against the budget of 31.

### 3.10 Relation to the other nulls and gated WPs

- Shuffled-expression (RGG) null, `module-engine-redesign.md` 5.1: it
  asks whether a species has module structure at all. The compact null
  does not touch that question. It conditions on the reference
  network's real cohesion and scores in the test species, where noise
  geometry is independent. A cohesive set made of reference noise maps
  to near-background density in the test species, so noise blobs pull
  the compact null down: the compact null errs liberal, never
  circular. On wood (65-107 samples) noise geometry is weaker than on
  the n = 20 leaf data.
- Complement, not prerequisite, for WP15 (joint modules), WP16 (edge
  history), WP17 (MUNK). None of them reads `module_preservation()`. A
  jointly detected module cannot be tested by `module_preservation()`
  at all (the test species shaped it), so the compact null applies to
  per-species modules from `detect_modules()` or `as_modules()` only.
- Hog level: WP-d is the same idea (match the abundance of what you
  shuffle) one level down. The clique level keeps its matched edge-set
  null.

## 4. Probe gate (WP-k0)

Script `prepare_data/evotree/e8_compact_null.R` (gitignored, beside
`e5_rcomplex_0.4.R`, whose loading it copies). No package change. Slow
R is fine.

- Data: wood, density 0.03, `compute_network()` per species with
  `gene_mean`. Modules: `detect_modules(net, seed = 1)` per species,
  default objective. Report module sizes; if the median module holds
  more than 7 % of the mappable genes, also run at a resolution that
  gives a median under 5 % and report both.
- Directions (6): Asp-Birch, Birch-Asp, Nor-Scots, Scots-Nor, Asp-Nor,
  Nor-Asp. Ortholog table as in e5, resolved by `resolve_ortholog_map()`
  without `edges`.
- Arms per direction and module:
  1. `genes`: `module_preservation(n_perm = 1000, seed = 1)`.
  2. `compact0`: sampler 3.3 at `beta = 0` (strata and ABin only).
  3. `compact`: calibrated `beta`.
  4. `compact_noexpr`: as 3 without the quintile stratum.
  Arms 2-4: 200 draws, a direct R port of 3.2-3.4 on the dense
  matrices (`mclapply` over modules; nothing inside the workers enters
  OpenMP).
- Per module: `size_mapped`; hubness = mean reference degree of the top
  5 % of its genes (at least one gene); mean `gene_mean`; raw
  `avg.weight` and `cor.degree`; per arm `Zsummary`, null mean and null
  sd of both statistics; for 3 and 4 `beta` and `cohesion_ratio`.
- Tables: (T1) Spearman of each score with size, hubness and mean
  expression, per direction and the median over directions, the
  MoBC-Net Table 1 layout; (T2) the same for the null mean and sd, the
  Table 2 layout; (T3) Spearman between `Zsummary_std` and each compact
  arm per direction, and the overlap of their top quartiles; (T4)
  `cohesion_ratio` and `n_bins` distributions; (T5) wall time per
  module per draw.

Stop/go, decided on the median over directions:

- Bias present: `|rho(Zsummary_std, size)|` or `|rho(Zsummary_std,
  hubness)|` at least 0.3. If neither, stop: nothing to fix.
- Removed: for arm 3, both `|rho|` at most 0.3 and at most half the
  `genes` value.
- Ranking changed: Spearman(`Zsummary_std`, arm 3) at most 0.9, or the
  top quartiles share at most 75 % of their modules.
- Feasible: `cohesion_ratio >= 0.9` for at least 75 % of modules.

Go only if all four hold. If arm 2 alone already meets "removed" and
"ranking changed", build `beta = 0` only: no calibration, no
`cohesion_*` columns. If arm 4 matches arm 3 within 0.1 on every
`rho`, drop the expression stratum.

## 5. Work packages

Each is one PR off main, after the gate.

### WP-k1: sampler kernel (internal; about 300 C++ lines, 80 R, 150 test)

Files: `src/compact_null.cpp`, `R/compact_null.R` (internal
`.compact_bins()`, `.compact_draws()`), `R/RcppExports.R`,
`tests/testthat/test-compact-null.R`. `src/module_preservation.cpp`
only to factor out the one-module scorer.

Tests:

1. ABin on three fixture degree vectors: every final bin has module
   share under `kappa` or is the only bin of its stratum; bins are
   contiguous degree ranges; each holds a module gene; equal to a
   20-line R reference port of `abinning_estimate_deg_bag_prob()`.
2. Every draw has the module's exact count per (stratum, bin), no
   module gene, no duplicate.
3. `beta = 0`: selection frequencies uniform within bins (chi-square,
   5,000 draws, p > 0.001).
4. Planted block fixture: median cohesion rises over `beta` = 0, 4,
   16; at 64 the draw recovers the block (Jaccard > 0.8).
5. `c(u)` at the end equals the intramodular degree computed from the
   matrix.
6. Same seed, same draws; `n_cores` 1 and 2 identical; dense and sparse
   networks identical.

### WP-k2: `module_preservation(compact = TRUE)` (about 150 R lines)

Files: `R/module_preservation.R` (formal, roxygen `@section Null
model`, `.pres_assemble()` columns), `R/compact_null.R` (calibration),
`tests/testthat/test-module-preservation.R`,
`tests/testthat/helper-rng-contract.R` (variant `[compact]`),
`tests/testthat/_snaps/surface.md`, `NEWS.md`, `CLAUDE.md` (one line
under "Modules: connectivity, not membership"), `man/`.

Tests:

1. `compact = FALSE`: every existing column identical to main on the
   preservation fixtures; new columns `NA`; `compact` table absent.
2. Fixture where every cohesive reference block is equally preserved:
   `Zsummary_std` above 10, `|Zsummary_compact| < 2`. Fixture with one
   block preserved and equally cohesive blocks not:
   `Zsummary_compact > 3` for that block.
3. Calibration: a compact draw used as the module, 50 seeds:
   `Zsummary_compact` mean within 0.3 of 0, sd in 0.7-1.3.
4. RNG contract passes with the new variant; `seed` reproduces.
5. `preservation_paired(..., compact = TRUE)` carries the columns.
6. Surface: 31 exports; `module_preservation` 12 formals.

Acceptance: `R CMD check` OK on the built tarball, lint clean,
`Rscript dev/sharpen-prose.R` clean. Wood, the six probe directions at
`n_perm = 10000`: per direction modules tested, median `Zsummary_std`
and `Zsummary_compact`, T1 within 0.1 of the probe, fraction with
`cohesion_ratio >= 0.9`, wall time against `compact = FALSE`.

### WP-k3: a call on the compact null (gated G-k3; about 40 lines)

Only if WP-k2 on wood shows modules high on `Zsummary_std` and near 0
on `Zsummary_compact`. `classify_preservation()` adds `compact_call`:
"beyond topology" when `q_compact < 0.1`, else "topology only";
`NA` without the columns. The existing `classification` does not
change.

## 6. NEWS (development version, WP-k2)

- `module_preservation()` gains `compact`. If `TRUE`, each module is
  also tested against random reference gene sets of the same size,
  degrees, copy number, expression level and cohesion, mapped and
  scored like the module. New columns `Zsummary_compact`, `p_compact`,
  `q_compact`, and a `compact` table with the sampler settings. The
  sampler borrows ABin and RandSDM from MoBC-Net (Lee et al. 2026) and
  replaces their path distance by the share of the drawn set a gene is
  wired to.

## 7. Out of scope

- Inter-module distance, crosstalk, MoBC link genes, S2B, FCS.
- Shortest paths of any kind; network sparsification.
- Within-species module significance (`module-engine-redesign.md` 5.1).
- A test-side compact null (3.1).
- Changes to `detect_modules()`, the hog-level null (WP-d) or the
  clique nulls.
- Batch Monte Carlo stopping.

## 8. Left open (gates for Martin)

- G-k0: run the probe (Orion, one array task per direction)?
- G-k1: after the probe: calibrated `beta`, `beta = 0` only, or
  nothing.
- G-k2: a 12th formal on `module_preservation()` (`compact = FALSE`),
  or the compact null on by default once it has run on real data.
- G-k3: does `classify_preservation()` get `compact_call`, and should
  `preservation_matrix_test()` ever read `Zsummary_compact`?
- Large modules. With the default modularity objective, wood modules
  may exceed 7 % of the pool; ABin then matches size and strata but not
  degree. Accept that, or test the compact null only at a finer
  resolution?
- One `n_perm` for both nulls (3.7). If the compact null at 10,000
  draws costs more than twice the gene null on wood, cap it at 1,000
  with its own floor, or speed up the sampler first? Default: speed up.
