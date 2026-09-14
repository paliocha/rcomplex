# MDL engine — research, design notes, and handoff

Status as of 2026-09-14: **not started.** This document exists so a
different agent (human or AI) can pick the work up without re-doing the
literature review. Everything below is research and planning; no R or
C++ code for an MDL engine exists anywhere in this repository yet. Grep
confirms it:

```sh
grep -rli "mdl\|logchoose\|logmultiset\|backbon" R/ src/
```

turns up nothing but an unrelated compiled `.o` file.

## 1. What problem this is meant to solve

Across several conversations the user raised the same underlying worry
about `rcomplex`'s density/threshold machinery (`compute_network()`'s
density cutoff, `density_sweep()`'s multiplier grid, the new
`coexpressolog_strength()` density grid, `suggest_reference_density()`'s
WGCNA-style scale-free stopping rule): **there are several free
parameters (which densities to test, how many, how to aggregate across
them) and no principled, parameter-free way to pick a single
"best" network sparsification level.** The user explicitly asked
whether an MDL (minimum description length) approach could replace
some of this guesswork, calling out Kirkley's backbone-extraction work
by name.

MDL backbone extraction answers a **different but complementary**
question to what `coexpressolog_strength()` already does:

- `coexpressolog_strength()` (implemented, committed) asks: *given a
  small grid of matched densities, how consistently does this specific
  edge's evidence hold up across the grid?* It is a robustness/
  sensitivity measure, deliberately keeping the whole grid rather than
  collapsing it to one answer.
- An **MDL engine** would ask: *of all the ways to keep only some
  edges from a weighted network, which single subset best compresses
  the information in the full weighted graph, with no threshold or
  density chosen by the analyst at all?* It is a **backbone
  selection** method — it would produce a candidate for a principled,
  data-driven `reference_density` (or a full replacement for choosing
  one), not a robustness curve.

These are not competing implementations. If both existed, the natural
integration point would be: use the MDL engine to propose (or sanity
check) `reference_density`, and keep `coexpressolog_strength()`'s
density grid as the separate robustness sweep around it.

## 2. Literature verified as legitimate

All of the following were checked directly against primary sources
(Crossref API, PubMed/PMC, journal/publisher pages, and the GitHub API
for license metadata) — not assumed from memory. DOIs resolve at
`https://doi.org/<DOI>`.

### MDL backbone extraction (the core method to port)

- **Kirkley, A.** "Fast Nonparametric Inference of Network Backbones
  for Weighted Graph Sparsification." *Physical Review X* 15, 031013
  (2025). **DOI: `10.1103/4pg6-mtmt`**
  (<https://doi.org/10.1103/4pg6-mtmt>). Introduces the MDL backbone
  objective (`fglobal`/`flocal`) ported below: a parameter-free,
  log-linear-time greedy method that trades off backbone size against
  how well it reconstructs the full weighted graph. Kirkley is faculty
  at HKU (Complex Networks Lab), PhD under Mark Newman — mainstream
  network-science lineage, not a fringe source.
- **Kirkley, A. & He, B.** "PANINIpy: Package of Algorithms for
  Nonparametric Inference with Networks In Python." *Journal of Open
  Source Software* 9(103), 7312 (2024). **DOI:
  `10.21105/joss.07312`** (<https://doi.org/10.21105/joss.07312>).
  The JOSS-reviewed, packaged Python implementation, including the MDL
  backboning module used as the porting source (Section 3).
  - GitHub: <https://github.com/HKU-Complex-Networks-Lab/PANINIpy>
  - PyPI: <https://pypi.org/project/paninipy/>
  - **License: MIT**, confirmed via the GitHub API
    (`GET /repos/HKU-Complex-Networks-Lab/PANINIpy` →
    `"license": {"key": "mit", ...}`). Reuse with attribution is
    legally clear.
- **Kirkley, A.** `mdl-network-backbones`
  (<https://github.com/aleckirkley/mdl-network-backbones>) — Kirkley's
  original, un-packaged repository, the one cited in the PRX paper's
  code-availability statement. **Confirmed via the GitHub API to carry
  no license** (`"license": null`, and `LICENSE`/`LICENSE.md` return
  HTTP 404 on both the `main` and `master` branches). **Do not port
  code from this repository.** Use PANINIpy instead (same author, same
  algorithm, JOSS-reviewed, MIT licensed).

### Adjacent MDL / statistically-validated-network lineage (future extensions, supplied by the user, verified but not yet acted on)

- **Betti, L., Musciotto, F., Battiston, F. & Mantegna, R. N.**
  "Identifying maximal sets of significantly interacting nodes in
  higher-order networks." *Physical Review E* (2026). **DOI:
  `10.1103/t86l-31d5`** (<https://doi.org/10.1103/t86l-31d5>); arXiv
  preprint: `arXiv:2209.12712`. Extends statistically-validated-network
  filtering (the Tumminello et al. 2011 lineage below) to higher-order
  (hypergraph) interactions. Relevant if the package ever validates
  cliques the way `coexpressolog_strength()` validates edges — not
  needed for a first static-graph MDL port.
- **Weng & Kirkley.** Spatiotemporal MDL network simplification,
  `arXiv:2605.05008` (<https://arxiv.org/abs/2605.05008>). Same
  author's extension of the MDL simplification idea to time-varying
  networks. Not needed for a first static-network port; noted for
  completeness since the user supplied it directly.
- **Baltakienė, M., Baltakys, K., Cardamone, D., Parisi, F.,
  Radicioni, T., Torricelli, M., van Lidth de Jeude, J. A. & Saracco,
  F.** "Maximum entropy approach to link prediction in bipartite
  networks." `arXiv:1805.04307` (2018)
  (<https://arxiv.org/abs/1805.04307>). **Preprint only — confirmed
  via OpenAlex (`is_published: false`) never appeared in a
  peer-reviewed venue; no separate DOI beyond the arXiv one
  (`10.48550/arXiv.1805.04307`).** A Complexity72h workshop product
  (IMT Lucca). Uses the bipartite configuration model (BiCM) maximum-
  entropy null to score/predict links in a bipartite graph — the same
  entropy-null-model lineage as the Tumminello statistically-
  validated-network approach below, applied to the bipartite case.
  Potentially relevant to this package's species x HOG structure
  (which is itself bipartite before per-species network projection),
  but orthogonal to the MDL backbone port in Section 3: BiCM scores
  edges against an entropy null, it does not compress/select a
  backbone by description length. Not needed for a first static-graph
  MDL port; noted for completeness since the user supplied it
  directly.

### Supporting literature already cited in `R/coexpressolog-strength.R` (cross-referenced here for a future MDL design-notes header, not the porting source)

- Tumminello, M., Miccichè, S., Lillo, F., Piilo, J. & Mantegna, R. N.
  "Statistically Validated Networks in Bipartite Complex Systems."
  *PLoS ONE* 6(3): e17994 (2011). **DOI: `10.1371/journal.pone.0017994`,
  PMID: 21483858.** This is the exact-hypergeometric-p + FDR method
  already reimplemented (independently, not via this paper's own code)
  in `coexpressolog_strength()`.
- Netotea, S., Sundell, D., Street, N. R. & Hvidsten, T. R.** "ComPlEx:
  conservation and divergence of co-expression networks in A.
  thaliana, Populus and O. sativa." *BMC Genomics* 15:106 (2014).
  **DOI: `10.1186/1471-2164-15-106`.** The base method this package
  extends.
- Curci, P. L. et al. "Identification of growth regulators using
  cross-species network analysis in plants." *Plant Physiology*
  190(4):2350–2365 (2022). **DOI: `10.1093/plphys/kiac374`.** The
  DS1–DS5 density-subnetwork design that motivated
  `coexpressolog_strength()`'s density grid.
- Zhang, B. & Horvath, S. "A General Framework for Weighted Gene
  Co-Expression Network Analysis." *Statistical Applications in
  Genetics and Molecular Biology* 4(1), Article 17 (2005). **DOI:
  `10.2202/1544-6115.1128`.**
- Langfelder, P. & Horvath, S. "WGCNA: an R package for weighted
  correlation network analysis." *BMC Bioinformatics* 9:559 (2008).
  **DOI: `10.1186/1471-2105-9-559`.**
- Barabási, A.-L. & Albert, R. "Emergence of scaling in random
  networks." *Science* 286:509–512 (1999). **DOI:
  `10.1126/science.286.5439.509`.**
- Albert, R., Jeong, H. & Barabási, A.-L. "Error and attack tolerance
  of complex networks." *Nature* 406:378–382 (2000). **DOI:
  `10.1038/35019019`.**
- Schneider, C. M., Moreira, A. A., Andrade, J. S., Havlin, S. &
  Herrmann, H. J. "Mitigation of malicious attacks on networks."
  *PNAS* 108:3838–3841 (2011). **DOI: `10.1073/pnas.1009440108`.**
- Borate, B. R., Chesler, E. J., Langston, M. A., Saxton, A. M. & Voy,
  B. H. "Comparison of threshold selection methods for microarray gene
  co-expression matrices." *BMC Research Notes* 2:240 (2009). **DOI:
  `10.1186/1756-0500-2-240`.**
- Bleker, C., Grady, S. K. & Langston, M. A. "A Comparative Study of
  Gene Co-Expression Thresholding Algorithms." *Journal of
  Computational Biology* 31(6):539–548 (2024). **DOI:
  `10.1089/cmb.2024.0509`**; PMID: 38781420.
- Su, Z., Liu, Y., Kurths, J. & Meyerhenke, H. "Generic network
  sparsification via degree- and subgraph-based edge sampling."
  *Information Sciences* (2024). **DOI: `10.1016/j.ins.2024.121096`.**
  ⚠️ **Correction flag:** `R/coexpressolog-strength.R`'s current
  header cites this as "Su, Meyer, Kurths, Marwan & Meyer" — the
  Crossref-verified author list is **Su, Liu, Kurths & Meyerhenke**.
  The existing in-repo citation has the wrong author names and should
  be fixed independently of the MDL work.
- Hassan, R. & Arifuzzaman, S. "Learning on Incomplete Graphs:
  Benchmarking GNN Robustness to Sparsification," presented at IEEE
  ICMLA. ⚠️ **Year flag, not fully resolved this session:**
  `R/coexpressolog-strength.R`'s header cites this as IEEE ICMLA 2021.
  The only accessible copy of the paper found this session is hosted
  at a filename suggesting `icmla26-gnn.pdf` (i.e. a 2026 event), and
  no DOI could be confirmed by Crossref lookup (a `10.1109/ICMLA...`
  DOI guess resolved to an unrelated paper). Authors and title are
  correct; **the venue year needs independent verification before
  citing further** — do not propagate "2021" without re-checking IEEE
  Xplore directly.

## 3. The actual algorithm to port (read in full this session)

Source: `paninipy/mdl_backboning/functions.py` in
`HKU-Complex-Networks-Lab/PANINIpy` (MIT licensed — attribution
required, but reuse is legally clear).
Repository: <https://github.com/HKU-Complex-Networks-Lab/PANINIpy>
File path: `paninipy/mdl_backboning/functions.py`

Input: a **directed, weighted edge list** `[(i, j, w_ij), ...]`
(undirected input is handled by duplicating each edge in both
directions before running the algorithm).

Core idea: choose the subset of edges ("the backbone") that minimizes
a two-part description length — the cost of transmitting *which* edges
are in the backbone, plus the cost of transmitting the full weighted
graph *given* that backbone — with no free threshold parameter. Edges
are ranked by weight and a single greedy prefix (highest weight first)
is chosen; the DL objective is evaluated at every prefix length and the
minimizer is the backbone size. No cutoff is chosen by the analyst;
the "cutoff" falls out of the DL minimization itself.

Two variants, both implemented in PANINIpy and both worth porting:

- **Global backbone** (`fglobal`): one DL objective over the whole
  network.
  ```
  fglobal(W, E, Wb, Eb) =
      log(E + 1) + log(W - E + 1)
      + logchoose(E, Eb)
      + logchoose(Wb - 1, Eb - 1)
      + logchoose(W - Wb - 1, E - Eb - 1)
  ```
  where `W` = total edge weight, `E` = total edge count, and `Wb`/`Eb`
  are the backbone's weight/count at a candidate prefix length.
- **Local backbone** (`flocal`): the analogous per-node objective,
  `flocal(si, ki, sbi, kbi)`, applied independently at each node using
  that node's own strength/degree and backbone strength/degree — this
  is what makes the local variant adapt to heterogeneous degree
  distributions instead of applying one global cutoff everywhere.

Both use `logchoose(n, k)` (log binomial coefficient, implemented via
`scipy.special.loggamma` in Python — R has `lchoose()` / `lgamma()`
built in, so this needs no external dependency) and a `logmultiset()`
helper for the multiset-counting term.

Algorithm shape (both variants): sort edges by decreasing weight once
(`O(E log E)`), do a single forward greedy scan accumulating running
`Wb`/`Eb` (or per-node running sums for the local variant), evaluate
the DL objective at every prefix, and keep the argmin. This is linear
in the number of edges after the sort, so it is cheap even at
whole-transcriptome density.

Other parameters PANINIpy exposes and their meaning if ported:

- `directed`: whether the input edge list is already directed or needs
  duplicating.
- `out_edges`: for the local variant, whether "local" means per
  out-edges or per in-edges of a node (irrelevant for an undirected
  co-expression network, where in/out are symmetric — can likely be
  dropped or hard-coded for rcomplex's use case).
- `allow_empty`: the DL objective is symmetric under "keep nothing" vs.
  "keep everything" in some degenerate cases; this flag forces a
  non-trivial (non-empty) backbone when set.
- `CR_type` (`"Relative"`, `"Max"`, `"Min"`): how the "inverse
  compression ratio" diagnostic (how much the backbone compresses
  relative to the full graph) is normalized for reporting — a
  diagnostic output, not something that changes which edges are kept.

Return value in PANINIpy: the backbone edge list(s) (global and/or
local) plus the compression-ratio diagnostic. A port for rcomplex
would presumably want to return an edge/density recommendation rather
than a full alternate network object, to slot into
`suggest_reference_density()`'s role.

## 4. Open questions that need answering before implementation starts

These were raised but never resolved in conversation — whoever picks
this up should get explicit direction from the user on them first:

1. **Integration point.** Four options were discussed, none chosen:
   - (a) An alternative/complementary path inside or alongside
     `suggest_reference_density()` in `R/coexpressolog-strength.R` —
     i.e., the MDL engine becomes a second, parameter-free way to
     propose a `reference_density`, to be compared against (not
     necessarily replace) the existing WGCNA-style scale-free-fit
     diagnostic.
   - (b) A general-purpose, standalone backbone-extraction utility
     (e.g. `mdl_backbone()`) usable on any rcomplex network
     independent of `coexpressolog_strength()` — closer to a direct,
     general port of PANINIpy's public API.
   - (c) A C++ kernel (consistent with the package's general
     convention of C++ kernels + R wrappers for anything
     performance-sensitive — see `CLAUDE.md`'s "Integer indices in
     C++, string mapping in R" convention). The algorithm's core loop
     (sort + single greedy scan + running sums) is simple enough to be
     a natural `src/` kernel if performance on whole-transcriptome
     edge counts matters.
   - (d) A thin R wrapper around Python's PANINIpy via `reticulate`,
     avoiding a from-scratch port entirely at the cost of adding a
     Python runtime dependency (a significant departure from the
     package's current all-R/C++ dependency profile — `DESCRIPTION`
     currently has no Python/`reticulate` dependency at all).
2. **Global vs. local variant.** Does rcomplex want both DL objectives
   ported, or just one? The local variant is likely the more relevant
   one for gene co-expression networks (highly heterogeneous degree
   distributions are the norm), but this needs confirming rather than
   assuming.
3. **Directed vs. undirected input.** rcomplex's co-expression networks
   are undirected (symmetric MR/CLR matrices); the port needs to use
   the "duplicate each edge both ways" convention PANINIpy documents
   for undirected input, and this should be tested against a known
   PANINIpy output on a toy graph to confirm the port is faithful
   before trusting it on real data.
4. **How this interacts with the sparse-store contract.** Per
   `CLAUDE.md`'s "Sparse network object" section, `compute_network()`
   stores a `dgCMatrix` above `store_density` and refuses analysis
   below `store_threshold`. If the MDL-selected backbone density falls
   below a network's `store_threshold`, the engine needs to either
   error clearly (matching `.net_check()`'s existing contract) or the
   calling convention needs to guarantee `store_density` is set high
   enough in advance — this should be designed, not discovered by a
   test failure.
5. **Attribution requirements.** MIT license permits reuse but still
   requires retaining the copyright/license notice; any ported
   function's roxygen header and/or `R/coexpressolog-strength.R`'s
   design-notes preamble should credit PANINIpy and Kirkley (2025)
   explicitly, following the same citation style already used
   throughout that file.

## 5. What is explicitly *not* done yet

- No R function, C++ kernel, roxygen doc, NAMESPACE export, or test
  exists for MDL backbone extraction anywhere in this repository.
- No design-notes header for an MDL engine exists in
  `R/coexpressolog-strength.R` (that file's current header covers only
  the SVN/mid-p density-profile design already implemented).
- No decision has been made on any of the five open questions in
  Section 4.
- No code has been copied, ported, or adapted from either PANINIpy or
  `aleckirkley/mdl-network-backbones` into this repository. The
  `functions.py` source was read in full in a prior session turn
  (for research purposes, to write this document accurately) but no
  lines from it were pasted or transcribed into `rcomplex`.

## 6. Suggested first implementation step (not yet started)

Once the open questions in Section 4 are resolved, the smallest
useful first step would be:

1. Implement `logchoose()` (trivial: `lchoose()` is already in base R)
   and the global-only `fglobal()` DL objective as a pure-R helper.
2. Implement the sort + greedy-scan + argmin loop as a pure-R function
   operating on a weight vector, validated against a hand-computed toy
   example (a handful of edges, computed by hand or against a PANINIpy
   run) before trusting it on real data.
3. Wire it into `suggest_reference_density()` (or a new function, per
   whichever integration point is chosen in Section 4.1) as an
   additional candidate density, reported alongside the existing
   scale-free-fit diagnostic rather than silently replacing it.
4. Add tests mirroring the rigor already used for
   `coexpressolog_strength()`: dense/sparse agreement if applicable,
   determinism (MDL backboning has no random component, so this should
   be exact reproducibility, not seed-scoped), and a toy-graph
   correctness check against a manually verified expected backbone.
5. Only after the R version is correct and tested should a C++ port
   be considered, and only if profiling shows it's needed at realistic
   gene-count edge lists.

## 7. Adjacent audit finding (not part of the MDL task, noted for continuity)

While checking the CRAN `SVN` package (the statistically-validated-
networks R package audited against `coexpressolog_strength()` in a
separate part of this session) its authorship was re-confirmed here:
the package is maintained by **Damien Challet**
(<https://CRAN.R-project.org/package=SVN>, package DOI
`10.32614/CRAN.package.SVN`), citing Tumminello et al. (2011, DOI
`10.1371/journal.pone.0017994`) as the methodological reference —
**Tumminello did not write the R package himself.** This corrects an
earlier loose phrasing in conversation ("the SVN package written by
Tumminello"). No code was borrowed from the `SVN` package; see the
session's separate SVN-audit notes for the full comparison.
