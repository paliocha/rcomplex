# Clique classifier consolidation

Date: 2026-09-22. Status: design note, nothing implemented. Martin:
"it is awkward to have two clique classifiers in this package."

## The two backends and who depends on them

| | species graph | gene graph |
|---|---|---|
| builder | `find_cliques()` (C++ Bron-Kerbosch, one best gene assignment per species clique) | `gene_clique_graph()` (every maximal clique, one per paralog combination) |
| classifier | `classify_cliques()`: complete / partial / differentiated / trait_specific / unclassified + `underpowered` flag (#24) | `classify_gene_cliques()`: complete_conserved / lineage_specific / partial_significant / partial_present / differentiated / underpowered (tier) / unclassified, with `missing_reason` |
| other consumers | `clique_stability()`, `clique_persistence()`, `clique_threshold_sweep()`, `clique_perturbation_test()`, `clique_intensity_test()`, the `.rcomplex` container | the `.rcomplex` container only |
| outside the package | nf-rcomplex calls `find_cliques()`, `clique_stability()`, `clique_hubs()` and hand-rolls its own taxonomy; it calls neither classifier | — |

So the species graph is load-bearing for every resampling test, but the
classifier on it has no dependants; the gene-graph classifier is the
published taxonomy (Rodriguez et al. 2026) and has no dependants either.

## Where the two classifiers disagree, and what the data said

1. **`underpowered`**: a flag beside the kept call in `classify_cliques()`
   (#24), a tier that replaces the call in `classify_gene_cliques()`.
2. **Missing species**: `classify_cliques()`'s `trait_specific` needs
   only a within-group clique; `classify_gene_cliques()`'s
   `lineage_specific` needs every outside species absent or untested,
   and `partial_present` refuses a species that was tested and failed
   (`missing_reason = tested_ns`) unless the failure was underpowered.
3. **Unit**: one best assignment per HOG versus every paralog
   combination. Only the gene graph can say which combinations form
   cliques, which the taxonomy needs.

Evidence (plan doc §10 and §11):

- Wood (EVOTREE): the gene-graph classifier reproduces every published
  tier exactly except partial_present (226 of 620; the 333 refused ones
  are `tested_ns` gaps — E3 says how many are genuinely well powered:
  [E3 pending]). Power gating moved 114 of 16,892 HOGs. The species-graph
  classifier flagged 319 of 416 `differentiated` and 648 of 7,817
  `trait_specific` HOGs as underpowered.
- Pooideae, analytical path with fold-enrichment power (P2): [P2 pending:
  species-graph flag table and gene-graph transitions, root and leaf].

## Options

**A. Keep both, align semantics.** Make `underpowered` a flag in the
gene graph too; give `trait_specific` the same absent/untested rule as
`lineage_specific`. Cheapest; leaves two names and two waterfalls.

**B. One classifier, on the gene graph.** `classify_gene_cliques()` is
the classifier; `classify_cliques()` is deprecated for one release with
a message, not a shim (its output shape differs). `find_cliques()` and
the resampling tests stay as they are: they answer "what happens to *a*
clique under resampling", which is a different question from
classification. Stability and persistence annotate a HOG, so they can
be joined onto the gene-graph result by `hog` instead of being an
argument of the classifier. `underpowered` becomes a flag here as well
(#24's argument applies: a replaced call is a lost call).

**C. One classifier, on the species graph.** Loses paralog
combinations, so it cannot produce the published taxonomy. Rejected.

## Recommendation

B. The gene-graph classifier already carries what the species-graph one
lacks (`missing_reason`, paralog combinations, the published tiers) and
nothing downstream needs the other. Sequence:

1. `underpowered` → flag column in `classify_gene_cliques()` (mirrors
   #24; breaking for code that filters on the tier).
2. Decide the partial_present rule from E3, and document it as the one
   rule.
3. Deprecate `classify_cliques()`; move its `stability` / `persistence`
   annotation to a join on `hog`.
4. Update CLAUDE.md "Two clique backends" to "one classifier, two
   builders", README, vignette.

Not before 0.4.0: 0.3.0 shipped today with both.
