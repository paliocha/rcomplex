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
  are `tested_ns` gaps). E3 (`prepare_data/evotree/e3_partial_present.R`)
  looked at the sixth species' tested edges to the five members in all
  8,944 such cliques: the maximum power among them is 1.0 in three
  quarters of cliques (minimum 0.33), the smallest q is 0.82 at the
  median, and `min_power` rescues 12 HOGs at 0.8, 17 at 0.9, 45 at 0.99.
  So the refusal is negative evidence, not a power artefact: the sixth
  species was there, was tested at full power, and its co-expression was
  not conserved with the members. The rule stands. (Aside: the missing
  species is aspen in 3,391 of the 8,944 cliques, birch and cherry about
  2,200 each, the conifers 187-604; aspen's data come from a different
  study than the other five.). Power gating moved 114 of 16,892 HOGs. The species-graph
  classifier flagged 319 of 416 `differentiated` and 648 of 7,817
  `trait_specific` HOGs as underpowered.
- Pooideae, analytical path with fold-enrichment power (P2, plan doc
  §12): the species-graph classifier calls 297 (root) / 384 (leaf)
  `trait_specific` HOGs and flags 44% / 31% of them at `min_power = 0.8`;
  no HOG is `differentiated`. The gene-graph classifier finds 0 / 4
  `lineage_specific` and 0 `differentiated` HOGs, and power gating moves
  nothing: the other trait group is present and tested, so
  `lineage_specific`'s absent-or-untested rule never fires. The gene
  graph has no tier for the study's actual question, one trait group
  conserved while the other is present but not co-conserved.

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

B, with one addition the Pooideae run makes non-negotiable: the gene
graph must first gain a `trait_specific` tier, or `classify_cliques()`'s
only unique signal disappears. Definition, mirroring the species-graph
rule but with the gene graph's evidence types: a complete clique over
one trait group, every outside species present and `tested_ns` (an
absent or untested outside species is `lineage_specific`, as now), no
within-group clique in the other group (else `differentiated`), and the
`underpowered` flag when any deciding outside edge is below `min_power`.
Then the gene-graph classifier carries everything the species-graph one
does, plus `missing_reason` and paralog combinations. Sequence:

0. Add `trait_specific` to `classify_gene_cliques()` as defined above;
   check it reproduces the 297 / 384 Pooideae HOGs and the flag shares.
1. `underpowered` → flag column in `classify_gene_cliques()` (mirrors
   #24; breaking for code that filters on the tier).
2. Keep the partial_present rule (E3: the refused gaps are well-powered
   negative evidence) and document it as the one rule.
3. Deprecate `classify_cliques()`; move its `stability` / `persistence`
   annotation to a join on `hog`.
4. Update CLAUDE.md "Two clique backends" to "one classifier, two
   builders", README, vignette.

Not before 0.4.0: 0.3.0 shipped today with both.
