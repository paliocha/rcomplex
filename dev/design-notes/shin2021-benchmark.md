# Shin 2021 reference benchmark

Probe, 2026-10-16. Scripts in `dev/probes/peo/`: `shin_tables.R`,
`shin_prep.py`, `shin_fastoma.sh`, `shin_nextflow.config`,
`stage_shin.R`, and the `fetch`, `prep`, `fastoma`, `shin_hogs`,
`stage`, `run` and `score` steps of `jobs.sh`. Orion directory
`$ROOT/complex_Shin2021/`.

## The paper

Shin et al. 2021, Nucleic Acids Res 49:e3 (doi:10.1093/nar/gkaa1041).
Six plants with proteomes measured on matched tissues.
- Arboretum-proteome clusters protein levels into seven modules per
  species, tied across the species tree (module 1 lowest, 7 highest).
- Clade-specific gene sets: orthogroups whose module pattern changes
  on one branch. Rule-based (23 sets) and de novo (286 sets).
- Muscari: joint spectral clustering of per-species co-expression
  networks with a phylogenetic prior.
Their validation runs from protein-level sets to co-expression. Ours
runs the other way: co-expression edges to tiers, then a check against
their protein-level sets.

## Data and ids

Six species, Ensembl Plants release 38 gene ids: Arabidopsis_thaliana
(`AT1G01020`), Medicago_truncatula (`MTR_5g009740`), Oryza_sativa
(RAP, `Os03g0309000`), Physcomitrella_patens (`PP1S7_216V6`),
Solanum_tuberosum (`PGSC0003DMG400022701`), Zea_mays (AGPv4,
`Zm00001d046467`). Tree, rooted on moss:
`((((mtr,ath)Anc4,stu)Anc3,(osa,zma)Anc5)Anc2,ppa)Anc1;`; it matches
their `SpeciesTree.txt`. Expression: the six Mendeley compendia
(doi:10.17632/bwzdv623xn, version 2), genes x SRA experiments, used
as published: no log, no level filter. Only genes with zero variance
or absent from the proteome are dropped. Reference tables: supplements
S4 (Arboretum modules), S5 (gene sets and their patterns), S6
(Muscari modules) and S8 (curated Medicago symbiosis genes), tidied
by `shin_tables.R` into `inputs/ref/`.

## Orthology: FastOMA on release 38

Their orthogroup file is not published, and our hog ids must match
their gene ids. So FastOMA v0.5.1 runs on the six release-38
proteomes: one protein per gene (the longest), header = gene id.
`hogs_at_level.py` cuts at the root, Anc1. A root hog born below the
root becomes one young hog (`Anc1:y<n>`), so lineage-specific families
stay in. `hog_age.tsv` gives each hog the TaxRange of its root hog.

## Design

One run (`all`): networks built with `block_size = 512` and passed as
`networks =`, `clades = clades_from_tree(tree)`, `null = TRUE`,
density 0.03, alpha 0.1, seed 1, 32 cores, 100 GB, 12 h. Then
`edge_history()`. No trait pass. `score.R Shin2021` writes the Part B
tables of every experiment (tier by age, signatures, tier by
expression quintile, tier by copies) and these:
- `ref_modules.tsv`: `hog_class` x Arboretum module (modal module of a
  hog's members, ties to the highest), with coverage.
- `ref_genesets.tsv`: tier mix for hogs in conserved, loss and change
  sets against all tested hogs, with the log2 ratio.
- `ref_genesets_detail.tsv`: per set, hogs, cliques, tier mix and the
  modal `edge_history()` signature.
- `curated_nodulation.tsv`: the S8 list. Medicago is the only legume,
  so nodulation hogs should sit outside complete cliques.

Expectation from the paper: loss sets go to `lineage_specific` or
partial tiers, change sets to `differentiated`, conserved sets to
`complete_conserved`; conserved tiers sit in high expression
quintiles.

## Results (pending)

| species | experiments | genes kept | genes with hog | in proteome |
|---|---|---|---|---|
| pending | | | | |

| pattern type | hogs | complete | lineage | partial | differentiated |
|---|---|---|---|---|---|
| pending | | | | | |

| Arboretum module | hogs | share conserved |
|---|---|---|
| pending | | |

## Related

Li, Nakayama & Sinha 2023 (Curr Opin Plant Biol 76:102474) review
comparative co-expression across plants and place this work among the
network-alignment approaches. Heterochrony is a user-side step: when
the same organ is sampled at stages that do not line up across
species, align the stages first. `block` centres a designed factor
within each species; it does not align stages across species.
