# Probes

Sources of functions cut from rcomplex 0.4.0 (tier C), kept for reference.
Each file starts with a header naming what it was.
Not built, not tested, not part of the package; no support.
Recover history with `git show main:<path>` (main @ b245ac0).
- `encodings-frc.R`: Forman-Ricci curvature per gene clique (note section 9.2).
- `encodings-twins.R`: false-twin classes and clique reduction (section 9.3).
- `encodings-common.R`: shared 12k setup for the two probes below (section 9.6).
- `encodings-rw.R`: equal-edge random walk on the mixed hypergraph (section 9.6).
- `encodings-orc.R`: hog-projected Ollivier-Ricci curvature and OT copy map (section 9.6).

- `copy-cap.R`, `copy-cap-exact.R`: is the 10-copy cap needed? Per-hog
  clique bounds and exact counts, capped against uncapped, on wood and
  leaf edge tables (`dev/design-notes/clique-budget.md`).

- `joint_louvain.cpp`, `joint_louvain.R`: probe E2 kernel; serial Louvain
  on per-species modularity plus a hog term with a species-blocked null.
- `joint_louvain_check.R`: E2 kernel acceptance (C++ against plain R,
  igraph at lambda 0, planted fixture, one wood run).
- `joint_planted.R`: E2 against two-step + match and E1 on a planted
  six-species model with defecting copies.
- `joint_halves.R`: the same engines on split halves of wood or leaf
  (split-half ARI, hog agreement, held-out AUROC z); per-run cache and an
  engine registry with `hleiden` / `bethe` hooks. Replaces
  `joint_wood_halves.R` (section 11 runs migrate into its cache).
- `bethe.R`, `bethe_check.R`: probe E3, species-blocked Bethe Hessian
  (K from negative eigenvalues); `joint_halves.R` hook. Fails check 1:
  K stays per species (78 = 6 x 13) or explodes past the cap.
- `hleiden.cpp`, `hleiden.R`, `hleiden_check.R`, `hleiden_disc.R`: WP15b
  probe, species-blocked h-Leiden on q(alpha) (hl0, hlS);
  `joint_halves.R` hook; disconnected-module count per cached run.
  Withdrawn (note section 12b).

- `structify.R`, `structify_check.R`, `structify_run.R`: Structify-Net
  toy networks on an 8-species tree with planted gains, losses, switches
  and duplications; scores `rcomplex()` tiers and `edge_history()`
  branches over a structure x noise grid (`dev/design-notes/structify-toy.md`).

- `peo/`: rcomplex 0.4.2 on the Plant Expression Omnibus (Koh et al.
  2023) on Orion. `sets.tsv` (species sets and traits),
  `hogs_at_level.py` (FastOMA hogs cut at the MRCA of a set),
  `stage.R`, `run.R`, `score.R`, and `jobs.sh` (push, submit, steps).
  Three experiments under `$ROOT/complex_<exp>/`: LeafvRoot (20 species;
  leaf, root, organ_block), C3vC4 (11 grasses, leaf), nodulation
  (12 species, root).

### PEO on Orion, 2026-10-16

Lib `$ROOT/peo/lib-0.4.2` (install job 1452230). Hogs job 1452229
(streaming orthoXML cut; pyham is in no env, and a new env would write
outside `$ROOT`). Stage jobs 1452231-1452233. Run jobs, 32 cores,
300 GB, 24 h each:

| job | experiment | run | state |
|---|---|---|---|
| 1452235 | LeafvRoot | leaf | cancelled |
| 1452236 | LeafvRoot | root | segfault (exit 139) |
| 1452237 | LeafvRoot | organ_block | segfault (exit 139) |
| 1452238 | C3vC4 | leaf | cancelled |
| 1452239 | nodulation | root | cancelled |
| 1452241 | LeafvRoot | leaf | blockwise rerun |
| 1452242 | LeafvRoot | root | blockwise rerun |
| 1452243 | LeafvRoot | organ_block | blockwise rerun |
| 1452244 | C3vC4 | leaf | blockwise rerun |
| 1452245 | nodulation | root | blockwise rerun |

The first round crashed at wheat root (71,330 genes): `rcomplex(expr)`
builds dense networks. `run.R` now builds each network with
`block_size = 512` and passes `networks =`. Diag job 1452240
(`diag.R`) runs the dense steps on the first n wheat genes.

Score after they finish: `jobs.sh submit score <exp>`.
