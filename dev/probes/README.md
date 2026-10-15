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
