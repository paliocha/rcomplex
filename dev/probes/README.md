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
