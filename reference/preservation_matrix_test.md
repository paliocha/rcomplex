# Trait relabelling test on the all-pairs preservation matrix

Tests whether modules are less well preserved between species that
differ in a trait than between species that share it. Uses every
contrast in an all-pairs
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
run.

## Usage

``` r
preservation_matrix_test(classification, group, block = NULL, seed = NULL)
```

## Arguments

- classification:

  The `classification` data frame from
  [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md),
  ideally over an all-pairs `pairs` table (see
  [`all_species_pairs()`](https://paliocha.github.io/rcomplex/reference/all_species_pairs.md)).
  Must carry `reference`, `test` and the effect column `Zsummary_std`;
  `q_value` is used only for the `$saturation` diagnostic. Rows whose
  effect is `NA` — `"untested"` modules, and any module whose null
  correlation was unavailable — are dropped, since nothing was measured
  for them.

- group:

  Named vector mapping species to trait values. Any number of levels is
  allowed. Species not appearing in `classification` are ignored, so the
  permuted label multiset is the one the matrix actually holds.

- block:

  Optional named list of species vectors, one per clade (a genus, say),
  as `clades` elsewhere. Only the top-level clades are read. A species
  in no clade forms its own clade. Supplying it enables the restricted
  within-block null and the within-block exclusion. With `block = NULL`
  only the free null is run.

- seed:

  Integer seed for the sampled branch, or `NULL` (default) to draw from
  the ambient stream and leave it advanced. A seed draws from a private
  stream and restores the caller's on exit, the package-wide contract
  described under
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).
  An enumerated null draws nothing, so a seed changes nothing there.

## Value

A list with components:

- observed:

  The statistic under the true labelling.

- form:

  `"difference"` for two trait levels, `"dispersion"` otherwise.

- class_means:

  Data frame with `class`, `n` and `mean_z`: the observed class means
  the statistic is built from.

- rows_per_pair:

  Data frame with `species1`, `species2`, `concordant`, `same_block` and
  `n`, one row per species pair contributing.

- free:

  List for the free null: `null_distribution`, `p_value`, `p_min`,
  `p_attainable`, `n_tied_max`, `n_distinct`, `rank`, `exact`,
  `n_labellings`, `n_scored`. `p_attainable` is the tie-aware floor —
  labellings sharing a statistic cannot be separated, so nothing scores
  below the tail of the maximum — and `n_tied_max` says how many
  labellings that tail holds. When `exact` is `FALSE` both are counted
  over the `n_scored` draws rather than over the `n_labellings` of the
  space.

- blocked:

  The same list for the within-block null, or `NULL` when `block` was
  not supplied.

- p_free,p_blocked:

  The two p-values, lifted out for convenience. `p_blocked` is `NA`
  without a `block`.

- saturation:

  List describing the resolution of the supplied q-values: `n_tests`,
  `n_distinct`, `q_floor`, `n_at_floor` and `n_at_one`. `n_distinct`,
  `n_at_floor` and `n_at_one` are counted with the tolerance
  [`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)
  uses, so the two agree on one q-vector — relative for `n_distinct` and
  `n_at_floor`, and a one-sided window below 1 for `n_at_one`, which
  counts any q-value at or above 1 so that invalid input is reported
  rather than dropped: values within a few last bits of each other are
  one value wobbled by arithmetic, not two resolvable ones. `n_tests`
  counts every row of the supplied `classification` carrying a finite
  q-value — the multiplicity population is the whole matrix, so this is
  deliberately *not* `n_rows` below, which is the subset this test
  averages over after dropping unmeasured rows and within-block pairs.
  Reported next to the p-values so the reader can see how much of the
  input's significance scale is collapsed.

- n_rows,n_pairs,n_excluded:

  Rows used, distinct species pairs behind them, and rows dropped for
  sharing a block.

- species,group,block:

  The design actually tested. `block` is the top-level clade of each
  species.

## Details

The null relabels species. It holds the label counts and the whole
preservation matrix fixed.

## Why the effect size and not the q-value

The statistic ranks contrasts by `Zsummary_std`, never by `p_value` or
`q_value`. A permutation p-value saturates: on the eight-species
Pooideae set with 2000 permutations and 511 module-directions, `q_value`
takes 172 distinct values, 35 modules sit tied at the floor of 0.00071
and 20 at exactly 1. Among those 35 floored modules `Zsummary_std` spans
6.4 to 66.7 — a tenfold range of effect that the p-value cannot see at
all. Anything that ranks, weights or orders modules must therefore read
the standardised effect, which is continuous; p and q are for the
significance *call* only. `$saturation` reports how far the supplied
q-values are collapsed, so a caller who reaches for them anyway can see
what resolution is left.

## The statistic

Each row of `classification` is one module in one direction of one
contrast, and its two species are trait-*concordant* when they share a
trait value and *discordant* otherwise. With two trait levels the
statistic is
`mean(Zsummary_std | concordant) - mean(Zsummary_std | discordant)`, so
a **positive** value means discordant pairs are the less preserved ones,
which is the biological hypothesis, and the test is one-sided upward.

With more than two levels a single difference no longer describes the
design, because discordance is no longer one thing. The classes are then
the concordant rows (pooled across levels) plus one class per unordered
pair of distinct levels, and the statistic is the row-weighted standard
deviation of the class means — their between-class dispersion, on the
`Zsummary_std` scale. Larger is still more structure, so the test stays
one-sided upward, but it is no longer signed and does not say which
direction the difference runs; read `$class_means` for that. Because the
concordant classes are pooled, a two-level trait has exactly two classes
and the dispersion form is never used for it: the binary case is the
difference above, exactly.

Every module-direction counts once, so a contrast with many modules
contributes more than one with few. That is deliberate — the unit is the
module, and a species pair with more modules carries more evidence — but
it means a single hugely partitioned species can dominate; check
`$rows_per_pair`.

## The two nulls

`p_free` permutes trait labels over all species, preserving the label
counts (`choose(8, 4) = 70` labellings for a balanced binary trait on
eight species). `p_blocked` permutes only *within* each `block`, which
holds the phylogeny fixed and is the conservative, phylogenetically
controlled null (16 labellings for four genera of two). Both are
reported because their agreement is the diagnostic: a free p-value far
smaller than the blocked one means the larger label space bought its
resolution by breaking phylogenetic control rather than by measuring the
trait.

Neither floor is `1 / n_labellings`, because relabellings that leave
every pair's concordance alone reproduce the statistic exactly. Swapping
the two trait values everywhere is one such, so a binary design always
carries at least a twofold tie at the maximum: 2/70 = 0.029 free on
eight species split 4/4, and 2/16 = 0.125 within four blocks of two —
which is why the restricted null of a four-genus paired design cannot
reach alpha on its own. With more levels the dispersion cannot tell the
level *names* apart either, so every renaming of the levels that
preserves the label counts ties as well. That is `g!` renamings only
when the design is balanced; in general it is
`prod(factorial(table(table(labels))))`, since a renaming has to map
each level onto one of the same size. Balanced three levels on six
species tie 6 of 90 labellings, a floor of 0.067; an unbalanced 2/2/1
split on five species ties only 2 of 30 — the same 0.067, but from two
ties rather than six. `p_attainable` reports whichever floor actually
binds and `n_tied_max` how many labellings share it; the function warns
when that floor exceeds 0.05, since no signal can then reach
significance.

Both floors describe an *enumerated* null. When a space is too large to
enumerate, `p_min`, `p_attainable` and `n_tied_max` are counted over the
10000 draws instead: `p_min` is then the sampling floor `1 / 10001`,
which can sit far below the labelling-space floor the design really
imposes, and `n_tied_max` counts drawn labellings with repeats rather
than distinct ones. `exact` says which regime a null is in, and only an
exact one bounds the design.

## Multiplicity scope of the q-values

[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
corrects within each contrast: `q_value` is Benjamini-Hochberg over that
contrast's own modules, a dozen or two tests. An all-pairs matrix
reframes those contrasts as one analysis — every module-direction tests
the same hypothesis — so the per-contrast correction under-corrects, by
roughly the number of contrasts. The relabelling test itself is
untouched — it never reads a q-value — but any per-module significance
call taken off the same matrix is.

## Why within-block pairs are excluded

In a paired design such as the Pooideae one — each genus contributing an
annual and a perennial — every within-genus pair is trait-discordant
while every trait-concordant pair is between-genus. Trait status and
phylogenetic distance are then perfectly confounded within that subset,
and the confound runs *against* the hypothesis: close relatives are
better preserved, so the within-genus discordant pairs are pulled upward
and mask a real effect. Those rows are dropped. The exclusion is by
block membership, which no relabelling changes, so the same rows are
excluded under every labelling and the null stays valid.

## See also

[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md),
[`all_species_pairs()`](https://paliocha.github.io/rcomplex/reference/all_species_pairs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pairs <- rcomplex:::all_species_pairs(names(group))
res <- preservation_paired(mods, nets, ortho, pairs, group = group)
pmt <- preservation_matrix_test(res$classification, group,
  block = split(names(genus), genus)
)
c(pmt$observed, pmt$p_free, pmt$p_blocked)
} # }
```
