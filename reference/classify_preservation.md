# Classify modules as conserved, moderately preserved, or diverged

Turns
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
output into a call per module. The combined permutation q-value drives
the call.

## Usage

``` r
classify_preservation(pres, species = NA_character_, pair_name = NA_character_)
```

## Arguments

- pres:

  Output of
  [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md).

- species:

  Optional species label recorded in the `species` column.

- pair_name:

  Optional contrast label recorded in the `pair_name` column.

## Value

A data frame with `module`, `species`, `pair_name`, `classification`,
`Zsummary`, `Zsummary_std`, `medianRank`, `q_value`, `size` and
`size_mapped`. `Zsummary_std` is the column the default criterion is
read against, so reproducing the call from `Zsummary` alone will not
match. `medianRank` is carried through from
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
and takes no part in the call; see *Ordering modules* for when to read
it instead.

## Details

The `Zsummary` thresholds of Langfelder et al. (2011) split the
significant modules a second time, so the familiar 2 / 10 cut points
still appear.

## Criteria

- conserved:

  `q_value < 0.1` and `Zsummary_std` is at or above 10

- moderate:

  `q_value < 0.1` and it is below 10

- diverged:

  `q_value >= 0.1`

- untested:

  `q_value` is `NA` – a statistic could not be computed, so neither
  preservation nor divergence was measured

A module whose degree correlation could not be computed – `cor.degree`
is undefined when intramodular connectivity is constant in either
network – gets `NA` for that statistic, so the `pmax` combination and
hence `q_value` are `NA` too. Such a module is reported `"untested"`
rather than `"diverged"`: nothing was measured, so divergence would be a
positive claim the data does not support. It is never called preserved
on the density statistic alone.

## Ordering modules

`Zsummary_std` decides the call, but it is not the only ordering worth
reading. It is a Z-score, so its denominator is a permutation standard
deviation, and that shrinks as a module grows: over 110 modules from 12
Pooideae contrasts the null sd of `cor.degree` tracks mapped module size
at Spearman -0.996 and that of `avg.weight` at -0.873. A large, weakly
preserved module can therefore outrank a small, strongly preserved one.

`medianRank` is the size-independent complement Langfelder et al. (2011)
report next to `Zsummary` for exactly that reason. It ranks the observed
`avg.weight` and `cor.degree` across the tested modules, 1 = strongest,
and averages the two ranks; no permutation moment enters it, so module
size cannot set its scale. Prefer `medianRank` when asking which modules
are the best preserved *relative to each other*, and `Zsummary_std` when
asking how far from its own null any one module sits. The two are not
redundant: on those same 12 contrasts they ordered 11% of module pairs
differently.

`medianRank` ranks across the modules of *one run*. Its scale is set by
how many modules that run tested, so it is not comparable between runs –
between the two directions of a contrast, or between contrasts in a
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
table – unless those runs tested equally many modules. It takes no part
in the classification.

## Examples

``` r
if (FALSE) { # \dontrun{
cls <- classify_preservation(pres)
table(cls$classification)
} # }
```
