# All unordered species pairs as a `pairs` data frame

Builds the `pairs` argument of
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md)
for every `choose(n, 2)` contrast. Saves writing an all-pairs table by
hand.

## Usage

``` r
all_species_pairs(species, sep = ".")
```

## Arguments

- species:

  Character vector of species identifiers. A *named* vector — a `group`
  map, say — is read as its names, since a bare species vector is not
  normally named.

- sep:

  Separator used to build `pair_name` (default `"."`), matching the
  `paste(species1, species2, sep = ".")` default of
  [`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md),
  so the two agree when `pair_name` is left implicit.

## Value

A data frame with columns `species1`, `species2` and `pair_name`, one
row per unordered pair. Both directions of each contrast are run by
[`preservation_paired()`](https://paliocha.github.io/rcomplex/reference/preservation_paired.md),
so each pair appears once.

## Details

Running all pairs rather than a designated few is what makes
[`preservation_matrix_test()`](https://paliocha.github.io/rcomplex/reference/preservation_matrix_test.md)
worth doing. Every extra contrast adds its modules to the class means
the statistic is built from. Any species the designated few left out
widens the relabelling null as well. Eight species split 4/4 give
`choose(8, 4) = 70` free labellings, against the 16 that a within-genus
null over four genera of two can reach.

## Examples

``` r
rcomplex:::all_species_pairs(c("BDIS", "BSYL", "HVUL"))
#>   species1 species2 pair_name
#> 1     BDIS     BSYL BDIS.BSYL
#> 2     BDIS     HVUL BDIS.HVUL
#> 3     BSYL     HVUL BSYL.HVUL
```
