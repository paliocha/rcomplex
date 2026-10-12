# Run module preservation across many species pairs

Applies
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
to each contrast in `pairs`. Runs both directions and reports them
separately.

## Usage

``` r
preservation_paired(modules, ...)

# Default S3 method
preservation_paired(
  modules,
  networks,
  orthologs,
  pairs,
  group = NULL,
  edges = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- modules:

  Named list of module assignments
  ([`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md)
  or
  [`as_modules()`](https://paliocha.github.io/rcomplex/reference/as_modules.md)),
  keyed by species.

- ...:

  Additional arguments passed to the default method.

- networks:

  Named list of
  [`compute_network()`](https://paliocha.github.io/rcomplex/reference/compute_network.md)
  results, keyed by species.

- orthologs:

  Data frame with columns `gene1`, `gene2`, `hog`.

- pairs:

  Data frame with columns `species1`, `species2` and optionally
  `pair_name`.

- group:

  Optional named vector mapping species to a clade. When supplied, a
  `group` column records `"conserved"` for preserved modules and the
  owning species' group for diverged ones. Modules reported `"untested"`
  are counted under `"untested"` rather than a clade: nothing was
  measured, so attributing them to one would overstate the evidence.

- edges:

  Optional
  [`find_coexpressologs()`](https://paliocha.github.io/rcomplex/reference/find_coexpressologs.md)
  results, used for paralog resolution.

- seed:

  Integer seed for the whole run, or `NULL` (default) to draw from the
  ambient stream and leave it advanced. The seed is applied once here
  and the per-direction
  [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
  calls are left unseeded, so the directions draw in sequence from one
  stream instead of every one of them reusing the same permutations. A
  seeded call restores the caller's stream on exit, the package-wide
  contract described under
  [`detect_modules()`](https://paliocha.github.io/rcomplex/reference/detect_modules.md).

  Before 0.3.0 `seed` reached
  [`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
  through `...` and handed every direction the identical seed.

## Value

A list with `classification` (one row per module per direction, carrying
`pair_name`, `module`, `species`, `reference`, `test`, `classification`
and both ordering columns, `Zsummary_std` and `medianRank` – the latter
ranks within a direction, so it is not comparable across the rows of
this table), `summary` (counts per contrast and direction) and `raw`
(the
[`module_preservation()`](https://paliocha.github.io/rcomplex/reference/module_preservation.md)
results, keyed by `"<reference>.<test>"`).

## Details

Preservation is directional. Whether the modules of species A survive in
B differs from the reverse.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- preservation_paired(mods, nets, ortho,
  pairs = data.frame(species1 = "BDIS", species2 = "BSYL"),
  group = c(BDIS = "annual", BSYL = "perennial")
)
res$summary
} # }
```
