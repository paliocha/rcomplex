# Write an rcomplex result as tab-separated files

`write_rcomplex()` writes each table of the result to its own file. Each
file takes the name of its element: `edges.tsv`, `cliques.tsv`,
`members.tsv`, `twins.tsv`, `capped.tsv`, `classification.tsv` and, with
`null = TRUE`, `edges_null.tsv` and `classification_null.tsv`.

## Usage

``` r
write_rcomplex(x, dir)
```

## Arguments

- x:

  A result of
  [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md).

- dir:

  Output directory. It is created when it does not exist.

## Value

The file paths, invisibly.

## See also

The example of
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md).
