# In-place mutual rank transformation

Overwrites `sim` with its mutual rank transform without allocating any n
x n temporaries. Pass 1 clamps each column to \[-1, 1\], optionally
negates it, and replaces the column by its average ranks. Pass 2
replaces each pair (i, j) by sqrt(R_ij \* R_ji) (log-normalized when
`log_transform`) and sets the diagonal to 0. Same formulas and tie
handling as
[`mutual_rank_transform_cached_cpp()`](https://paliocha.github.io/rcomplex/reference/mutual_rank_transform_cached_cpp.md),
which is kept as the reference implementation.

## Usage

``` r
mutual_rank_inplace_cpp(sim, log_transform, negate, n_cores)
```

## Arguments

- sim:

  Symmetric correlation/similarity matrix (n x n). Modified in place:
  must be a fresh allocation not shared with any other R object, and
  must be stored as double (`REALSXP`); any other storage type is an
  error, since coercion would silently operate on a copy. Must not
  contain `NaN`/`NA` (ranking would be undefined); an error is raised
  and the contents of `sim` are then unspecified.

- log_transform:

  If FALSE, raw mutual rank with ascending ranks (original Rmd formula).
  If TRUE, Obayashi & Kinoshita (2009) log-normalized formula with
  descending ranks (values in 0 to 1 range).

- negate:

  If TRUE, negate every value but the diagonal before ranking, so the
  strongest anticorrelation ranks first.

- n_cores:

  Number of OpenMP threads

## Value

Invisible `NULL`; `sim` is modified in place.
