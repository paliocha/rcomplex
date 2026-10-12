# Extract dgCMatrix slots from a dense matrix at a store threshold

Keeps the off-diagonal entries `>= thr` of both triangles and returns
the column-compressed slots of the corresponding `dgCMatrix`: 0-based
row indices `i` (ascending within each column), column pointers `p`
(length n + 1) and values `x`.

## Usage

``` r
extract_sparse_cpp(m, thr, n_cores = 1L)
```

## Arguments

- m:

  Dense numeric matrix (n x n)

- thr:

  Store threshold; entries below it are discarded

- n_cores:

  Number of threads for parallel computation (default 1)

## Value

List with integer vectors `i` and `p` and numeric vector `x`
