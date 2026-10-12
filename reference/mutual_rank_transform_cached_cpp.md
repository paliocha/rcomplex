# Cached mutual rank transformation

Transforms a correlation matrix using mutual rank normalization.

## Usage

``` r
mutual_rank_transform_cached_cpp(sim, log_transform = FALSE, n_cores = 1L)
```

## Arguments

- sim:

  Symmetric correlation/similarity matrix (n x n)

- log_transform:

  If FALSE (default), uses raw mutual rank with ascending ranks
  (original Rmd formula). If TRUE, applies Obayashi & Kinoshita (2009)
  log-normalized formula with descending ranks (values in 0 to 1 range).

- n_cores:

  Number of OpenMP threads (default: 1)

## Value

Mutual rank normalized matrix (n x n)

## Details

When `log_transform = TRUE`: \$\$S\_{ij} = 1 - \log(\sqrt{R\_{ij} \cdot
R\_{ji}}) / \log(n)\$\$ where ranks are descending (highest correlation
= rank 1).

When `log_transform = FALSE`: \$\$MR\_{ij} = \sqrt{R\_{ij} \cdot
R\_{ji}}\$\$ where ranks are ascending (lowest correlation = rank 1),
matching the original RComPlEx R Markdown formula.
