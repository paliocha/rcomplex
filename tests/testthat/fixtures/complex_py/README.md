# complex_py fixture

Cross-implementation equivalence fixture for `tests/testthat/test-equivalence.R`:
rcomplex must reproduce the co-expressolog calls of canonical ComPlEx
([natstreet/ComPlEx_python](https://github.com/natstreet/ComPlEx_python))
on this data.

## Files

| File | What |
|------|------|
| `make_fixture.R` | Generator (reads no external input). Seed 20260901; 8 co-expression modules x 30 samples. sp1 = 320 module genes + 60 noise genes; sp2 = 352 module genes (4 extra paralogs per module) + 60 noise genes. Writes `sp1_expr.tsv`, `sp2_expr.tsv`, `ortho_pairs.tsv` and `ortho_long.tsv` (ComPlEx_python ortholog format; not committed). |
| `sp1_expr.tsv`, `sp2_expr.tsv` | Expression matrices, genes x samples, first column `Genes`. |
| `ortho_pairs.tsv` | Ortholog pairs (`Species1`, `Species2`, `hog`) in `compare_neighborhoods()` format: 320 one-to-one ortholog groups, plus 32 sp2 paralogs joining the first 4 groups of each module; 64 sp2 partners randomly permuted among themselves (most, not all, land in a different module; fixed points stay conserved). |
| `expected_calls.tsv` | ComPlEx_python output `pyout/RData/comparison_tables/comparison_sp1_sp2.tsv`: the 149 pairs called conserved (`Max.p.val` = pairwise max of BH-adjusted p-values, < 0.05). |

ComPlEx_python builds each network on the ortholog-restricted expression
table, so the test restricts each species to the genes present in
`ortho_pairs.tsv` before calling `compute_network()`.

## Generation environment

The committed TSVs were written with R version 4.6.1 (2026-06-24),
`RNGkind()` = `Mersenne-Twister`, `Inversion`, `Rejection` (R defaults).
`test-equivalence.R` re-runs `make_fixture.R` in a temporary directory and
checks `ortho_pairs.tsv` byte-for-byte (strings only) and the two
expression TSVs as parsed numbers (`read.delim`, tolerance 1e-12).
Byte-identical expression TSVs are only guaranteed on the reference
platform (arm64/clang, md5s below); x86/gcc agrees to within 1e-12
relative (last-ULP differences in some printed doubles).
`tools::md5sum()` of the committed files:

| File | md5 |
|------|-----|
| `sp1_expr.tsv` | `24573a93d4fd715646ec41aff99aa8bd` |
| `sp2_expr.tsv` | `a3964101db876045981d274fcac0f416` |
| `ortho_pairs.tsv` | `355266f7e104898971c4d626b2bc273c` |
| `expected_calls.tsv` | `8e61c790aecbbd358471adbe2626473c` |

## Regenerate

Run `make_fixture.R`, then the Python command:

1. In an empty directory: `Rscript make_fixture.R` (writes the three input
   TSVs and `ortho_long.tsv`).
2. Run ComPlEx_python at commit `1838a5c`. `statsmodels` is not needed: put
   the BH shim below on `PYTHONPATH` (identical to R `p.adjust(method = "BH")`).

```bash
git clone https://github.com/natstreet/ComPlEx_python
git -C ComPlEx_python checkout 1838a5c
mkdir -p shim/statsmodels/stats
touch shim/statsmodels/__init__.py shim/statsmodels/stats/__init__.py
# write shim/statsmodels/stats/multitest.py (below), then:
PYTHONPATH=shim python3 ComPlEx_python/complex_py.py \
  --s1-expr sp1_expr.tsv --s2-expr sp2_expr.tsv --orthologs ortho_long.tsv \
  --s1-name sp1 --s2-name sp2 --out-dir pyout --density 0.03 --min-expr=-1000
cp pyout/RData/comparison_tables/comparison_sp1_sp2.tsv expected_calls.tsv
```

```python
# shim/statsmodels/stats/multitest.py  -- BH identical to R p.adjust(method = "BH")
import numpy as np
def multipletests(pvals, alpha=0.05, method="fdr_bh", **kw):
    p = np.asarray(pvals, dtype=float); n = p.size
    order = np.argsort(p)[::-1]
    ranked = p[order] * n / np.arange(n, 0, -1)
    adj = np.minimum(1.0, np.minimum.accumulate(ranked))
    out = np.empty(n); out[order] = adj
    return out < alpha, out, None, None
```

3. Copy `sp1_expr.tsv`, `sp2_expr.tsv`, `ortho_pairs.tsv` and
   `expected_calls.tsv` into this directory.
