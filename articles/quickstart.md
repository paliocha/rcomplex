# Quickstart

This quickstart runs the whole analysis with one call,
[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md).
It answers one question: which genes keep their co-expression partners
in other species?

## Data

The package ships leaf and root RNA-seq from eight Pooideae grasses.
Each species has 20 leaf samples. We use four species: two annual and
two perennial.

``` r
library(rcomplex)
library(SummarizedExperiment)

se <- readRDS(system.file("extdata", "pooideae_vignette.rds",
  package = "rcomplex"
))
leaf <- lapply(
  se[c("BDIS", "BSYL", "HVUL", "HJUB")],
  function(s) s[, s$tissue == "leaf"]
)
```

[`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md)
needs two inputs. The first is one expression object per species, with
genes in rows and samples in columns. The second is a table of hogs:
ortholog groups with the columns `species`, `gene` and `hog`. For an
OrthoFinder file, call
[`read_orthologs()`](https://paliocha.github.io/rcomplex/reference/read_orthologs.md)
instead.

``` r
hogs <- do.call(rbind, lapply(names(leaf), function(s) {
  data.frame(
    species = s, gene = rownames(leaf[[s]]),
    hog = rowData(leaf[[s]])$hog
  )
}))
head(hogs, 3)
#>   species                 gene         hog
#> 1    BDIS Bradi0012s00100.v3.2 HOG:0033476
#> 2    BDIS    Bradi1g00480.v3.2 HOG:0031154
#> 3    BDIS    Bradi1g00552.v3.2 HOG:0027601
```

## Run

`clades` names the species sets. `null = TRUE` repeats the run with hog
labels shuffled among genes of the same expression level. The run takes
about 10 seconds.

``` r
res <- rcomplex(leaf, hogs,
  clades = list(
    annual = c("BDIS", "HVUL"),
    perennial = c("BSYL", "HJUB")
  ),
  null = TRUE, seed = 1
)
#> gene_clique_graph: 1245 genes collapsed into 577 twin classes; cliques carry multiplicity
res
#> rcomplex: 4 species, 1,500 hogs, 3 tiers
#>   BDIS  1,935 genes  20 samples  density 0.03  r >= 0.35
#>   BSYL  1,967 genes  20 samples  density 0.03  r >= 0.33
#>   HVUL  1,925 genes  20 samples  density 0.03  r >= 0.29
#>   HJUB  3,983 genes  20 samples  density 0.03  r >= 0.28
#>   edges 33,575 tested, 10,387 called at q < 0.1
#>   cliques 3,125 distinct (17,653 with paralog copies): complete_conserved 1,977, partial_present    68, unclassified 1,080
```

The header lists, per species, the genes, the samples and the density.
The value after “r \>=” is the weakest correlation that became an edge.
At 20 samples, a correlation of 0.3 is weak evidence. Read this value
first.

The last line counts the cliques by tier. A clique is a set of genes,
one per species, that are all co-expressed. A hog with paralogs gives
one clique per gene combination. Paralogs with the same partners are
twins. One clique stands for all their combinations. The line counts
distinct cliques, and the count in brackets includes the combinations.

Gene combinations can explode in a large hog. The clique search keeps at
most `max_cliques` cliques per hog (default 100,000). A hog over it
drops the genes with the weakest edges. `res$dropped` lists them.

## Read the results

[`summary()`](https://rdrr.io/r/base/summary.html) gives the tier table.
With `null = TRUE` it also gives the calls and cliques of the shuffled
hog labels.

``` r
summary(res)
#>                tier cliques cliques_null multiplicity hogs patterns
#>  complete_conserved    1977            0        16201  290      290
#>     partial_present      68            0          114   57       65
#>        unclassified    1080            0         1338  427      605
#> null: calls at q < 0.1 
#>  species1 species2 calls calls_null false_call_rate
#>      BDIS     BSYL  1340          0     0.000000000
#>      BDIS     HVUL  1239          2     0.001614205
#>      BDIS     HJUB  2090          0     0.000000000
#>      BSYL     HVUL  1134          8     0.007054674
#>      BSYL     HJUB  2056          0     0.000000000
#>      HVUL     HJUB  2528         22     0.008702532
```

`calls_null` and `cliques_null` count what abundance and network
structure alone produce. Compare them with `calls` and `cliques`.

`res$edges` has one row per ortholog pair and species pair.

``` r
head(res$edges[, c(
  "gene1", "gene2", "hog", "score", "evalue",
  "q_value", "effect_size", "power"
)], 4)
#>               gene1                gene2         hog     score      evalue
#> 1 Bradi5g10580.v3.2 Brasyl.9G105800.v2.1 HOG:0000016 11.550557   1.2461643
#> 2 Bradi2g26807.v3.2 Brasyl.8G139000.v2.1 HOG:0000024 12.106031   0.8479316
#> 3 Bradi3g42030.v3.2 Brasyl.3G295700.v2.1 HOG:0000036  8.322596  11.6758414
#> 4 Bradi3g13345.v3.2 Brasyl.3G006800.v2.1 HOG:0000038  2.866293 512.6243391
#>        q_value effect_size     power
#> 1 0.0011109275    3.170334 0.9889071
#> 2 0.0008152895    4.648396 0.8530924
#> 3 0.0090246912    3.401209 0.9738610
#> 4 0.2328586241    1.709270 0.9992404
```

- `score` is `-log2(p)`. Higher is stronger.
- `evalue` is the expected number of pairs this strong by chance.
- `q_value` is the false discovery rate.
- `effect_size` is the fold enrichment of shared partners.
- `power` is the chance to call the pair if it were conserved.

A low `power` means a non-significant pair is weak evidence of
divergence.

`res$classification` has one row per clique. `as.data.frame(res)`
returns the same table.

``` r
table(as.data.frame(res)$classification)
#> 
#> complete_conserved    partial_present       unclassified 
#>               1977                 68               1080
```

## Save

[`write_rcomplex()`](https://paliocha.github.io/rcomplex/reference/write_rcomplex.md)
writes `edges.tsv`, `cliques.tsv`, `members.tsv`, `twins.tsv` and
`classification.tsv`. After a null run it adds `edges_null.tsv` and
`classification_null.tsv`.

``` r
out <- file.path(tempdir(), "rcomplex-out")
write_rcomplex(res, out)
list.files(out)
#> [1] "classification_null.tsv" "classification.tsv"     
#> [3] "cliques.tsv"             "dropped.tsv"            
#> [5] "edges_null.tsv"          "edges.tsv"              
#> [7] "members.tsv"             "twins.tsv"
```

## Check the p-value floor

Run
[`pvalue_resolution()`](https://paliocha.github.io/rcomplex/reference/pvalue_resolution.md)
on any p-values you report. It shows how many values tie at the floor.

``` r
pvalue_resolution(res$edges$p_value)
#> p-value resolution
#>   values          33575
#>   distinct        26617 (79% of the values)
#>   minimum         1.44e-60 (unique)
#>   at 1            6297
#> 
#>   The minimum is untied, so the values still order the tests. Check again
#>   after any change to n_perm or to the number of tests.
#>   6297 values sit at 1. That is a real result -- no draw was less extreme
#>   than the observed statistic -- not a truncation, but those tests are
#>   unordered too.
```

## Next steps

- Designed experiment: pass `block` to remove a designed factor.
- Networks built elsewhere: pass `networks =` to
  [`rcomplex()`](https://paliocha.github.io/rcomplex/reference/rcomplex.md).
- Modules: pass `modules = TRUE`.
- Each step alone: see the walkthrough article.
