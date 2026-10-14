# Read an ortholog group file

Reads an ortholog group file into the long table. The table has one row
per gene: `species`, `gene` and `hog`.

## Usage

``` r
read_orthologs(
  file,
  species = NULL,
  format = c("auto", "orthofinder", "plaza", "long")
)
```

## Arguments

- file:

  Path to a tab-delimited ortholog file. It can be gzipped.

- species:

  Character vector of species to keep, or `NULL` for all species in the
  file.

- format:

  One of `"auto"`, `"orthofinder"`, `"plaza"` or `"long"`. `"auto"`
  reads the format from the header.

## Value

A data frame with columns `species`, `gene` and `hog`.

## Details

Three formats are read:

- orthofinder:

  OrthoFinder `N0.tsv` (HOGs, column `HOG`) or `Orthogroups.tsv` (column
  `Orthogroup`). The column name may use any case, and `N1.HOG` is read
  too. There is one column per species, with genes separated by commas.
  The columns `OG` and `Gene Tree Parent Clade` are skipped. Empty and
  `NA` cells hold no gene. The reader keeps gene ids as written; map
  suffixes such as `.p1` to expression ids yourself.

- plaza:

  PLAZA orthologs: columns `species`, `gene_id` and `gene_content`,
  where `gene_content` lists the orthologs of the anchor gene as
  `code:gene1,gene2;code:gene3`. PLAZA lists pairwise orthologs; a hog
  is a connected component of those lists over the whole file, so every
  gene is in exactly one hog. Hogs are integers, numbered by the
  smallest `species` and gene id they contain.

- long:

  The output shape itself: columns `species`, `gene` and `hog`. Other
  columns are dropped.

`"auto"` picks `long` when the header has `species`, `gene` and `hog`,
`orthofinder` when it has `HOG` or `Orthogroup`, and `plaza` when it has
`gene_content`. FastOMA and eggNOG output have no reader: write them as
the long table, one row per gene and group.

## Examples

``` r
f <- system.file("extdata", "N0.tsv", package = "rcomplex")
long <- read_orthologs(f, species = c("SpA", "SpB"))
head(long)
#>   species      gene           hog
#> 1     SpA GeneA_001 N0.HOG0000001
#> 2     SpA GeneA_002 N0.HOG0000002
#> 3     SpA GeneA_012 N0.HOG0000002
#> 4     SpA GeneA_003 N0.HOG0000003
#> 5     SpA GeneA_004 N0.HOG0000004
#> 6     SpA GeneA_005 N0.HOG0000005
prepare_orthologs(long)
#>        gene1     gene2           hog species1 species2
#> 1  GeneA_001 GeneB_001 N0.HOG0000001      SpA      SpB
#> 2  GeneA_001 GeneB_011 N0.HOG0000001      SpA      SpB
#> 3  GeneA_002 GeneB_002 N0.HOG0000002      SpA      SpB
#> 4  GeneA_012 GeneB_002 N0.HOG0000002      SpA      SpB
#> 5  GeneA_003 GeneB_003 N0.HOG0000003      SpA      SpB
#> 6  GeneA_004 GeneB_004 N0.HOG0000004      SpA      SpB
#> 7  GeneA_005 GeneB_005 N0.HOG0000005      SpA      SpB
#> 8  GeneA_006 GeneB_006 N0.HOG0000006      SpA      SpB
#> 9  GeneA_007 GeneB_007 N0.HOG0000007      SpA      SpB
#> 10 GeneA_008 GeneB_008 N0.HOG0000008      SpA      SpB
#> 11 GeneA_009 GeneB_009 N0.HOG0000009      SpA      SpB
#> 12 GeneA_010 GeneB_010 N0.HOG0000010      SpA      SpB
```
