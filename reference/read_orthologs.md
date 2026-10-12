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
  `Orthogroup`). There is one column per species, with genes separated
  by commas. The columns `OG` and `Gene Tree Parent Clade` are skipped.

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

`"auto"` picks `orthofinder` when the header has `HOG` or `Orthogroup`,
`plaza` when it has `gene_content`, and `long` when it has `species`,
`gene` and `hog`. FastOMA and eggNOG output have no reader: write them
as the long table, one row per gene and group.

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
#>        gene1     gene2           hog
#> 1  GeneA_001 GeneB_001 N0.HOG0000001
#> 2  GeneA_001 GeneB_011 N0.HOG0000001
#> 3  GeneA_002 GeneB_002 N0.HOG0000002
#> 4  GeneA_012 GeneB_002 N0.HOG0000002
#> 5  GeneA_003 GeneB_003 N0.HOG0000003
#> 6  GeneA_004 GeneB_004 N0.HOG0000004
#> 7  GeneA_005 GeneB_005 N0.HOG0000005
#> 8  GeneA_006 GeneB_006 N0.HOG0000006
#> 9  GeneA_007 GeneB_007 N0.HOG0000007
#> 10 GeneA_008 GeneB_008 N0.HOG0000008
#> 11 GeneA_009 GeneB_009 N0.HOG0000009
#> 12 GeneA_010 GeneB_010 N0.HOG0000010
```
