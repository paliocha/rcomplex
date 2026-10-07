# Tests for read_orthologs() and prepare_orthologs() on the long table.

ext <- function(f) system.file("extdata", f, package = "rcomplex")

test_that("N0 and long fixtures give the same table for the shared species", {
  n0 <- read_orthologs(ext("N0.tsv"), species = c("SpA", "SpB"))
  long <- read_orthologs(ext("orthologs_small.txt"))
  expect_identical(n0, long)
  expect_identical(
    read_orthologs(ext("N0.tsv"), c("SpA", "SpB"), format = "orthofinder"),
    n0
  )
  expect_identical(
    read_orthologs(ext("orthologs_small.txt"), format = "long"),
    long
  )
})

test_that("OrthoFinder input keeps every species and skips empty cells", {
  all <- read_orthologs(ext("N0.tsv"))
  expect_named(all, c("species", "gene", "hog"))
  expect_setequal(unique(all$species), c("SpA", "SpB", "SpC"))
  sp_c <- all[all$species == "SpC", ]
  expect_identical(nrow(sp_c), 10L)
  expect_identical(sp_c$hog[sp_c$gene == "GeneC_013"], "N0.HOG0000012")
  expect_false(any(all$gene == ""))
})

test_that("Orthogroups.tsv is read through its Orthogroup column", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c(
    "Orthogroup\tSpA\tSpB",
    "OG1\ta1, a2\tb1",
    "OG2\t\tb2"
  ), f)
  expect_identical(
    read_orthologs(f),
    data.frame(
      species = c("SpA", "SpA", "SpB", "SpB"),
      gene = c("a1", "a2", "b1", "b2"),
      hog = c("OG1", "OG1", "OG1", "OG2")
    )
  )
})

test_that("PLAZA input gives anchor and member genes per group", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c(
    "species\tgene_id\tgene_content",
    "spa\ta1\tspb:b1,b2;spc:c1",
    "spa\ta2\tspb:b1,b2;spc:c1",
    "spa\ta3\tspb:b3"
  ), f)
  out <- read_orthologs(f, species = c("spa", "spb"))
  expect_identical(out, data.frame(
    species = c("spa", "spa", "spa", "spb", "spb", "spb"),
    gene = c("a1", "a2", "a3", "b1", "b2", "b3"),
    hog = c(1L, 1L, 2L, 1L, 1L, 2L)
  ))
  pairs <- prepare_orthologs(out)
  expect_identical(nrow(pairs), 5L)
  expect_setequal(
    paste(pairs$gene1, pairs$gene2),
    c("a1 b1", "a1 b2", "a2 b1", "a2 b2", "a3 b3")
  )
})

test_that("read_orthologs refuses unknown species, formats and files", {
  expect_error(read_orthologs(ext("N0.tsv"), "SpX"), "SpX")
  expect_error(read_orthologs("no_such_file.tsv"), "not found")
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("a\tb", "1\t2"), f)
  expect_error(read_orthologs(f), "format")
  expect_error(read_orthologs(f, format = "long"), "species, gene, hog")
  expect_error(read_orthologs(f, format = "plaza"), "gene_content")
  expect_error(read_orthologs(f, format = "orthofinder"), "HOG")
})

test_that("prepare_orthologs pairs genes within each HOG", {
  long <- read_orthologs(ext("N0.tsv"))
  pairs <- prepare_orthologs(long)
  expect_named(pairs, c("gene1", "gene2", "hog"))
  is_ab <- startsWith(pairs$gene1, "GeneA") & startsWith(pairs$gene2, "GeneB")
  ab <- pairs[is_ab, ]
  # 10 HOGs, one extra paralog in each of the first two
  expect_identical(nrow(ab), 12L)
  expect_setequal(
    ab$gene2[ab$hog == "N0.HOG0000001"],
    c("GeneB_001", "GeneB_011")
  )
  # A-C and B-C pairs: 8 each (HOG1 and HOG2 carry paralogs)
  expect_identical(sum(startsWith(pairs$gene2, "GeneC")), 16L)
})

test_that("prepare_orthologs maps genes through reductions", {
  long <- data.frame(
    species = c("SP_A", "SP_A", "SP_A", "SP_B", "SP_B"),
    gene = c("A1", "A2", "A3", "B1", "B2"),
    hog = c("HOG1", "HOG1", "HOG2", "HOG1", "HOG2")
  )
  red <- list(
    SP_A = list(gene_map = data.frame(
      original = c("A1", "A2"), representative = c("A1", "A1")
    )),
    SP_B = list(gene_map = data.frame(original = "B1", representative = "B1"))
  )
  pairs <- prepare_orthologs(long, red)
  # A1 and A2 merge to A1; A3 is missing from the map and keeps its name
  expect_identical(nrow(pairs), 2L)
  expect_identical(pairs$gene1[pairs$hog == "HOG1"], "A1")
  expect_identical(pairs$gene1[pairs$hog == "HOG2"], "A3")
  # Without reductions, every paralog pair survives
  expect_identical(nrow(prepare_orthologs(long)), 3L)
})

test_that("prepare_orthologs validates inputs", {
  long <- data.frame(
    species = c("A", "B"), gene = c("a", "b"),
    hog = c("H", "H")
  )
  expect_error(prepare_orthologs(long[, 1:2]), "species, gene, hog")
  expect_error(prepare_orthologs(long[1, ]), "at least two species")
  expect_error(prepare_orthologs(long, list(1, 2)), "named list")
  expect_error(
    prepare_orthologs(long, list(A = list(gene_map = data.frame()))),
    "missing species"
  )
  expect_error(
    prepare_orthologs(long, list(A = list(), B = list())), "\\$gene_map"
  )
})

test_that("prepare_orthologs ignores genes without a HOG", {
  long <- data.frame(
    species = c("A", "B", "A", "B"),
    gene = c("a1", "b1", "a2", "b2"),
    hog = c("H", "H", NA, NA)
  )
  expect_identical(nrow(prepare_orthologs(long)), 1L)
})
