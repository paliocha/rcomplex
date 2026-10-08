# Tests for resolve_ortholog_map()

# Fixture: HOG "H1" is multi-copy (3 species-1 genes x 2 species-2 genes),
# "H2" and "H3" are single-copy.  Candidate species-2 genes are
# A1, A2, B1, C1 - the set every resolution path must preserve.
map_fixture <- function() {
  list(
    genes1 = c("a1", "a2", "a3", "b1", "c1"),
    genes2 = c("A1", "A2", "A3", "B1", "C1"),
    ortho = data.frame(
      gene1 = c("a1", "a1", "a2", "a2", "a3", "a3", "b1", "c1"),
      gene2 = c("A1", "A2", "A1", "A2", "A1", "A2", "B1", "C1"),
      hog = c("H1", "H1", "H1", "H1", "H1", "H1", "H2", "H3"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(
      gene1 = c("a2", "a3", "a2"),
      gene2 = c("A1", "A1", "A2"),
      species1 = "SP_A", species2 = "SP_B",
      hog = "H1",
      q_value = c(0.001, 0.01, 0.02),
      effect_size = c(0.9, 0.5, 0.2),
      jaccard = c(0.8, 0.4, 0.1),
      type = "conserved",
      stringsAsFactors = FALSE
    )
  )
}

candidate_sp2 <- function(fx) {
  cand <- fx$ortho[fx$ortho$gene1 %in% fx$genes1 &
                     fx$ortho$gene2 %in% fx$genes2, ]
  sort(unique(cand$gene2))
}


# ---- structure ----

test_that("resolve_ortholog_map returns the documented structure", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("gene1", "gene2", "hog", "source"))
  expect_true(all(res$source %in%
                    c("coexpressolog", "unresolved")))
})

test_that("without evidence every candidate pair is unresolved", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2)

  expect_true(all(res$source == "unresolved"))
  expect_equal(nrow(res), nrow(fx$ortho))
})

test_that("genes outside the universes are dropped", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, c("a1", "b1"), c("A1", "B1"))

  expect_setequal(res$gene1, c("a1", "b1"))
  expect_setequal(res$gene2, c("A1", "B1"))
})


# ---- coexpressolog layer ----

test_that("mutual-best coexpressologs resolve remaining copies", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = fx$edges
  )

  co <- res[res$source == "coexpressolog", ]
  expect_equal(nrow(co), 1L)
  expect_equal(co$gene1, "a2")
  expect_equal(co$gene2, "A1")
})

test_that("one-sided best pairs are not accepted", {
  fx <- map_fixture()
  # a3's best is A1, but A1's best is a2 - so a3-A1 is not mutual.
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = fx$edges
  )

  expect_false(any(res$source == "coexpressolog" & res$gene1 == "a3"))
})

test_that("reversed edge orientation is handled", {
  fx <- map_fixture()
  rev_edges <- fx$edges
  names(rev_edges)[names(rev_edges) == "gene1"] <- "tmp"
  names(rev_edges)[names(rev_edges) == "gene2"] <- "gene1"
  names(rev_edges)[names(rev_edges) == "tmp"] <- "gene2"
  rev_edges$species1 <- "SP_B"
  rev_edges$species2 <- "SP_A"

  fwd <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = fx$edges
  )
  rev <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = rev_edges
  )

  expect_equal(fwd, rev)
})

test_that("non-significant edges do not resolve copies", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$type <- "ns"
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B", edges = edges
  )

  expect_false(any(res$source == "coexpressolog"))
})

# ---- the preserved-gene-set invariant ----

test_that("resolution never changes the set of mappable species-2 genes", {
  fx <- map_fixture()
  expected <- candidate_sp2(fx)

  variants <- list(
    none = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2),
    coexpr = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      species1 = "SP_A", species2 = "SP_B",
      edges = fx$edges
    )
  )

  for (nm in names(variants)) {
    expect_setequal(sort(unique(variants[[nm]]$gene2)), expected)
  }
})

test_that("each resolved species-2 gene has exactly one species-1 partner", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = fx$edges
  )

  resolved <- res[res$source != "unresolved", ]
  expect_false(anyDuplicated(resolved$gene2) > 0L)
})

test_that("resolution reduces the number of candidate pairs", {
  fx <- map_fixture()
  none <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2)
  both <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    species1 = "SP_A", species2 = "SP_B",
    edges = fx$edges
  )

  expect_lt(nrow(both), nrow(none))
})


# ---- validation ----

test_that("resolve_ortholog_map validates its inputs", {
  fx <- map_fixture()

  expect_error(
    resolve_ortholog_map("nope", fx$genes1, fx$genes2),
    "must be a data.frame"
  )
  expect_error(
    resolve_ortholog_map(data.frame(a = 1), fx$genes1, fx$genes2),
    "must have columns"
  )
  expect_error(
    resolve_ortholog_map(fx$ortho, 1:3, fx$genes2),
    "must be character vectors"
  )
  expect_error(
    resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      edges = fx$edges
    ),
    "species1 and species2 are required"
  )
  expect_error(
    resolve_ortholog_map(fx$ortho, "zzz", "ZZZ"),
    "no gene pair with both genes"
  )
})


# ---- defensive paths ----

test_that("edges missing hog are reported by name", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$hog <- NULL

  expect_error(
    resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      species1 = "SP_A", species2 = "SP_B", edges = edges
    ),
    "edges missing columns: hog"
  )
})

test_that("tied copies resolve the same way under any collation", {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old))
  # Two species-1 copies per HOG whose names C and en_US order differently
  # (B2 before b1 in C, after it in en_US), tied in the coexpressolog layer.
  ortho <- data.frame(
    gene1 = c("b1", "B2", "c3", "C4"),
    gene2 = c("t1", "t1", "t2", "t2"),
    hog = c("H1", "H1", "H2", "H2"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    gene1 = c("c3", "C4"), gene2 = c("t2", "t2"),
    species1 = "SP_A", species2 = "SP_B", hog = "H2",
    q_value = 0.01, effect_size = 0.7, jaccard = 0.5, type = "conserved",
    stringsAsFactors = FALSE
  )
  run <- function() {
    resolve_ortholog_map(ortho, c("b1", "B2", "c3", "C4"), c("t1", "t2"),
      species1 = "SP_A", species2 = "SP_B", edges = edges
    )
  }
  # the C-locale (radix) order decides, uppercase first; needs no other locale
  Sys.setlocale("LC_COLLATE", "C")
  r_c <- run()
  expect_identical(r_c$gene1[r_c$source == "coexpressolog"], "C4")
  # factor ID columns whose levels are in another collation's order
  # (b1 before B2) resolve as the character columns do
  ortho_f <- ortho
  ortho_f$gene1 <- factor(ortho$gene1, levels = c("b1", "B2", "c3", "C4"))
  edges_f <- edges
  edges_f$gene1 <- factor(edges$gene1, levels = c("c3", "C4"))
  r_f <- resolve_ortholog_map(ortho_f, c("b1", "B2", "c3", "C4"), c("t1", "t2"),
    species1 = "SP_A", species2 = "SP_B", edges = edges_f
  )
  expect_identical(r_f, r_c)
  en <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if(!nzchar(en), "en_US.UTF-8 collation not available")
  skip_if(
    identical(sort(c("b1", "B2")), c("B2", "b1")),
    "en_US collates like C here, so the comparison would prove nothing"
  )
  expect_identical(run(), r_c)
})
