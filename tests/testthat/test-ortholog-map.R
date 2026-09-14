# Tests for resolve_ortholog_map()

# Fixture: HOG "H1" is multi-copy (3 species-1 genes x 2 species-2 genes),
# "H2" and "H3" are single-copy.  Candidate species-2 genes are
# A1, A2, B1, C1 - the set every resolution path must preserve.
map_fixture <- function() {
  list(
    genes1 = c("a1", "a2", "a3", "b1", "c1"),
    genes2 = c("A1", "A2", "A3", "B1", "C1"),
    ortho = data.frame(
      Species1 = c("a1", "a1", "a2", "a2", "a3", "a3", "b1", "c1"),
      Species2 = c("A1", "A2", "A1", "A2", "A1", "A2", "B1", "C1"),
      hog = c("H1", "H1", "H1", "H1", "H1", "H1", "H2", "H3"),
      stringsAsFactors = FALSE
    ),
    cliques = data.frame(
      hog = "H1", SP_A = "a1", SP_B = "A2",
      n_species = 2L, mean_q = 0.01,
      stringsAsFactors = FALSE
    ),
    edges = data.frame(
      gene1 = c("a2", "a3", "a2"),
      gene2 = c("A1", "A1", "A2"),
      species1 = "SP_A", species2 = "SP_B",
      hog = "H1",
      q.value = c(0.001, 0.01, 0.02),
      effect_size = c(0.9, 0.5, 0.2),
      jaccard = c(0.8, 0.4, 0.1),
      type = "conserved",
      stringsAsFactors = FALSE
    )
  )
}

candidate_sp2 <- function(fx) {
  cand <- fx$ortho[fx$ortho$Species1 %in% fx$genes1 &
                     fx$ortho$Species2 %in% fx$genes2, ]
  sort(unique(cand$Species2))
}


# ---- structure ----

test_that("resolve_ortholog_map returns the documented structure", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("gene1", "gene2", "hog", "source"))
  expect_true(all(res$source %in%
                    c("clique", "coexpressolog", "unresolved")))
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


# ---- clique layer ----

test_that("a clique resolves the copy pair it names", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    cliques = fx$cliques
  )

  clique_rows <- res[res$source == "clique", ]
  expect_equal(nrow(clique_rows), 1L)
  expect_equal(clique_rows$gene1, "a1")
  expect_equal(clique_rows$gene2, "A2")
})

test_that("the best clique per HOG wins on n_species then mean_q", {
  fx <- map_fixture()
  cl <- data.frame(
    hog = c("H1", "H1", "H1"),
    SP_A = c("a1", "a2", "a3"),
    SP_B = c("A2", "A1", "A1"),
    n_species = c(2L, 3L, 3L),
    mean_q = c(0.001, 0.05, 0.01),
    stringsAsFactors = FALSE
  )
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B", cliques = cl
  )

  # n_species = 3 beats the lower mean_q at n_species = 2; among the two
  # three-species cliques the lower mean_q (a3) wins.
  clique_rows <- res[res$source == "clique", ]
  expect_equal(clique_rows$gene1, "a3")
  expect_equal(clique_rows$gene2, "A1")
})

test_that("clique pairs absent from the ortholog table are ignored", {
  fx <- map_fixture()
  cl <- data.frame(
    hog = "H1", SP_A = "a1", SP_B = "A3",
    n_species = 2L, mean_q = 0.01,
    stringsAsFactors = FALSE
  )
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B", cliques = cl
  )

  expect_false(any(res$source == "clique"))
})

test_that("clique columns must exist for both species", {
  fx <- map_fixture()
  expect_error(
    resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "NOPE", cliques = fx$cliques
    ),
    "no column for species"
  )
})


# ---- coexpressolog layer ----

test_that("mutual-best coexpressologs resolve remaining copies", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
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
    sp1 = "SP_A", sp2 = "SP_B",
    edges = fx$edges
  )

  expect_false(any(res$source == "coexpressolog" & res$gene1 == "a3"))
})

test_that("cliques take precedence over coexpressologs", {
  fx <- map_fixture()
  # Give a1 a strong coexpressolog to A1; the clique already put a1 on A2.
  edges <- rbind(fx$edges, data.frame(
    gene1 = "a1", gene2 = "A1", species1 = "SP_A", species2 = "SP_B",
    hog = "H1", q.value = 1e-6, effect_size = 0.99, jaccard = 0.99,
    type = "conserved", stringsAsFactors = FALSE
  ))
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = edges, cliques = fx$cliques
  )

  a1 <- res[res$gene1 == "a1", ]
  expect_equal(nrow(a1), 1L)
  expect_equal(a1$gene2, "A2")
  expect_equal(a1$source, "clique")
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
    sp1 = "SP_A", sp2 = "SP_B",
    edges = fx$edges
  )
  rev <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = rev_edges
  )

  expect_equal(fwd, rev)
})

test_that("non-significant edges do not resolve copies", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$type <- "ns"
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B", edges = edges
  )

  expect_false(any(res$source == "coexpressolog"))
})

test_that("alpha filters edges when the table has no type column", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$type <- NULL
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = edges, alpha = 1e-6
  )

  expect_false(any(res$source == "coexpressolog"))
})


# ---- q.value guard ----

test_that("rank_by = q.value errors on HOG-level (constant) q-values", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$q.value <- 0.001 # what method = "permutation" produces

  expect_error(
    resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "SP_B",
      edges = edges, rank_by = "q.value"
    ),
    "constant within every multi-copy HOG"
  )
})

test_that("rank_by = q.value works when q-values vary within a HOG", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = fx$edges, rank_by = "q.value"
  )

  co <- res[res$source == "coexpressolog", ]
  expect_equal(co$gene1, "a2")
  expect_equal(co$gene2, "A1")
})


# ---- the preserved-gene-set invariant ----

test_that("resolution never changes the set of mappable species-2 genes", {
  fx <- map_fixture()
  expected <- candidate_sp2(fx)

  variants <- list(
    none = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2),
    clique = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "SP_B",
      cliques = fx$cliques
    ),
    coexpr = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "SP_B",
      edges = fx$edges
    ),
    both = resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "SP_B",
      edges = fx$edges, cliques = fx$cliques
    )
  )

  for (nm in names(variants)) {
    expect_setequal(sort(unique(variants[[nm]]$gene2)), expected)
  }
})

test_that("each resolved species-2 gene has exactly one species-1 partner", {
  fx <- map_fixture()
  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = fx$edges, cliques = fx$cliques
  )

  resolved <- res[res$source != "unresolved", ]
  expect_false(anyDuplicated(resolved$gene2) > 0L)
})

test_that("resolution reduces the number of candidate pairs", {
  fx <- map_fixture()
  none <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2)
  both <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B",
    edges = fx$edges, cliques = fx$cliques
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
      cliques = fx$cliques
    ),
    "sp1 and sp2 are required"
  )
  expect_error(
    resolve_ortholog_map(fx$ortho, "zzz", "ZZZ"),
    "No orthologs found"
  )
})


# ---- cross-layer collisions ----

test_that("a coexpressolog cannot re-claim a gene the clique layer resolved", {
  fx <- map_fixture()
  # The clique puts a1 on A2. These edges make a3-A2 the mutual best, which
  # would give A2 a second resolved partner; .pres_project() then majority-
  # votes over the two labels and drops the gene on a tie, removing it from
  # the mappable set - exactly what resolution must never do.
  edges <- data.frame(
    gene1 = c("a3", "a3", "a2"),
    gene2 = c("A2", "A1", "A2"),
    species1 = "SP_A", species2 = "SP_B", hog = "H1",
    q.value = 0.001,
    effect_size = c(0.99, 0.10, 0.50),
    jaccard = c(0.9, 0.1, 0.5),
    type = "conserved",
    stringsAsFactors = FALSE
  )

  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B", edges = edges, cliques = fx$cliques
  )

  resolved <- res[res$source != "unresolved", ]
  expect_false(anyDuplicated(resolved$gene2) > 0L)
  expect_equal(resolved$gene1[resolved$gene2 == "A2"], "a1")
  expect_setequal(sort(unique(res$gene2)), candidate_sp2(fx))
})


# ---- defensive paths ----

test_that("cliques without n_species or mean_q still resolve", {
  fx <- map_fixture()
  cl <- data.frame(
    hog = c("H1", "H1"),
    SP_A = c("a1", "a2"),
    SP_B = c("A2", "A1"),
    stringsAsFactors = FALSE
  )

  res <- resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
    sp1 = "SP_A", sp2 = "SP_B", cliques = cl
  )
  expect_true(any(res$source == "clique"))
})

test_that("edges missing hog are reported by name", {
  fx <- map_fixture()
  edges <- fx$edges
  edges$hog <- NULL

  expect_error(
    resolve_ortholog_map(fx$ortho, fx$genes1, fx$genes2,
      sp1 = "SP_A", sp2 = "SP_B", edges = edges, rank_by = "q.value"
    ),
    "edges missing columns: hog"
  )
})
