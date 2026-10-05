# Tests for as_modules(): any gene partition into the preservation tests.
# Fixture: pres_fixture() / true_modules() in helper-preservation.R.

# nolint start: object_usage_linter. (fixture functions from helper files)
as_mod_fixture <- function() {
  fx <- pres_fixture()
  hand_a <- true_modules(fx$netA, fx$mods)
  hand_b <- true_modules(fx$netB, fx$mods)
  c(fx, list(
    handA = hand_a, handB = hand_b,
    labA = hand_a$modules, labB = hand_b$modules
  ))
}
# nolint end

test_that("vector and list inputs give the same assignment", {
  v <- c(g1 = "a", g2 = "a", g3 = "b", g4 = NA)
  l <- list(a = c("g1", "g2"), b = "g3")
  expect_identical(as_modules(v), as_modules(l))
  m <- as_modules(v)
  expect_identical(m$modules, c(g1 = "a", g2 = "a", g3 = "b"))
  expect_identical(m$module_genes, list(a = c("g1", "g2"), b = "g3"))
  expect_identical(m$n_modules, 2L)
})

test_that("as_modules() is idempotent and passes module objects through", {
  m <- as_modules(list(a = c("g1", "g2"), b = "g3"))
  expect_identical(as_modules(m), m)
  fake <- list(modules = c(x = 1L), module_genes = list(`1` = "x"))
  expect_identical(as_modules(fake), fake)
})

test_that("min_size unassigns small modules", {
  m <- as_modules(list(a = c("g1", "g2"), b = "g3"), min_size = 2)
  expect_identical(names(m$module_genes), "a")
  expect_false("g3" %in% names(m$modules))
})

test_that("bad input fails with the offending value", {
  expect_error(
    as_modules(list(a = c("g1", "g2"), b = c("g2", "g3"))),
    "g2 \\(a, b\\)"
  )
  expect_error(as_modules(c(g1 = "a", g1 = "b")), "assigned twice: g1")
  expect_error(as_modules(c(g1 = "a", g2 = "")), "empty module label")
  expect_error(as_modules(list(a = 1:3)), "character vector")
  expect_error(as_modules(list(c("g1"))), "named by module label")
  expect_error(as_modules(c("a", "b")), "named by gene")
  expect_error(as_modules(c(g1 = "a"), min_size = 0), "min_size")
  expect_error(as_modules(list(a = c("g1", NA))), "NA or empty")
  expect_error(as_modules(list(a = c("g1", ""))), "NA or empty")
  expect_error(as_modules(list(a = "g1", a = "g2")), "label used twice: a")
  m <- as_modules(list(a = c("g1", "g2")))
  expect_error(as_modules(m, min_size = 2), "not to an existing module")
  expect_error(
    as_modules(data.frame(gene = c("g1", "g2"), module = c("a", "b"))),
    "not a data frame"
  )
})

test_that("factor labels and empty input are handled", {
  v <- c(g1 = "b", g2 = "a", g3 = "b")
  f <- factor(v)
  names(f) <- names(v)
  expect_identical(as_modules(f), as_modules(v))
  empty <- list(
    as_modules(c(g1 = NA, g2 = NA)),
    as_modules(list(a = character(0)))
  )
  for (z in empty) {
    expect_identical(z$n_modules, 0L)
    expect_length(z$module_genes, 0L)
    expect_length(z$modules, 0L)
  }
})

test_that("module order follows the labels, numerically for numbers", {
  # numeric vector labels sort numerically; a list keeps its own order
  expect_identical(names(as_modules(c(g1 = 2, g2 = 1))$module_genes),
                   c("1", "2"))
  expect_identical(
    names(as_modules(list("2" = "g1", "1" = "g2"))$module_genes),
    c("2", "1")
  )
  v <- stats::setNames(c(1:12, 10L), sprintf("g%02d", 1:13))
  expect_identical(names(as_modules(v)$module_genes), as.character(1:12))
  expect_identical(
    names(as_modules(list(b = "x", a = "y"))$module_genes), c("b", "a")
  )
})

test_that("preservation results are the same through as_modules()", {
  fx <- as_mod_fixture()
  a <- module_preservation(fx$handA, fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, seed = 11
  )
  b <- module_preservation(as_modules(fx$labA), fx$netA, fx$netB, fx$ortho,
    n_perm = 50L, seed = 11
  )
  expect_identical(b$preservation, a$preservation)

  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  expect_identical(
    module_correspondence(as_modules(fx$labA), as_modules(fx$labB), map)$pairs,
    module_correspondence(fx$handA, fx$handB, map)$pairs
  )

  nets <- list(A = fx$netA, B = fx$netB)
  pairs <- data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE)
  p1 <- preservation_paired(list(A = fx$handA, B = fx$handB), nets, fx$ortho,
    pairs,
    n_perm = 50L, seed = 1
  )
  am <- list(A = as_modules(fx$labA), B = as_modules(fx$labB))
  p2 <- preservation_paired(am,
    nets, fx$ortho, pairs,
    n_perm = 50L, seed = 1
  )
  expect_identical(p2$classification, p1$classification)
})

test_that("twelve integer-labelled modules give the whole result unchanged", {
  fx <- as_mod_fixture()
  g <- rownames(fx$netA$network)[1:156]
  lab <- stats::setNames(rep(1:12, each = 13), g)
  hand <- list(modules = lab, module_genes = split(names(lab), lab),
               n_modules = 12L)
  a <- module_preservation(hand, fx$netA, fx$netB, fx$ortho,
                           n_perm = 30L, seed = 3)
  b <- module_preservation(as_modules(lab), fx$netA, fx$netB, fx$ortho,
                           n_perm = 30L, seed = 3)
  expect_identical(b, a)
})

test_that("named gene sets in any order give the integer-labelled result", {
  fx <- as_mod_fixture()
  lab <- fx$labA
  names_of <- c("ribosome", "auxin", "zinc", "calvin") # not alphabetical
  sets <- stats::setNames(lapply(1:4, function(k) names(lab)[lab == k]),
                          names_of)
  a <- module_preservation(fx$handA, fx$netA, fx$netB, fx$ortho,
                           n_perm = 50L, seed = 11)
  b <- module_preservation(as_modules(sets), fx$netA, fx$netB, fx$ortho,
                           n_perm = 50L, seed = 11)
  expect_setequal(b$coverage$module, names_of)
  pa <- a$preservation
  pb <- b$preservation
  pb <- pb[match(names_of[as.integer(pa$module)], pb$module), ]
  # observed statistics match; the permutation columns do not, because
  # module_preservation() draws its null per module in label order, so
  # renaming the modules reorders the draws
  obs <- c("size", "size_mapped", "avg.weight", "cor.degree")
  expect_equal(unname(as.list(pb[obs])), unname(as.list(pa[obs])))
})

test_that("doubles that print alike share one label", {
  m <- as_modules(c(g1 = 0.1 + 0.2, g2 = 0.3, g3 = 1))
  expect_identical(m$n_modules, 2L)
  expect_identical(names(m$module_genes), c("0.3", "1"))
})

test_that("label and gene-ID order do not depend on the collation locale", {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old))
  fx <- as_mod_fixture()
  lab <- fx$labA
  sets <- stats::setNames(
    lapply(1:4, function(k) names(lab)[lab == k]),
    c("Zinc", "auxin", "Beta", "calvin") # case-mixed: C and en_US differ
  )
  run <- function() {
    module_preservation(as_modules(sets), fx$netA, fx$netB, fx$ortho,
                        n_perm = 30L, seed = 5)
  }
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  lab_b <- fx$labB
  sets_b <- stats::setNames(
    lapply(1:4, function(k) names(lab_b)[lab_b == k]),
    c("yew", "Ash", "oak", "Birch")
  )
  corr <- function() {
    module_correspondence(as_modules(sets), as_modules(sets_b), map,
                          seed = 1)$pairs
  }
  # the deterministic pick among tied reference genes, mixed-case IDs
  tie <- data.frame(
    gene1 = c("b1", "B2"), gene2 = "t1", module = "m",
    source = "unresolved", stringsAsFactors = FALSE
  )
  # C-locale (radix) outcomes, checked without any other locale
  Sys.setlocale("LC_COLLATE", "C")
  r_c <- run()
  c_c <- corr()
  c_order <- c("Beta", "Zinc", "auxin", "calvin")
  expect_identical(r_c$preservation$module, c_order)
  expect_identical(unique(c_c$module_sp1), c_order)
  expect_identical(rcomplex:::.pres_project(tie)$gene1, "B2")
  # and the same under en_US, which collates these differently
  en <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if(!nzchar(en), "en_US.UTF-8 collation not available")
  skip_if(
    identical(sort(c("b1", "B2")), c("B2", "b1")),
    "en_US collates like C here, so the comparison would prove nothing"
  )
  expect_identical(run(), r_c)
  expect_identical(corr(), c_c)
  expect_identical(rcomplex:::.pres_project(tie)$gene1, "B2")
})

test_that("a detect_modules() result and its membership give one result", {
  fx <- as_mod_fixture()
  m <- detect_modules(fx$netA, seed = 1L)
  expect_s3_class(m$modules, "membership")
  expect_identical(as_modules(m), m)
  a <- module_preservation(m, fx$netA, fx$netB, fx$ortho,
                           n_perm = 30L, seed = 2)
  b <- module_preservation(as_modules(m$modules), fx$netA, fx$netB, fx$ortho,
                           n_perm = 30L, seed = 2)
  expect_identical(b, a)
})

test_that("named gene sets reach tag_permutation() via preservation_paired()", {
  fx <- as_mod_fixture()
  # an unstructured partner makes A's modules diverge, so the HOG pool
  # tag_permutation() reads is non-empty (as in the hand-built test)
  set.seed(99)
  e <- matrix(stats::rnorm(500 * 40), 500, 40)
  rownames(e) <- paste0("B", sprintf("%04d", seq_len(500)))
  net_b <- compute_network(e, density = 0.03, sparse = FALSE)
  hand <- list(A = fx$handA, B = true_modules(net_b, fx$mods))
  named_of <- function(m, nm) {
    lab <- m$modules
    as_modules(stats::setNames(
      lapply(seq_along(nm), function(k) names(lab)[lab == k]), nm
    ))
  }
  named <- list(
    A = named_of(hand$A, c("ribosome", "auxin", "zinc", "calvin")),
    B = named_of(hand$B, c("yew", "Ash", "oak", "Birch"))
  )
  nets <- list(A = fx$netA, B = net_b)
  grp <- c(A = "annual", B = "perennial")
  pairs <- data.frame(sp1 = "A", sp2 = "B", pair_name = "AB",
                      stringsAsFactors = FALSE)
  pool <- function(mods) {
    res <- suppressWarnings(preservation_paired(mods, nets, fx$ortho, pairs,
      group = grp, n_perm = 50L, min_module_size = 3L, seed = 1
    ))
    tp <- suppressMessages(suppressWarnings(tag_permutation(
      res$classification, mods, fx$ortho, pairs,
      group = grp, target_group = "annual",
      n_perm = 50L, min_recurrence = 1L
    )))
    list(observed = tp$observed, hogs = sort(tp$recurrence_table$hog))
  }
  p_hand <- pool(hand)
  p_named <- pool(named)
  expect_gt(p_hand$observed, 0L)
  expect_identical(p_named, p_hand)
})

test_that("a map with factor gene columns gives the character result", {
  fx <- as_mod_fixture()
  map <- resolve_ortholog_map(
    fx$ortho, rownames(fx$netA$network), rownames(fx$netB$network)
  )
  map_f <- map
  map_f$gene1 <- factor(map$gene1, levels = rev(unique(map$gene1)))
  map_f$gene2 <- factor(map$gene2, levels = rev(unique(map$gene2)))
  pres <- function(m) {
    module_preservation(fx$handA, fx$netA, fx$netB, map = m,
                        n_perm = 30L, seed = 4)
  }
  expect_identical(pres(map_f), pres(map))
  expect_identical(
    module_correspondence(fx$handA, fx$handB, map_f, seed = 1),
    module_correspondence(fx$handA, fx$handB, map, seed = 1)
  )
})
