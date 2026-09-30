# Tests for as_modules(): any gene partition into the preservation tests.
# pres_fixture() / true_modules() live in test-module-preservation.R; the
# fixture is rebuilt here from the same generator so this file runs alone.

as_mod_fixture <- function() {
  set.seed(77)
  per <- 40L
  loadings <- lapply(1:4, function(k) {
    l <- stats::rlnorm(per, 0, 0.9)
    l / max(l)
  })
  expr <- function(seed, n, prefix) {
    set.seed(seed)
    e <- matrix(stats::rnorm(n * 40), n, 40)
    for (k in 1:4) {
      f <- stats::rnorm(40)
      rows <- ((k - 1) * per + 1):(k * per)
      for (j in seq_along(rows)) {
        lam <- loadings[[k]][j]
        sd_j <- sqrt(max(1e-6, 1 - lam^2))
        e[rows[j], ] <- lam * f + stats::rnorm(40, sd = sd_j)
      }
    }
    rownames(e) <- paste0(prefix, sprintf("%04d", seq_len(n)))
    e
  }
  netA <- compute_network(expr(31, 300, "A"), density = 0.03, sparse = FALSE) # nolint
  netB <- compute_network(expr(32, 500, "B"), density = 0.03, sparse = FALSE) # nolint
  label <- function(net) {
    g <- rownames(net$network)[1:160]
    stats::setNames(rep(1:4, each = per), g)
  }
  hand <- function(net) {
    m <- label(net)
    list(modules = m, module_genes = split(names(m), m), n_modules = 4L)
  }
  list(
    netA = netA, netB = netB, labA = label(netA), labB = label(netB),
    handA = hand(netA), handB = hand(netB),
    ortho = data.frame(
      Species1 = paste0("A", sprintf("%04d", 1:300)),
      Species2 = paste0("B", sprintf("%04d", 1:300)),
      hog = paste0("H", 1:300), stringsAsFactors = FALSE
    )
  )
}

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
