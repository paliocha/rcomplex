# The driver on the shipped two-species fixture.

drv_expr <- function() {
  rd <- function(f) {
    as.matrix(read.delim(system.file("extdata", f, package = "rcomplex"),
      row.names = 1L
    ))
  }
  list(SpA = rd("expr_sp1_small.txt"), SpB = rd("expr_sp2_small.txt"))
}
drv_ortho <- system.file("extdata", "orthologs_small.txt",
  package = "rcomplex"
)

# Three species with three shared modules of 30 genes and 1:1 orthologs.
# The shipped fixture gives each gene one partner, so no pair passes the
# hypergeometric overlap gate there and it has no cliques.
drv_syn <- function() {
  sp <- c("SpA", "SpB", "SpC")
  ld <- rep(list(seq(0.9, 0.6, length.out = 30L)), 3L)
  expr <- withr::with_preserve_seed(lapply(seq_along(sp), function(i) {
    pres_expr(i, 120L, sp[i], ld, 30L, n_samp = 20L) # nolint
  }))
  names(expr) <- sp
  ortho <- data.frame(
    species = rep(sp, each = 120L),
    gene = unlist(lapply(expr, rownames), use.names = FALSE),
    hog = rep(sprintf("H%03d", 1:120), 3L)
  )
  list(expr = expr, ortho = ortho)
}


test_that("rcomplex() runs the pipeline and counts tiers", {
  res <- rcomplex(drv_expr(), drv_ortho, density = 0.1, seed = 1L)
  expect_s3_class(res, "rcomplex")
  expect_named(res, c(
    "networks", "edges", "cliques", "classification", "call"
  ))
  expect_identical(names(res$networks), c("SpA", "SpB"))
  expect_identical(
    names(res$edges)[1:8],
    c(
      "gene1", "gene2", "hog", "score", "evalue", "q_value",
      "effect_size", "power"
    )
  )
  expect_identical(nrow(res$classification), 0L)
  expect_identical(as.data.frame(res), res$classification)

  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, seed = 1L)
  tiers <- table(res$classification$classification)
  expect_identical(c(tiers), c(complete_conserved = 3L))
})


test_that("expr and networks give identical edges", {
  expr <- drv_expr()
  a <- rcomplex(expr, drv_ortho, density = 0.1, seed = 1L)
  nets <- lapply(expr, compute_network, density = 0.1)
  b <- rcomplex(networks = nets, orthologs = drv_ortho, seed = 1L)
  expect_identical(a$edges, b$edges)
  expect_identical(unique(a$edges$sign), "positive")
  neg <- compute_network(expr$SpB, density = 0.1, sign = "negative")
  expect_error(
    rcomplex(networks = list(SpA = nets$SpA, SpB = neg), orthologs = drv_ortho),
    "networks differ in sign: SpA positive, SpB negative"
  )
  expect_identical(a$classification, b$classification)
})


test_that("print() shows species, edges and tiers", {
  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, seed = 1L)
  expect_snapshot(print(res))
})


test_that("null = TRUE reports calls beside null calls", {
  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, null = TRUE, seed = 1L)
  expect_true(is.data.frame(res$edges_null))
  expect_identical(names(res$edges_null), names(res$edges))
  expect_identical(tail(names(res$edges), 2L), c("type", "sign"))
  expect_output(s <- summary(res), "calls_null")
  expect_named(s$null, c(
    "species1", "species2", "calls", "calls_null", "false_call_rate"
  ))
  expect_identical(nrow(s$null), 3L)
  expect_true(all(s$null$calls > 0L))
  expect_identical(s$tiers$cliques, 3L)
  off <- rcomplex(drv_expr(), drv_ortho, density = 0.1, seed = 1L)
  expect_output(s_off <- summary(off), "null: not run")
  expect_named(s_off$tiers, c("tier", "cliques"))
})


test_that("block runs each species on its wiring layer", {
  expr <- drv_expr()
  block <- list(
    SpA = rep(c("t1", "t2", "t3"), length.out = 10L),
    SpB = rep(c("t1", "t2"), each = 5L)
  )
  expect_message(
    res <- rcomplex(expr, drv_ortho,
      block = block, density = 0.1,
      null = TRUE, seed = 1L
    ),
    "SpB lacks block level: t3"
  )
  wiring <- split_layers(expr$SpB, block$SpB)$wiring
  expect_identical(
    res$networks$SpB$threshold,
    compute_network(wiring, density = 0.1)$threshold
  )
})


test_that("modules = TRUE adds modules and preservation", {
  syn <- drv_syn()
  res <- suppressWarnings(suppressMessages(
    rcomplex(syn$expr, syn$ortho, modules = TRUE, seed = 1L)
  ))
  expect_named(res$modules, c("SpA", "SpB", "SpC"))
  expect_true(is.data.frame(res$preservation$classification))
})


test_that("write_rcomplex() writes one TSV per table", {
  res <- rcomplex(drv_expr(), drv_ortho, density = 0.1, seed = 1L)
  dir <- withr::local_tempdir()
  paths <- write_rcomplex(res, dir)
  expect_identical(
    basename(paths),
    c("edges.tsv", "cliques.tsv", "classification.tsv")
  )
  back <- read.delim(paths[[3L]])
  expect_identical(nrow(back), nrow(res$classification))
})


test_that("rcomplex() refuses bad input", {
  expr <- drv_expr()
  nets <- lapply(expr, compute_network, density = 0.1)
  expect_error(rcomplex(orthologs = drv_ortho), "one of expr and networks")
  expect_error(
    rcomplex(expr, drv_ortho, networks = nets),
    "one of expr and networks"
  )
  names(expr) <- c("SpA", "SpX")
  expect_error(rcomplex(expr, drv_ortho), "SpX")
  expect_error(
    rcomplex(networks = nets, orthologs = drv_ortho, null = TRUE),
    "needs expr"
  )
})
