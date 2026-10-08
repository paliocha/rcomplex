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
    "networks", "edges", "cliques", "members", "classification", "call"
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


test_that("r_threshold is the weakest correlation that passes", {
  set.seed(5)
  x <- matrix(rnorm(40 * 12), 40, 12, dimnames = list(sprintf("g%02d", 1:40)))
  cm <- cor(t(x))
  weakest <- function(net, val) {
    m <- as.matrix(net$network)
    pass <- m >= net$threshold & upper.tri(m)
    val[pass]
  }
  dense <- compute_network(x, density = 0.1, sparse = FALSE)
  expect_equal(dense$params$r_threshold, min(weakest(dense, cm)))
  sp <- compute_network(x, density = 0.1)
  expect_identical(sp$params$r_threshold, dense$params$r_threshold)
  blk <- compute_network(x, density = 0.1, block_size = 8L)
  expect_identical(blk$params$r_threshold, sp$params$r_threshold)

  neg <- compute_network(x, density = 0.1, sign = "negative")
  expect_equal(neg$params$r_threshold, max(weakest(neg, cm)))
  expect_lt(neg$params$r_threshold, 0)
  expect_identical(
    compute_network(x,
      density = 0.1, sign = "negative",
      block_size = 8L
    )$params$r_threshold,
    neg$params$r_threshold
  )

  f <- rep(c("a", "b"), each = 6L)
  rect <- (pmax(cor(t(x[, 1:6])), 0) + pmax(cor(t(x[, 7:12])), 0)) / 2
  pt <- compute_network(x, density = 0.1, partition = f)
  expect_equal(pt$params$r_threshold, min(weakest(pt, rect)))
  rect_neg <- (pmax(-cor(t(x[, 1:6])), 0) + pmax(-cor(t(x[, 7:12])), 0)) / 2
  pt_neg <- compute_network(x, density = 0.1, partition = f, sign = "negative")
  expect_equal(pt_neg$params$r_threshold, -min(weakest(pt_neg, rect_neg)))

  expect_null(as_network(sp$network)$params$r_threshold)
})


test_that("print() shows species, edges and tiers", {
  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, seed = 1L)
  expect_snapshot(print(res))
  nets <- lapply(syn$expr, function(x) as_network(compute_network(x)$network))
  expect_snapshot(print(rcomplex(networks = nets, orthologs = syn$ortho)))
})


# One generator at 6, 20 and 200 samples: same genes, same orthologs.
# Modules hold the first 75 genes; their orthologs are the planted pairs.
test_that("the driver runs and stays calibrated at every sample size", {
  sp <- c("SpA", "SpB")
  ld <- rep(list(seq(0.9, 0.6, length.out = 25L)), 3L)
  ortho <- data.frame(
    species = rep(sp, each = 200L),
    gene = c(sprintf("SpA%04d", 1:200), sprintf("SpB%04d", 1:200)),
    hog = rep(sprintf("H%03d", 1:200), 2L)
  )
  planted <- sprintf("H%03d", 1:75)
  run <- function(n) {
    expr <- withr::with_preserve_seed(lapply(1:2, function(i) {
      pres_expr(i, 200L, sp[i], ld, 25L, n_samp = n) # nolint
    }))
    names(expr) <- sp
    res <- rcomplex(expr, ortho, density = 0.05, null = TRUE, seed = 1L)
    expect_output(s <- summary(res), "calls_null")
    hit <- res$edges$q_value[res$edges$hog %in% planted] < 0.1
    list(res = res, s = s, power = mean(hit %in% TRUE))
  }
  runs <- lapply(c(6L, 20L, 200L), run)
  ref <- runs[[1L]]
  for (r in runs) {
    expect_identical(names(r$res), names(ref$res))
    expect_identical(names(r$res$edges), names(ref$res$edges))
    expect_identical(
      names(r$res$classification), names(ref$res$classification)
    )
    expect_identical(names(r$s$null), names(ref$s$null))
    expect_lte(r$s$null$false_call_rate, 0.1)
  }
  power <- vapply(runs, `[[`, 1, "power")
  expect_true(all(diff(power) >= 0))
  expect_gt(power[[1L]], 0)
  r_thr <- vapply(runs, function(r) {
    r$res$networks$SpA$params$r_threshold
  }, 1)
  expect_true(all(diff(r_thr) < 0))
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
  # the fixture makes no calls, so the false-call rate is undefined
  expect_output(s <- summary(res), "calls_null")
  expect_identical(s$null$calls, 0L)
  expect_identical(s$null$false_call_rate, NA_real_)
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
  dir <- file.path(withr::local_tempdir(), "new")
  paths <- write_rcomplex(res, dir)
  expect_true(dir.exists(dir))
  expect_identical(
    basename(paths),
    c("edges.tsv", "cliques.tsv", "members.tsv", "classification.tsv")
  )
  back <- read.delim(paths[[4L]])
  expect_identical(nrow(back), nrow(res$classification))
})


test_that("rcomplex() refuses bad input", {
  expr <- drv_expr()
  nets <- lapply(expr, compute_network, density = 0.1)
  expect_error(rcomplex(orthologs = drv_ortho), "give exactly one")
  expect_error(
    rcomplex(expr, drv_ortho, networks = nets),
    "give exactly one"
  )
  names(expr) <- c("SpA", "SpX")
  expect_error(rcomplex(expr, drv_ortho), "SpX")
  expect_error(
    rcomplex(networks = nets, orthologs = drv_ortho, null = TRUE),
    "needs expr"
  )
})
