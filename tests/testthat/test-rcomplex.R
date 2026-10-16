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
    "networks", "edges", "cliques", "members", "twins", "capped",
    "classification", "call"
  ))
  expect_identical(nrow(res$capped), 0L)
  expect_identical(names(res$networks), c("SpA", "SpB"))
  expect_identical(
    names(res$edges)[1:8],
    c(
      "gene1", "gene2", "hog", "score", "evalue", "q_value",
      "effect_size", "power"
    )
  )
  # No call and no clades: no graph gives a clique.
  expect_identical(nrow(res$classification), 0L)
  expect_identical(as.data.frame(res), res$classification)

  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, seed = 1L)
  cls <- res$classification
  hogs <- table(cls$hog_class[!duplicated(cls$hog)])
  expect_identical(hogs[["complete_conserved"]], 3L)
  expect_identical(
    sort(unique(cls$run[cls$classification == "complete_conserved"])), 1:2
  )
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
  expect_false(any(res$classification$capped))
  res$capped <- data.frame(hog = c("H001", "H002"))
  expect_output(
    print(res),
    "cliques [^\n]*; 2 hogs over the clique budget, reported as blocks"
  )
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
  expect_true(is.data.frame(res$classification_null))
  expect_identical(
    names(res$classification_null), names(res$classification)
  )
  real <- s$tiers$tier == "complete_conserved"
  # Each complete clique is maximal on the strict and the loose graph.
  expect_identical(s$tiers$cliques[real], 6L)
  expect_identical(s$tiers$multiplicity[real], 6)
  expect_identical(s$tiers$hogs[real], 3L)
  expect_identical(s$tiers$patterns[real], 3L)
  expect_identical(
    s$tiers$cliques_null,
    as.integer(table(factor(
      res$classification_null$classification, levels = s$tiers$tier
    )))
  )
  again <- rcomplex(syn$expr, syn$ortho, null = TRUE, seed = 1L)
  expect_identical(again$edges_null, res$edges_null)
  other <- rcomplex(syn$expr, syn$ortho, null = TRUE, seed = 2L)
  expect_false(identical(other$edges_null, res$edges_null))
  off <- rcomplex(drv_expr(), drv_ortho, density = 0.1, seed = 1L)
  expect_output(s_off <- summary(off), "null: not run")
  expect_named(s_off$tiers, c(
    "tier", "cliques", "cliques_null", "multiplicity", "hogs", "patterns"
  ))
})


test_that("the null shuffle keeps each gene's expression decile", {
  g <- sprintf("g%02d", 1:40)
  net <- function(p) {
    m <- matrix(0, 40L, 40L, dimnames = list(paste0(p, g), paste0(p, g)))
    list(network = m, gene_mean = stats::setNames(1:40, paste0(p, g)))
  }
  nets <- list(A = net("A"), B = net("B"))
  hog <- rep(sprintf("H%02d", 1:20), each = 2L)
  long <- data.frame(
    species = rep(c("A", "B", "B"), c(40L, 40L, 1L)),
    gene = c(paste0("A", g), paste0("B", g), "Bout"),
    hog = c(hog, hog, "H01")
  )
  shuf <- withr::with_seed(1L, rcomplex:::.driver_shuffle(long, nets))
  expect_identical(shuf[, c("species", "hog")], long[, c("species", "hog")])
  expect_identical(shuf$gene[81L], "Bout")
  decile <- function(x) ceiling(as.integer(substring(x, 3L)) / 4)
  expect_identical(decile(shuf$gene[1:80]), decile(long$gene[1:80]))
  expect_false(identical(shuf$gene, long$gene))
  expect_setequal(shuf$gene, long$gene)

  nets$B$gene_mean <- NULL
  expect_message(
    withr::with_seed(1L, rcomplex:::.driver_shuffle(long, nets)),
    "null: B has no gene_mean, so the shuffle ignores expression"
  )
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
  kept <- rownames(res$networks$SpB$network)
  expect_identical(res$networks$SpB$gene_mean, rowMeans(expr$SpB)[kept])
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
    c(
      "edges.tsv", "cliques.tsv", "members.tsv", "twins.tsv",
      "capped.tsv", "classification.tsv"
    )
  )
  back <- read.delim(paths[[6L]])
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
    rcomplex(networks = nets, orthologs = drv_ortho, method = "rank"),
    "need expr"
  )
  expect_error(
    rcomplex(networks = nets, orthologs = drv_ortho, block_size = 8L),
    "block_size"
  )
})


test_that("block_size gives the edges and tiers of the dense build", {
  syn <- drv_syn()
  key <- c("hog", "species1", "species2", "gene1", "gene2")
  srt <- function(e) {
    e <- e[do.call(order, e[key]), , drop = FALSE]
    rownames(e) <- NULL
    e
  }
  two <- lapply(syn$expr, function(x) rep(c("t1", "t2"), ncol(x) / 2L))
  for (block in list(NULL, two)) {
    a <- rcomplex(syn$expr, syn$ortho, block = block, seed = 1L)
    b <- rcomplex(syn$expr, syn$ortho,
      block = block, block_size = 32L, seed = 1L
    )
    ea <- srt(a$edges)
    eb <- srt(b$edges)
    expect_identical(ea[key], eb[key])
    expect_equal(ea$q_value, eb$q_value, tolerance = 1e-8)
    expect_gt(nrow(a$classification), 0L)
    expect_identical(a$classification, b$classification)
  }
})


test_that("a species over 40,000 genes gets a block_size message", {
  note <- rcomplex:::.driver_size_note
  expect_message(
    note(list(A = matrix(0, 40001L, 2L)), NULL),
    "40001.*block_size = 512"
  )
  expect_no_message(note(list(A = matrix(0, 40000L, 2L)), NULL))
  expect_no_message(note(list(A = matrix(0, 40001L, 2L)), 512L))
})


test_that("null = TRUE runs on networks, with or without gene_mean", {
  syn <- drv_syn()
  nets <- lapply(syn$expr, compute_network)
  expect_named(nets$SpA$gene_mean, rownames(nets$SpA$network))
  expect_equal(nets$SpA$gene_mean, rowMeans(syn$expr$SpA))
  a <- rcomplex(networks = nets, orthologs = syn$ortho, null = TRUE, seed = 1L)
  b <- rcomplex(syn$expr, syn$ortho, null = TRUE, seed = 1L)
  expect_identical(a$edges_null, b$edges_null)
  nets$SpC <- as_network(nets$SpC$network)
  expect_message(
    rcomplex(networks = nets, orthologs = syn$ortho, null = TRUE, seed = 1L),
    "null: SpC has no gene_mean"
  )
})


# Two species on one reference genome share gene ids.
test_that("shared gene ids stay in their own species pair", {
  syn <- drv_syn()
  ids <- rownames(syn$expr$SpB)
  expr <- syn$expr
  rownames(expr$SpC) <- ids
  ortho <- syn$ortho
  ortho$gene[ortho$species == "SpC"] <- ids
  nets <- lapply(expr, compute_network)
  nets_d <- nets
  nets_d$SpC <- compute_network(syn$expr$SpC)
  # SpC genes in other hogs than their SpB namesakes
  rev_c <- function(o) {
    o$hog[o$species == "SpC"] <- rev(o$hog[o$species == "SpC"])
    prepare_orthologs(o)
  }
  shared <- find_coexpressologs(nets, rev_c(ortho))
  distinct <- find_coexpressologs(nets_d, rev_c(syn$ortho))
  n <- function(e) c(table(paste(e$species1, e$species2)))
  expect_identical(n(shared), n(distinct))
  ac <- function(e) as.list(e[e$species2 == "SpC", c("gene1", "hog")])
  expect_identical(ac(shared), ac(distinct))

  res <- rcomplex(expr, ortho, null = TRUE, seed = 1L)
  real <- n(res$edges)
  null <- n(res$edges_null)[names(real)]
  expect_true(all(null <= real * 1.05))
  expect_true(all(null >= real * 0.95))
})


# The driver's three gene graphs: strict, loose, and all tested pairs.

# One gene per species in each hog; q holds one value per species pair,
# in combn() order.
tier_edges <- function(hog, sp, q) {
  p <- utils::combn(sp, 2L)
  data.frame(
    gene1 = paste0(p[1L, ], "_", hog), gene2 = paste0(p[2L, ], "_", hog),
    species1 = p[1L, ], species2 = p[2L, ], hog = hog, q_value = q
  )
}

tier_cliques <- function(e, sp, clades = NULL) {
  rcomplex:::.driver_cliques(e, sp, clades, alpha = 0.1, min_power = 0.8)
}


test_that("the loose graph gives partial_significant cliques", {
  sp <- c("A", "B", "C", "D")
  # AB AC AD BC BD CD: D joins A at q < 0.1 and B, C at q = 0.5
  e <- tier_edges("H1", sp, c(0.01, 0.01, 0.01, 0.01, 0.5, 0.5))
  g <- tier_cliques(e, sp)
  cls <- g$classification
  expect_identical(cls$run, 1:2)
  expect_identical(cls$classification, c("unclassified", "partial_significant"))
  expect_identical(cls$n_species, 3:4)
  expect_identical(unique(cls$hog_class), "partial_significant")
  expect_identical(g$cliques$run, 1:2)
  expect_identical(g$cliques$alpha_graph, c(0.1, 0.9))
})


# Two clades, each complete within, every cross pair tested at q = 1
# (H1); and a hog complete over all six species (H2).
tier_two_clades <- function() {
  sp <- c("A", "B", "C", "D", "E", "F")
  p <- utils::combn(sp, 2L)
  within <- (p[1L, ] %in% sp[1:3]) == (p[2L, ] %in% sp[1:3])
  e <- rbind(
    tier_edges("H1", sp, ifelse(within, 0.01, 1)),
    tier_edges("H2", sp, 0.01)
  )
  clades <- list(one = sp[1:3], two = sp[4:6])
  list(sp = sp, g = tier_cliques(e, sp, clades))
}


test_that("the all-tested graph gives differentiated cliques", {
  cls <- tier_two_clades()$g$classification
  h1 <- cls[cls$hog == "H1", ]
  expect_identical(unique(h1$hog_class), "differentiated")
  expect_identical(h1$classification[h1$run == 3L], "differentiated")
  expect_false(any(h1$n_species[h1$run < 3L] == 6L))
})


test_that("a claimed hog gets no clique on the all-tested graph", {
  g <- tier_two_clades()$g
  cls <- g$classification
  expect_identical(
    unique(cls$hog_class[cls$hog == "H2"]), "complete_conserved"
  )
  expect_false(any(cls$run == 3L & cls$hog == "H2"))
  run3 <- g$cliques$run == 3L
  expect_identical(unique(g$cliques$hog[run3]), "H1")
  expect_identical(unique(g$cliques$alpha_graph[run3]), Inf)
  expect_identical(nrow(g$capped), 0L)
})


test_that("run and clique_id join the result tables", {
  syn <- drv_syn()
  res <- rcomplex(syn$expr, syn$ortho, seed = 1L)
  cls <- res$classification
  key <- function(x) paste(x$run, x$clique_id)
  expect_true(all(key(cls) %in% key(res$cliques)))
  expect_true(all(key(cls) %in% key(res$members)))
  expect_identical(
    names(res$cliques)[1:3], c("run", "alpha_graph", "clique_id")
  )
  expect_identical(names(res$members)[1:2], c("run", "alpha_graph"))
  # Without clades there is no all-tested graph.
  expect_identical(sort(unique(res$cliques$alpha_graph)), c(0.1, 0.9))
  expect_false(any(cls$run == 3L))
  res3 <- suppressMessages(rcomplex(syn$expr, syn$ortho,
    clades = list(ab = c("SpA", "SpB")), seed = 1L
  ))
  expect_true(any(res3$cliques$run == 3L))
  expect_true(any(res3$classification$run == 3L))
  expect_output(print(res3), "at q < 0.9, [0-9,]+ on all tested pairs")
  cl <- NULL
  res2 <- rcomplex(syn$expr, syn$ortho, clades = cl, seed = 1L)
  expect_output(print(res2), "graphs: [^\n]*at q < 0.9\n")
  paths <- write_rcomplex(res, withr::local_tempdir())
  back <- read.delim(paths[basename(paths) == "classification.tsv"])
  expect_identical(back$run, cls$run)
  back <- read.delim(paths[basename(paths) == "cliques.tsv"])
  expect_identical(back$run, res$cliques$run)
})


test_that("the null run classifies three graphs", {
  res <- rcomplex(drv_expr(), drv_ortho,
    density = 0.1, null = TRUE, seed = 1L
  )
  expect_identical(names(res$classification_null)[1L], "run")
})


test_that("rcomplex() forwards r_min and prints the floor", {
  res <- rcomplex(drv_expr(), drv_ortho, density = 0.1, r_min = 0.5,
    seed = 1L
  )
  expect_identical(res$networks$SpA$params$r_min, 0.5)
  expect_gte(res$networks$SpA$params$r_threshold, 0.5)
  expect_gt(res$networks$SpA$params$n_floored, 0L)
  expect_output(print(res), "floor r >= 0.50: [0-9]+ edges removed")
  nets <- lapply(drv_expr(), compute_network, density = 0.1)
  expect_error(
    rcomplex(networks = nets, orthologs = drv_ortho, r_min = 0.3),
    "r_min needs expr; you gave networks"
  )
})
