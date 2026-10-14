# joint_modules(): star expansion over the species networks and the hogs.

ari <- function(x, y) {
  tab <- table(x, y)
  pairs <- function(n) sum(n * (n - 1) / 2)
  a <- pairs(tab)
  b <- pairs(rowSums(tab))
  d <- pairs(colSums(tab))
  e <- b * d / pairs(sum(tab))
  (a - e) / ((b + d) / 2 - e)
}

long_ortho <- function(o) {
  data.frame(
    species = rep(c("A", "B"), each = nrow(o)),
    gene = c(o$gene1, o$gene2), hog = rep(o$hog, 2L),
    stringsAsFactors = FALSE
  )
}

# ARI against the planted modules, on the planted genes only
planted_ari <- function(res, s, net, mods) {
  truth <- true_modules(net, mods)$modules # nolint
  ari(truth, res$modules[[s]]$modules[names(truth)])
}

# pres_fixture() shape with loading 0.7 on every planted gene. Its own
# heavy-tailed loadings (median near 0.14) leave most planted genes
# without a within-module edge at density 0.03, so no method reaches the
# planted modules there; detect_modules() alone gets an ARI near 0.
jm_fixture <- function() {
  per <- 40L
  ld <- rep(list(rep(0.7, per)), 4L)
  eB <- pres_expr(32, 500, "B", ld, per) # nolint
  list(
    netA = compute_network(pres_expr(31, 300, "A", ld, per), # nolint
      density = 0.03, sparse = FALSE
    ),
    netB = compute_network(eB, density = 0.03, sparse = FALSE),
    netB_sparse = compute_network(eB,
      density = 0.03, sparse = TRUE, store_density = 0.03
    ),
    ortho = data.frame(
      gene1 = sprintf("A%04d", 1:300), gene2 = sprintf("B%04d", 1:300),
      hog = paste0("H", 1:300), stringsAsFactors = FALSE
    ),
    mods = lapply(1:4, function(k) ((k - 1) * per + 1):(k * per))
  )
}

# Hog weight 4 couples the species. The default (median scaled weight, near
# 1 for mutual-rank stores) lies below the coupling transition between 1
# and 2: on this fixture it recovers each species but shares no label.
jm_fx <- local({
  f <- jm_fixture()
  long <- long_ortho(f$ortho)
  list(
    f = f, long = long,
    res = joint_modules(list(A = f$netA, B = f$netB), long,
      weight = 4, seed = 1
    )
  )
})

test_that("the default weight recovers each species on its own", {
  f <- jm_fx$f
  res <- joint_modules(list(A = f$netA, B = f$netB), jm_fx$long, seed = 1)
  expect_gt(res$params$weight, 0)
  expect_lte(res$params$weight, 1)
  expect_gt(planted_ari(res, "A", f$netA, f$mods), 0.9)
  expect_gt(planted_ari(res, "B", f$netB, f$mods), 0.9)
})

test_that("joint modules recover the planted modules in both species", {
  f <- jm_fx$f
  res <- jm_fx$res
  expect_gt(planted_ari(res, "A", f$netA, f$mods), 0.9)
  expect_gt(planted_ari(res, "B", f$netB, f$mods), 0.9)
  majority <- function(s, net, k) {
    g <- rownames(net$network)[f$mods[[k]]]
    lab <- res$modules[[s]]$modules[g]
    as.integer(names(which.max(table(lab))))
  }
  for (k in seq_along(f$mods)) {
    expect_identical(majority("A", f$netA, k), majority("B", f$netB, k))
  }
  expect_named(res, c("modules", "hogs", "params", "modularity"))
  expect_named(res$modules, c("A", "B"))
  expect_identical(res$params$weight, 4)
})

test_that("a non-conserved copy falls outside module_main", {
  f <- jm_fx$f
  long <- jm_fx$long
  # B0041 (planted module 2) moves from H41 into H1 (module 1)
  long$hog[long$species == "B" & long$gene == "B0041"] <- "H1"
  res <- joint_modules(list(A = f$netA, B = f$netB), long,
    weight = 4, seed = 1
  )
  h1 <- res$hogs[res$hogs$hog == "H1", ]
  expect_identical(h1$n_copies, 3L)
  expect_identical(h1$n_modules, 2L)
  expect_identical(h1$copies_main, 2L)
  expect_false(res$modules$B$modules[["B0041"]] == h1$module_main)
  rest <- res$hogs[res$hogs$hog %in% paste0("H", 2:160), ]
  expect_gte(mean(rest$n_modules == 1L), 0.9)
  expect_false(is.unsorted(res$hogs$hog))
})

test_that("a sparse store gives the same modules", {
  f <- jm_fx$f
  res <- joint_modules(list(A = f$netA, B = f$netB_sparse), jm_fx$long,
    weight = 4, seed = 1
  )
  expect_gt(planted_ari(res, "B", f$netB, f$mods), 0.9)
})

test_that("per-species modules feed module_preservation()", {
  f <- jm_fx$f
  mods <- jm_fx$res$modules$A
  expect_identical(as_modules(mods), mods)
  expect_identical(mods$method, "joint")
  pres <- suppressMessages(module_preservation(mods, f$netA, f$netB,
    f$ortho,
    n_perm = 20L, seed = 1
  ))
  expect_s3_class(pres$preservation, "data.frame")
})

test_that("bad input names the argument", {
  f <- jm_fx$f
  long <- jm_fx$long
  expect_error(joint_modules(list(f$netA, f$netB), long), "nets")
  expect_error(joint_modules(list(A = f$netA), long), "nets")
  expect_error(
    joint_modules(list(A = f$netA, C = f$netB), long), "orthologs"
  )
  expect_error(
    joint_modules(list(A = f$netA, B = f$netB), long, engine = "x"),
    "engine"
  )
  expect_error(
    joint_modules(list(A = f$netA, B = f$netB), long, weight = -1),
    "weight"
  )
})
