# joint_modules(): per-species modules lined up by hog overlap.

# Three species, two planted modules of 30 genes each. Hogs H1-H20 and
# H31-H50 are one-to-one inside their module. The 10 defector hogs D21-D30
# hold genes 21-30 of A and B (module 1) and genes 51-60 of C (module 2).
jm_fixture <- function() {
  sp <- c("A", "B", "C")
  ld <- rep(list(rep(0.8, 30L)), 2L)
  nets <- lapply(seq_along(sp), function(k) {
    compute_network(pres_expr(40 + k, 60, sp[k], ld, 30L), # nolint
      density = 0.25, sparse = TRUE, store_density = 0.25
    )
  })
  names(nets) <- sp
  one <- c(1:20, 31:50)
  i <- c(one, 21:30)
  gene <- sprintf("%s%04d", rep(sp, each = 50L), c(i, i, one, 51:60))
  hog <- rep(c(paste0("H", one), paste0("D", 21:30)), 3L)
  ortho <- data.frame(species = rep(sp, each = 50L), gene = gene, hog = hog)
  planted <- lapply(nets, function(n) {
    stats::setNames(rep(1:2, each = 30L), rownames(n$network))
  })
  list(nets = nets, ortho = ortho, planted = planted)
}

jm_fx <- local({
  f <- jm_fixture()
  f$res <- joint_modules(f$nets, f$ortho, modules = f$planted, seed = 1)
  f
})

test_that("planted modules line up into two shared labels", {
  res <- jm_fx$res
  expect_named(res, c("modules", "members", "hogs", "params"))
  expect_named(res$modules, c("A", "B", "C"))
  m <- res$members
  expect_identical(sort(unique(m$label)), 1:2)
  for (l in 1:2) {
    expect_setequal(m$species[m$label == l], c("A", "B", "C"))
    expect_identical(sum(m$label == l), 3L)
  }
  expect_identical(m$n_genes, rep(30L, 6L))
  lab_a <- res$modules$A$modules
  expect_identical(lab_a[["A0001"]], m$label[m$species == "A"][1L])
  expect_identical(res$modules$A$method, "joint")
  expect_identical(as_modules(res$modules$A), res$modules$A)

  h <- res$hogs
  expect_false(is.unsorted(h$hog))
  expect_identical(nrow(h), 50L)
  expect_true(all(h$n_modules[startsWith(h$hog, "D")] == 2L))
  expect_true(all(h$n_modules[startsWith(h$hog, "H")] == 1L))
  expect_identical(h$copies_main[h$hog == "D21"], 2L)

  p <- res$params
  expect_identical(p$lambda, 1)
  expect_identical(p$gamma, 1)
  expect_true(p$modules_given)
  expect_true(is.finite(p$quality) && p$quality > 0)
})

test_that("lambda = 0 lines up nothing across species", {
  res <- joint_modules(jm_fx$nets, jm_fx$ortho,
    modules = jm_fx$planted, lambda = 0, seed = 1
  )
  n_sp <- tapply(res$members$species, res$members$label, function(x) {
    length(unique(x))
  })
  expect_true(all(n_sp == 1L))
  expect_lt(res$params$quality, jm_fx$res$params$quality)
})

test_that("modules = NULL runs detect_modules() per species", {
  f <- jm_fx
  res <- joint_modules(f$nets, f$ortho, seed = 7)
  given <- lapply(seq_along(f$nets), function(k) {
    detect_modules(f$nets[[k]], seed = rcomplex:::.task_seed(7, 1L, k))
  })
  names(given) <- names(f$nets)
  res2 <- joint_modules(f$nets, f$ortho, modules = given, seed = 7)
  expect_false(res$params$modules_given)
  expect_identical(res$members, res2$members)
  expect_identical(
    lapply(res$modules, `[[`, "modules"),
    lapply(res2$modules, `[[`, "modules")
  )
})

test_that("the same seed gives the same result", {
  res <- joint_modules(jm_fx$nets, jm_fx$ortho,
    modules = jm_fx$planted, seed = 1
  )
  expect_identical(res, jm_fx$res)
})

test_that("modules may be vectors, in any species order", {
  f <- jm_fx
  mods <- rev(lapply(f$planted, function(m) {
    stats::setNames(c("x", "y")[m], names(m))
  }))
  res <- joint_modules(f$nets, f$ortho, modules = mods, seed = 1)
  expect_identical(res$members$label, f$res$members$label)
  expect_type(res$members$module, "character")
})

test_that("bad input names the argument", {
  f <- jm_fx
  n <- f$nets
  pl <- f$planted
  expect_error(joint_modules(unname(n), f$ortho, pl), "`nets`")
  expect_error(joint_modules(n["A"], f$ortho, pl["A"]), "`nets`")
  expect_error(joint_modules(n, f$ortho, pl[c("A", "B")]), "`modules`")
  pl$C <- pl$C[-1L]
  expect_error(joint_modules(n, f$ortho, pl), "modules\\$C.*C0001")
  pl <- f$planted
  expect_error(joint_modules(n, f$ortho, pl, lambda = -1), "`lambda`")
  expect_error(joint_modules(n, f$ortho, pl, lambda = "a"), "`lambda`")
  expect_error(joint_modules(n, f$ortho, pl, lambda = 1:2), "`lambda`")
})
