# Tests for edge_history(): Fitch gains and losses on a species tree

eh_table <- function() {
  data.frame(
    run = 1:4, clique_id = as.character(1:4), hog = "H1",
    missing_species = c("C,D", "D", "D", "C,D"),
    missing_reason = c(
      "tested_ns,tested_ns", "tested_ns", "untested", "absent,underpowered"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("edge_history places gains and losses on a four-tip tree", {
  skip_if_not_installed("ape")
  phy <- ape::read.tree(text = "((A,B),(C,D));")
  res <- edge_history(eh_table(), phy, c("A", "B", "C", "D"))
  cl <- res$cliques
  expect_identical(cl$run, 1:4)
  expect_identical(cl$n_gain, c(1L, 0L, 0L, 0L))
  expect_identical(cl$n_loss, c(0L, 1L, 0L, 0L))
  expect_identical(cl$gains, c("node6", "", "", ""))
  expect_identical(cl$losses, c("", "D", "", ""))
  expect_identical(cl$state_root, c(0L, 1L, 1L, 1L))
  expect_identical(cl$root_tied, c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(cl$n_mpr, c(2, 1, 1, 1))
  expect_identical(cl$n_present, c(2L, 3L, 3L, 2L))
  expect_identical(cl$n_absent, c(2L, 1L, 0L, 0L))
  expect_identical(cl$n_unknown, c(0L, 0L, 1L, 2L))

  br <- res$branches
  expect_identical(br$branch, c("A", "B", "C", "D", "node6", "node7"))
  expect_identical(br$n_species, c(1L, 1L, 1L, 1L, 2L, 2L))
  expect_identical(br$n_gain, c(0L, 0L, 0L, 0L, 1L, 0L))
  expect_identical(br$n_loss, c(0L, 0L, 0L, 1L, 0L, 0L))
})

test_that("edge_history treats a tip outside species as unknown", {
  skip_if_not_installed("ape")
  phy <- ape::read.tree(text = "((A,B),(C,D));")
  tab <- eh_table()[2, -1]
  res <- edge_history(tab, phy, c("A", "B", "D"))
  expect_false("run" %in% names(res$cliques))
  expect_identical(res$cliques$n_unknown, 1L)
  # C unknown: one gain on (A,B) ties one loss on D at the root
  expect_identical(res$cliques$gains, "node6")
  expect_true(res$cliques$root_tied)
  expect_identical(res$cliques$n_mpr, 3)
})

test_that("edge_history checks its inputs", {
  skip_if_not_installed("ape")
  phy <- ape::read.tree(text = "((A,B),(C,D));")
  sp <- c("A", "B", "C", "D")
  expect_error(edge_history(eh_table(), phy, c(sp, "E")), "E")
  expect_error(
    edge_history(eh_table(), ape::unroot(phy), sp), "ape::root",
    fixed = TRUE
  )
  expect_error(
    edge_history(eh_table()[, -5], phy, sp), "missing_reason"
  )
  none <- edge_history(eh_table()[0, ], phy, sp)
  expect_identical(nrow(none$cliques), 0L)
  expect_identical(nrow(none$branches), 6L)
  expect_identical(sum(none$branches$n_gain + none$branches$n_loss), 0L)
})

test_that("edge_history matches phangorn::fitch() on random trees", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  sp <- paste0("t", 1:8)
  lab <- c("0" = "tested_ns", "?" = "untested")
  withr::with_seed(7L, {
    for (r in 1:20) {
      tree <- ape::rtree(8, tip.label = sp)
      mat <- replicate(50, {
        x <- sample(c("1", "0", "?"), 8, replace = TRUE)
        x[sample(8, 2)] <- "1"
        x
      })
      rownames(mat) <- sp
      tab <- data.frame(
        clique_id = as.character(1:50), hog = "H",
        missing_species = apply(mat, 2, function(x) {
          paste(sp[x != "1"], collapse = ",")
        }),
        missing_reason = apply(mat, 2, function(x) {
          paste(lab[x[x != "1"]], collapse = ",")
        })
      )
      out <- edge_history(tab, tree, sp)$cliques
      dat <- phangorn::phyDat(
        mat,
        type = "USER", levels = c("0", "1"), ambiguity = "?"
      )
      ref <- phangorn::fitch(tree, dat, site = "site")[attr(dat, "index")]
      expect_equal(out$n_gain + out$n_loss, as.vector(ref))
      expect_true(all(out$n_mpr >= 1))
      expect_true(all(!out$root_tied[out$n_mpr == 1]))
    }
  })
})
