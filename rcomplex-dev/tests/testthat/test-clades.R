# Tests for clades: .check_clades() and clades_from_tree()

test_that(".check_clades refuses malformed and crossing clades", {
  sp <- c("A", "B", "C", "D")
  expect_error(rcomplex:::.check_clades(c(A = "x"), sp), "named list")
  expect_error(
    rcomplex:::.check_clades(list(c("A", "B")), sp),
    "named list"
  )
  expect_error(
    rcomplex:::.check_clades(list(x = c("A", NA)), sp),
    "clade x"
  )
  expect_error(
    rcomplex:::.check_clades(
      list(x = c("A", "B"), y = c("B", "C")), sp
    ),
    "x and y overlap"
  )
  expect_error(
    rcomplex:::.check_clades(list(D = c("A", "B")), sp),
    "species outside every clade: D"
  )
})


test_that(".check_clades restricts clades to the species analysed", {
  cl <- rcomplex:::.check_clades(
    list(x = c("A", "B", "Z"), y = c("A", "B"), z = "Z"),
    c("A", "B", "C")
  )
  # z holds no analysed species; y duplicates x once Z is gone.
  expect_equal(cl, list(x = c("A", "B")))
})


test_that("clade helpers find the smallest clade and the top groups", {
  cl <- list(
    outer = c("A", "B", "C", "D"), mid = c("A", "B", "C"),
    inner = c("A", "B"), other = c("E", "F")
  )
  expect_equal(rcomplex:::.clade_home(cl, c("A", "B")), "inner")
  expect_equal(rcomplex:::.clade_home(cl, c("A", "C")), "mid")
  expect_equal(rcomplex:::.clade_home(cl, c("A", "D")), "outer")
  expect_true(is.na(rcomplex:::.clade_home(cl, c("A", "E"))))
  expect_equal(
    rcomplex:::.clade_groups(cl, c("A", "D", "E", "G")),
    c(A = "outer", D = "outer", E = "other", G = "G")
  )
})


test_that("clades_from_tree round-trips a four-tip tree", {
  skip_if_not_installed("ape")
  phy <- ape::read.tree(text = "((A,B)AB,(C,D)CD)root;")
  expect_equal(
    clades_from_tree(phy),
    list(AB = c("A", "B"), CD = c("C", "D"))
  )
  phy$node.label <- NULL
  expect_equal(
    clades_from_tree(phy),
    list(node6 = c("A", "B"), node7 = c("C", "D"))
  )
  one <- clades_from_tree(phy, min_size = 1L)
  expect_equal(one[c("A", "D")], list(A = "A", D = "D"))
  expect_length(one, 6L)
  # The list passes the clade check unchanged.
  expect_equal(
    rcomplex:::.check_clades(clades_from_tree(phy), phy$tip.label),
    clades_from_tree(phy)
  )
})


test_that("clades_from_tree asks for ape when it is missing", {
  local_mocked_bindings(
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
  expect_error(clades_from_tree(list()), "Install the ape package")
})
