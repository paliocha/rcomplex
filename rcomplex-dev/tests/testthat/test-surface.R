test_that("exported surface is stable", {
  ns <- asNamespace("rcomplex")
  exports <- sort(grep("^\\.__", getNamespaceExports("rcomplex"),
    invert = TRUE, value = TRUE
  ))
  n_formals <- vapply(exports, function(f) {
    sum(setdiff(names(formals(get(f, envir = ns))), "...") != "")
  }, integer(1))
  expect_snapshot({
    exports
    n_formals
  })
})
