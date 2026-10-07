test_that("exported surface is stable", {
  ns <- asNamespace("rcomplex")
  exports <- sort(grep("^\\.__", getNamespaceExports("rcomplex"),
    invert = TRUE, value = TRUE
  ))
  # An S3 generic only forwards: count the formals of its .default method
  n_formals <- vapply(exports, function(f) {
    fn <- get(f, envir = ns)
    dflt <- paste0(f, ".default")
    if (any(grepl("UseMethod", deparse(body(fn)))) &&
          exists(dflt, envir = ns, inherits = FALSE)) {
      fn <- get(dflt, envir = ns)
    }
    sum(setdiff(names(formals(fn)), "...") != "")
  }, integer(1))
  expect_snapshot({
    exports
    n_formals
  })
})
