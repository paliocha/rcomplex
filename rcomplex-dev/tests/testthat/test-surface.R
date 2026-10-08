# Wide-signature exceptions, by name (WP5/WP6 parent decisions)
wide_ok <- c(
  "rcomplex", "module_preservation", "find_coexpressologs", "density_sweep"
)

surface <- function() {
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
  list(exports = exports, n_formals = n_formals)
}

test_that("exported surface is stable", {
  s <- surface()
  exports <- s$exports
  n_formals <- s$n_formals
  expect_snapshot({
    exports
    n_formals
  })
})

test_that("exports stay within budget", {
  s <- surface()
  expect_lte(length(s$exports), 31)
  over <- s$n_formals[s$n_formals > 8]
  expect_setequal(setdiff(names(over), wide_ok), character())
})

test_that("source files and README stay within budget", {
  # Skipped when the source tree is absent (installed-package check:
  # R/ is not installed and README is Rbuildignored)
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "source tree not available")
  n <- vapply(
    list.files(rdir, "[.]R$", full.names = TRUE),
    function(f) length(readLines(f)),
    integer(1)
  )
  expect_true(all(n <= 1500), info = names(n)[n > 1500])
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not available")
  expect_lte(length(readLines(readme)), 150)
})
