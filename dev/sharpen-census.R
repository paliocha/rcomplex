ex <- grep("^export\\(", readLines("NAMESPACE"), value = TRUE)
ex <- sub("^export\\((.*)\\)$", "\\1", ex)
rd <- function(p) paste(unlist(lapply(p, readLines)), collapse = "\n")
rsrc <- rd(list.files("R", full.names = TRUE))
readme <- rd("README.md")
tut <- rd("vignettes/rcomplex-tutorial.Rmd")
meth <- rd("vignettes/articles/methods.Rmd")
tests <- rd(list.files("tests/testthat", pattern = "[.]R$", full.names = TRUE))
cnt <- function(f, txt) {
  m <- gregexpr(paste0("(?<![A-Za-z_.])", f, "\\("), txt, perl = TRUE)[[1]]
  sum(m > 0)
}
d <- data.frame(
  fn = ex,
  internal = sapply(ex, function(f) {
    defs <- sum(gregexpr(paste0("\n", f, "(\\.[a-z]+)? <- function"), rsrc,
      perl = TRUE
    )[[1]] > 0)
    cnt(f, rsrc) - defs
  }),
  readme = sapply(ex, cnt, readme),
  tutorial = sapply(ex, cnt, tut),
  methods = sapply(ex, cnt, meth),
  tests = sapply(ex, cnt, tests)
)
d <- d[order(d$readme + d$tutorial, d$internal), ]
print(d, row.names = FALSE)
