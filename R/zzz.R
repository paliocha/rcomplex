# The package's only .onLoad(); a second one elsewhere would silently
# replace it (collation order), so add to this one instead.

# Environment captured when the package loads (see .blas_fork_safe() in
# R/rng.R, which reads the VECLIB_MAXIMUM_THREADS snapshot).
.load_env <- new.env(parent = emptyenv())
.load_env$veclib_threads <- ""

.onLoad <- function(libname, pkgname) {
  .load_env$veclib_threads <- Sys.getenv("VECLIB_MAXIMUM_THREADS")
}
