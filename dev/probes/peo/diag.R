# Where does the dense network build crash on wheat root (71,330 genes)?
# Usage: Rscript diag.R <n> <mode>; jobs.sh diag runs one call per n and
# mode, so a segfault at one n does not hide the others.
# mode cora:  Rfast::cora(t(x)) alone.
# mode steps: the dense steps of compute_network() one at a time:
#   cor_rfast, mutual_rank_inplace_cpp, density_threshold_cpp (density and
#   store_density), extract_sparse_cpp.
# mode full:  compute_network(x, density = 0.03, n_cores = 8).
suppressPackageStartupMessages(library(rcomplex))
args <- commandArgs(TRUE)
n <- as.integer(args[1])
mode <- args[2]
ROOT <- "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
for (p in c("Rcpp", "RcppArmadillo", "Rfast", "rcomplex")) {
  cat(p, as.character(packageVersion(p)), "\n")
}
x <- readRDS(file.path(
  ROOT, "complex_LeafvRoot/inputs/root",
  "Triticum_aestivum.rds"
))
x <- x[seq_len(min(n, nrow(x))), , drop = FALSE]
cat(sprintf(
  "n = %d genes x %d samples; n^2 = %.4g (2^31 = %.4g)\n",
  nrow(x), ncol(x), as.numeric(nrow(x))^2, 2^31
))
mem <- function() {
  g <- gc()
  sprintf("%.0f MB max", sum(g[, "max used"] * c(56, 8)) / 2^20)
}
run <- function(name, expr) {
  cat("before", name, "\n")
  t <- system.time(v <- force(expr))[["elapsed"]]
  cat(sprintf("after %s: %.0f s, %s\n", name, t, mem()))
  v
}
if (mode == "cora") {
  r <- run("cora", Rfast::cora(t(x)))
} else if (mode == "steps") {
  net <- run("cor_rfast", rcomplex:::cor_rfast(x))
  run(
    "mutual_rank_inplace_cpp",
    rcomplex:::mutual_rank_inplace_cpp(net, FALSE, FALSE, N_CORES)
  )
  thr <- run(
    "density_threshold_cpp",
    rcomplex:::density_threshold_cpp(net, 0.03)
  )
  sthr <- run(
    "density_threshold_cpp(store)",
    rcomplex:::density_threshold_cpp(net, 0.05)
  )
  s <- run(
    "extract_sparse_cpp",
    rcomplex:::extract_sparse_cpp(net, sthr, N_CORES)
  )
  cat("thr", thr, "store_thr", sthr, "nnz", length(s$x), "\n")
} else {
  net <- run(
    "compute_network",
    compute_network(x, density = 0.03, n_cores = N_CORES)
  )
}
cat("PASS", mode, n, "\n")
