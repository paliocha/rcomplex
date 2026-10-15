# WP15c dev check: joint_modules() reproduces the probe's match step.
# Wood half A, seed 1, lambda 1. Reads the probe E2 cache: halves.rds
# (per-half networks and the hog table), the two-step labels and the
# match labels of joint_louvain(init, freeze = TRUE).
#   R_LIBS=<lib> Rscript dev/probes/joint_match_check.R [cache_dir]
# Pass: ARI >= 0.95, or quality >= probe quality - 1e-8. Wall time of
# the lining-up step (modules given) under 5 s.
library(rcomplex)
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) {
  args[1]
} else {
  "/Users/martin/Documents/R/rcomplex/prepare_data/cache/probe_e2"
}

ari <- function(x, y) {
  tab <- table(x, y)
  pairs <- function(n) sum(n * (n - 1) / 2)
  a <- pairs(tab)
  b <- pairs(rowSums(tab))
  d <- pairs(colSums(tab))
  e <- b * d / pairs(sum(tab))
  (a - e) / ((b + d) / 2 - e)
}

halves <- readRDS(file.path(dir, "halves.rds"))
nets <- halves$nets$A
two <- readRDS(file.path(dir, "wood_runs_two_step_0_1_A.rds"))
ref <- readRDS(file.path(dir, "wood_runs_match_1_1_A.rds"))

# Probe labels are named "species<TAB>gene", species in network order
key <- function(s) paste(s, rownames(nets[[s]]$network), sep = "\t")
mods <- lapply(stats::setNames(names(nets), names(nets)), function(s) {
  as_modules(stats::setNames(two$lab[key(s)], rownames(nets[[s]]$network)))
})

t0 <- Sys.time()
res <- joint_modules(nets, halves$ortho, modules = mods, lambda = 1, seed = 1)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

lab <- unlist(lapply(names(nets), function(s) res$modules[[s]]$modules),
  use.names = FALSE
)
ref_lab <- ref$lab[unlist(lapply(names(nets), key))]
a <- ari(lab, ref_lab)
cat(sprintf("species modules: %d\n", nrow(res$members)))
cat(sprintf(
  "shared labels: ours %d, probe %d\n",
  length(unique(lab)), length(unique(ref_lab))
))
cat(sprintf("ARI: %.6f\n", a))
cat(sprintf(
  "quality: ours %.10f, probe %.10f\n",
  res$params$quality, ref$quality
))
cat(sprintf("wall time (modules given): %.2f s\n", secs))
pass <- (a >= 0.95 || res$params$quality >= ref$quality - 1e-8) && secs < 5
cat(if (pass) "PASS\n" else "FAIL\n")
