# Do detect_modules() modules replicate across independent samples?
# Pooideae leaf, per species: samples split into two halves (two replicates of
# each time point each), 2,000 genes of highest variance over all samples,
# compute_network() and detect_modules() on each half, adjusted Rand index
# between the two partitions. Floors: the same with every gene's samples
# shuffled (independent noise), and the seed ARI (same half, two seeds) as the
# ceiling of what the algorithm itself reproduces. Package defaults (CPM) and
# the modularity objective. Writes dev/bench/modules_split_half.tsv.
# Run from the package root: Rscript dev/bench/modules_split_half.R [n_cores]
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(SummarizedExperiment)
})
n_cores <- as.integer(commandArgs(TRUE)[1] %||% 4L)
RES <- c(0.5, 1, 2)
ari <- function(a, b) {
  g <- intersect(names(a), names(b))
  igraph::compare(as.integer(factor(a[g])), as.integer(factor(b[g])), method = "adjusted.rand")
}
mods <- function(x, obj, seed) {
  net <- compute_network(x, n_cores = n_cores)
  m <- detect_modules(net,
    resolution = RES, objective_function = obj, seed = seed,
    n_cores = n_cores, test_k1 = FALSE
  )
  m$modules
}
rows <- list()
for (f in list.files("prepare_data/data", "_se\\.rds$", full.names = TRUE)) {
  sp <- sub("_se\\.rds$", "", basename(f))
  se <- readRDS(f)
  se <- se[, colData(se)$tissue == "leaf"]
  x <- as.matrix(assay(se))
  x <- x[apply(x, 1L, var) > 0, , drop = FALSE]
  x <- x[order(-apply(x, 1L, var))[1:2000], , drop = FALSE]
  tp <- as.character(colData(se)$time_point)
  rk <- stats::ave(seq_along(tp), tp, FUN = seq_along)
  h1 <- which(rk %% 2L == 1L)
  h2 <- which(rk %% 2L == 0L)
  set.seed(1L)
  xs <- t(apply(x, 1L, sample))
  for (obj in c("CPM", "modularity")) {
    a <- mods(x[, h1], obj, 1L)
    b <- mods(x[, h2], obj, 1L)
    a2 <- mods(x[, h1], obj, 2L)
    sa <- mods(xs[, h1], obj, 1L)
    sb <- mods(xs[, h2], obj, 1L)
    r <- data.frame(
      species = sp, objective = obj,
      n_mod_half1 = length(unique(a)), n_mod_half2 = length(unique(b)),
      split_half_ari = round(ari(a, b), 3), seed_ari = round(ari(a, a2), 3),
      shuffled_ari = round(ari(sa, sb), 3),
      n_mod_shuffled = length(unique(sa))
    )
    print(r, row.names = FALSE)
    rows[[length(rows) + 1L]] <- r
  }
}
res <- do.call(rbind, rows)
utils::write.table(res, "dev/bench/modules_split_half.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(stats::aggregate(cbind(split_half_ari, seed_ari, shuffled_ari) ~ objective, res, median))
