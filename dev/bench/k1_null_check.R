# Does detect_modules()'s K = 1 test reject "no community structure" on null
# data? Its null rewires the graph keeping degrees, which removes the geometric
# clustering that a correlation graph from few samples has even on noise.
#
# Null data: (a) seeded rnorm expression; (b) Pooideae leaf per species with
# every gene's samples shuffled independently (skipped if prepare_data/data/
# is absent). Reference: (c) the same Pooideae data unshuffled.
# 2,000 top-variance genes, 20 samples, package defaults, the documented
# multi-resolution call. Writes dev/bench/k1_null_check.tsv.
#
# Run from the package root: Rscript dev/bench/k1_null_check.R [n_cores]
# On macOS with R's Accelerate BLAS set VECLIB_MAXIMUM_THREADS=1: forked K = 1
# workers otherwise segfault in arma::eigs_sym (Accelerate is not fork-safe).
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
n_cores <- as.integer(commandArgs(TRUE)[1] %||% 4L)
N_GENES <- 2000L
top_var <- function(x) x[order(-apply(x, 1L, var))[seq_len(min(N_GENES, nrow(x)))], , drop = FALSE]
shuffle_rows <- function(x) t(apply(x, 1L, sample))
run <- function(x, data, id, seed) {
  set.seed(seed)
  net <- compute_network(x, n_cores = n_cores)
  t0 <- proc.time()[["elapsed"]]
  m <- detect_modules(net,
    resolution = c(0.5, 1, 2), seed = seed, n_cores = n_cores,
    test_k1 = TRUE, n_perm_k1 = 100L
  )
  k1 <- m$k1_test
  r <- data.frame(
    data = data, id = id, seed = seed, n_genes = nrow(x),
    n_modules = m$n_modules, has_structure = isTRUE(k1$has_structure),
    p_value = k1$p_value %||% NA_real_, n_perm = k1$n_perm_completed %||% NA_integer_,
    seconds = round(proc.time()[["elapsed"]] - t0, 1)
  )
  print(r, row.names = FALSE)
  r
}
rows <- list()
for (s in 1:10) {
  set.seed(s)
  x <- matrix(rnorm(N_GENES * 20L), N_GENES, 20L, dimnames = list(sprintf("g%04d", seq_len(N_GENES)), NULL))
  rows[[length(rows) + 1L]] <- run(x, "rnorm", paste0("rnorm", s), s)
}
se_files <- list.files("prepare_data/data", "_se\\.rds$", full.names = TRUE)
if (length(se_files)) {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  for (f in se_files) {
    sp <- sub("_se\\.rds$", "", basename(f))
    se <- readRDS(f)
    x <- as.matrix(assay(se[, colData(se)$tissue == "leaf"]))
    x <- top_var(x[apply(x, 1L, var) > 0, , drop = FALSE])
    set.seed(1L)
    rows[[length(rows) + 1L]] <- run(shuffle_rows(x), "pooideae_shuffled", sp, 1L)
    rows[[length(rows) + 1L]] <- run(x, "pooideae_real", sp, 1L)
  }
}
res <- do.call(rbind, rows)
utils::write.table(res, "dev/bench/k1_null_check.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nK = 1 rejection rate (has_structure) by data set:\n")
print(stats::aggregate(cbind(rejected = has_structure, n_modules) ~ data, res, mean))
