# Does detect_modules()'s K = 1 test reject "no community structure" on null
# data? Its null rewires the graph keeping degrees, which removes the geometric
# clustering that a correlation graph from few samples has even on noise.
#
# Null data: (a) seeded rnorm expression; (b) Pooideae leaf per species with
# every gene's samples shuffled independently (skipped if prepare_data/data/
# is absent). Reference: (c) the same Pooideae data unshuffled.
# 2,000 top-variance genes, 20 samples, package defaults and the documented
# multi-resolution call, under both Leiden objectives; the observed
# statistic and the null range are kept so a tie shows apart from a margin. Writes dev/bench/k1_null_check.tsv.
#
# Run from the package root: Rscript dev/bench/k1_null_check.R [n_cores]
# On macOS with R's Accelerate BLAS set VECLIB_MAXIMUM_THREADS=1: forked K = 1
# workers otherwise segfault in arma::eigs_sym (Accelerate is not fork-safe).
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
n_cores <- if (length(args)) as.integer(args[1]) else 4L
N_GENES <- 2000L
top_var <- function(x) x[order(-apply(x, 1L, var))[seq_len(min(N_GENES, nrow(x)))], , drop = FALSE]
shuffle_rows <- function(x) t(apply(x, 1L, sample))
OBJ <- c("CPM", "modularity")
run <- function(x, data, id, seed) {
  set.seed(seed)
  net <- compute_network(x, n_cores = n_cores)
  do.call(rbind, lapply(OBJ, function(obj) {
    t0 <- proc.time()[["elapsed"]]
    m <- detect_modules(net,
      resolution = c(0.5, 1, 2), objective_function = obj, seed = seed,
      n_cores = n_cores, test_k1 = TRUE, n_perm_k1 = 100L
    )
    k1 <- m$k1_test
    r <- data.frame(
      data = data, id = id, seed = seed, objective = obj, n_genes = nrow(x),
      n_modules = m$n_modules, has_structure = isTRUE(k1$has_structure),
      p_value = k1$p_value %||% NA_real_,
      n_perm = k1$n_perm_completed %||% NA_integer_,
      lambda_obs = signif(k1$lambda_obs %||% NA_real_, 4),
      null_min = signif(min(k1$lambda_null), 4),
      null_max = signif(max(k1$lambda_null), 4),
      seconds = round(proc.time()[["elapsed"]] - t0, 1)
    )
    print(r, row.names = FALSE)
    r
  }))
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
    # materialise the shuffle here: as a lazy argument it would be drawn
    # inside run(), after run()'s own set.seed(), coupling the shuffle to
    # the detection seed (same data as before, since both seeds are 1)
    set.seed(1L)
    xs <- shuffle_rows(x)
    rows[[length(rows) + 1L]] <- run(xs, "pooideae_shuffled", sp, 1L)
    rows[[length(rows) + 1L]] <- run(x, "pooideae_real", sp, 1L)
  }
}
res <- do.call(rbind, rows)
# only a full run (with the Pooideae data) replaces the committed table
if (length(se_files)) {
  utils::write.table(res, "dev/bench/k1_null_check.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  message("prepare_data/data/ not found: random-data rows only, committed table left as is")
}
cat("\nK = 1 rejection rate (has_structure) by data set:\n")
print(stats::aggregate(cbind(rejected = has_structure, n_modules) ~ data + objective, res, mean))
