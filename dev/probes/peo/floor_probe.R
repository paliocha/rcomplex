# WP22 probe: where does the correlation floor sit, per species?
# Usage: Rscript floor_probe.R <exp> <run>
# Per species, from the staged matrix:
# (a) Aoki density valley: for r cut-offs 0.00 ... 0.99, pairs with
#     r >= cut-off, genes with at least one such pair, and density =
#     pairs / choose(genes linked, 2); r_valley is the cut-off at the
#     minimum density. Correlations are computed in row blocks, never
#     the full n x n matrix.
# (b) Null floor: each gene's samples permuted (.shuffle_genes) on a
#     random 2,000-gene subset; r_null_floor is the r above which the
#     shuffled tail supplies 1 % of the real edge count k = 0.03 *
#     choose(n, 2), i.e. the 1 - 0.0003 quantile of the shuffled r.
# (c) r_threshold of the run, and the raw r of the network's edges
#     (rebuilt with block_size = 512, read as .r_threshold() does).
# Writes score/<run>_floor.tsv and score/<run>_floor_curve.tsv.
suppressPackageStartupMessages(library(rcomplex))
args <- commandArgs(TRUE)
EXP <- args[1]
RUN <- args[2]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
E <- file.path(ROOT, paste0("complex_", EXP))
D <- file.path(E, "inputs")
S <- file.path(E, "score")
info <- read.delim(file.path(D, "species.tsv"))
info <- info[info$run == RUN, ]
nets <- read.delim(file.path(E, "runs", RUN, "networks.tsv"))
CUT <- 0:99 / 100
BLOCK <- 1000L
set.seed(1)

floor_one <- function(s) {
  t0 <- Sys.time()
  x <- readRDS(file.path(D, RUN, paste0(s, ".rds")))
  n <- nrow(x)
  zt <- rcomplex:::.standardise_for_cor(x, "pearson")
  # (a) histogram of pair r in 0.01 bins and each gene's best partner
  blocks <- split(seq_len(n), (seq_len(n) - 1L) %/% BLOCK)
  parts <- parallel::mclapply(blocks, function(b) {
    r <- crossprod(zt[, b, drop = FALSE], zt)
    r[cbind(seq_along(b), b)] <- -Inf
    up <- r[, -seq_len(max(b)), drop = FALSE]
    lo <- r[, seq_len(max(b)), drop = FALSE]
    lo[outer(seq_along(b), seq_len(max(b)), function(k, j) j <= b[k])] <- -Inf
    v <- c(up[up >= 0], lo[lo >= 0])
    list(
      h = tabulate(pmin(floor(v * 100) + 1L, 100L), 100L),
      max = apply(r, 1L, max)
    )
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  bad <- vapply(parts, inherits, NA, "try-error")
  if (any(bad)) stop(s, ": a block failed: ", parts[bad][[1L]])
  h <- Reduce(`+`, lapply(parts, `[[`, "h"))
  mx <- unlist(lapply(parts, `[[`, "max"))
  pairs <- rev(cumsum(rev(h)))
  linked <- vapply(CUT, function(cc) sum(mx >= cc), 1)
  dens <- ifelse(linked >= 2, pairs / choose(linked, 2), NA)
  curve <- data.frame(
    species = s, cutoff = CUT, pairs = pairs,
    genes_linked = linked, density = dens
  )
  r_valley <- CUT[which.min(dens)]
  # (b) shuffled tail on 2,000 genes
  idx <- sample.int(n, min(2000L, n))
  r0 <- stats::cor(t(rcomplex:::.shuffle_genes(x[idx, , drop = FALSE])))
  r0 <- r0[upper.tri(r0)]
  r_null <- unname(stats::quantile(r0, 1 - 0.01 * 0.03))
  # (c) raw r of the network's edges
  net <- compute_network(x,
    density = 0.03, sign = "positive", block_size = 512,
    n_cores = N_CORES
  )
  e <- rcomplex:::.adj_edges(net$network, net$threshold)
  up <- e$rows < e$cols
  i <- match(rownames(net$network)[e$rows[up]], rownames(x))
  j <- match(rownames(net$network)[e$cols[up]], rownames(x))
  rm(net, e)
  ch <- split(seq_along(i), (seq_along(i) - 1L) %/% (2^24 %/% nrow(zt)))
  re <- unlist(parallel::mclapply(ch, function(k) {
    colSums(zt[, i[k], drop = FALSE] * zt[, j[k], drop = FALSE])
  }, mc.cores = N_CORES), use.names = FALSE)
  cat(sprintf(
    "%s: %d genes, %.0f s\n", s, n,
    as.numeric(Sys.time() - t0, units = "secs")
  ))
  list(curve = curve, row = data.frame(
    species = s, n_runs = ncol(x), n_genes = n,
    r_threshold = nets$r_threshold[nets$species == s],
    r_threshold_rebuilt = min(re), r_valley = r_valley,
    r_null_floor = round(r_null, 3), n_edges = length(re),
    edges_below_valley = sum(re < r_valley),
    share_below_valley = round(mean(re < r_valley), 4),
    edges_below_null = sum(re < r_null),
    share_below_null = round(mean(re < r_null), 4),
    agree_0.1 = abs(r_valley - r_null) <= 0.1
  ))
}

out <- lapply(info$species, floor_one)
tab <- do.call(rbind, lapply(out, `[[`, "row"))
data.table::fwrite(tab, file.path(S, paste0(RUN, "_floor.tsv")),
  sep = "\t", quote = FALSE
)
data.table::fwrite(do.call(rbind, lapply(out, `[[`, "curve")),
  file.path(S, paste0(RUN, "_floor_curve.tsv")),
  sep = "\t", quote = FALSE
)
print(tab, row.names = FALSE)
