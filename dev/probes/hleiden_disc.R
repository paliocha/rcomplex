# Probe WP15b: modules disconnected in the 2-section (co-expression plus
# hog edges) for every cached joint_halves.R run of a data set, and
# hleiden's n_split / stalls. Writes <ds>_disconnected.tsv.
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/hleiden_disc.R wood
suppressPackageStartupMessages(library(rcomplex))
DS <- commandArgs(TRUE)[1]
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
jl <- new.env(parent = asNamespace("rcomplex"))
jl$JL_DIR <- "dev/probes"
sys.source("dev/probes/joint_louvain.R", jl)
hv <- readRDS(file.path(
  OUT, if (DS == "wood") "halves.rds" else "leaf_halves.rds"
))
ENG <- c("E2", "match", "hleiden0", "hleidenS")
rows <- list()
for (half in c("A", "B")) {
  d <- jl$jl_data(hv$nets[[half]], hv$ortho)
  key <- paste(rep(d$species, lengths(d$genes)), unlist(d$genes), sep = "\t")
  ij <- function(m) {
    cbind(m@i + 1L, rep.int(seq_len(ncol(m)), diff(m@p)))
  }
  e <- rbind(ij(d$co), ij(d$hog))
  e <- e[e[, 1] < e[, 2], ]
  rm(d)
  invisible(gc())
  fs <- list.files(OUT, sprintf(
    "^%s_runs_(%s)_.*_%s\\.rds$", DS, paste(ENG, collapse = "|"), half
  ), full.names = TRUE)
  for (f in fs) {
    r <- readRDS(f)
    lab <- match(r$lab[key], unique(r$lab[key]))
    k <- lab[e[, 1]] == lab[e[, 2]]
    g <- igraph::make_graph(as.vector(t(e[k, ])),
      n = length(lab),
      directed = FALSE
    )
    x <- unique(data.frame(lab, comp = igraph::components(g)$membership))
    big <- tabulate(lab) >= 20L
    bad <- unique(x$lab[duplicated(x$lab)])
    rows[[length(rows) + 1L]] <- data.frame(
      ds = DS, engine = r$engine, lambda = r$lambda, seed = r$seed,
      half = half, n_modules = length(unique(lab)),
      disconnected = length(bad), disconnected_20 = sum(big[bad]),
      n_split = if (is.null(r$params$n_split)) NA else r$params$n_split,
      n_moved = if (is.null(r$params$n_moved)) NA else r$params$n_moved,
      stalls = if (is.null(r$params$stalls)) NA else r$params$stalls,
      converged = if (is.null(r$params$converged)) {
        NA
      } else {
        r$params$converged
      },
      secs = r$secs
    )
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, file.path(OUT, paste0(DS, "_disconnected.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(stats::aggregate(
  cbind(disconnected, disconnected_20, n_split, stalls, secs) ~ engine,
  tab, function(x) c(mean = mean(x), max = max(x)),
  na.action = na.pass
))
