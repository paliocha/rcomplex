# Run rcomplex() on one staged PEO experiment and run, then
# edge_history() and, for a trait experiment, a trait pass that
# classifies the same edges with the two trait groups as clades.
# Usage: Rscript run.R <exp> <run>
# Networks are built blockwise and passed as networks = (see below).
# Writes runs/<run>/: write_rcomplex() tables, res.rds (no networks),
# networks.tsv, edge_history_{cliques,branches}.tsv, trait_pass{,_null}.tsv,
# timing.tsv. max_mb is the R heap of the main process from gc(); forked
# workers are not in it, so read MaxRSS from sacct as well.
cat("rcomplex", as.character(packageVersion("rcomplex")), "\n")
print(sessionInfo())
suppressPackageStartupMessages(library(rcomplex))
args <- commandArgs(TRUE)
EXP <- args[1]
RUN <- args[2]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
D <- file.path(ROOT, paste0("complex_", EXP), "inputs")
OUT <- file.path(ROOT, paste0("complex_", EXP), "runs", RUN)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

info <- read.delim(file.path(D, "species.tsv"))
info <- info[info$run == RUN, ]
if (nrow(info) == 0L) stop("run: no staged species for ", EXP, " ", RUN)
sp <- info$species
tree <- ape::read.tree(file.path(D, "tree.nwk"))

timing <- list()
step <- function(name, expr) {
  gc(reset = TRUE)
  t <- system.time(v <- force(expr))[["elapsed"]]
  # max used cells: Ncells are 56 bytes, Vcells 8
  mb <- sum(gc()[, "max used"] * c(56, 8)) / 2^20
  timing[[name]] <<- data.frame(
    step = name, seconds = round(t, 1), max_mb = round(mb)
  )
  cat(sprintf("step %s: %.0f s\n", name, t))
  v
}

big <- which.max(info$n_genes_kept)
cat("largest species:", sp[big], info$n_genes_kept[big], "genes\n")
ortho <- read_orthologs(file.path(D, "orthologs.tsv"), species = sp)
block <- NULL
if (RUN == "organ_block") block <- readRDS(file.path(D, RUN, "block.rds"))[sp]
# WP24 partition runs: per-sample factor in inputs/<run>/part.rds; the
# matrix may come from another run (species.tsv column source_run)
part <- NULL
pf <- file.path(D, RUN, "part.rds")
if (file.exists(pf)) part <- readRDS(pf)[sp]
src <- if (is.null(info$source_run)) RUN else info$source_run
src[is.na(src)] <- RUN
names(src) <- sp

# The driver builds dense networks from expr (a 71k-gene wheat matrix
# segfaults), so build each blockwise here and pass networks =. This
# copies the driver's expr path: the wiring layer under block, and
# n_samples and gene_mean from the raw matrix (the null's deciles).
nets <- list()
for (s in sp) {
  x <- readRDS(file.path(D, src[[s]], paste0(s, ".rds")))
  xw <- if (is.null(block)) x else split_layers(x, block[[s]])$wiring
  # partition needs the dense build (64-bit Armadillo library)
  nets[[s]] <- step(paste0("network_", s), if (is.null(part)) {
    compute_network(xw,
      density = 0.03, sign = "positive", block_size = 512,
      n_cores = N_CORES
    )
  } else {
    compute_network(xw,
      density = 0.03, sign = "positive", partition = part[[s]],
      n_cores = N_CORES
    )
  })
  nets[[s]]$params$n_samples <- ncol(x)
  nets[[s]]$gene_mean <- rowMeans(x)[rownames(nets[[s]]$network)]
  rm(x, xw)
}

res <- step("rcomplex", rcomplex(
  networks = nets, orthologs = ortho, clades = clades_from_tree(tree),
  alpha = 0.1, modules = FALSE, null = TRUE, n_cores = N_CORES, seed = 1
))
rm(nets)
print(res)
summary(res)
num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
nets <- data.frame(
  species = sp,
  n_genes = vapply(res$networks, function(n) num(n$n_genes), 1),
  n_samples = vapply(res$networks, function(n) num(n$params$n_samples), 1),
  density = vapply(res$networks, function(n) num(n$params$density), 1),
  r_threshold = vapply(res$networks, function(n) {
    num(n$params$r_threshold)
  }, 1)
)
res$networks <- NULL
tsv <- function(x, f) {
  data.table::fwrite(x, file.path(OUT, f),
    sep = "\t", na = "NA", quote = FALSE
  )
}
tsv(nets, "networks.tsv")
step("write", {
  write_rcomplex(res, OUT)
  saveRDS(res, file.path(OUT, "res.rds"))
})

eh <- step("edge_history", edge_history(res$classification, tree, sp))
tsv(eh$cliques, "edge_history_cliques.tsv")
tsv(eh$branches, "edge_history_branches.tsv")

if (!all(info$trait == "-")) {
  groups <- split(info$species, info$trait)
  tp <- step("trait_pass", rcomplex:::.driver_cliques(res$edges, sp,
    clades = groups, alpha = 0.1, min_power = 0.8
  ))
  tsv(tp$classification, "trait_pass.tsv")
  print(table(tp$classification$hog_class[!duplicated(tp$classification$hog)]))
  # the same pass on the null edges: the chance count of trait_specific
  tp0 <- step("trait_pass_null", rcomplex:::.driver_cliques(res$edges_null,
    sp, clades = groups, alpha = 0.1, min_power = 0.8
  ))
  tsv(tp0$classification, "trait_pass_null.tsv")
}

timing <- do.call(rbind, timing)
timing$job_id <- Sys.getenv("SLURM_JOB_ID", NA)
tsv(timing, "timing.tsv")
print(timing, row.names = FALSE)
