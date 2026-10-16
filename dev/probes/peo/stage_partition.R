# WP24 probe: stage a partition run (per-sample factor per species).
# Usage: Rscript stage_partition.R <exp> <run>
#   LeafvRoot organ_partition: all PO-labelled runs; organs with >= 50
#     runs per species; partition = organ. Same gene filter as stage.R
#     on the pooled kept runs. Writes inputs/organ_partition/<sp>.rds.
#   Shin2021 / C3vC4 kmeans_partition: the staged matrix of the pooled
#     run (all / leaf); k = max(2, min(20, floor(n_runs / 100))) and
#     kmeans(t(x), k, nstart = 10) after set.seed(1). No matrix copy:
#     species.tsv names the source run.
# Writes inputs/<run>/part.rds (named list of factors in column order),
# inputs/<run>/species.tsv (levels and run counts) and appends the run's
# rows to inputs/species.tsv (source_run, levels).
args <- commandArgs(TRUE)
EXP <- args[1]
RUN <- args[2]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
PEO <- "/mnt/users/martpali/AnnualPerennial/coexpr/PEO"
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
E <- file.path(ROOT, paste0("complex_", EXP))
D <- file.path(E, "inputs")
dir.create(file.path(D, RUN), showWarnings = FALSE)
info0 <- read.delim(file.path(D, "species.tsv"))
ortho <- data.table::fread(file.path(D, "orthologs.tsv"), data.table = FALSE)
lv_string <- function(f) {
  t <- sort(table(f), decreasing = TRUE)
  paste0(names(t), ":", t, collapse = ",")
}

organ_one <- function(sp) {
  map <- read.delim(file.path(PEO, "taxid_species_mapping.tsv"))
  taxid <- map$taxid[match(sp, map$species)]
  ed <- file.path(PEO, "expression_data")
  po <- read.delim(file.path(ed, sprintf("taxid%s_PO.tsv", taxid)))
  f <- file.path(ed, sprintf("taxid%s_tpm.tsv.gz", taxid))
  cols <- names(data.table::fread(f, nrows = 0L))
  po <- po[po$SampleID %in% cols[-1L] & !duplicated(po$SampleID) &
    grepl("^PO:", po$PO), ]
  n <- table(po$PO)
  keep <- names(n)[n >= 50L]
  po <- po[po$PO %in% keep, ]
  tpm <- data.table::fread(f,
    select = c(1L, match(po$SampleID, cols)),
    data.table = FALSE
  )
  x <- as.matrix(tpm[-1L])
  rownames(x) <- tpm[[1L]]
  rm(tpm)
  x <- x[rowMeans(x >= 1) >= 0.05, , drop = FALSE]
  y <- log2(x + 1)
  y <- y[rowSums((y - rowMeans(y))^2) > 0, , drop = FALSE]
  saveRDS(y, file.path(D, RUN, paste0(sp, ".rds")))
  part <- factor(po$PO)
  list(part = part, row = data.frame(
    species = sp, n_samples = ncol(y), n_genes_kept = nrow(y),
    n_genes_with_hog = sum(rownames(y) %in% ortho$gene[ortho$species == sp]),
    source_run = RUN, n_levels = nlevels(part), levels = lv_string(part)
  ))
}

kmeans_one <- function(sp, src) {
  x <- readRDS(file.path(D, src, paste0(sp, ".rds")))
  k <- max(2L, min(20L, floor(ncol(x) / 100)))
  set.seed(1)
  km <- stats::kmeans(t(x), k, nstart = 10)
  part <- factor(paste0("k", km$cluster))
  r <- info0[info0$species == sp & info0$run == src, ][1L, ]
  list(part = part, row = data.frame(
    species = sp, n_samples = ncol(x), n_genes_kept = nrow(x),
    n_genes_with_hog = r$n_genes_with_hog, source_run = src,
    n_levels = k, levels = lv_string(part)
  ))
}

if (RUN == "organ_partition") {
  sp <- unique(info0$species[info0$run == "leaf"])
  res <- parallel::mclapply(sp, organ_one,
    mc.cores = min(N_CORES, 4L), mc.preschedule = FALSE
  )
} else {
  src <- if (EXP == "Shin2021") "all" else "leaf"
  sp <- unique(info0$species[info0$run == src])
  res <- lapply(sp, kmeans_one, src = src)
}
bad <- vapply(res, inherits, NA, "try-error")
if (any(bad)) stop("staging failed: ", res[bad][[1L]])
names(res) <- sp
part <- lapply(res, `[[`, "part")
saveRDS(part, file.path(D, RUN, "part.rds"))
tab <- do.call(rbind, lapply(res, `[[`, "row"))

# Shin: the most frequent tissue annotation per k-means cluster
if (EXP == "Shin2021") {
  meta <- file.path(
    D, "mendeley",
    "Transcriptome_expression_compendia_metadata.xlsx"
  )
  sheet <- c(
    Arabidopsis_thaliana = "A.thaliana", Oryza_sativa = "O.sativa",
    Zea_mays = "Z.mays", Medicago_truncatula = "M.truncatula",
    Solanum_tuberosum = "S.tuberosum", Physcomitrella_patens = "P.patens"
  )
  tab$top_annot <- vapply(sp, function(s) {
    m <- tryCatch(readxl::read_excel(meta, sheet[[s]]),
      error = function(e) NULL
    )
    if (is.null(m)) return(NA_character_)
    x <- readRDS(file.path(D, "all", paste0(s, ".rds")))
    an <- m$Annot[match(colnames(x), m$experimentAccession)]
    paste(vapply(levels(part[[s]]), function(l) {
      t <- sort(table(an[part[[s]] == l]), decreasing = TRUE)
      paste0(
        l, "=", if (length(t)) names(t)[1L] else "NA",
        "(", if (length(t)) t[[1L]] else 0L, "/", sum(part[[s]] == l), ")"
      )
    }, ""), collapse = ",")
  }, "")
}
rownames(tab) <- NULL
write.table(tab, file.path(D, RUN, "species.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(tab, row.names = FALSE)

# rows for run.R and score.R in the experiment's species.tsv
new <- data.frame(
  species = tab$species,
  taxid = info0$taxid[match(tab$species, info0$species)],
  trait = info0$trait[match(tab$species, info0$species)], run = RUN,
  n_samples = tab$n_samples, n_genes_kept = tab$n_genes_kept,
  n_genes_with_hog = tab$n_genes_with_hog, source_run = tab$source_run
)
old <- info0[info0$run != RUN, ]
if (is.null(old$source_run)) old$source_run <- old$run
all_cols <- union(names(old), names(new))
fill <- function(d) {
  for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA
  d[all_cols]
}
write.table(rbind(fill(old), fill(new)), file.path(D, "species.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
