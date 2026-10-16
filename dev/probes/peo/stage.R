# Stage the expression inputs of one PEO experiment for rcomplex().
# Usage: Rscript stage.R <exp>     (exp: LeafvRoot, C3vC4 or nodulation)
# Needs inputs/orthologs.tsv from hogs_at_level.py. Per species and run:
# keep the runs of the organ, log2(TPM + 1), drop genes with TPM < 1 in
# more than 95 % of the runs and genes with zero variance. Writes
# inputs/<run>/<species>.rds (genes x runs), inputs/organ_block/block.rds,
# inputs/species.tsv and inputs/tree.nwk.
args <- commandArgs(TRUE)
EXP <- args[1]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
PEO <- Sys.getenv("PEO_DATA", "/mnt/users/martpali/AnnualPerennial/coexpr/PEO")
TREE <- Sys.getenv("PEO_TREE", paste0(
  "/mnt/project/FjellheimLab/martpali/AnnualPerennial/peo-fastoma/",
  "output/species_tree_checked.nwk"
))
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
RUNS <- list(
  LeafvRoot = c("leaf", "root", "organ_block"), C3vC4 = "leaf",
  nodulation = "root"
)
ORGANS <- list(leaf = "leaf", root = "root", organ_block = c("leaf", "root"))
PO <- c(leaf = "PO:0025034", root = "PO:0009005")

D <- file.path(ROOT, paste0("complex_", EXP), "inputs")
sets <- read.delim(file.path(ROOT, "peo", "scripts", "sets.tsv"))
sets <- sets[sets$experiment == EXP, ]
if (nrow(sets) == 0L) stop("exp: no species for ", EXP, " in sets.tsv")
map <- read.delim(file.path(PEO, "taxid_species_mapping.tsv"))
ortho <- data.table::fread(file.path(D, "orthologs.tsv"), data.table = FALSE)
for (r in RUNS[[EXP]]) dir.create(file.path(D, r), showWarnings = FALSE)

stage_one <- function(sp) {
  taxid <- map$taxid[match(sp, map$species)]
  ed <- file.path(PEO, "expression_data")
  po <- read.delim(file.path(ed, sprintf("taxid%s_PO.tsv", taxid)))
  f <- file.path(ed, sprintf("taxid%s_tpm.tsv.gz", taxid))
  cols <- names(data.table::fread(f, nrows = 0L))
  organ <- names(PO)[match(po$PO, PO)]
  ok <- !is.na(organ) & po$SampleID %in% cols[-1L] &
    !duplicated(po$SampleID)
  ids <- po$SampleID[ok]
  organ <- organ[ok]
  tpm <- data.table::fread(f,
    select = c(1L, match(ids, cols)),
    data.table = FALSE
  )
  m <- as.matrix(tpm[-1L])
  rownames(m) <- tpm[[1L]]
  rm(tpm)
  if (anyDuplicated(rownames(m))) stop(sp, ": duplicated gene ids")
  with_hog <- ortho$gene[ortho$species == sp]
  info <- list()
  block <- NULL
  for (r in RUNS[[EXP]]) {
    j <- organ %in% ORGANS[[r]]
    x <- m[, j, drop = FALSE]
    x <- x[rowMeans(x >= 1) >= 0.05, , drop = FALSE]
    y <- log2(x + 1)
    y <- y[rowSums((y - rowMeans(y))^2) > 0, , drop = FALSE]
    saveRDS(y, file.path(D, r, paste0(sp, ".rds")))
    if (r == "organ_block") block <- factor(organ[j], c("leaf", "root"))
    info[[r]] <- data.frame(
      species = sp, taxid = taxid,
      trait = sets$trait[sets$species == sp], run = r,
      n_leaf = sum(organ == "leaf"), n_root = sum(organ == "root"),
      n_samples = ncol(y), n_genes_kept = nrow(y),
      n_genes_with_hog = sum(rownames(y) %in% with_hog)
    )
  }
  list(info = do.call(rbind, info), block = block)
}

res <- parallel::mclapply(sets$species, stage_one,
  mc.cores = N_CORES, mc.preschedule = FALSE
)
bad <- vapply(res, inherits, NA, "try-error")
if (any(bad)) {
  stop(
    "staging failed for ", paste(sets$species[bad], collapse = ", "),
    ": ", res[bad][[1L]]
  )
}
names(res) <- sets$species
if ("organ_block" %in% RUNS[[EXP]]) {
  saveRDS(lapply(res, `[[`, "block"), file.path(D, "organ_block", "block.rds"))
}
info <- do.call(rbind, lapply(res, `[[`, "info"))
rownames(info) <- NULL
write.table(info, file.path(D, "species.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
tree <- ape::keep.tip(ape::read.tree(TREE), sets$species)
ape::write.tree(tree, file.path(D, "tree.nwk"))
print(info, row.names = FALSE)
cat("tree:", ape::Ntip(tree), "tips, rooted", ape::is.rooted(tree), "\n")
