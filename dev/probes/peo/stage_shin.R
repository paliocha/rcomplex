# Stage the Shin et al. 2021 compendia for rcomplex().
# Usage: Rscript stage_shin.R
# The compendia are used as published (quantile-normalised by the
# authors): no log, no expression-level filter. Drops genes with zero
# variance and genes absent from the FastOMA proteome. Writes
# inputs/all/<species>.rds (genes x experiments), inputs/species.tsv
# (the columns run.R reads) and inputs/tree.nwk.
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
E <- file.path(ROOT, "complex_Shin2021")
D <- file.path(E, "inputs")
dir.create(file.path(D, "all"), showWarnings = FALSE)
SP <- c(
  ath = "Arabidopsis_thaliana", mtr = "Medicago_truncatula",
  osa = "Oryza_sativa", ppa = "Physcomitrella_patens",
  stu = "Solanum_tuberosum", zma = "Zea_mays"
)
ortho <- data.table::fread(file.path(D, "orthologs.tsv"), data.table = FALSE)
info <- do.call(rbind, lapply(names(SP), function(code) {
  sp <- SP[[code]]
  x <- data.table::fread(
    file.path(D, "mendeley", paste0(code, "_compendium_matrix.txt")),
    data.table = FALSE
  )
  m <- as.matrix(x[-1L])
  rownames(m) <- x[[1L]]
  fa <- readLines(file.path(
    E, "fastoma", "in_folder", "proteome",
    paste0(sp, ".fa")
  ))
  prot <- sub("^>", "", fa[startsWith(fa, ">")])
  n_in <- nrow(m)
  m <- m[rownames(m) %in% prot, , drop = FALSE]
  m <- m[rowSums((m - rowMeans(m))^2) > 0, , drop = FALSE]
  saveRDS(m, file.path(D, "all", paste0(sp, ".rds")))
  data.frame(
    species = sp, taxid = NA, trait = "-", run = "all", n_leaf = NA,
    n_root = NA, n_samples = ncol(m), n_genes_in = n_in,
    n_genes_kept = nrow(m),
    n_genes_with_hog = sum(rownames(m) %in% ortho$gene[ortho$species == sp]),
    value_range = paste(signif(range(m), 3), collapse = " to ")
  )
}))
write.table(info, file.path(D, "species.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
file.copy(file.path(E, "fastoma", "in_folder", "species_tree.nwk"),
  file.path(D, "tree.nwk"),
  overwrite = TRUE
)
tree <- ape::read.tree(file.path(D, "tree.nwk"))
print(info, row.names = FALSE)
cat("tree:", ape::Ntip(tree), "tips, rooted", ape::is.rooted(tree), "\n")
