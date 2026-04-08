#!/usr/bin/env Rscript
# One-time conversion: long-format vst_hog.RDS -> per-species SummarizedExperiment .rds
#
# Input:  nf-rcomplex/vst_hog.RDS (14.6M-row tibble, long format)
# Output: prepare_data/data/{species}_se.rds (one SE per species)
#
# The input tibble has columns:
#   abbrev, species, life_cycle, pair, tissue, time.point, day, real_day,
#   gene_id, protein_id, HOG, vst.count, vst.zscore, sample_id, running.no, R1, R2
#
# Each SE contains:
#   assay:   gene x sample matrix of vst.count values
#   rowData: gene_id, protein_id, hog (HOG renamed for rcomplex convention)
#   colData: sample_id, tissue, time.point, day, real_day
#
# Usage:
#   Rscript prepare_data/convert_to_se.R [species_csv]
#
# If species_csv is omitted, converts the 8 species in the plan:
#   BDIS, BSYL, HVUL, HJUB, BMAX, BMED, VBRO, FPRA
#
# Additional species in the data not converted by default:
#   PANN (annual, Poa pair), PSUP (perennial, Poa pair), BMED2BMAX (BMED re-annotated on BMAX genome)

args <- commandArgs(trailingOnly = TRUE)

default_species <- c("BDIS", "BSYL", "HVUL", "HJUB", "BMAX", "BMED", "VBRO", "FPRA")
target_species <- if (length(args) >= 1) {
  strsplit(args[1], ",")[[1]]
} else {
  default_species
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
    !requireNamespace("S4Vectors", quietly = TRUE)) {
  stop("Install Bioconductor packages:\n",
       "  BiocManager::install(c('SummarizedExperiment', 'S4Vectors'))")
}

input_path <- "nf-rcomplex/vst_hog.RDS"
output_dir <- "prepare_data/data"

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path,
       "\nRun from the rcomplex repo root.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading ", input_path, " ...")
dat <- readRDS(input_path)
message("  ", nrow(dat), " rows, ", length(unique(dat$abbrev)), " species")

available <- unique(dat$abbrev)
missing <- setdiff(target_species, available)
if (length(missing) > 0) {
  stop("Species not found in data: ", paste(missing, collapse = ", "),
       "\nAvailable: ", paste(sort(available), collapse = ", "))
}

for (sp in target_species) {
  message("\n--- ", sp, " ---")
  sp_data <- dat[dat$abbrev == sp, ]

  genes   <- unique(sp_data$gene_id)
  samples <- unique(sp_data$sample_id)
  message("  ", length(genes), " genes x ", length(samples), " samples")

  # Pivot to gene x sample matrix
  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))
  idx <- match(sp_data$gene_id, genes)
  jdx <- match(sp_data$sample_id, samples)
  mat[cbind(idx, jdx)] <- sp_data$vst.count

  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message("  WARNING: ", n_na, " NA values in expression matrix (",
            round(100 * n_na / length(mat), 2), "%)")
  }

  # Build rowData (gene-level metadata, deduplicated)
  gene_info <- sp_data[!duplicated(sp_data$gene_id), ]
  gene_info <- gene_info[match(genes, gene_info$gene_id), ]
  rd <- S4Vectors::DataFrame(
    gene_id    = gene_info$gene_id,
    protein_id = gene_info$protein_id,
    hog        = gene_info$HOG,          # renamed from HOG to hog
    row.names  = genes
  )
  n_hog <- sum(!is.na(rd$hog))
  message("  ", n_hog, "/", length(genes), " genes with HOG assignment (",
          round(100 * n_hog / length(genes), 1), "%)")

  # Build colData (sample-level metadata, deduplicated)
  sample_info <- sp_data[!duplicated(sp_data$sample_id), ]
  sample_info <- sample_info[match(samples, sample_info$sample_id), ]
  cd <- S4Vectors::DataFrame(
    sample_id  = sample_info$sample_id,
    tissue     = sample_info$tissue,
    time_point = as.character(sample_info$time.point),
    day        = sample_info$day,
    real_day   = sample_info$real_day,
    row.names  = samples
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(vst = mat),
    rowData = rd,
    colData = cd
  )

  out_path <- file.path(output_dir, paste0(sp, "_se.rds"))
  saveRDS(se, out_path)
  message("  Saved: ", out_path, " (", round(file.size(out_path) / 1e6, 1), " MB)")
}

message("\nDone. ", length(target_species), " SummarizedExperiment files in ", output_dir, "/")
message("Samplesheet species: ", paste(target_species, collapse = ", "))
