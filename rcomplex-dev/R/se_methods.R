#' Build per-species SummarizedExperiment from long-format data
#'
#' Pivots a long-format expression table (one row per gene-sample
#' observation) into a SummarizedExperiment for a single species.
#' This is a convenience helper, not core rcomplex functionality.
#'
#' @param data Data frame in long format with columns for species,
#'   gene ID, sample ID, and expression value.
#' @param species Character, species abbreviation to extract.
#' @param assay_col Column containing expression values.
#' @param gene_col Column with gene identifiers.
#' @param sample_col Column with sample identifiers.
#' @param species_col Column with species abbreviations.
#' @param hog_col Column with HOG identifiers, or \code{NULL} to omit.
#' @param gene_metadata Additional columns for \code{rowData}.
#' @param sample_metadata Additional columns for \code{colData}.
#' @return A SummarizedExperiment.
#' @noRd
build_se <- function(data, species,
                     assay_col = "vst.count",
                     gene_col = "gene_id",
                     sample_col = "sample_id",
                     species_col = "abbrev",
                     hog_col = "HOG",
                     gene_metadata = NULL,
                     sample_metadata = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop(
      "SummarizedExperiment and S4Vectors packages are required. ",
      "Install with: ",
      "BiocManager::install(c('SummarizedExperiment', 'S4Vectors'))"
    )
  }
  required <- c(species_col, gene_col, sample_col, assay_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "data missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Filter to species
  sp_data <- data[data[[species_col]] %in% species, , drop = FALSE]
  if (nrow(sp_data) == 0) {
    stop(
      "No rows found for species '", species, "' in column '",
      species_col, "'"
    )
  }

  # Pivot to genes x samples matrix
  genes <- unique(sp_data[[gene_col]])
  samples <- unique(sp_data[[sample_col]])
  mat <- matrix(NA_real_,
    nrow = length(genes), ncol = length(samples),
    dimnames = list(genes, samples)
  )
  idx <- match(sp_data[[gene_col]], genes)
  jdx <- match(sp_data[[sample_col]], samples)
  mat[cbind(idx, jdx)] <- sp_data[[assay_col]]

  # Build rowData
  gene_info <- sp_data[!duplicated(sp_data[[gene_col]]), , drop = FALSE]
  gene_info <- gene_info[match(genes, gene_info[[gene_col]]), , drop = FALSE]
  rd_cols <- gene_col
  if (!is.null(hog_col) && hog_col %in% names(gene_info)) {
    rd_cols <- c(rd_cols, hog_col)
  }
  if (!is.null(gene_metadata)) {
    rd_cols <- c(rd_cols, intersect(gene_metadata, names(gene_info)))
  }
  rd <- S4Vectors::DataFrame(gene_info[, rd_cols, drop = FALSE],
    row.names = genes
  )
  # Rename hog_col to "hog" for consistency with rcomplex conventions
  if (!is.null(hog_col) && hog_col %in% names(rd)) {
    names(rd)[names(rd) == hog_col] <- "hog"
  }

  # Build colData
  sample_info <- sp_data[!duplicated(sp_data[[sample_col]]), , drop = FALSE]
  sample_info <- sample_info[match(samples, sample_info[[sample_col]]), ,
    drop = FALSE
  ]
  cd_cols <- sample_col
  if (!is.null(sample_metadata)) {
    cd_cols <- c(cd_cols, intersect(sample_metadata, names(sample_info)))
  }
  cd <- S4Vectors::DataFrame(sample_info[, cd_cols, drop = FALSE],
    row.names = samples
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(mat), assay_col),
    rowData = rd,
    colData = cd,
    metadata = list(species = species)
  )
}
