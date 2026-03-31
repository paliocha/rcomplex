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
    stop("SummarizedExperiment and S4Vectors packages are required. ",
         "Install with: BiocManager::install(c('SummarizedExperiment', 'S4Vectors'))")
  }
  required <- c(species_col, gene_col, sample_col, assay_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("data missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Filter to species
  sp_data <- data[data[[species_col]] %in% species, , drop = FALSE]
  if (nrow(sp_data) == 0) {
    stop("No rows found for species '", species, "' in column '",
         species_col, "'")
  }

  # Pivot to genes x samples matrix
  genes <- unique(sp_data[[gene_col]])
  samples <- unique(sp_data[[sample_col]])
  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))
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
                              row.names = genes)
  # Rename hog_col to "hog" for consistency with rcomplex conventions
  if (!is.null(hog_col) && hog_col %in% names(rd)) {
    names(rd)[names(rd) == hog_col] <- "hog"
  }

  # Build colData
  sample_info <- sp_data[!duplicated(sp_data[[sample_col]]), , drop = FALSE]
  sample_info <- sample_info[match(samples, sample_info[[sample_col]]), ,
                              drop = FALSE]
  cd_cols <- sample_col
  if (!is.null(sample_metadata)) {
    cd_cols <- c(cd_cols, intersect(sample_metadata, names(sample_info)))
  }
  cd <- S4Vectors::DataFrame(sample_info[, cd_cols, drop = FALSE],
                              row.names = samples)

  SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(mat), assay_col),
    rowData = rd,
    colData = cd,
    metadata = list(species = species)
  )
}


#' Extract pairwise orthologs from two SummarizedExperiment objects
#'
#' Derives ortholog pairs by matching HOG identifiers in
#' \code{rowData()} of two species' SummarizedExperiment objects.
#' Returns a data frame in the same format as
#' \code{\link{parse_orthologs}}.
#'
#' @param se1,se2 \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   objects with a \code{hog} column in \code{rowData}.
#' @param hog_col Column name in \code{rowData} containing HOG
#'   identifiers (default \code{"hog"}).
#'
#' @return Data frame with columns \code{Species1} (gene IDs from
#'   \code{se1}), \code{Species2} (gene IDs from \code{se2}), and
#'   \code{hog}.
#'
#' @examples
#' \dontrun{
#' ortho <- extract_orthologs(se_bdis, se_hvul)
#' head(ortho)
#' # Can be passed directly to compare_neighborhoods()
#' comp <- compare_neighborhoods(net1, net2, ortho)
#' }
#'
#' @export
extract_orthologs <- function(se1, se2, hog_col = "hog") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required. ",
         "Install with: BiocManager::install('SummarizedExperiment')")
  }
  rd1 <- SummarizedExperiment::rowData(se1)
  rd2 <- SummarizedExperiment::rowData(se2)

  if (!hog_col %in% names(rd1)) {
    stop("rowData(se1) missing '", hog_col, "' column")
  }
  if (!hog_col %in% names(rd2)) {
    stop("rowData(se2) missing '", hog_col, "' column")
  }

  hogs1 <- as.character(rd1[[hog_col]])
  hogs2 <- as.character(rd2[[hog_col]])
  genes1 <- rownames(se1)
  genes2 <- rownames(se2)

  # Find shared HOGs
  shared <- intersect(hogs1[!is.na(hogs1)], hogs2[!is.na(hogs2)])
  if (length(shared) == 0) {
    return(data.frame(Species1 = character(0),
                      Species2 = character(0),
                      hog = character(0)))
  }

  # Build all pairwise combinations within each shared HOG
  sp1_by_hog <- split(genes1, hogs1)
  sp2_by_hog <- split(genes2, hogs2)

  pairs <- do.call(rbind, lapply(shared, function(h) {
    g1 <- sp1_by_hog[[h]]
    g2 <- sp2_by_hog[[h]]
    expand.grid(Species1 = g1, Species2 = g2, stringsAsFactors = FALSE)
  }))
  pairs$hog <- rep(shared, vapply(shared, function(h) {
    length(sp1_by_hog[[h]]) * length(sp2_by_hog[[h]])
  }, integer(1)))

  pairs
}
