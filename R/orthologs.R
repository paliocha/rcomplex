#' Parse PLAZA ortholog group file
#'
#' Reads a PLAZA ortholog group file and extracts pairwise ortholog
#' relationships between two species.
#'
#' @param file Path to the PLAZA ortholog group file (tab-delimited,
#'   optionally gzipped).
#' @param species1 PLAZA species code for species 1 (e.g., `"potri"`).
#' @param species2 PLAZA species code for species 2 (e.g., `"piabi"`).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Species1}{Gene identifier for species 1}
#'     \item{Species2}{Gene identifier for species 2}
#'     \item{hog}{Integer ortholog group identifier}
#'   }
#'
#' @examples
#' \dontrun{
#' ortho <- parse_orthologs("orthologs.tsv", "species_A", "species_B")
#' head(ortho)
#' }
#'
#' @export
parse_orthologs <- function(file, species1, species2) {
  if (!file.exists(file)) {
    stop("Ortholog file not found: ", file)
  }

  ortho <- utils::read.delim(file)

  ortho <- ortho |>
    dplyr::filter(.data$species == .env$species1) |>
    dplyr::group_by(.data$gene_content) |>
    dplyr::mutate(hog = dplyr::cur_group_id()) |>
    dplyr::ungroup() |>
    tidyr::separate_longer_delim("gene_content", delim = ";") |>
    dplyr::filter(grepl(.env$species2, .data$gene_content, fixed = TRUE)) |>
    tidyr::separate_wider_delim("gene_content",
                                delim = ":",
                                names = c("prefix", "gene_content")) |>
    tidyr::separate_longer_delim("gene_content", delim = ",") |>
    dplyr::mutate(
      Species1 = .data$gene_id,
      Species2 = .data$gene_content
    ) |>
    dplyr::select("Species1", "Species2", "hog")

  as.data.frame(ortho)
}


#' Reduce orthogroups by merging correlated paralogs
#'
#' Within each ortholog group (HOG), paralogs with Pearson correlation above
#' \code{cor_threshold} are merged into a single representative gene via
#' Ward.D2 agglomerative clustering. Merged genes are replaced by their
#' averaged expression profile.
#'
#' This is an optional preprocessing step before \code{\link{compute_network}}.
#' It reduces redundancy from recent duplications where paralogs retain nearly
#' identical expression patterns, shrinking the expression matrix and avoiding
#' combinatorial blowup in downstream clique detection.
#'
#' Subfunctionalized paralogs (distinct expression programs) are preserved
#' as separate clusters. Zero-variance genes and genes not assigned to any
#' HOG are kept as-is.
#'
#' @param expr_matrix Numeric matrix (genes x samples) with gene identifiers
#'   as row names.
#' @param orthologs Data frame with columns \code{Species1} (or the column
#'   matching gene row names), \code{Species2}, and \code{hog}, as
#'   returned by \code{\link{parse_orthologs}}.
#' @param gene_col Character: which column of \code{orthologs} contains gene
#'   IDs matching row names of \code{expr_matrix} (default \code{"Species1"}).
#' @param cor_threshold Pearson correlation threshold for merging paralogs
#'   within a HOG (default 0.7). Higher values are more conservative (fewer
#'   merges).
#'
#' @return A list with components:
#'   \describe{
#'     \item{expr_matrix}{Reduced expression matrix (genes x samples) with
#'       row names. Merged genes have averaged expression.}
#'     \item{gene_map}{Data frame with columns \code{original} and
#'       \code{representative}, mapping every original gene to its
#'       representative in the reduced matrix.}
#'     \item{n_original}{Number of genes before reduction.}
#'     \item{n_reduced}{Number of genes after reduction.}
#'     \item{n_merged}{Number of genes absorbed into representatives.}
#'   }
#'
#' @examples
#' \dontrun{
#' reduced <- reduce_orthogroups(expr_matrix, orthologs)
#' reduced$expr_matrix  # reduced expression matrix
#' reduced$gene_map     # original -> representative mapping
#' }
#'
#' @export
reduce_orthogroups <- function(expr_matrix, orthologs,
                                gene_col = "Species1",
                                cor_threshold = 0.7) {
  if (!is.matrix(expr_matrix) || !is.numeric(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix")
  }
  if (is.null(rownames(expr_matrix))) {
    stop("expr_matrix must have row names (gene identifiers)")
  }
  if (!gene_col %in% names(orthologs)) {
    stop("orthologs must have column '", gene_col, "'")
  }
  if (!"hog" %in% names(orthologs)) {
    stop("orthologs must have column 'hog'")
  }
  if (cor_threshold < 0 || cor_threshold > 1) {
    stop("cor_threshold must be between 0 and 1")
  }

  gene_names <- rownames(expr_matrix)
  n_genes <- nrow(expr_matrix)

  # Build HOG membership: list of row indices (1-based) per HOG
  ortho_sub <- orthologs[orthologs[[gene_col]] %in% gene_names, , drop = FALSE]
  gene_to_row <- stats::setNames(seq_len(n_genes), gene_names)

  hog_genes <- split(ortho_sub[[gene_col]], ortho_sub$hog)
  hog_members <- lapply(hog_genes, function(genes) {
    as.integer(gene_to_row[genes[genes %in% gene_names]])
  })
  hog_members <- hog_members[lengths(hog_members) > 0]

  # Genes not in any HOG
  in_hog <- unique(unlist(hog_members))
  non_hog_idx <- as.integer(setdiff(seq_len(n_genes), in_hog))

  # Call C++
  result <- reduce_orthogroups_cpp(
    expr_matrix, hog_members, non_hog_idx, cor_threshold
  )

  # Attach row names to reduced matrix
  reduced_mat <- result$expr_matrix
  rownames(reduced_mat) <- gene_names[result$out_row_source]
  colnames(reduced_mat) <- colnames(expr_matrix)

  # Build gene mapping
  gene_map <- data.frame(
    original = gene_names[result$map_from],
    representative = gene_names[result$map_to]
  )

  list(
    expr_matrix = reduced_mat,
    gene_map = gene_map,
    n_original = result$n_original,
    n_reduced = result$n_reduced,
    n_merged = result$n_merged
  )
}
