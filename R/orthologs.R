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
#'     \item{OrthoGroup}{Integer ortholog group identifier}
#'   }
#'
#' @export
parse_orthologs <- function(file, species1, species2) {
  if (!file.exists(file)) {
    stop("Ortholog file not found: ", file)
  }

  ortho <- utils::read.delim(file, stringsAsFactors = FALSE)

  ortho <- ortho |>
    dplyr::filter(.data$species == .env$species1) |>
    dplyr::group_by(.data$gene_content) |>
    dplyr::mutate(OrthoGroup = dplyr::cur_group_id()) |>
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
    dplyr::select("Species1", "Species2", "OrthoGroup")

  as.data.frame(ortho)
}
