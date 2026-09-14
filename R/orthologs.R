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
#' @details
#' Implemented with \pkg{data.table} (\code{\link[data.table]{fread}} only)
#' rather than \pkg{dplyr}/\pkg{tidyr} (used through rcomplex 0.3.0). This
#' was the only site in the package using either package, so the switch
#' drops both from Imports; it is a dependency-economy change
#' (data.table has zero hard dependencies), not a hot-path optimization --
#' this parser runs once per ortholog file, not inside any permutation or
#' per-replicate loop. `hog` ids are assigned by the sorted order of each
#' row's full `gene_content` string, matching the previous
#' \code{dplyr::cur_group_id()} convention (which numbers groups by sorted
#' key, not order of appearance) so pipelines that persist HOG ids across
#' a rerun see the same numbering.
#'
#' @export
parse_orthologs <- function(file, species1, species2) {
  if (!file.exists(file)) {
    stop("Ortholog file not found: ", file)
  }

  dt <- data.table::fread(file, sep = "\t", header = TRUE,
                          showProgress = FALSE, data.table = FALSE)
  dt <- dt[dt$species == species1, , drop = FALSE]

  # hog id: sorted-key group numbering, matching dplyr::cur_group_id()
  hog <- match(dt$gene_content, sort(unique(dt$gene_content)))

  # one row per ";"-delimited per-species chunk of gene_content
  chunks <- strsplit(dt$gene_content, ";", fixed = TRUE)
  n_chunks <- lengths(chunks)
  gene_id <- rep(dt$gene_id, n_chunks)
  hog <- rep(hog, n_chunks)
  chunk <- unlist(chunks, use.names = FALSE)

  keep <- grepl(species2, chunk, fixed = TRUE)
  gene_id <- gene_id[keep]
  hog <- hog[keep]
  chunk <- chunk[keep]

  # split "prefix:gene[,gene...]" on the first ":", then explode on ","
  rest <- sub("^[^:]*:", "", chunk)
  gene_content_list <- strsplit(rest, ",", fixed = TRUE)
  n_genes <- lengths(gene_content_list)

  data.frame(
    Species1 = rep(gene_id, n_genes),
    Species2 = unlist(gene_content_list, use.names = FALSE),
    hog = rep(hog, n_genes)
  )
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


#' Extract pairwise orthologs, optionally with paralog-reduced gene names
#'
#' Convenience wrapper that calls \code{\link{extract_orthologs}} for every
#' species pair and, when \code{reductions} is supplied, maps gene names
#' through the gene maps produced by \code{\link{reduce_orthogroups}}.  This
#' replaces ~20 lines of boilerplate when combining SummarizedExperiment
#' objects with paralog reduction outputs.
#'
#' Paralog reduction is a lossy step: correlated paralogs are averaged into
#' one representative gene, so per-paralog identity is gone by the time
#' \code{orthologs} reaches downstream consumers. That trade-off pays for
#' itself for consumers that need one counterpart per gene (module
#' preservation, the species-graph clique backend), but the gene-graph
#' clique backend (\code{\link{gene_clique_graph}} /
#' \code{\link{classify_gene_cliques}}) is built to resolve which paralog
#' copy is conserved, and needs every original paralog as its own node to do
#' that. Pass \code{reductions = NULL} (the default) to skip paralog
#' reduction entirely and keep every gene at its original identity.
#'
#' @param se_list Named list of
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} objects, keyed
#'   by species code.
#' @param reductions Named list of \code{\link{reduce_orthogroups}} outputs,
#'   keyed by the same species codes as \code{se_list}, or \code{NULL}
#'   (the default) to skip paralog reduction and keep original gene
#'   identities. When supplied, each element must contain a \code{$gene_map}
#'   data frame with columns \code{original} and \code{representative}.
#' @param hog_col Column name in \code{rowData} containing HOG identifiers
#'   (default \code{"hog"}).
#'
#' @return A data frame with columns \code{Species1}, \code{Species2}, and
#'   \code{hog}. When \code{reductions} is supplied, gene names have been
#'   replaced by their reduced representatives and duplicate rows (arising
#'   when multiple original genes map to the same representative) are
#'   removed; when \code{reductions} is \code{NULL}, every original paralog
#'   pair within a shared HOG is retained as its own row.
#'
#' @examples
#' \dontrun{
#' # With paralog reduction
#' ortho <- prepare_orthologs(se_list, reductions)
#'
#' # Without paralog reduction -- keeps every paralog, for gene_clique_graph()
#' ortho_full <- prepare_orthologs(se_list)
#' head(ortho_full)
#' }
#'
#' @export
prepare_orthologs <- function(se_list, reductions = NULL, hog_col = "hog") {
  # --- validation ---
  if (!is.list(se_list) || is.null(names(se_list))) {
    stop("se_list must be a named list")
  }
  sp_names <- names(se_list)
  if (length(sp_names) < 2) {
    stop("se_list must contain at least two species")
  }

  if (!is.null(reductions)) {
    if (!is.list(reductions) || is.null(names(reductions))) {
      stop("reductions must be a named list")
    }
    missing_sp <- setdiff(names(se_list), names(reductions))
    if (length(missing_sp) > 0) {
      stop("reductions missing species present in se_list: ",
           paste(missing_sp, collapse = ", "))
    }
    for (sp in names(se_list)) {
      if (is.null(reductions[[sp]]$gene_map)) {
        stop("reductions[['", sp, "']] must have a $gene_map element")
      }
    }
  }

  pairs <- utils::combn(sp_names, 2, simplify = FALSE)

  result_list <- lapply(pairs, function(pair) {
    sp1 <- pair[1]
    sp2 <- pair[2]

    ortho <- extract_orthologs(se_list[[sp1]], se_list[[sp2]],
                               hog_col = hog_col)
    if (nrow(ortho) == 0 || is.null(reductions)) return(ortho)

    # Map Species1 through sp1 gene_map
    gm1 <- reductions[[sp1]]$gene_map
    idx1 <- match(ortho$Species1, gm1$original)
    mapped1 <- gm1$representative[idx1]
    ortho$Species1 <- ifelse(is.na(mapped1), ortho$Species1, mapped1)

    # Map Species2 through sp2 gene_map
    gm2 <- reductions[[sp2]]$gene_map
    idx2 <- match(ortho$Species2, gm2$original)
    mapped2 <- gm2$representative[idx2]
    ortho$Species2 <- ifelse(is.na(mapped2), ortho$Species2, mapped2)

    ortho
  })

  result <- do.call(rbind, result_list)
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame(Species1 = character(0),
                      Species2 = character(0),
                      hog = character(0)))
  }

  unique(result)
}
