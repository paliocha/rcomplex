#' Read an ortholog group file
#'
#' Reads an ortholog group file into the long table. The table has one
#' row per gene: `species`, `gene` and `hog`.
#'
#' @param file Path to a tab-delimited ortholog file. It can be gzipped.
#' @param species Character vector of species to keep, or `NULL` for all
#'   species in the file.
#' @param format One of `"auto"`, `"orthofinder"`, `"plaza"` or `"long"`.
#'   `"auto"` reads the format from the header.
#'
#' @return A data frame with columns `species`, `gene` and `hog`.
#'
#' @details
#' Three formats are read:
#' \describe{
#'   \item{orthofinder}{OrthoFinder `N0.tsv` (HOGs, column `HOG`) or
#'     `Orthogroups.tsv` (column `Orthogroup`). There is one column per
#'     species, with genes separated by commas. The columns `OG` and
#'     `Gene Tree Parent Clade` are skipped.}
#'   \item{plaza}{PLAZA ortholog groups: columns `species`, `gene_id` and
#'     `gene_content`, where `gene_content` lists the members as
#'     `code:gene1,gene2;code:gene3`. Each row is the group of its anchor
#'     gene. Rows with the same `gene_content` share one integer `hog`,
#'     numbered by the sorted `gene_content` strings of the whole file.}
#'   \item{long}{The output shape itself: columns `species`, `gene` and
#'     `hog`. Other columns are dropped.}
#' }
#' `"auto"` picks `orthofinder` when the header has `HOG` or `Orthogroup`,
#' `plaza` when it has `gene_content`, and `long` when it has `species`,
#' `gene` and `hog`. FastOMA and eggNOG output have no reader: write them
#' as the long table, one row per gene and group.
#'
#' @examples
#' f <- system.file("extdata", "N0.tsv", package = "rcomplex")
#' long <- read_orthologs(f, species = c("SpA", "SpB"))
#' head(long)
#' prepare_orthologs(long)
#'
#' @export
read_orthologs <- function(file, species = NULL,
                           format = c(
                             "auto", "orthofinder", "plaza",
                             "long"
                           )) {
  format <- match.arg(format)
  if (!file.exists(file)) stop("ortholog file not found: ", file)
  dt <- data.table::fread(file,
    sep = "\t", header = TRUE,
    showProgress = FALSE, data.table = FALSE
  )
  cols <- names(dt)
  is_of <- any(c("HOG", "Orthogroup") %in% cols)
  is_long <- all(c("species", "gene", "hog") %in% cols)
  if (format == "auto") {
    format <- if (is_of) {
      "orthofinder"
    } else if ("gene_content" %in% cols) {
      "plaza"
    } else if (is_long) {
      "long"
    } else {
      stop("cannot read the format from the header; set format")
    }
  }
  out <- switch(format,
    orthofinder = {
      if (!is_of) stop("an OrthoFinder file needs a HOG or Orthogroup column")
      .orthofinder_long(dt)
    },
    plaza = {
      if (!all(c("species", "gene_id", "gene_content") %in% cols)) {
        stop("a PLAZA file needs columns species, gene_id, gene_content")
      }
      .plaza_long(dt)
    },
    long = {
      if (!is_long) stop("a long file needs columns species, gene, hog")
      data.frame(
        species = as.character(dt$species),
        gene = as.character(dt$gene), hog = dt$hog
      )
    }
  )
  if (!is.null(species)) {
    miss <- setdiff(species, out$species)
    if (length(miss) > 0L) {
      stop("species not in file: ", paste(miss, collapse = ", "))
    }
    out <- out[out$species %in% species, , drop = FALSE]
  }
  out <- unique(out)
  rownames(out) <- NULL
  out
}


#' Character cells of a column, with NA as ""
#' @noRd
.cells <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}


#' OrthoFinder N0.tsv / Orthogroups.tsv to the long table
#' @noRd
.orthofinder_long <- function(dt) {
  hog <- dt[[if ("HOG" %in% names(dt)) "HOG" else "Orthogroup"]]
  sp <- setdiff(
    names(dt), c("HOG", "OG", "Gene Tree Parent Clade", "Orthogroup")
  )
  parts <- lapply(sp, function(s) {
    g <- strsplit(.cells(dt[[s]]), ",", fixed = TRUE)
    gene <- trimws(unlist(g, use.names = FALSE))
    h <- rep(hog, lengths(g))
    keep <- gene != ""
    data.frame(species = rep(s, sum(keep)), gene = gene[keep], hog = h[keep])
  })
  do.call(rbind, parts)
}


#' PLAZA species / gene_id / gene_content to the long table
#'
#' The anchor gene and every member of its `gene_content` share the hog
#' of that `gene_content` string (sorted-key numbering, as in the
#' pre-0.4.0 `parse_orthologs()`).
#' @noRd
.plaza_long <- function(dt) {
  content <- .cells(dt$gene_content)
  hog <- match(content, sort(unique(content)))
  chunks <- strsplit(content, ";", fixed = TRUE)
  chunk <- unlist(chunks, use.names = FALSE)
  genes <- strsplit(sub("^[^:]*:", "", chunk), ",", fixed = TRUE)
  n <- lengths(genes)
  data.frame(
    species = c(as.character(dt$species), rep(sub(":.*$", "", chunk), n)),
    gene = c(as.character(dt$gene_id), unlist(genes, use.names = FALSE)),
    hog = c(hog, rep(rep(hog, lengths(chunks)), n))
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
#' @param orthologs Data frame with columns \code{gene1} (or the column
#'   matching gene row names), \code{gene2}, and \code{hog}, as
#'   returned by \code{\link{prepare_orthologs}}.
#' @param gene_col Character: which column of \code{orthologs} contains gene
#'   IDs matching row names of \code{expr_matrix} (default \code{"gene1"}).
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
#' reduced$expr_matrix # reduced expression matrix
#' reduced$gene_map # original -> representative mapping
#' }
#'
#' @export
reduce_orthogroups <- function(expr_matrix, orthologs,
                               gene_col = "gene1",
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


#' Pair orthologs for every species pair
#'
#' Pairs the genes of each HOG for every pair of species in the long
#' table. Genes can first be renamed through [reduce_orthogroups()]
#' output.
#'
#' @details
#' Species pairs follow the order in which species first appear in
#' `orthologs`. Rows with `NA` hog are skipped. To build the long table
#' from `SummarizedExperiment` objects, bind one
#' `data.frame(species = , gene = rownames(se), hog = rowData(se)$hog)`
#' per species.
#'
#' Paralog reduction is a lossy step: correlated paralogs are averaged into
#' one representative gene, so per-paralog identity is gone by the time
#' the pairs reach downstream consumers. That trade-off pays for
#' itself for consumers that need one counterpart per gene (module
#' preservation, the species-graph clique backend), but the gene-graph
#' clique backend (\code{gene_clique_graph} /
#' \code{classify_gene_cliques}) is built to resolve which paralog
#' copy is conserved, and needs every original paralog as its own node to do
#' that. Pass \code{reductions = NULL} (the default) to skip paralog
#' reduction entirely and keep every gene at its original identity.
#'
#' @param orthologs The long table from [read_orthologs()]: columns
#'   `species`, `gene` and `hog`.
#' @param reductions Named list of [reduce_orthogroups()] outputs, one per
#'   species in `orthologs`, or `NULL` to keep every gene. Each element
#'   needs a `gene_map` data frame with columns `original` and
#'   `representative`.
#'
#' @return A data frame with columns `gene1`, `gene2` and `hog`, one row
#'   per gene pair in a shared HOG. With `reductions`, gene names are the
#'   representatives and duplicate rows are removed.
#'
#' @examples
#' f <- system.file("extdata", "N0.tsv", package = "rcomplex")
#' ortho <- prepare_orthologs(read_orthologs(f))
#' head(ortho)
#'
#' @export
prepare_orthologs <- function(orthologs, reductions = NULL) {
  if (!is.data.frame(orthologs) ||
        !all(c("species", "gene", "hog") %in% names(orthologs))) {
    stop("orthologs must be a data frame with columns species, gene, hog")
  }
  sp_names <- unique(as.character(orthologs$species))
  if (length(sp_names) < 2) {
    stop("orthologs must contain at least two species")
  }

  if (!is.null(reductions)) {
    if (!is.list(reductions) || is.null(names(reductions))) {
      stop("reductions must be a named list")
    }
    missing_sp <- setdiff(sp_names, names(reductions))
    if (length(missing_sp) > 0) {
      stop(
        "reductions missing species present in orthologs: ",
        paste(missing_sp, collapse = ", ")
      )
    }
    for (sp in sp_names) {
      if (is.null(reductions[[sp]]$gene_map)) {
        stop("reductions[['", sp, "']] must have a $gene_map element")
      }
    }
  }

  pairs <- utils::combn(sp_names, 2, simplify = FALSE)

  result_list <- lapply(pairs, function(pair) {
    species1 <- pair[1]
    species2 <- pair[2]

    m <- merge(
      orthologs[orthologs$species == species1, c("hog", "gene")],
      orthologs[orthologs$species == species2, c("hog", "gene")],
      by = "hog", incomparables = NA
    )
    ortho <- data.frame(gene1 = m$gene.x, gene2 = m$gene.y, hog = m$hog)
    if (nrow(ortho) == 0 || is.null(reductions)) {
      return(ortho)
    }

    # Map gene1 through species1 gene_map
    gm1 <- reductions[[species1]]$gene_map
    idx1 <- match(ortho$gene1, gm1$original)
    mapped1 <- gm1$representative[idx1]
    ortho$gene1 <- ifelse(is.na(mapped1), ortho$gene1, mapped1)

    # Map gene2 through species2 gene_map
    gm2 <- reductions[[species2]]$gene_map
    idx2 <- match(ortho$gene2, gm2$original)
    mapped2 <- gm2$representative[idx2]
    ortho$gene2 <- ifelse(is.na(mapped2), ortho$gene2, mapped2)

    ortho
  })

  result <- do.call(rbind, result_list)
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame(
      gene1 = character(0),
      gene2 = character(0),
      hog = character(0)
    ))
  }

  unique(result)
}
