#' Compare co-expression neighborhoods across species
#'
#' For each ortholog pair, tests whether co-expression neighborhoods are
#' conserved between species using hypergeometric tests in both directions.
#'
#' @section Ortholog preparation for de novo transcriptomes:
#' When working with de novo transcriptome assemblies, it is **critical** to
#' filter ortholog tables to expressed genes **before** applying any
#' maximum-paralogs filter. De novo assemblies can produce 10--50x more gene
#' models than are actually expressed (e.g. 205K models but only 14K expressed).
#' Unexpressed fragments inflate paralog counts per ortholog group, causing the
#' paralogs filter to discard groups that would pass if only expressed members
#' were counted. Failing to do this can reduce the usable gene set by 40--70%
#' and eliminate all statistical power from the hypergeometric test.
#'
#' Recommended workflow:
#' 1. Build expression matrices for both species
#' 2. Filter ortholog table to genes present in expression data
#' 3. Apply maximum-paralogs filter on the reduced ortholog table
#' 4. Compute networks and run comparison
#'
#' @param net1 Network object for species 1 (output of [compute_network()]).
#' @param net2 Network object for species 2 (output of [compute_network()]).
#' @param orthologs Data frame with columns `Species1`, `Species2`, and
#'   `hog` (output of [parse_orthologs()]).
#' @param n_cores Number of threads for parallel computation (default 1).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Species1}{Gene identifier for species 1}
#'     \item{Species2}{Gene identifier for species 2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{Species1.neigh}{Number of neighbors of Species1 gene in net1}
#'     \item{Species1.ortho.neigh}{Number of ortholog-mapped neighbors
#'       from net2}
#'     \item{Species1.neigh.overlap}{Intersection size}
#'     \item{Species1.p.val.con}{Upper-tail hypergeometric p-value for
#'       conservation (direction 1)}
#'     \item{Species1.p.val.div}{Lower-tail hypergeometric p-value for
#'       divergence (direction 1)}
#'     \item{Species1.effect.size}{Fold enrichment (direction 1). Values > 1
#'       indicate conservation, < 1 indicate divergence.}
#'     \item{Species2.neigh}{Number of neighbors of Species2 gene in net2}
#'     \item{Species2.ortho.neigh}{Number of ortholog-mapped neighbors
#'       from net1}
#'     \item{Species2.neigh.overlap}{Intersection size}
#'     \item{Species2.p.val.con}{Upper-tail hypergeometric p-value for
#'       conservation (direction 2)}
#'     \item{Species2.p.val.div}{Lower-tail hypergeometric p-value for
#'       divergence (direction 2)}
#'     \item{Species2.effect.size}{Fold enrichment (direction 2). Values > 1
#'       indicate conservation, < 1 indicate divergence.}
#'   }
#'
#' @examples
#' \dontrun{
#' comparison <- compare_neighborhoods(net_A, net_B, orthologs)
#' head(comparison[, c("Species1", "Species2", "hog",
#'                      "Species1.effect.size")])
#' }
#'
#' @export
compare_neighborhoods <- function(net1, net2, orthologs, n_cores = 1L) {
  # Validate inputs
  if (!is.list(net1) || is.null(net1$network)) {
    stop("net1 must be a network object from compute_network()")
  }
  if (!is.list(net2) || is.null(net2$network)) {
    stop("net2 must be a network object from compute_network()")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }

  net1_mat <- net1$network
  net2_mat <- net2$network
  thr1 <- net1$threshold
  thr2 <- net2$threshold

  net1_genes <- rownames(net1_mat)
  net2_genes <- rownames(net2_mat)

  # Filter orthologs to genes present in both networks and deduplicate
  orthologs <- orthologs[orthologs$Species1 %in% net1_genes &
                           orthologs$Species2 %in% net2_genes, ,
                         drop = FALSE]
  orthologs <- unique(orthologs[, c("Species1", "Species2", "hog"),
                                drop = FALSE])

  if (nrow(orthologs) == 0) {
    stop("No orthologs found in both networks")
  }

  # Build gene name -> 0-based index maps
  idx1 <- stats::setNames(seq_along(net1_genes) - 1L, net1_genes)
  idx2 <- stats::setNames(seq_along(net2_genes) - 1L, net2_genes)

  sp1_idx <- as.integer(idx1[orthologs$Species1])
  sp2_idx <- as.integer(idx2[orthologs$Species2])

  # Call C++
  result <- compare_neighborhoods_cpp(
    net1 = net1_mat,
    net2 = net2_mat,
    thr1 = thr1,
    thr2 = thr2,
    pair_sp1_idx = sp1_idx,
    pair_sp2_idx = sp2_idx,
    ortho_sp1_idx = sp1_idx,
    ortho_sp2_idx = sp2_idx,
    n_cores = n_cores
  )

  # Combine with ortholog info
  cbind(
    orthologs[, c("Species1", "Species2", "hog"), drop = FALSE],
    result
  )
}


#' Build clique edges from pairwise comparison results
#'
#' Converts output from \code{\link{summarize_comparison}} into the edge
#' format expected by \code{\link{find_cliques}} and
#' \code{\link{clique_persistence}}. Renames columns, injects species
#' identity, and computes per-pair effect sizes and classification.
#'
#' @param comparison Data frame from \code{\link{summarize_comparison}}
#'   (the \code{$results} element). Must contain columns \code{Species1},
#'   \code{Species2}, \code{hog}, plus q-value and effect-size columns
#'   from both directions.
#' @param sp1 Species abbreviation for \code{Species1} genes (e.g.,
#'   \code{"BDIS"}).
#' @param sp2 Species abbreviation for \code{Species2} genes.
#' @param alternative Which test direction to use for q-values and
#'   classification: \code{"greater"} (conservation, default) or
#'   \code{"less"} (divergence).
#' @param alpha Significance threshold for the \code{type} column
#'   (default 0.05).
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (Species1)}
#'     \item{gene2}{Gene identifier (Species2)}
#'     \item{species1}{Species abbreviation for gene1 (\code{sp1})}
#'     \item{species2}{Species abbreviation for gene2 (\code{sp2})}
#'     \item{hog}{Ortholog group identifier}
#'     \item{q.value}{Minimum of the two directional q-values}
#'     \item{effect_size}{Geometric mean of directional effect sizes}
#'     \item{type}{\code{"conserved"} or \code{"diverged"} if
#'       \code{q.value < alpha}; \code{"ns"} otherwise}
#'   }
#'
#' @examples
#' \dontrun{
#' summary_AB <- summarize_comparison(comparison_AB)
#' edges_AB <- comparison_to_edges(summary_AB$results, "SP_A", "SP_B")
#'
#' # Combine multiple species pairs for find_cliques()
#' edges <- rbind(edges_AB, edges_AC, edges_BC)
#' }
#'
#' @export
comparison_to_edges <- function(comparison, sp1, sp2,
                                alternative = c("greater", "less"),
                                alpha = 0.05) {
  alternative <- match.arg(alternative)

  suffix <- if (alternative == "greater") "con" else "div"
  q1_col <- paste0("Species1.q.val.", suffix)
  q2_col <- paste0("Species2.q.val.", suffix)

  required <- c("Species1", "Species2", "hog",
                 "Species1.effect.size", "Species2.effect.size",
                 q1_col, q2_col)
  missing_cols <- setdiff(required, names(comparison))
  if (length(missing_cols) > 0) {
    stop("comparison missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Did you pass summarize_comparison()$results?")
  }

  q_min <- pmin(comparison[[q1_col]], comparison[[q2_col]], na.rm = TRUE)
  q_min[is.infinite(q_min)] <- NA_real_
  eff_geo <- sqrt(comparison$Species1.effect.size *
                  comparison$Species2.effect.size)

  type_label <- if (alternative == "greater") "conserved" else "diverged"
  type <- ifelse(q_min < alpha, type_label, "ns")

  data.frame(
    gene1 = comparison$Species1,
    gene2 = comparison$Species2,
    species1 = sp1,
    species2 = sp2,
    hog = comparison$hog,
    q.value = q_min,
    effect_size = eff_geo,
    type = type
  )
}
