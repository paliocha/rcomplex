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
#'   `OrthoGroup` (output of [parse_orthologs()]).
#' @param n_cores Number of threads for parallel computation (default 1).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Species1}{Gene identifier for species 1}
#'     \item{Species2}{Gene identifier for species 2}
#'     \item{OrthoGroup}{Ortholog group identifier}
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
#' @export
compare_neighborhoods <- function(net1, net2, orthologs, n_cores = 1L) {
  n_cores <- as.integer(n_cores)

  # Validate inputs
  if (!is.list(net1) || is.null(net1$network)) {
    stop("net1 must be a network object from compute_network()")
  }
  if (!is.list(net2) || is.null(net2$network)) {
    stop("net2 must be a network object from compute_network()")
  }
  if (!all(c("Species1", "Species2", "OrthoGroup") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, OrthoGroup")
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
  orthologs <- unique(orthologs[, c("Species1", "Species2", "OrthoGroup"),
                                drop = FALSE])

  if (nrow(orthologs) == 0) {
    stop("No orthologs found in both networks")
  }

  # Build gene name -> 0-based index maps
  idx1 <- stats::setNames(seq_along(net1_genes) - 1L, net1_genes)
  idx2 <- stats::setNames(seq_along(net2_genes) - 1L, net2_genes)

  # Map ortholog pair genes to indices
  pair_sp1_idx <- idx1[orthologs$Species1]
  pair_sp2_idx <- idx2[orthologs$Species2]

  # Map full ortholog table to indices
  ortho_sp1_idx <- idx1[orthologs$Species1]
  ortho_sp2_idx <- idx2[orthologs$Species2]

  # Call C++
  result <- compare_neighborhoods_cpp(
    net1 = net1_mat,
    net2 = net2_mat,
    thr1 = thr1,
    thr2 = thr2,
    pair_sp1_idx = as.integer(pair_sp1_idx),
    pair_sp2_idx = as.integer(pair_sp2_idx),
    ortho_sp1_idx = as.integer(ortho_sp1_idx),
    ortho_sp2_idx = as.integer(ortho_sp2_idx),
    n_cores = n_cores
  )

  # Combine with ortholog info
  cbind(
    orthologs[, c("Species1", "Species2", "OrthoGroup"), drop = FALSE],
    result
  )
}
