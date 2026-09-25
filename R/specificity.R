#' Neighbourhood specificity of ortholog pairs
#'
#' @description
#' Pair-level co-expressolog test that ranks each ortholog pair against
#' every other gene of the partner species. The anchor's co-expression
#' partners are translated to the partner species, and the paired
#' ortholog is scored by how well its own co-expression ranking
#' recognises that translated list compared with how well every other
#' partner-species gene recognises it. It follows the co-expression
#' conservation score of Suresh et al. (2023). It is the engine of
#' `method = "rank"` in [find_coexpressologs()] and
#' [density_sweep()]; [summarize_specificity()] turns its p-values into
#' q-values.
#'
#' @details
#' For anchor gene \eqn{i} in species 1 and direction 1 to 2:
#' \enumerate{
#'   \item \eqn{T} is the set of species-2 genes orthologous, through
#'     any copy, to a network neighbour of \eqn{i} (an entry at or above
#'     the analysis threshold), minus every species-2 gene in \eqn{i}'s
#'     own ortholog group. Removing the own group keeps a pair from
#'     scoring itself through its paralogs.
#'   \item For every species-2 gene \eqn{j}, the AUROC of
#'     \eqn{T} (without \eqn{j}) among the other \eqn{n_2 - 1} genes,
#'     ranked by their co-expression with \eqn{j}. Stored entries of
#'     column \eqn{j} are ranked exactly; unstored entries tie at the
#'     bottom with their mid-rank. A sparse network therefore ranks only
#'     its stored top entries, so `store_density` sets the resolution.
#'   \item The p-value of the ortholog pair \eqn{(i, o)} is
#'     \eqn{(1 + g) / n_2}, where \eqn{g} counts the genes
#'     \eqn{j \ne o} whose AUROC is at least that of \eqn{o}: the rank
#'     of the ortholog among all species-2 genes. A translated list that
#'     thousands of genes recognise as well as \eqn{o} gives a large
#'     p-value however high its AUROC.
#' }
#' Direction 2 to 1 swaps the roles. All p-values of one direction lie on
#' the same \eqn{1 / n_2} grid. They are not calibrated on their own:
#' pass them with the p-values of a comparison against [null_network()]
#' to [summarize_specificity()].
#'
#' Runs on the sparse representation; dense networks are converted with
#' every nonzero entry stored.
#'
#' @inheritParams compare_neighborhoods
#' @param directions `"both"` (default), `"1to2"` (anchors in species 1)
#'   or `"2to1"`.
#'
#' @return A data frame with `Species1`, `Species2`, `hog` and, per
#'   requested direction (`Species1.` for 1 to 2, `Species2.` for 2 to
#'   1):
#'   \describe{
#'     \item{neigh}{Anchor degree.}
#'     \item{mapped}{Size of the translated set \eqn{T}.}
#'     \item{auroc}{AUROC of \eqn{T} in the paired ortholog's ranking.}
#'     \item{p.val}{Rank p-value of the paired ortholog, as above.}
#'     \item{jaccard}{As in [compare_neighborhoods()].}
#'     \item{effect.size}{Equal to `auroc`.}
#'   }
#'   `auroc` and `p.val` are `NA` when \eqn{T} is empty or spans every
#'   other partner gene.
#'
#' @references
#' Suresh, H., Crow, M., Jorstad, N., Hodge, R., Lein, E., Dobin, A.,
#' Bakken, T. & Gillis, J. (2023). Comparative single-cell transcriptomic
#' analysis of primate brains highlights human-specific regulatory
#' evolution. \emph{Nature Ecology & Evolution}, 7(11), 1930--1943.
#' \doi{10.1038/s41559-023-02186-7}
#'
#' @export
compare_specificity <- function(net1, net2, orthologs, n_cores = 1L,
                                directions = c("both", "1to2", "2to1")) {
  directions <- match.arg(directions)
  op <- .ortholog_pair_index(net1, net2, orthologs)
  net1 <- .net_as_sparse(net1, n_cores)
  net2 <- .net_as_sparse(net2, n_cores)
  a1 <- .net_cpp_args(net1, net1$threshold)
  a2 <- .net_cpp_args(net2, net2$threshold)
  do_12 <- directions != "2to1"
  do_21 <- directions != "1to2"

  res <- specificity_sparse_cpp(
    p1 = a1$p, i1 = a1$i, x1 = a1$x, thr1 = a1$thr,
    p2 = a2$p, i2 = a2$i, x2 = a2$x, thr2 = a2$thr,
    pair_sp1_idx = op$sp1_idx, pair_sp2_idx = op$sp2_idx,
    ortho_sp1_idx = op$sp1_idx, ortho_sp2_idx = op$sp2_idx,
    do_12 = do_12, do_21 = do_21, n_cores = n_cores
  )
  for (s in c("Species1", "Species2")[c(do_12, do_21)]) {
    res[[paste0(s, ".effect.size")]] <- res[[paste0(s, ".auroc")]]
  }
  cbind(op$orthologs, res)
}


#' Sparse copy of a network object with every nonzero entry stored
#'
#' Dense networks are repacked as a `dgCMatrix` at the smaller of the
#' analysis threshold and the smallest nonzero value, so the unstored
#' entries are exactly the zeros (or nothing) and tie at the bottom of
#' every column, as the specificity ranks assume. Sparse networks pass
#' through.
#' @noRd
.net_as_sparse <- function(net, n_cores = 1L) {
  if (.net_is_sparse(net)) {
    return(net)
  }
  m <- .net_check(net, net$threshold)
  # every nonzero is stored, so the store threshold is exact even when
  # the analysis threshold sits below the smallest nonzero entry
  net$store_threshold <- min(net$threshold, m[m != 0])
  slots <- extract_sparse_cpp(m, net$store_threshold, n_cores)
  net$network <- methods::new(
    "dgCMatrix",
    i = slots$i, p = slots$p, x = slots$x,
    Dim = dim(m), Dimnames = dimnames(m)
  )
  net
}


#' Argument checks shared by the specificity path of find_coexpressologs()
#' and density_sweep()
#' @noRd
.check_specificity_args <- function(method, alternative, null_networks) {
  if (method != "rank" && !is.null(null_networks)) {
    stop("null_networks is only used with method = \"rank\"")
  }
  if (method == "rank") {
    if (is.null(null_networks)) {
      stop(
        "method = \"rank\" needs null_networks; build one per ",
        "species with null_network()"
      )
    }
    if (alternative == "less") {
      stop("method = \"rank\" supports alternative = \"greater\" only")
    }
  }
}


#' Specificity edges for one species pair against pooled null draws
#'
#' Direction 1 -> 2 is calibrated by comparing `net1` with each null of
#' species 2, direction 2 -> 1 by each null of species 1 against `net2`.
#' @noRd
.specificity_pair_edges <- function(net1, net2, nulls1, nulls2, orthologs,
                                    sp1, sp2, alpha, n_cores, pi0_method,
                                    pval_combine) {
  cmp <- compare_specificity(net1, net2, orthologs, n_cores)
  null_p <- list(
    sp1 = unlist(lapply(nulls2, function(nb) {
      compare_specificity(net1, nb, orthologs, n_cores,
        directions = "1to2"
      )$Species1.p.val
    })),
    sp2 = unlist(lapply(nulls1, function(na) {
      compare_specificity(na, net2, orthologs, n_cores,
        directions = "2to1"
      )$Species2.p.val
    }))
  )
  # the specificity p-values are continuous-ish ranks, so the randomized
  # pi0 of the hypergeometric path reduces to plain Storey
  if (pi0_method == "randomized") pi0_method <- "storey"
  summarize_specificity(cmp, null_p, alpha,
    pi0_method = pi0_method, sp1 = sp1, sp2 = sp2,
    pval_combine = pval_combine
  )$edges
}


#' Accept the pre-0.3.1 name of the hypergeometric arm
#' @noRd
.method_alias <- function(method) {
  if (identical(method, "analytical")) "hypergeometric" else method
}
