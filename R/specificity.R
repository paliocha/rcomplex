#' Neighbourhood specificity of ortholog pairs
#'
#' For each ortholog pair and direction, the anchor's co-expression
#' neighbourhood is mapped through the ortholog table to the partner
#' species (orthologs of the anchor's own group removed) and every gene of
#' the partner species is scored by the AUROC of that mapped set against
#' the gene's co-expression column. The reported p-value is the rank of
#' the paired gene among all partner-species genes, on the `1 / n` grid.
#' Runs on the sparse representation; dense networks are converted with
#' every nonzero entry stored.
#'
#' @inheritParams compare_neighborhoods
#' @param directions `"both"` (default), `"1to2"` (anchors in species 1)
#'   or `"2to1"`.
#'
#' @return A data frame with `Species1`, `Species2`, `hog` and, per
#'   requested direction (`Species1.` for 1 to 2, `Species2.` for 2 to
#'   1): `neigh` (anchor degree), `mapped` (size of the mapped set),
#'   `auroc`, `p.val`, `jaccard` (as in [compare_neighborhoods()]) and
#'   `effect.size` (equal to `auroc`). `auroc` and `p.val` are `NA` when
#'   the mapped set is empty or spans every other partner gene.
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
  slots <- extract_sparse_cpp(m, min(net$threshold, m[m != 0]), n_cores)
  net$network <- methods::new(
    "dgCMatrix",
    i = slots$i, p = slots$p, x = slots$x,
    Dim = dim(m), Dimnames = dimnames(m)
  )
  net
}
