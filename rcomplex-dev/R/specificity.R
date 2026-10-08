#' Neighbourhood specificity of ortholog pairs
#'
#' @description
#' Ranks each ortholog pair against every other gene of the partner
#' species. This is the engine of `method = "rank"` in
#' [find_coexpressologs()] and [density_sweep()].
#' [summarize_specificity()] turns its p-values into q-values.
#'
#' @details
#' The anchor's co-expression partners are translated to the partner
#' species. The paired ortholog is scored by how well its own
#' co-expression ranking recognises that translated list, compared with how
#' well every other partner-species gene recognises it. The score follows
#' the co-expression conservation score of Suresh et al. (2023).
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
#' @return A data frame with `gene1`, `gene2`, `hog` and, per
#'   requested direction (`species1.` for 1 to 2, `species2.` for 2 to
#'   1):
#'   \describe{
#'     \item{neigh}{Anchor degree.}
#'     \item{mapped}{Size of the translated set \eqn{T}.}
#'     \item{auroc}{AUROC of \eqn{T} in the paired ortholog's ranking.}
#'     \item{p_value}{Rank p-value of the paired ortholog, as above.}
#'     \item{jaccard}{As in [compare_neighborhoods()].}
#'     \item{n.cand}{Number of candidate genes in the partner network.}
#'     \item{effect_size}{Equal to `auroc`.}
#'     \item{auroc.grid}{Matrix, one row per pair and one column per
#'       raw-p fraction f (1e-5 to 1, the column names): the
#'       `ceiling(f * n.cand)`-th largest candidate AUROC for the anchor,
#'       i.e. the AUROC the ortholog needs to reach raw p of about f.
#'       [summarize_specificity()] reads it for the edge `power`. It is
#'       computed on every call; a comparison against a [null_network()]
#'       partner, which only needs `p_value`, pays for it too (the package's
#'       own null runs skip it).}
#'   }
#'   `auroc` and `p_value` are `NA` when \eqn{T} is empty or spans every
#'   other partner gene.
#'
#' @references
#' Suresh, H., Crow, M., Jorstad, N., Hodge, R., Lein, E., Dobin, A.,
#' Bakken, T. & Gillis, J. (2023). Comparative single-cell transcriptomic
#' analysis of primate brains highlights human-specific regulatory
#' evolution. \emph{Nature Ecology & Evolution}, 7(11), 1930--1943.
#' \doi{10.1038/s41559-023-02186-7}
#'
#' @keywords internal
compare_specificity <- function(net1, net2, orthologs, n_cores = 1L,
                                directions = c("both", "1to2", "2to1")) {
  .specificity_run(net1, net2, orthologs, n_cores, match.arg(directions),
    grid_frac = .rank_grid_frac
  )
}


# Raw-p fractions at which each anchor's AUROC grid is recorded; the power
# of a rank-test edge interpolates the call threshold on this grid. From
# 1e-5 (below 1 / n for up to 100,000 candidates) to 1, so neither the
# reference rank nor a call threshold is clamped in practice.
.rank_grid_frac <- c(
  1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05,
  0.1, 0.2, 0.3, 0.5, 0.7, 1
)


# compare_specificity() without the argument checks; the null runs pass
# grid_frac = numeric(0), since only the observed comparison needs a grid.
.specificity_run <- function(net1, net2, orthologs, n_cores, directions,
                             grid_frac) {
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
    do_12 = do_12, do_21 = do_21, n_cores = n_cores,
    grid_frac = as.numeric(grid_frac)
  )
  grid_cols <- grepl("\\.auroc\\.grid$", names(res))
  out <- cbind(op$orthologs, as.data.frame(res[!grid_cols]))
  for (s in c("species1", "species2")[c(do_12, do_21)]) {
    out[[paste0(s, ".effect_size")]] <- out[[paste0(s, ".auroc")]]
    if (length(grid_frac)) {
      g <- res[[paste0(s, ".auroc.grid")]]
      colnames(g) <- format(grid_frac,
        scientific = FALSE,
        drop0trailing = TRUE
      )
      out[[paste0(s, ".auroc.grid")]] <- g
    }
  }
  out
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


#' Argument checks for the specificity path of find_coexpressologs()
#' @noRd
.check_specificity_args <- function(method, null_networks, rho0 = NULL) {
  if (method != "rank" && !is.null(null_networks)) {
    stop("null_networks is only used with method = \"rank\"")
  }
  if (method == "rank" && !is.null(rho0)) {
    stop("rho0 is only used with method = \"hypergeometric\"")
  }
  if (method == "rank") {
    if (is.null(null_networks)) {
      stop(
        "method = \"rank\" needs null_networks; build one per ",
        "species with null_network()"
      )
    }
  }
}


#' Specificity edges for one species pair against pooled null draws
#'
#' Direction 1 -> 2 is calibrated by comparing `net1` with each null of
#' species 2, direction 2 -> 1 by each null of species 1 against `net2`.
#' @noRd
.specificity_pair_edges <- function(net1, net2, nulls1, nulls2, orthologs,
                                    species1, species2, n_cores,
                                    pval_combine) {
  cmp <- compare_specificity(net1, net2, orthologs, n_cores)
  null_p <- list(
    species1 = unlist(lapply(nulls2, function(nb) {
      .specificity_run(net1, nb, orthologs, n_cores, "1to2",
        grid_frac = numeric(0)
      )$species1.p_value
    })),
    species2 = unlist(lapply(nulls1, function(na) {
      .specificity_run(na, net2, orthologs, n_cores, "2to1",
        grid_frac = numeric(0)
      )$species2.p_value
    }))
  )
  summarize_specificity(cmp, null_p,
    species1 = species1, species2 = species2, pval_combine = pval_combine
  )$edges
}


# Detection power of each rank-test edge: the probability that the pair
# would have been called had the ortholog ranked at the reference raw p
# `p0` among the anchor's candidates (default: the median raw p of the
# called pairs, per direction). Per direction, the call threshold on the
# raw-p scale is the largest raw p called. The anchor's AUROC grid
# (log-linear interpolation in the raw-p fraction) turns both into AUROCs,
# G(p0) and G(p_cut), and the ortholog's AUROC over its t translated genes
# against n - 1 - t others is taken as normal around G(p0) with the
# Hanley-McNeil (1982) standard error. A fixed reference AUROC does not
# work here: anchors whose candidates all score high need a high AUROC and
# also have orthologs that reach one, so a fixed alternative inverted the
# power (validation in dev/design-notes/module-engine-redesign.md, 11.15).
# With p0 below the threshold the power is at least 0.5. The grid ranks
# include the ortholog itself, a shift of at most one rank. Directions
# combine like pval_combine, as in the hypergeometric .edge_power().
.rank_power <- function(res, alpha, pval_combine = c("max", "min"),
                        p0 = NULL) {
  pval_combine <- match.arg(pval_combine)
  .check_p0(p0)
  na_out <- rep(NA_real_, nrow(res))
  dirs <- c("species1", "species2")
  need <- as.vector(outer(dirs, c(
    ".p_value", ".q_value_con", ".mapped", ".n.cand", ".auroc.grid"
  ), paste0))
  if (nrow(res) == 0L) {
    return(na_out)
  }
  if (!all(need %in% names(res))) {
    .warn_rank_power_na(paste(setdiff(need, names(res)), collapse = ", "))
    return(na_out)
  }
  combine <- if (pval_combine == "min") pmin else pmax
  q_comb <- combine(res$species1.q_value_con, res$species2.q_value_con,
    na.rm = TRUE
  )
  called <- is.finite(q_comb) & q_comb < alpha
  pw <- lapply(dirs, function(d) {
    col <- function(s) res[[paste0(d, s)]]
    q <- col(".q_value_con")
    p <- col(".p_value")
    sig <- !is.na(q) & q < alpha
    grid <- col(".auroc.grid")
    lf <- suppressWarnings(log(as.numeric(colnames(grid))))
    # the fractions travel as column names; a frame that lost them cannot
    # be read, so it gets no power rather than a wrong one
    if (length(lf) < 2L || length(lf) != ncol(grid) || anyNA(lf) ||
          any(diff(lf) <= 0)) {
      .warn_rank_power_na(paste0(d, ".auroc.grid fraction names"))
      return(NULL) # a damaged direction voids the edge, whatever "min" does
    }
    # nothing called in this direction: nothing was detectable, so a miss
    # here is no evidence -- power 0, not NA (which the classifiers would
    # read as "every miss is a rejection")
    if (!any(sig)) {
      return(rep(0, nrow(res)))
    }
    f0 <- p0
    if (is.null(f0)) {
      # pairs called and significant in this direction: under "min" a pair
      # called through the other direction alone sits above this
      # direction's threshold and would pull the reference past it
      use <- called & sig & is.finite(p)
      # no pair called overall (e.g. under "max" the other direction never
      # calls): fall back on this direction's own calls
      if (!any(use)) use <- sig & is.finite(p)
      f0 <- stats::median(p[use])
    }
    at <- function(f) {
      x <- min(max(log(f), lf[1L]), lf[length(lf)])
      k <- min(findInterval(x, lf), length(lf) - 1L)
      w <- (x - lf[k]) / (lf[k + 1L] - lf[k])
      grid[, k] + w * (grid[, k + 1L] - grid[, k])
    }
    a_crit <- at(max(p[sig]))
    a_ref <- at(f0)
    a <- pmin(pmax(a_ref, 0.5 + 1e-6), 1 - 1e-6)
    t <- as.numeric(col(".mapped"))
    m <- as.numeric(col(".n.cand")) - 1 - t
    q1 <- a / (2 - a)
    q2 <- 2 * a^2 / (1 + a)
    se <- sqrt((a * (1 - a) + (t - 1) * (q1 - a^2) +
                  (m - 1) * (q2 - a^2)) / (t * m))
    out <- stats::pnorm((a_ref - a_crit) / se)
    out[t <= 0 | m <= 0] <- 0
    out
  })
  if (any(vapply(pw, is.null, logical(1)))) {
    return(na_out)
  }
  if (pval_combine == "max") {
    pmin(pw[[1L]], pw[[2L]])
  } else {
    pmax(pw[[1L]], pw[[2L]], na.rm = TRUE)
  }
}


.check_p0 <- function(p0) {
  ok <- is.null(p0) || (is.numeric(p0) && length(p0) == 1L &&
                          !is.na(p0) && p0 > 0 && p0 < 1)
  if (!ok) stop("p0 must be NULL or a single number in (0, 1)")
}


# Rank power that cannot be computed for a structural reason (a damaged
# frame), as opposed to "no pair called": say so, because NA power makes
# the clique classifiers read every miss as a rejection.
.warn_rank_power_na <- function(what) {
  warning(
    "rank-test power is NA: the frame is missing or has unreadable ", what,
    "; keep rank frames with saveRDS()",
    call. = FALSE
  )
}
