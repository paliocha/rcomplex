#' Gains and losses of each clique on the species tree
#'
#' Reads clique presence over species as a binary character. Fitch
#' parsimony on a rooted species tree places each gain and loss on a
#' branch. A second table sums the changes per branch.
#'
#' @details
#' A member species has state 1. A species with reason `tested_ns` has
#' state 0: the pair was tested and failed. Every other reason is unknown
#' (`?`), and so is a tip outside `species`. An unknown tip never counts
#' as a change. An `absent` species is unknown because gene presence is
#' its own character, and a gene loss is not a co-expression loss.
#'
#' The pass is unit-cost Fitch, run once per distinct pattern of missing
#' species. The tree may hold polytomies. A clique restricted to one clade
#' costs the same as one gain on its stem or one loss outside it. Such a
#' root tie resolves to absent, so the clique reads as one gain, as
#' `lineage_specific` does. `root_tied` flags the tie. Below the root, a
#' tie takes the change as early as possible (ACCTRAN). `n_mpr` counts
#' the equally parsimonious reconstructions.
#'
#' @param cliques A [classify_gene_cliques()] table. It needs `clique_id`,
#'   `hog`, `missing_species` and `missing_reason`. `run` is kept.
#' @param tree A rooted `phylo` from the ape package, with tips named by
#'   species.
#' @param species Every species of the analysis. Each must be a tip.
#' @return A list of two data frames. `cliques` has one row per input row,
#'   in input order: `n_present`, `n_absent` (state 0), `n_unknown`,
#'   `state_root`, `root_tied`, `n_gain`, `n_loss`, `n_mpr`, and the
#'   branches in `gains` and `losses`. `branches` has one row per branch,
#'   named by the node below it as in [clades_from_tree()]: `n_species`,
#'   and `n_gain` and `n_loss` summed over the cliques.
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   tree <- ape::read.tree(text = "((A,B),(C,D));")
#'   cl <- data.frame(
#'     clique_id = c("1", "2"), hog = "H1",
#'     missing_species = c("C,D", "D"),
#'     missing_reason = c("tested_ns,tested_ns", "tested_ns")
#'   )
#'   edge_history(cl, tree, c("A", "B", "C", "D"))
#' }
#' @export
edge_history <- function(cliques, tree, species) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Install the ape package to read a tree.")
  }
  need <- c("clique_id", "hog", "missing_species", "missing_reason")
  gone <- setdiff(need, names(cliques))
  if (length(gone) > 0L) stop("cliques needs the column ", gone[1L], ".")
  if (!inherits(tree, "phylo")) stop("tree must be an ape phylo object.")
  if (!ape::is.rooted(tree)) stop("tree is unrooted. Root it with ape::root().")
  tips <- tree$tip.label
  off <- setdiff(species, tips)
  if (length(off) > 0L) stop("species ", off[1L], " is not a tip of tree.")

  n_tip <- length(tips)
  n_node <- n_tip + tree$Nnode
  edge <- tree$edge[ape::postorder(tree), , drop = FALSE]
  base <- ifelse(tips %in% species, 1L, NA_integer_)
  key <- paste(cliques$missing_species, cliques$missing_reason, sep = "|")
  u <- unique(key)
  g <- match(key, u)
  res <- lapply(match(u, key), function(i) {
    ms <- strsplit(cliques$missing_species[i], ",", fixed = TRUE)[[1L]]
    mr <- strsplit(cliques$missing_reason[i], ",", fixed = TRUE)[[1L]]
    bad <- setdiff(ms, species)
    if (length(bad) > 0L) {
      stop("missing_species holds ", bad[1L], ", which is not in species.")
    }
    st <- base
    st[match(ms, tips)] <- ifelse(mr == "tested_ns", 0L, NA_integer_)
    .fitch_binary(st, edge, n_tip, n_node)
  })
  pick <- function(f, v) vapply(res, `[[`, v, f)
  chg <- matrix(pick("chg", integer(n_node)), nrow = n_node)
  st <- matrix(pick("tip", integer(n_tip)), nrow = n_tip)
  nm <- .node_names(tree)
  lst <- function(k) {
    vapply(
      seq_along(u), function(j) paste(nm[chg[, j] == k], collapse = ","),
      ""
    )
  }
  n <- function(x) as.integer(colSums(x, na.rm = TRUE))[g]
  out <- data.frame(
    clique_id = cliques$clique_id, hog = cliques$hog,
    n_present = n(st == 1L), n_absent = n(st == 0L),
    n_unknown = n(is.na(st)),
    state_root = pick("root", integer(1))[g],
    root_tied = pick("tied", logical(1))[g],
    n_gain = n(chg == 1L), n_loss = n(chg == -1L),
    n_mpr = pick("mpr", numeric(1))[g],
    gains = lst(1L)[g], losses = lst(-1L)[g], stringsAsFactors = FALSE
  )
  if ("run" %in% names(cliques)) out <- cbind(run = cliques$run, out)
  w <- tabulate(g, length(u))
  ns <- c(rep(1L, n_tip), lengths(ape::prop.part(tree)))
  keep <- seq_len(n_node)[-(n_tip + 1L)]
  branches <- data.frame(
    branch = nm[keep], n_species = ns[keep],
    n_gain = as.integer((chg == 1L) %*% w)[keep],
    n_loss = as.integer((chg == -1L) %*% w)[keep],
    stringsAsFactors = FALSE
  )
  list(cliques = out, branches = branches)
}

#' Unit-cost Fitch on one binary character, as a Sankoff pass
#'
#' `st` holds 0, 1 or `NA` per tip; `edge` is in post-order. Returns the
#' tip states, the root state, the root tie, the count of equally
#' parsimonious reconstructions, and the change per node (1 gain, -1 loss).
#' @noRd
.fitch_binary <- function(st, edge, n_tip, n_node) {
  cost <- matrix(0, n_node, 2L)
  cnt <- matrix(1, n_node, 2L)
  known <- which(!is.na(st))
  cost[cbind(known, 2L - st[known])] <- Inf
  flip <- matrix(c(0, 1, 1, 0), 2L)
  for (k in seq_len(nrow(edge))) {
    p <- edge[k, 1L]
    ch <- edge[k, 2L]
    for (s in 1:2) {
      tot <- cost[ch, ] + flip[, s]
      m <- min(tot)
      cost[p, s] <- cost[p, s] + m
      cnt[p, s] <- cnt[p, s] * sum(cnt[ch, tot == m])
    }
  }
  root <- n_tip + 1L
  rc <- cost[root, ]
  state <- integer(n_node)
  state[root] <- which.min(rc)
  for (k in rev(seq_len(nrow(edge)))) {
    p <- state[edge[k, 1L]]
    tot <- cost[edge[k, 2L], ] + flip[, p]
    state[edge[k, 2L]] <- if (tot[1L] == tot[2L]) 3L - p else which.min(tot)
  }
  chg <- integer(n_node)
  chg[edge[, 2L]] <- state[edge[, 2L]] - state[edge[, 1L]]
  list(
    tip = st, root = state[root] - 1L, tied = rc[1L] == rc[2L],
    mpr = sum(cnt[root, rc == min(rc)]), chg = chg
  )
}
