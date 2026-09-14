#' Exchangeability blocks for the trait-relabelling null (internal)
#'
#' The null for \code{\link{tag_permutation}} holds each contrast's pair of
#' trait labels fixed and varies only which species in the contrast carries
#' which. That is a constraint per contrast, not per species, so contrasts
#' sharing a species are coupled: relabelling one changes the other. The
#' unit of independence is therefore the connected component of the graph
#' whose nodes are species and whose edges are contrasts, and the null is
#' the product over components of each component's admissible labellings.
#'
#' A disjoint pairing --- every species in exactly one contrast --- gives
#' one component per contrast and two labellings each (swap or not), so the
#' familiar \code{2^k} space is the special case, not the assumption.
#'
#' @param pairs Data frame with \code{sp1}, \code{sp2}.
#' @param group Named trait vector covering every species in \code{pairs}.
#' @return List with \code{membership} (component index per species, named),
#'   \code{labellings} (list, one entry per component, each a list of named
#'   trait vectors over that component's species) and \code{n} (labellings
#'   per component).
#' @noRd
.tp_blocks <- function(pairs, group) {
  species <- unique(c(pairs$sp1, pairs$sp2))
  g <- igraph::graph_from_data_frame(
    pairs[, c("sp1", "sp2"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(name = species, stringsAsFactors = FALSE)
  )
  membership <- igraph::components(g)$membership
  names(membership) <- igraph::V(g)$name

  labellings <- lapply(sort(unique(membership)), function(cid) {
    members <- names(membership)[membership == cid]
    .tp_block_labellings(members, pairs, group)
  })
  list(
    membership = membership, labellings = labellings,
    n = vapply(labellings, length, integer(1))
  )
}


#' Admissible labellings of one exchangeability block (internal)
#'
#' Each contrast constrains its two species to carry that contrast's two
#' labels, in either order. Fixing one species therefore forces its
#' partners, and forces theirs in turn, so a component's whole labelling is
#' determined by the label given to one seed species. Enumerating the
#' alphabet for the seed and propagating is linear in the contrasts, and
#' needs no assumption that the trait is binary or that the design is
#' balanced --- a contrast whose two species share a label pins them both,
#' and an inconsistent propagation is simply discarded.
#'
#' @noRd
.tp_block_labellings <- function(members, pairs, group) {
  in_block <- pairs$sp1 %in% members
  edges <- pairs[in_block, c("sp1", "sp2"), drop = FALSE]
  alphabet <- unique(unname(group[members]))
  seed <- members[1L]

  out <- list()
  for (start in alphabet) {
    lab <- stats::setNames(rep(NA_character_, length(members)), members)
    lab[seed] <- start
    queue <- seed
    ok <- TRUE
    while (length(queue) > 0L && ok) {
      u <- queue[1L]
      queue <- queue[-1L]
      touching <- which(edges$sp1 == u | edges$sp2 == u)
      for (ei in touching) {
        v <- if (edges$sp1[ei] == u) edges$sp2[ei] else edges$sp1[ei]
        # The contrast's own two labels, as a multiset. Removing the one
        # u carries leaves exactly what v must carry.
        multiset <- c(group[[edges$sp1[ei]]], group[[edges$sp2[ei]]])
        pos <- match(lab[[u]], multiset)
        if (is.na(pos)) {
          ok <- FALSE
          break
        }
        want <- multiset[-pos]
        if (is.na(lab[[v]])) {
          lab[v] <- want
          queue <- c(queue, v)
        } else if (!identical(lab[[v]], want)) {
          ok <- FALSE
          break
        }
      }
    }
    if (ok && !anyNA(lab)) {
      out[[length(out) + 1L]] <- lab
    }
  }
  # Two seed labels can propagate to the same assignment (a component
  # pinned by a same-label contrast), so collapse duplicates.
  keys <- vapply(
    out, function(l) paste(l[members], collapse = "\r"),
    character(1)
  )
  out <- out[!duplicated(keys)]
  if (length(out) == 0L) {
    # Cannot happen for a labelling derived from `group` itself, which is
    # admissible by construction, but never return an empty null.
    out <- list(stats::setNames(unname(group[members]), members))
  }
  out
}


#' Expected recurrence under independent sides of the observed sizes
#'
#' The raw recurrence count scales with how many HOGs the chosen sides
#' happen to hold, so a contrast whose two sides differ greatly in size
#' dominates the statistic and the relabelling null then measures set size
#' rather than shared identity. Subtracting the count expected from
#' independent uniform sides of exactly those sizes removes that term
#' without assuming anything about the trait.
#'
#' With \code{p_j} the chance one HOG falls in side \code{j}, the number of
#' sides containing it is Poisson-binomial, so the expected number of HOGs
#' reaching \code{min_recurrence} sides is \code{n_universe} times that
#' distribution's upper tail. The convolution is exact and costs
#' \code{O(k^2)}.
#'
#' @noRd
.tp_expected <- function(sizes, n_universe, min_recurrence) {
  if (n_universe <= 0L || length(sizes) == 0L) {
    return(0)
  }
  p <- pmin(1, sizes / n_universe)
  dist <- 1
  for (pj in p) {
    dist <- c(dist * (1 - pj), 0) + c(0, dist * pj)
  }
  if (min_recurrence > length(dist) - 1L) {
    return(0)
  }
  n_universe * sum(dist[(min_recurrence + 1L):length(dist)])
}
