#' Encode clique edge data as 0-based integer vectors for C++
#'
#' Shared helper for [find_cliques()] and [clique_stability()]. Builds
#' string-to-int maps, converts edge columns to 0-based integer vectors,
#' and filters out edges where either species is not in `target_species`.
#'
#' @param edges Data frame with columns gene1, gene2, species1, species2,
#'   hog, q_value, effect_size (already type-filtered).
#' @param target_species Character vector of species abbreviations.
#' @return A list with components: sp_map, gene_map, hog_map, all_genes,
#'   unique_hogs, edge_hog, edge_g1, edge_g2, edge_species1, edge_species2,
#'   edge_qval, edge_effect, and a logical `any_valid` flag.
#' @noRd
encode_clique_edges <- function(edges, target_species) {
  sp_map <- stats::setNames(seq_along(target_species) - 1L, target_species)

  # Gene identity is scoped by (species, gene), not gene name alone: gene
  # identifiers are not guaranteed unique across species, so keying only on
  # the raw name would collapse two different species' genes that happen to
  # share a string into one graph vertex. `all_genes` stays a plain name
  # vector for display (`enc$all_genes[idx + 1]`); `gene_key_map` is the
  # species-scoped lookup used to build the integer index space itself.
  gene_key1 <- paste(edges$species1, edges$gene1, sep = "\x02")
  gene_key2 <- paste(edges$species2, edges$gene2, sep = "\x02")
  all_gene_keys <- unique(c(gene_key1, gene_key2))
  gene_key_map <- stats::setNames(
    seq_along(all_gene_keys) - 1L, all_gene_keys
  )
  all_genes <- sub("^[^\x02]*\x02", "", all_gene_keys)

  unique_hogs <- unique(edges$hog)
  hog_map <- stats::setNames(
    seq_along(unique_hogs) - 1L,
    as.character(unique_hogs)
  )

  # Convert to 0-based integer vectors
  edge_hog <- as.integer(hog_map[as.character(edges$hog)])
  edge_g1 <- as.integer(gene_key_map[gene_key1])
  edge_g2 <- as.integer(gene_key_map[gene_key2])
  edge_species1 <- as.integer(sp_map[edges$species1])
  edge_species2 <- as.integer(sp_map[edges$species2])

  # Filter out edges where either species is not in target_species
  valid <- !is.na(edge_species1) & !is.na(edge_species2)
  any_valid <- any(valid)

  edge_hog <- edge_hog[valid]
  edge_g1 <- edge_g1[valid]
  edge_g2 <- edge_g2[valid]
  edge_species1 <- edge_species1[valid]
  edge_species2 <- edge_species2[valid]
  edge_qval <- as.numeric(edges$q_value[valid])
  edge_effect <- as.numeric(edges$effect_size[valid])

  list(
    sp_map = sp_map, gene_map = gene_key_map, hog_map = hog_map,
    all_genes = all_genes, unique_hogs = unique_hogs,
    edge_hog = edge_hog, edge_g1 = edge_g1, edge_g2 = edge_g2,
    edge_species1 = edge_species1, edge_species2 = edge_species2,
    edge_qval = edge_qval, edge_effect = edge_effect,
    any_valid = any_valid
  )
}


#' Entropy-maximising scale for the ensemble weight map
#'
#' Solves `S'(z) = 0` for the connection-probability map
#' `p(w, z) = z w / (1 + z w)`, where `S` is the Shannon entropy of the
#' induced binary ensemble (Garlaschelli, Ahnert, Fink & Caldarelli
#' 2013). `S'` is positive below `1 / max(w)` and negative above
#' `1 / min(w)`, so the root is bracketed and unique. A collapsed
#' bracket -- one weight, or all weights equal -- has nothing to
#' separate and returns `1 / mean(w)`, giving `p = 1/2`.
#'
#' @param w Numeric vector of non-negative weights (one species pair).
#' @return The scale `z`, or `NA_real_` when `w` holds no usable value.
#' @noRd
.maxent_scale <- function(w) {
  # Entropy-maximising scale z* for p(w, z) = z w / (1 + z w)
  # (Garlaschelli, Ahnert, Fink & Caldarelli 2013). S'(z) is positive
  # below 1/max(w) and negative above 1/min(w), so the root is bracketed
  # and unique; uniroot finds it in one pass over the weights.
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) {
    return(NA_real_)
  }
  lo <- 1 / max(w)
  hi <- 1 / min(w)
  # One edge, or all weights equal: the bracket collapses and there is
  # nothing to tell the weights apart. Maximum entropy for a single
  # probability is 1/2, which z = 1/w delivers -- the least biased answer
  # rather than an NA that would silently void the clique's intensity.
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    return(1 / mean(w))
  }
  ds <- function(z) sum(w / (1 + z * w)^2 * log(1 / (z * w)))
  # Endpoints are the bracket by construction; guard against the
  # degenerate case where floating point puts both on the same side.
  if (ds(lo) <= 0 || ds(hi) >= 0) {
    return(NA_real_)
  }
  stats::uniroot(ds, lower = lo, upper = hi, tol = .Machine$double.eps^0.5)$root
}


#' Ensemble connection probability of every edge row
#'
#' The Onnela weight used by [find_cliques()]. Association strength
#' (`effect_size`, the observed overlap over its expectation) is mapped to
#' a connection probability `p = z w / (1 + z w)`, with the scale `z` fixed
#' per species pair by maximum entropy. Intensity is then the geometric
#' mean of those probabilities, i.e. the per-edge probability that the
#' whole clique exists in the binary ensemble the weighted graph induces.
#'
#' Association strength is used rather than the Jaccard index because
#' `E[Jaccard]` grows with neighbourhood size, so a Jaccard weight ranks
#' hub genes above equally conserved low-degree ones (van Eck & Waltman
#' 2009; measured here as Spearman +0.58 with clique mean member degree,
#' against -0.08 for this weight).
#'
#' @param edges Edge table with `effect_size`, `species1`, `species2`.
#' @return Numeric vector parallel to the rows of `edges`, in (0, 1).
#' @references
#' Onnela, Saramaki, Kertesz & Kaski (2005) Phys. Rev. E 71, 065103.
#' Garlaschelli, Ahnert, Fink & Caldarelli (2013) LNCS 8852, 107-118.
#' van Eck & Waltman (2009) J. Am. Soc. Inf. Sci. Technol. 60, 1635-1651.
#' @noRd
.onnela_weight <- function(edges) {
  n <- nrow(edges)
  out <- rep(NA_real_, n)
  if (n == 0L) {
    return(out)
  }
  as_w <- if ("effect_size" %in% names(edges)) {
    as.numeric(edges$effect_size)
  } else {
    rep(NA_real_, n)
  }
  as_w[!is.finite(as_w) | as_w <= 0] <- NA_real_
  if (all(is.na(as_w))) {
    rlang::warn(
      c(
        "`edges` has no usable `effect_size` column.",
        i = paste0(
          "Clique intensity maps each edge's association strength to ",
          "a connection probability and is NA without it; ",
          "find_coexpressologs() output carries it. This warning is ",
          "shown once per session."
        )
      ),
      class = "rcomplex_missing_effect_size",
      .frequency = "once",
      .frequency_id = "rcomplex_missing_effect_size"
    )
    return(out)
  }
  pair <- paste(pmin(edges$species1, edges$species2),
    pmax(edges$species1, edges$species2),
    sep = "\x01"
  )
  for (idx in split(seq_len(n), pair)) {
    ok <- idx[!is.na(as_w[idx])]
    if (length(ok) == 0L) next
    z <- .maxent_scale(as_w[ok])
    if (!is.finite(z)) next
    out[ok] <- z * as_w[ok] / (1 + z * as_w[ok])
  }
  out
}


#' Warn when an edge table looks already cut to `conserved`
#'
#' Clique intensity fits each edge's weight scale on every tested pair of
#' its species pair. A table holding only `conserved` rows fits that scale
#' on significant edges alone instead, which is a different quantity that
#' nothing downstream can tell apart. Tables without a
#' `type` column cannot be judged and pass silently. Warns once per
#' session: clique_stability(), clique_threshold_sweep() and
#' classify_cliques() call find_cliques()
#' many times on one table.
#'
#' @param edges Edge data frame.
#' @noRd
.warn_if_prefiltered <- function(edges) {
  if (!"type" %in% names(edges) || nrow(edges) == 0L ||
        !all(edges$type %in% "conserved")) {
    return(invisible(FALSE))
  }
  rlang::warn(
    c(
      paste0(
        "`edges` holds only `conserved` rows (",
        paste(unique(edges$type), collapse = ", "),
        "), so it looks pre-filtered."
      ),
      i = paste0(
        "Clique intensity fits the weight scale on every ",
        "tested pair of its species pair; here that population is only ",
        "the rows already kept."
      ),
      i = paste0(
        "Pass the unfiltered edge table, e.g. find_coexpressologs() ",
        "output. This warning is shown once per session."
      )
    ),
    class = "rcomplex_prefiltered_edges",
    .frequency = "once",
    .frequency_id = "rcomplex_prefiltered_edges"
  )
  invisible(TRUE)
}


# Which rows of `edges` belong to each clique. Keyed by hog + sorted
# (species, gene) endpoints: gene identifiers are not guaranteed unique
# across species, so a key built from gene names alone could match an
# edge from an unrelated species pair that happens to share both gene
# strings; folding species into each endpoint before sorting keeps the
# pair (and its species) intact. compute_clique_edge_stats() and the
# matched-edge null both route through this, so the observed intensity
# and its null are scored over provably the same edge set.
.clique_edge_rows <- function(cliques, edges, target_species) {
  ek1 <- paste(edges$species1, edges$gene1, sep = "\x02")
  ek2 <- paste(edges$species2, edges$gene2, sep = "\x02")
  edge_key <- paste(edges$hog, pmin(ek1, ek2), pmax(ek1, ek2), sep = "\x01")
  edge_idx <- stats::setNames(seq_len(nrow(edges)), edge_key)
  lapply(seq_len(nrow(cliques)), function(i) {
    row_vals <- cliques[i, target_species, drop = TRUE]
    present <- !is.na(row_vals)
    genes <- unlist(row_vals[present], use.names = FALSE)
    gene_sp <- target_species[present]
    if (length(genes) < 2L) {
      return(integer(0))
    }
    pairs <- utils::combn(seq_along(genes), 2L)
    k1 <- paste(gene_sp[pairs[1L, ]], genes[pairs[1L, ]], sep = "\x02")
    k2 <- paste(gene_sp[pairs[2L, ]], genes[pairs[2L, ]], sep = "\x02")
    keys <- paste(cliques$hog[i], pmin(k1, k2), pmax(k1, k2), sep = "\x01")
    matched <- edge_idx[keys]
    unname(matched[!is.na(matched)])
  })
}


#' Compute per-clique edge statistics (intensity, min effect size)
#'
#' Onnela weights are ensemble connection probabilities, not
#' `1 - q_value`: a clique only ever contains edges that passed alpha, so
#' `1 - q` sat above `1 - alpha` on every edge and left intensity flat
#' (#11). They are built from association strength rather than the
#' Jaccard index, whose null expectation grows with neighbourhood size
#' (#15).
#'
#' @param cliques Data frame from find_cliques (with hog + species columns).
#' @param edges Data frame with gene1, gene2, hog, q_value, effect_size.
#' @param target_species Character vector of species names.
#' @param weights Per-row edge weights parallel to `edges`; by default
#'   the ensemble connection probability of each row's association
#'   strength, scaled within its species pair. Callers that filter
#'   `edges` must weight first and subset with the rows, so the scale is
#'   still fitted on every tested pair.
#' @return Data frame with columns intensity, min_effect_size.
#' @noRd
compute_clique_edge_stats <- function(cliques, edges, target_species,
                                      weights = .onnela_weight(edges)) {
  n <- nrow(cliques)
  intensity <- rep(NA_real_, n)
  min_eff <- rep(NA_real_, n)
  rows <- .clique_edge_rows(cliques, edges, target_species)

  for (i in seq_len(n)) {
    matched <- rows[[i]]
    if (length(matched) == 0L) next

    effs <- edges$effect_size[matched]
    w <- weights[matched]
    # An edge without a finite weight leaves the clique's weight set
    # incomplete. Scoring the remaining edges would silently score a
    # different edge set, so an incomplete clique reports NA instead.
    if (!anyNA(w)) {
      intensity[i] <- exp(mean(log(w)))
    }
    min_eff[i] <- min(effs)
  }

  data.frame(intensity = intensity, min_effect_size = min_eff)
}


#' Find co-expression cliques using C++ two-level decomposition
#'
#' For each Hierarchical Ortholog Group (HOG), uses Bron-Kerbosch with
#' pivoting on the species-level adjacency graph to find maximal species
#' cliques, then backtracking to assign the best gene per species
#' (minimising mean q-value across all present edges).
#'
#' @param edges Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (first gene in pair)}
#'     \item{gene2}{Gene identifier (second gene in pair)}
#'     \item{species1}{Species for gene1}
#'     \item{species2}{Species for gene2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{q_value}{q-value for the edge (from pair-level testing)}
#'     \item{effect_size}{Numeric effect size}
#'   }
#'   Optionally includes a \code{type} column for filtering. Pass the
#'   unfiltered table (every tested pair, e.g. the output of
#'   \code{\link{find_coexpressologs}}): cliques are built from
#'   \code{"conserved"} rows only, but \code{intensity} fits each edge's
#'   weight scale on every row of its species pair, so a pre-filtered
#'   table changes what it measures.
#'   A table whose \code{type} column holds only \code{"conserved"} rows
#'   triggers a warning (class \code{rcomplex_prefiltered_edges}, shown
#'   once per session).
#' @param target_species Character vector of species abbreviations.
#' @param min_species Minimum number of species per clique
#'   (default: \code{length(target_species)}).
#' @param cost_weights Named numeric vector with elements \code{"q"} and
#'   \code{"effect"} controlling the composite cost used for gene-assignment
#'   ranking. The cost is \code{q * mean_q - effect * mean_effect}
#'   (lower is better). Default \code{c(q = 1, effect = 0)} reproduces
#'   the original mean-q-only ranking.
#'
#' @return A data frame with one row per clique:
#'   \describe{
#'     \item{hog}{Ortholog group identifier}
#'     \item{<species>}{One column per target species (gene ID or NA)}
#'     \item{n_species}{Number of species in the clique}
#'     \item{mean_q}{Mean q-value across present clique edges}
#'     \item{max_q}{Maximum q-value across present clique edges}
#'     \item{mean_effect_size}{Mean effect size across present edges}
#'     \item{n_edges}{Number of present edges}
#'     \item{n_missing}{Number of missing edges (always 0)}
#'     \item{intensity}{Onnela intensity: geometric mean, across present
#'       edges, of each edge's ensemble connection probability (in
#'       (0, 1); higher = stronger conservation). Each edge's
#'       association strength (\code{effect_size}, observed overlap over
#'       its expectation) is mapped to a probability
#'       \code{p = z w / (1 + z w)}, with \code{z} fitted per species
#'       pair by maximum entropy, so intensity reads as the per-edge
#'       probability that the whole clique exists. Not \code{1 - q_value}:
#'       every clique edge already passed alpha, so \code{1 - q} left no
#'       range (#11). Not the Jaccard index, whose null expectation grows
#'       with neighbourhood size, which ranked hub genes above equally
#'       conserved low-degree ones (#15). \code{NA} when \code{edges} has
#'       no usable \code{effect_size} column, or when any present clique
#'       edge lacks a finite value.}
#'     \item{min_effect_size}{Minimum effect size across present edges
#'       (bottleneck enrichment)}
#'   }
#'
#' @examples
#' \dontrun{
#' cliques <- find_cliques(edges, target_species = c("SP_A", "SP_B", "SP_C"))
#'
#' # Allow partial species cliques (2 of 3 species)
#' partial <- find_cliques(edges, target_species, min_species = 2L)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
find_cliques <- function(edges, ...) UseMethod("find_cliques")

#' @rdname find_cliques
#' @export
find_cliques.default <- function(edges, target_species,
                                 min_species = length(target_species),
                                 cost_weights = c(q = 1.0, effect = 0.0), ...) {
  # Validate inputs
  required_cols <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q_value", "effect_size"
  )
  missing <- setdiff(required_cols, names(edges))
  if (length(missing) > 0) {
    stop("edges missing required columns: ", paste(missing, collapse = ", "))
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  min_species <- as.integer(min_species)
  if (min_species < 2L) stop("min_species must be >= 2")
  if (min_species > length(target_species)) {
    stop("min_species must be <= length(target_species)")
  }

  # Validate cost_weights
  if (!is.numeric(cost_weights) || length(cost_weights) != 2L) {
    stop("cost_weights must be a named numeric vector of length 2")
  }
  if (is.null(names(cost_weights)) ||
        !all(c("q", "effect") %in% names(cost_weights))) {
    stop("cost_weights must have names 'q' and 'effect'")
  }
  if (any(cost_weights < 0)) {
    stop("cost_weights values must be >= 0")
  }
  w_q <- as.double(cost_weights[["q"]])
  w_eff <- as.double(cost_weights[["effect"]])

  # Empty result template
  empty_cols <- c(
    list(hog = character(0)),
    stats::setNames(lapply(target_species, \(x) character(0)), target_species),
    list(
      n_species = integer(0), mean_q = numeric(0), max_q = numeric(0),
      mean_effect_size = numeric(0), n_edges = integer(0),
      n_missing = integer(0),
      intensity = numeric(0), min_effect_size = numeric(0)
    )
  )
  empty_result <- as.data.frame(empty_cols)

  # Onnela weights fit their scale on every tested pair of its species
  # pair, so they are taken before the type filter; fitted among
  # conserved edges alone they would only describe edges that already
  # passed alpha.
  .warn_if_prefiltered(edges)
  weights <- .onnela_weight(edges)
  if ("type" %in% names(edges)) {
    keep <- edges$type %in% "conserved"
    edges <- edges[keep, , drop = FALSE]
    weights <- weights[keep]
  }
  if (nrow(edges) == 0) {
    return(empty_result)
  }

  # Encode edges as 0-based integer vectors
  enc <- encode_clique_edges(edges, target_species)
  if (!enc$any_valid) {
    return(empty_result)
  }

  # Call C++
  result <- find_cliques_cpp(
    enc$edge_hog, enc$edge_g1, enc$edge_g2,
    enc$edge_species1, enc$edge_species2,
    enc$edge_qval, enc$edge_effect,
    length(target_species), min_species,
    length(enc$unique_hogs), length(enc$all_genes),
    10L, 0L,
    w_q, w_eff
  )

  # Map back to strings
  if (length(result$hog_idx) == 0) {
    return(empty_result)
  }

  # HOG names
  hog_names <- enc$unique_hogs[result$hog_idx + 1L]

  # Gene matrix: map 0-based indices back to gene names
  gene_matrix <- result$genes # IntegerMatrix (n_cliques x n_species)
  gene_df <- as.data.frame(matrix(NA_character_,
    nrow = nrow(gene_matrix),
    ncol = ncol(gene_matrix)
  ))
  names(gene_df) <- target_species
  for (j in seq_len(ncol(gene_matrix))) {
    idx <- gene_matrix[, j]
    present <- !is.na(idx)
    gene_df[present, j] <- enc$all_genes[idx[present] + 1L]
  }

  # Build result data frame
  out <- data.frame(hog = hog_names)
  out <- cbind(out, gene_df)
  out$n_species <- as.integer(result$n_species)
  out$mean_q <- result$mean_q
  out$max_q <- result$max_q
  out$mean_effect_size <- result$mean_effect_size
  out$n_edges <- as.integer(result$n_edges)
  out$n_missing <- as.integer(result$n_missing)

  # Compute Onnela intensity and min effect size
  stats <- compute_clique_edge_stats(out, edges, target_species,
    weights = weights
  )
  out$intensity <- stats$intensity
  out$min_effect_size <- stats$min_effect_size
  out
}


#' Classify HOGs by clique conservation pattern
#'
#' Convenience wrapper that runs \code{\link{find_cliques}} internally
#' (once for all species, once per top-level clade) and applies a sequential
#' waterfall classification. For fine-grained control over the clique
#' detection parameters per step, call \code{find_cliques()} directly.
#'
#' Each HOG is assigned to exactly one category; earlier categories
#' take precedence.
#'
#' The pipeline:
#' \enumerate{
#'   \item \strong{complete}: all target species form a clique (all
#'     \code{C(N,2)} edges conserved).
#'   \item \strong{partial}: a clique with \code{min_species}
#'     to \code{N-1} species lies in no clade.
#'   \item \strong{differentiated}: cliques sit in two or more disjoint
#'     clades, and no conserved edge joins two top-level clades.
#'   \item \strong{trait_specific}: all cliques sit in one clade.
#'   \item \strong{unclassified}: none of the above.
#' }
#' The home clade of a clique is the smallest clade that holds all its
#' species. A HOG keeps the home clades that no other of its home clades
#' holds. Clades are laminar, so these are disjoint. Nested homes count as
#' one clade, so they are not \code{differentiated}. Cliques are searched
#' within each top-level clade. A species in no clade forms its own clade.
#'
#' @section Underpowered calls:
#' \code{differentiated} and \code{trait_specific} both rest on edges
#' that were \emph{not} conserved. An edge too weak to have been called
#' is no evidence either way, so such a call has to survive reading that
#' edge as conserved. Where it does not, the row is flagged
#' \code{underpowered = TRUE} and keeps its classification -- the call
#' is reported and qualified, not replaced, so nothing is lost and a
#' caller can filter on it. Needs a \code{power} column in
#' \code{edges}; without one, or where it is \code{NA}, the flag is
#' \code{FALSE}.
#'
#' What the flag marks: power is computed against the typical fold
#' enrichment of a called pair (\code{rho0}, about 2.5 on the Pooideae
#' data), and at that enrichment a low-degree gene's shared partners
#' amount to less than one, which no count test can call. So
#' \code{underpowered} marks genes whose neighbourhoods are too small
#' for ordinary conservation to be visible -- on the Pooideae data most
#' of the lowest degree decile and none of the highest.
#'
#' @section Choosing between the two clique classifiers:
#' This function and \code{\link{classify_gene_cliques}} answer
#' different questions; neither is deprecated in favour of the other.
#'
#' \code{classify_cliques()} works on the per-orthogroup \emph{species}
#' graph. It asks which species are joined by conserved co-expression
#' and whether that pattern respects the trait split, and it returns one
#' row per HOG. \code{\link{find_cliques}} commits to one best gene
#' assignment per species clique, so a multi-copy HOG still gets a
#' single answer, and that answer is what \code{\link{clique_stability}}
#' and \code{\link{clique_threshold_sweep}} consume -- the
#' \code{stability_class} / \code{robust} columns exist only on this
#' side.
#'
#' \code{\link{classify_gene_cliques}} works on the \emph{gene} graph
#' built by \code{\link{gene_clique_graph}}. It asks which individual
#' gene copies are mutually conserved, so one HOG can yield several
#' overlapping cliques and the answer names paralogs rather than
#' species. It applies the taxonomy of Rodriguez et al. (2026), plus
#' \code{trait_specific}, with two explicit tolerance tiers
#' (\code{partial_significant}
#' for weak wiring, \code{partial_present} for a missing gene).
#'
#' Reach for \code{classify_gene_cliques()} when which copy sits in the
#' conserved core matters, or when the published taxonomy is what has to
#' be reported. Reach for \code{classify_cliques()} for a
#' one-row-per-HOG trait summary wired into the stability machinery.
#'
#' @param edges Data frame with columns \code{gene1}, \code{gene2},
#'   \code{species1}, \code{species2}, \code{hog}, \code{q_value},
#'   \code{effect_size}, and \code{type}. Must contain ALL edges
#'   (conserved + ns + diverged), not pre-filtered, because the
#'   differentiated check needs to verify absence of cross-group
#'   conserved edges. An optional \code{power} column (from
#'   \code{\link{find_coexpressologs}}) enables the
#'   \code{"underpowered"} class; without it, or where it is \code{NA},
#'   the classification is unchanged.
#' @param target_species Character vector of all species.
#' @param clades Named list of species vectors, one per clade. Clades
#'   may nest but must not cross. A species in no clade forms its own
#'   clade. A flat trait is \code{split(names(trait), trait)}. A tree
#'   gives \code{\link{clades_from_tree}(phy)}.
#' @param min_species Minimum species for a partial or within-group
#'   clique (default 2).
#' @param stability Optional output of \code{\link{clique_stability}}.
#' @param min_power Detection power below which a non-conserved edge is
#'   read as uninformative rather than as evidence against conservation
#'   (default 0.8). Only used when \code{edges} has \code{power}. For
#'   rank-test edges (\code{find_coexpressologs(method = "rank")}) use 0.9:
#'   under the default reference rank \code{p0} their power is at least
#'   0.5 by construction (0 when a direction of the species pair has no
#'   call under \code{pval_combine = "max"}, or neither has under
#'   \code{"min"}), and it overstates detection.
#'
#' @return A data frame with one row per HOG:
#'   \describe{
#'     \item{hog}{HOG identifier}
#'     \item{classification}{One of \code{"complete"}, \code{"partial"},
#'       \code{"differentiated"}, \code{"trait_specific"},
#'       \code{"unclassified"}}
#'     \item{underpowered}{\code{TRUE} where a \code{differentiated}
#'       or \code{trait_specific} call rests on an edge whose
#'       \code{power} is below \code{min_power}. The classification is
#'       kept; this qualifies it. \code{FALSE} for every other tier and
#'       whenever \code{edges} carries no \code{power} column.}
#'     \item{n_species}{Species count in the best clique (NA for
#'       unclassified)}
#'     \item{best_mean_q}{Mean q-value of the best clique (NA for
#'       unclassified)}
#'     \item{trait_groups}{Comma-separated top-level clades that hold
#'       the cliques (NA for complete/partial/unclassified)}
#'     \item{clade}{For \code{trait_specific}, the home clade of the
#'       cliques. For \code{differentiated}, the smallest clade that holds
#'       all home clades. Otherwise, the smallest clade that holds the
#'       reported clique. \code{NA} when no clade holds them, and for
#'       unclassified HOGs.}
#'     \item{stability_class}{From stability results (NA if not
#'       provided)}
#'     \item{robust}{Logical: stability class is at least 0 (NA if no
#'       stability results provided)}
#'   }
#'
#' @examples
#' \dontrun{
#' clades <- list(annual = c("SP_A", "SP_B"), perennial = c("SP_C", "SP_D"))
#' result <- classify_cliques(edges, target_species, clades)
#' table(result$classification)
#' }
#'
#' @seealso \code{\link{find_cliques}} for the species-graph backend
#'   this wraps; \code{\link{classify_gene_cliques}} and
#'   \code{\link{gene_clique_graph}} for the copy-level alternative
#'   described above.
#' @references
#' Rodriguez E, Birkeland S, Chapple ED, et al. (2026).
#' Comparative regulomics of wood formation across dicot and
#' conifer trees. \emph{Nature Communications} 17(1).
#' \doi{10.1038/s41467-026-75624-2}
#'
#' @param ... Additional arguments passed to the default method.
#' @export
classify_cliques <- function(edges, ...) UseMethod("classify_cliques")

#' @rdname classify_cliques
#' @export
classify_cliques.default <- function(
  edges, target_species, clades,
  min_species = 2L,
  stability = NULL,
  min_power = 0.8, ...
) {
  # --- Validation ---
  required_cols <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q_value", "effect_size", "type"
  )
  missing_cols <- setdiff(required_cols, names(edges))
  if (length(missing_cols) > 0) {
    stop(
      "edges missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  clades <- .check_clades(clades, target_species)
  min_species <- as.integer(min_species)
  if (min_species < 2L) stop("min_species must be >= 2")
  ok_power <- is.numeric(min_power) && length(min_power) == 1L &&
    !is.na(min_power) && min_power >= 0 && min_power <= 1
  if (!ok_power) {
    stop("min_power must be a single number in [0, 1]")
  }

  if (!is.null(stability)) {
    if (!is.list(stability) || is.null(stability$stability)) {
      stop("stability must be output of clique_stability()")
    }
  }

  trait_char <- .clade_groups(clades, target_species)
  trait_levels <- unique(trait_char)
  n_sp <- length(target_species)
  all_hogs <- unique(edges$hog)

  # Empty result template
  empty <- data.frame(
    hog = character(0), classification = character(0),
    n_species = integer(0), best_mean_q = numeric(0),
    trait_groups = character(0), clade = character(0),
    underpowered = logical(0),
    stability_class = integer(0),
    robust = logical(0),
    stringsAsFactors = FALSE
  )
  if (length(all_hogs) == 0) {
    return(empty)
  }

  # --- Steps 1+2: Find all cliques (complete + partial in one pass) ---
  all_cliques <- .cc_with_clade(
    find_cliques(edges, target_species, min_species = min_species),
    clades, target_species
  )

  # Complete = all N species, no missing edges
  is_complete <- all_cliques$n_species == n_sp & all_cliques$n_missing == 0L
  complete_hogs <- unique(all_cliques$hog[is_complete])
  # Partial cliques must span at least 2 trait groups (cross-group)
  # to distinguish from trait_specific / differentiated patterns.
  partial_candidates <- setdiff(unique(all_cliques$hog), complete_hogs)
  partial_hogs <- character(0)
  for (h in partial_candidates) {
    hog_rows <- all_cliques[all_cliques$hog == h, , drop = FALSE]
    # Check if any clique row spans 2+ trait groups
    any_cross <- FALSE
    for (r in seq_len(nrow(hog_rows))) {
      spp <- target_species[!is.na(hog_rows[r, target_species])]
      traits_in <- unique(trait_char[spp])
      if (length(traits_in) >= 2L) {
        any_cross <- TRUE
        break
      }
    }
    if (any_cross) partial_hogs <- c(partial_hogs, h)
  }

  # --- Step 3: Within-group cliques per trait group ---
  within_group_cliques <- list()

  for (group in trait_levels) {
    group_sp <- names(trait_char[trait_char == group])
    if (length(group_sp) < 2L) {
      within_group_cliques[[group]] <- NULL
      next
    }
    wg <- .cc_with_clade(
      find_cliques(edges, group_sp, min_species = min_species),
      clades, group_sp
    )
    within_group_cliques[[group]] <- wg
  }

  # --- Step 4: Differentiated (2+ groups w/ cliques, no cross-group) ---
  remaining <- setdiff(all_hogs, c(complete_hogs, partial_hogs))

  # Identify cross-group conserved edges
  conserved <- edges[edges$type %in% "conserved", , drop = FALSE]
  if (nrow(conserved) > 0) {
    t1 <- trait_char[conserved$species1]
    t2 <- trait_char[conserved$species2]
    # Only keep edges where both species are in target_species
    valid <- !is.na(t1) & !is.na(t2)
    cross_conserved <- conserved[valid & t1 != t2, , drop = FALSE]
    hogs_with_cross <- unique(cross_conserved$hog)
  } else {
    hogs_with_cross <- character(0)
  }

  # Home clades of each HOG's within-group cliques, kept when no other
  # home holds them. Clades are laminar, so the kept homes are disjoint.
  wg_all <- do.call(rbind, lapply(within_group_cliques, function(df) {
    df[, c("hog", "clade"), drop = FALSE]
  }))
  homes_by_hog <- list()
  if (!is.null(wg_all)) homes_by_hog <- split(wg_all$clade, wg_all$hog)
  first_sp <- vapply(clades, function(v) min(match(v, target_species)), 1)
  top_homes <- function(h) {
    hm <- unique(homes_by_hog[[h]])
    held <- vapply(hm, function(a) {
      any(vapply(setdiff(hm, a), function(b) {
        all(clades[[a]] %in% clades[[b]])
      }, logical(1)))
    }, logical(1))
    hm <- hm[!held]
    hm[order(first_sp[hm])]
  }
  groups_of <- function(hm) {
    paste(unique(trait_char[vapply(clades[hm], `[`, "", 1L)]),
      collapse = ","
    )
  }

  # Differentiated: two or more disjoint homes, no cross-group edge.
  # Trait-specific: one home. Nested homes count as one.
  diff_hogs <- character(0)
  ts_hogs <- character(0)
  hog_homes <- list()
  for (h in remaining) {
    hm <- top_homes(h)
    if (length(hm) >= 2L && !h %in% hogs_with_cross) {
      diff_hogs <- c(diff_hogs, h)
    } else if (length(hm) == 1L) {
      ts_hogs <- c(ts_hogs, h)
    } else {
      next
    }
    hog_homes[[h]] <- hm
  }
  hog_groups <- vapply(hog_homes, groups_of, "")
  hog_clade <- vapply(hog_homes, function(hm) {
    if (length(hm) == 1L) hm else .clade_home(clades, unlist(clades[hm]))
  }, "")

  # --- Step 6: Unclassified ---
  classified <- c(complete_hogs, partial_hogs, diff_hogs, ts_hogs)
  unclass_hogs <- setdiff(all_hogs, classified)

  # --- Build output ---
  # Helper: best clique per HOG from a cliques df
  best_per_hog <- function(cliques_df, hogs) {
    sub <- cliques_df[cliques_df$hog %in% hogs, , drop = FALSE]
    if (nrow(sub) == 0) {
      return(data.frame(
        hog = character(0), n_species = integer(0),
        best_mean_q = numeric(0), clade = character(0)
      ))
    }
    sub <- sub[order(sub$mean_q), , drop = FALSE]
    sub <- sub[!duplicated(sub$hog), , drop = FALSE]
    data.frame(
      hog = sub$hog, n_species = sub$n_species,
      best_mean_q = sub$mean_q, clade = sub$clade,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()

  # Complete
  if (length(complete_hogs) > 0) {
    info <- best_per_hog(
      all_cliques[is_complete, , drop = FALSE], complete_hogs
    )
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "complete",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = NA_character_, clade = info$clade,
      stringsAsFactors = FALSE
    )
  }

  # Partial
  if (length(partial_hogs) > 0) {
    info <- best_per_hog(all_cliques, partial_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "partial",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = NA_character_, clade = info$clade,
      stringsAsFactors = FALSE
    )
  }

  # Differentiated
  if (length(diff_hogs) > 0) {
    # Best within-group clique for each differentiated HOG
    wg_summary <- do.call(rbind, lapply(
      within_group_cliques[trait_levels],
      function(df) {
        if (!is.null(df) && nrow(df) > 0) {
          df[, c("hog", "n_species", "mean_q", "clade"), drop = FALSE]
        } else {
          NULL
        }
      }
    ))
    info <- best_per_hog(wg_summary, diff_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog,
      classification = "differentiated",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = unname(hog_groups[info$hog]),
      clade = unname(hog_clade[info$hog]), stringsAsFactors = FALSE
    )
  }

  # Trait-specific
  if (length(ts_hogs) > 0) {
    wg_summary2 <- do.call(rbind, lapply(
      within_group_cliques[trait_levels],
      function(df) {
        if (!is.null(df) && nrow(df) > 0) {
          df[, c("hog", "n_species", "mean_q", "clade"), drop = FALSE]
        } else {
          NULL
        }
      }
    ))
    info <- best_per_hog(wg_summary2, ts_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog,
      classification = "trait_specific",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = unname(hog_groups[info$hog]),
      clade = unname(hog_clade[info$hog]), stringsAsFactors = FALSE
    )
  }

  # Unclassified
  if (length(unclass_hogs) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      hog = unclass_hogs, classification = "unclassified",
      n_species = NA_integer_, best_mean_q = NA_real_,
      trait_groups = NA_character_, clade = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(empty)
  }
  rownames(out) <- NULL

  # --- Underpowered annotation ---
  # A specificity or divergence call that rests on an edge too weak to
  # have been called is not evidence, but it is also not a different
  # kind of clique: overwriting the classification threw the call away
  # and left no way to recover it. Carry it as a flag instead, so the
  # call survives and a caller can filter on it. Each call is read
  # against the home clades of its cliques.
  up_hogs <- character(0)
  if ("power" %in% names(edges)) {
    up_hogs <- .cc_underpowered_hogs(
      edges, c(diff_hogs, ts_hogs), within_group_cliques, trait_char,
      min_power, hog_homes, clades
    )
  }
  out$underpowered <- out$hog %in% up_hogs

  # --- Stability annotation ---
  out$stability_class <- NA_integer_
  if (!is.null(stability$stability) && nrow(stability$stability) > 0) {
    stab_df <- stability$stability
    sc_vec <- stability$stability_class
    if ("hog" %in% names(stab_df) && length(sc_vec) > 0) {
      # stability_class is per-clique (same order as full_cliques).
      # Map to HOG via the stability DF.
      all_clique_idx <- unique(stab_df$clique_idx)
      all_clique_hogs <- stab_df$hog[match(all_clique_idx, stab_df$clique_idx)]
      # Best (max) stability_class per HOG: a HOG's classification is
      # determined by its best clique, so we take the most stable one.
      # Use max (optimistic) not min (conservative) — consistent with
      # best_per_hog() which picks the lowest-mean-q clique per HOG.
      best_sc <- tapply(sc_vec[all_clique_idx], all_clique_hogs, max)
      idx <- match(out$hog, names(best_sc))
      out$stability_class[!is.na(idx)] <- as.integer(best_sc[idx[!is.na(idx)]])
    }
  }

  # --- Robust flag ---
  has_stab <- !is.null(stability$stability) && nrow(stability$stability) > 0
  if (has_stab) {
    out$robust <- !is.na(out$stability_class) & out$stability_class >= 0L
  } else {
    out$robust <- NA
  }

  out
}


#' HOGs whose specificity or divergence call rests on underpowered edges
#'
#' A deciding edge is a non-conserved row of the HOG with `power` below
#' `min_power`, one endpoint a member of one of the HOG's within-group
#' cliques and the other across a boundary. The boundaries are the home
#' clades of the HOG's cliques, then the top-level clades. Such an edge
#' could not have been called, so the call cannot rule it out as
#' conserved.
#'
#' @param edges Full edge table carrying `power`.
#' @param hogs Candidate HOGs (differentiated and trait-specific).
#' @param wg_cliques Named list of within-group [find_cliques()] tables.
#' @param trait_char Named trait group of every target species.
#' @param min_power As in [classify_cliques()].
#' @param hog_homes Home clades of each candidate HOG, named by HOG.
#' @param clades Checked clade list.
#' @return Character vector of the HOGs to reclassify.
#' @noRd
.cc_underpowered_hogs <- function(edges, hogs, wg_cliques, trait_char,
                                  min_power, hog_homes, clades) {
  hogs <- as.character(hogs)
  e_hog <- as.character(edges$hog)
  pw <- as.numeric(edges$power)
  # Inside one of its HOG's home clades, a species takes that name.
  sep <- "\x01"
  home <- unlist(hog_homes, use.names = FALSE)
  n_in <- lengths(clades[home])
  key <- paste(rep(rep(names(hog_homes), lengths(hog_homes)), n_in),
    unlist(clades[home], use.names = FALSE),
    sep = sep
  )
  val <- rep(home, n_in)
  side <- function(sp) {
    v <- val[match(paste(e_hog, sp, sep = sep), key)]
    ifelse(is.na(v), unname(trait_char[sp]), v)
  }
  t1 <- side(as.character(edges$species1))
  t2 <- side(as.character(edges$species2))
  cand <- e_hog %in% hogs & !(edges$type %in% "conserved") &
    !is.na(pw) & pw < min_power & !is.na(t1) & !is.na(t2) & t1 != t2
  if (!any(cand)) {
    return(character(0))
  }
  members <- unlist(lapply(wg_cliques, function(df) {
    if (is.null(df) || nrow(df) == 0L) {
      return(NULL)
    }
    df <- df[as.character(df$hog) %in% hogs, , drop = FALSE]
    sp_cols <- intersect(names(trait_char), names(df))
    unlist(lapply(sp_cols, function(s) {
      g <- df[[s]]
      ok <- !is.na(g)
      paste(df$hog[ok], s, g[ok], sep = sep)
    }))
  }), use.names = FALSE)
  k1 <- paste(e_hog, edges$species1, edges$gene1, sep = sep)[cand]
  k2 <- paste(e_hog, edges$species2, edges$gene2, sep = sep)[cand]
  unique(e_hog[cand][k1 %in% members | k2 %in% members])
}


#' Add a `clade` column to a [find_cliques()] table
#'
#' The clade is the smallest one that holds every species of the row.
#' @noRd
.cc_with_clade <- function(cliques, clades, species) {
  sp <- intersect(species, names(cliques))
  present <- !is.na(as.matrix(cliques[, sp, drop = FALSE]))
  cliques$clade <- vapply(seq_len(nrow(cliques)), function(i) {
    .clade_home(clades, sp[present[i, ]])
  }, character(1))
  cliques
}
