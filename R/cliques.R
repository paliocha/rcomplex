#' Encode clique edge data as 0-based integer vectors for C++
#'
#' Shared helper for [find_cliques()] and [clique_stability()]. Builds
#' string-to-int maps, converts edge columns to 0-based integer vectors,
#' and filters out edges where either species is not in `target_species`.
#'
#' @param edges Data frame with columns gene1, gene2, species1, species2,
#'   hog, q.value, effect_size (already type-filtered).
#' @param target_species Character vector of species abbreviations.
#' @return A list with components: sp_map, gene_map, hog_map, all_genes,
#'   unique_hogs, edge_hog, edge_g1, edge_g2, edge_sp1, edge_sp2,
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
  edge_sp1 <- as.integer(sp_map[edges$species1])
  edge_sp2 <- as.integer(sp_map[edges$species2])

  # Filter out edges where either species is not in target_species
  valid <- !is.na(edge_sp1) & !is.na(edge_sp2)
  any_valid <- any(valid)

  edge_hog <- edge_hog[valid]
  edge_g1 <- edge_g1[valid]
  edge_g2 <- edge_g2[valid]
  edge_sp1 <- edge_sp1[valid]
  edge_sp2 <- edge_sp2[valid]
  edge_qval <- as.numeric(edges$q.value[valid])
  edge_effect <- as.numeric(edges$effect_size[valid])

  list(
    sp_map = sp_map, gene_map = gene_key_map, hog_map = hog_map,
    all_genes = all_genes, unique_hogs = unique_hogs,
    edge_hog = edge_hog, edge_g1 = edge_g1, edge_g2 = edge_g2,
    edge_sp1 = edge_sp1, edge_sp2 = edge_sp2,
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


#' Warn when an edge table looks already cut to `edge_type`
#'
#' Clique intensity fits each edge's weight scale on every tested pair of
#' its species pair. A table holding only `edge_type` rows fits that scale
#' on significant edges alone instead, which is a different quantity that
#' nothing downstream can tell apart. Tables without a
#' `type` column cannot be judged and pass silently. Warns once per
#' session: clique_stability(), clique_threshold_sweep(),
#' clique_perturbation_test() and classify_cliques() call find_cliques()
#' many times on one table.
#'
#' @param edges Edge data frame.
#' @param edge_type Edge types find_cliques() keeps.
#' @noRd
.warn_if_prefiltered <- function(edges, edge_type) {
  if (!"type" %in% names(edges) || nrow(edges) == 0L ||
        !all(edges$type %in% edge_type)) {
    return(invisible(FALSE))
  }
  rlang::warn(
    c(
      paste0(
        "`edges` holds only `edge_type` rows (",
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
#' `1 - q.value`: a clique only ever contains edges that passed alpha, so
#' `1 - q` sat above `1 - alpha` on every edge and left intensity flat
#' (#11). They are built from association strength rather than the
#' Jaccard index, whose null expectation grows with neighbourhood size
#' (#15).
#'
#' @param cliques Data frame from find_cliques (with hog + species columns).
#' @param edges Data frame with gene1, gene2, hog, q.value, effect_size.
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
#' When \code{max_missing_edges > 0}, tolerates up to that many missing
#' species-pair edges per clique. Instead of Bron-Kerbosch (which only
#' finds fully connected subgraphs), all species subsets with at most
#' \code{max_missing_edges} missing edges are enumerated. Assignments
#' prefer fewer missing edges, then lower mean q-value.
#'
#' @param edges Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (first gene in pair)}
#'     \item{gene2}{Gene identifier (second gene in pair)}
#'     \item{species1}{Species for gene1}
#'     \item{species2}{Species for gene2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{q.value}{q-value for the edge (from pair-level testing)}
#'     \item{effect_size}{Numeric effect size}
#'   }
#'   Optionally includes a \code{type} column for filtering. Pass the
#'   unfiltered table (every tested pair, e.g. the output of
#'   \code{\link{find_coexpressologs}}): cliques are built from
#'   \code{edge_type} rows only, but \code{intensity} fits each edge's
#'   weight scale on every row of its species pair, so a pre-filtered
#'   table changes what it measures.
#'   A table whose \code{type} column holds only \code{edge_type} rows
#'   triggers a warning (class \code{rcomplex_prefiltered_edges}, shown
#'   once per session).
#' @param target_species Character vector of species abbreviations.
#' @param min_species Minimum number of species per clique
#'   (default: \code{length(target_species)}).
#' @param max_genes_per_sp Maximum genes considered per species per HOG
#'   (default 10). Keeps the most-connected genes.
#' @param max_missing_edges Maximum number of missing species-pair edges
#'   tolerated per clique (default 0 = all edges required). When > 0,
#'   uses subset enumeration instead of Bron-Kerbosch. Practical limit:
#'   ~25 species; falls back to exact BK for larger species counts.
#' @param edge_type If \code{edges} has a \code{type} column, keep only
#'   edges with \code{type \%in\% edge_type} (default \code{"conserved"}).
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
#'     \item{n_missing}{Number of missing edges (0 when
#'       \code{max_missing_edges = 0})}
#'     \item{intensity}{Onnela intensity: geometric mean, across present
#'       edges, of each edge's ensemble connection probability (in
#'       (0, 1); higher = stronger conservation). Each edge's
#'       association strength (\code{effect_size}, observed overlap over
#'       its expectation) is mapped to a probability
#'       \code{p = z w / (1 + z w)}, with \code{z} fitted per species
#'       pair by maximum entropy, so intensity reads as the per-edge
#'       probability that the whole clique exists. Not \code{1 - q.value}:
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
#'
#' # Tolerate 1 missing edge (e.g., 5 of 6 edges in a 4-species clique)
#' tolerant <- find_cliques(edges, target_species, max_missing_edges = 1L)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
find_cliques <- function(edges, ...) UseMethod("find_cliques")

#' @rdname find_cliques
#' @export
find_cliques.default <- function(edges, target_species,
                                 min_species = length(target_species),
                                 max_genes_per_sp = 10L,
                                 max_missing_edges = 0L,
                                 edge_type = "conserved",
                                 cost_weights = c(q = 1.0, effect = 0.0), ...) {
  # Validate inputs
  required_cols <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value", "effect_size"
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
  # pair, so they are taken before the edge_type filter; fitted among
  # conserved edges alone they would only describe edges that already
  # passed alpha.
  .warn_if_prefiltered(edges, edge_type)
  weights <- .onnela_weight(edges)
  if ("type" %in% names(edges)) {
    keep <- edges$type %in% edge_type
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
    enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
    enc$edge_qval, enc$edge_effect,
    length(target_species), min_species,
    length(enc$unique_hogs), length(enc$all_genes),
    as.integer(max_genes_per_sp), as.integer(max_missing_edges),
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


#' Leave-k-out jackknife structural stability for cliques
#'
#' Tests how structurally robust each clique is to species removal.
#' The analysis removes k = 1, 2, \ldots, \code{max_k} species at a time,
#' re-runs full clique detection on each reduced species set, and checks
#' whether matching gene assignments are preserved (Jaccard similarity).
#' ALL cliques are tested, regardless of trait composition. Trait
#' annotations are added post-hoc if \code{species_trait} is provided.
#'
#' @param edges Data frame (same format as \code{\link{find_cliques}}).
#' @param target_species Character vector of species that define clique
#'   membership.
#' @param species_trait Optional named character or factor vector mapping
#'   species to trait groups. Names must include all \code{all_species}.
#'   If provided, trait annotations (\code{traits}, \code{sole_rep}) are
#'   added to the output. If \code{NULL} (default), trait columns are
#'   \code{NA}.
#' @param all_species Character vector of ALL species in the analysis
#'   universe (default: \code{target_species}). Leave-k-out subsets are
#'   drawn from \code{all_species}. Must be a superset of
#'   \code{target_species}. When larger than \code{target_species},
#'   removing a non-target species tests whether the clique signal is
#'   robust to changes in the broader phylogenetic context.
#' @param full_cliques Output of \code{\link{find_cliques}}, or \code{NULL} to
#'   compute internally (default).
#' @param min_species Minimum species per clique in the full dataset
#'   (default: \code{length(target_species)}). During leave-k-out, reduced
#'   cliques require only 2 species (the minimum meaningful clique size);
#'   a clique is testable if at least 2 of its species remain active.
#' @param max_k Maximum number of species to leave out
#'   (default: \code{length(all_species) - 2}, leaving at least 2 species).
#' @param max_genes_per_sp Maximum genes per species per HOG (default 10).
#' @param jaccard_threshold Minimum Jaccard similarity for matching
#'   reduced-dataset cliques to full-dataset cliques (default 0.8).
#' @param edge_type Edge type filter (default \code{"conserved"}).
#' @param n_cores Number of OpenMP threads (default 1).
#' @param cost_weights Cost weights for gene-assignment ranking
#'   (same as in \code{\link{find_cliques}}).
#'   Default \code{c(q = 1, effect = 0)}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{stability}{Data frame with columns: \code{clique_idx} (1-based),
#'       \code{hog}, \code{k}, \code{n_subsets}, \code{n_stable},
#'       \code{stability_score}, \code{species_present}, \code{traits},
#'       \code{sole_rep}. One row per (clique, k) pair.}
#'     \item{clique_disruption}{Data frame with columns: \code{species},
#'       \code{n_cliques_disrupted} (k=1 only), and optionally
#'       \code{trait_value} if \code{species_trait} is provided.
#'       One row per species in \code{all_species}.}
#'     \item{stability_class}{Integer vector (length = number of cliques):
#'       highest k at which each clique is structurally stable across ALL
#'       subsets (0 = unstable at k=1)}
#'     \item{novel_cliques}{Integer: total count of novel cliques across
#'       all subsets}
#'   }
#'
#' @details
#' ## How it works
#'
#' For each combination of k species removed from \code{all_species}
#' (k = 1, ..., max_k):
#' \enumerate{
#'   \item All edges involving the removed species are dropped
#'   \item Clique detection re-runs among remaining target species
#'   \item Reduced cliques are matched to full-dataset cliques by Jaccard
#'     similarity of gene assignments (considering only non-removed species)
#' }
#'
#' A clique is \emph{testable} in a subset if at least 2 of its species
#' remain active. It is \emph{stable at level k} if its gene assignments
#' are preserved across ALL C(N, k) subsets where it is testable.
#'
#' Stability is purely structural — trait annotations are orthogonal and
#' added post-hoc from \code{species_trait} if provided.
#'
#' ## sole_rep column
#'
#' The \code{sole_rep} column is \code{TRUE} if this clique has a single
#' trait value and that trait value has only one target species
#' representative. Only populated when \code{species_trait} is provided.
#'
#' ## full_cliques parameter
#'
#' When \code{full_cliques = NULL} (default), cliques are computed internally
#' via \code{\link{find_cliques}}. You can also precompute them:
#' \preformatted{
#' fc <- find_cliques(edges, target_species)
#' stab <- clique_stability(edges, target_species,
#'                          full_cliques = fc)
#' }
#'
#' @examples
#' \dontrun{
#' # Structural stability for all cliques
#' cliques <- find_cliques(edges, all_sp, min_species = 2L)
#' stab <- clique_stability(edges, all_sp,
#'   full_cliques = cliques
#' )
#'
#' # With trait annotation
#' trait <- setNames(rep(c("annual", "perennial"), each = 4), all_sp)
#' stab <- clique_stability(edges, all_sp,
#'   species_trait = trait,
#'   full_cliques = cliques
#' )
#'
#' # Cliques surviving any single species dropout
#' k1 <- stab$stability[stab$stability$k == 1, ]
#' stable_cliques <- cliques[k1$clique_idx[k1$stability_score == 1], ]
#'
#' # Multi-level: stability_class >= 2 survives any pair of dropouts
#' deeply_stable <- which(stab$stability_class >= 2)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
clique_stability <- function(edges, ...) UseMethod("clique_stability")

#' @rdname clique_stability
#' @export
clique_stability.default <- function(
  edges, target_species,
  species_trait = NULL,
  all_species = target_species,
  full_cliques = NULL,
  min_species = length(target_species),
  max_k = length(all_species) - 2L,
  max_genes_per_sp = 10L,
  jaccard_threshold = 0.8,
  edge_type = "conserved", n_cores = 1L,
  cost_weights = c(q = 1.0, effect = 0.0), ...
) {
  # Validate inputs
  required_cols <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value", "effect_size"
  )
  missing <- setdiff(required_cols, names(edges))
  if (length(missing) > 0) {
    stop("edges missing required columns: ", paste(missing, collapse = ", "))
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  if (!all(target_species %in% all_species)) {
    stop("target_species must be a subset of all_species")
  }
  if (!is.null(species_trait)) {
    if (!is.character(species_trait) && !is.factor(species_trait)) {
      stop("species_trait must be a named character or factor vector")
    }
    if (is.null(names(species_trait))) {
      stop("species_trait must be a named vector with species as names")
    }
    missing_sp <- setdiff(all_species, names(species_trait))
    if (length(missing_sp) > 0) {
      stop(
        "species_trait missing entries for: ",
        paste(missing_sp, collapse = ", ")
      )
    }
  }
  max_k <- as.integer(max_k)
  if (max_k < 1L) {
    stop("max_k must be >= 1")
  }
  if (max_k >= length(all_species)) {
    stop("max_k must be < length(all_species)")
  }

  # Build is_target: 1 for target species, 0 for non-target
  is_target <- as.integer(all_species %in% target_species)

  # Empty result template
  empty_stability <- data.frame(
    clique_idx = integer(0), hog = character(0),
    k = integer(0), n_subsets = integer(0), n_stable = integer(0),
    stability_score = numeric(0), species_present = character(0),
    traits = character(0), sole_rep = logical(0)
  )
  empty_disruption <- data.frame(
    species = character(0),
    n_cliques_disrupted = integer(0)
  )
  empty_result <- list(
    stability = empty_stability,
    clique_disruption = empty_disruption,
    stability_class = integer(0),
    novel_cliques = 0L
  )

  # Filter edges by type if applicable. find_cliques() below gets the
  # unfiltered table instead: it applies the same filter itself, ranks the
  # intensity weights over every tested pair, and would otherwise read the
  # filtered copy as a caller's pre-filtered table and warn.
  edges_all <- edges
  if ("type" %in% names(edges)) {
    edges <- edges[edges$type %in% edge_type, , drop = FALSE]
  }
  if (nrow(edges) == 0) {
    return(empty_result)
  }

  # Encode edges with ALL species (full universe)
  enc <- encode_clique_edges(edges, all_species)
  if (!enc$any_valid) {
    return(empty_result)
  }

  # Compute full cliques if not provided
  if (is.null(full_cliques)) {
    full_cliques <- find_cliques(edges_all, target_species,
      min_species = min_species,
      max_genes_per_sp = max_genes_per_sp,
      edge_type = edge_type,
      cost_weights = cost_weights
    )
  }
  if (nrow(full_cliques) == 0) {
    return(empty_result)
  }

  # Re-encode full_cliques into all_species index space
  fc_hog_idx <- as.integer(enc$hog_map[as.character(full_cliques$hog)])
  n_fc <- nrow(full_cliques)
  fc_genes <- matrix(NA_integer_, nrow = n_fc, ncol = length(all_species))
  for (j in seq_along(all_species)) {
    sp <- all_species[j]
    if (sp %in% names(full_cliques)) {
      gnames <- full_cliques[[sp]]
      present <- !is.na(gnames)
      if (any(present)) {
        # gene_map is keyed by "species\x02gene" (see encode_clique_edges()):
        # gene identifiers are not unique across species, so the lookup must
        # be scoped by `sp`, not by the raw gene name alone.
        fc_genes[present, j] <- as.integer(
          enc$gene_map[paste(sp, gnames[present], sep = "\x02")]
        )
      }
    }
  }
  raw_cliques <- list(
    hog_idx = fc_hog_idx,
    genes = fc_genes,
    n_species = as.integer(full_cliques$n_species),
    mean_q = full_cliques$mean_q,
    max_q = full_cliques$max_q,
    mean_effect_size = full_cliques$mean_effect_size,
    n_edges = as.integer(full_cliques$n_edges)
  )

  # Call C++ stability function (trait-agnostic)
  cpp_result <- find_cliques_stability_cpp(
    enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
    enc$edge_qval, enc$edge_effect,
    length(all_species),
    length(enc$unique_hogs), length(enc$all_genes),
    is_target, raw_cliques,
    as.integer(max_k), as.integer(max_genes_per_sp),
    jaccard_threshold, n_cores,
    as.double(cost_weights[["q"]]), as.double(cost_weights[["effect"]])
  )

  # Post-process: stability data frame
  stab <- cpp_result$stability
  if (nrow(stab) > 0) {
    stab$clique_idx <- stab$clique_idx + 1L
    stab$hog <- enc$unique_hogs[raw_cliques$hog_idx[stab$clique_idx] + 1L]
  } else {
    stab$hog <- character(0)
  }

  # Post-hoc trait annotation per clique (vectorized via presence matrix)
  sp_cols <- intersect(target_species, names(full_cliques))
  present_mat <- !is.na(full_cliques[, sp_cols, drop = FALSE])

  sp_present <- apply(present_mat, 1L, function(row) {
    paste(sp_cols[row], collapse = ",")
  })

  if (!is.null(species_trait)) {
    trait_char <- as.character(species_trait[sp_cols])
    trait_counts <- table(as.character(species_trait[target_species]))

    trait_annot <- vapply(seq_len(n_fc), \(i) {
      tv <- trait_char[present_mat[i, ]]
      if (length(tv) == 0L) {
        NA_character_
      } else {
        paste(sort(unique(tv)), collapse = ",")
      }
    }, character(1))

    sole_rep_vec <- vapply(seq_len(n_fc), \(i) {
      tv <- unique(trait_char[present_mat[i, ]])
      if (length(tv) == 1L) unname(trait_counts[tv]) == 1L else FALSE
    }, logical(1))
  } else {
    trait_annot <- rep(NA_character_, n_fc)
    sole_rep_vec <- rep(NA, n_fc)
  }

  # Merge annotations onto stability data frame by clique_idx
  if (nrow(stab) > 0) {
    stab$species_present <- sp_present[stab$clique_idx]
    stab$traits <- trait_annot[stab$clique_idx]
    stab$sole_rep <- sole_rep_vec[stab$clique_idx]
  } else {
    stab$species_present <- character(0)
    stab$traits <- character(0)
    stab$sole_rep <- logical(0)
  }

  # Post-process: clique_disruption (one row per all_species)
  disrupt <- cpp_result$clique_disruption
  if (nrow(disrupt) > 0) {
    disrupt$species <- all_species[disrupt$species_idx + 1L]
    cols <- c("species", "n_cliques_disrupted")
    if (!is.null(species_trait)) {
      disrupt$trait_value <- as.character(species_trait[disrupt$species])
      cols <- c("species", "trait_value", "n_cliques_disrupted")
    }
    disrupt <- disrupt[, cols, drop = FALSE]
  } else {
    disrupt <- empty_disruption
  }

  sc <- cpp_result$stability_class

  list(
    stability = stab,
    clique_disruption = disrupt,
    stability_class = sc,
    novel_cliques = cpp_result$novel_cliques
  )
}


#' Compute co-expressolog persistence scores for cliques
#'
#' For each clique, measures how robust the conservation signal is to
#' threshold tightening. Co-expressologs are genes that are co-expression
#' neighbours of the clique gene in both species (connected through the
#' ortholog mapping). The persistence score is the ratio of the weakest
#' co-expressolog edge's MR value to its species threshold.
#'
#' A persistence of 1.0 means the weakest co-expressolog edge is exactly
#' at threshold (marginal). Values above 1.0 indicate the conservation
#' signal survives at stricter density thresholds.
#'
#' @param cliques Output of \code{\link{find_cliques}} (data frame with
#'   \code{hog}, one column per species, and summary statistics).
#' @param target_species Character vector matching column names in
#'   \code{cliques} and names in \code{networks}.
#' @param networks Named list of \code{\link{compute_network}} outputs
#'   keyed by species abbreviation. Each element must have \code{$network}
#'   (named numeric matrix or \code{dgCMatrix}) and \code{$threshold}
#'   (scalar).
#' @param edges Data frame with columns \code{gene1}, \code{gene2},
#'   \code{species1}, \code{species2} (same format as
#'   \code{\link{find_cliques}} input). Used solely as the ortholog
#'   mapping to identify co-expressologs; \code{q.value},
#'   \code{effect_size}, and \code{type} columns are ignored. Should
#'   include all comparison pairs, not just conserved edges.
#'
#' @return The input \code{cliques} data frame with two appended columns:
#'   \describe{
#'     \item{persistence}{Minimum co-expressolog ratio across all
#'       species pairs. The weakest shared-neighbour edge.}
#'     \item{mean_persistence}{Mean co-expressolog ratio. Overall
#'       robustness of the conservation signal.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- clique_persistence(cliques, target_species, networks, edges)
#' result[result$persistence > 2.0, ] # robust to 2x threshold tightening
#' }
#'
#' @export
clique_persistence <- function(cliques, target_species, networks, edges) {
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop(
      "cliques must be a data frame from find_cliques(); got ",
      paste(class(cliques), collapse = "/")
    )
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop(
      "networks must be a named list keyed by species; got ",
      paste(class(networks), collapse = "/")
    )
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  missing_sp <- setdiff(target_species, names(cliques))
  if (length(missing_sp) > 0) {
    stop(
      "cliques missing columns for species: ",
      paste(missing_sp, collapse = ", ")
    )
  }
  missing_net <- setdiff(target_species, names(networks))
  if (length(missing_net) > 0) {
    stop(
      "networks missing entries for species: ",
      paste(missing_net, collapse = ", ")
    )
  }
  for (sp in target_species) {
    .net_check(networks[[sp]], networks[[sp]]$threshold)
  }
  ortho_cols <- c("gene1", "gene2", "species1", "species2")
  missing_cols <- setdiff(ortho_cols, names(edges))
  if (length(missing_cols) > 0) {
    stop(
      "edges missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Pre-build ortholog lookup: ortho[["sp_a.sp_b"]][[gene_a]] -> c(gene_b, ...)
  ortho <- list()
  for (pair in utils::combn(target_species, 2, simplify = FALSE)) {
    sp_a <- pair[1]
    sp_b <- pair[2]
    fwd <- edges[
      edges$species1 == sp_a & edges$species2 == sp_b, , drop = FALSE
    ]
    rev <- edges[
      edges$species1 == sp_b & edges$species2 == sp_a, , drop = FALSE
    ]
    gene_a <- c(fwd$gene1, rev$gene2)
    gene_b <- c(fwd$gene2, rev$gene1)
    if (length(gene_a) == 0L) next
    ortho[[paste(sp_a, sp_b, sep = ".")]] <-
      lapply(split(gene_b, gene_a), unique)
  }

  n <- nrow(cliques)
  persistence <- rep(NA_real_, n)
  mean_persistence <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    ratios <- numeric(0)

    # Active species and genes for this clique
    active_sp <- character(0)
    active_g <- character(0)
    for (sp in target_species) {
      g <- cliques[[sp]][i]
      if (!is.na(g)) {
        active_sp <- c(active_sp, sp)
        active_g <- c(active_g, g)
      }
    }
    names(active_g) <- active_sp
    if (length(active_sp) < 2L) next

    # For each species pair in the clique
    for (a_idx in seq_along(active_sp)[-length(active_sp)]) {
      for (b_idx in (a_idx + 1L):length(active_sp)) {
        sp_a <- active_sp[a_idx]
        sp_b <- active_sp[b_idx]
        g_a <- active_g[[sp_a]]
        g_b <- active_g[[sp_b]]

        net_a <- networks[[sp_a]]$network
        thr_a <- networks[[sp_a]]$threshold
        net_b <- networks[[sp_b]]$network
        thr_b <- networks[[sp_b]]$threshold
        if (!g_a %in% rownames(net_a) || !g_b %in% rownames(net_b)) next

        # Neighbours (excluding self); column access returns a named
        # numeric vector for both dense and dgCMatrix networks
        mr_a <- net_a[, g_a]
        mr_a <- mr_a[names(mr_a) != g_a]
        neigh_a <- names(mr_a[mr_a >= thr_a])
        if (length(neigh_a) == 0L) next

        mr_b <- net_b[, g_b]
        mr_b <- mr_b[names(mr_b) != g_b]
        neigh_b <- names(mr_b[mr_b >= thr_b])
        if (length(neigh_b) == 0L) next

        # Map neighbours of g_a through orthologs to sp_b
        lookup <- ortho[[paste(sp_a, sp_b, sep = ".")]]
        if (is.null(lookup)) next

        mapped <- lookup[neigh_a]
        mapped <- mapped[!vapply(mapped, is.null, logical(1))]
        if (length(mapped) == 0L) next

        gene_a_vec <- rep(names(mapped), lengths(mapped))
        gene_b_vec <- unlist(mapped, use.names = FALSE)

        # Filter to co-expressologs (ortholog is also a neighbour)
        is_coexpr <- gene_b_vec %in% neigh_b
        if (!any(is_coexpr)) next

        # Vectorised ratio: min(MR_a/thr_a, MR_b/thr_b) per pair
        pair_ratios <- pmin(
          mr_a[gene_a_vec[is_coexpr]] / thr_a,
          mr_b[gene_b_vec[is_coexpr]] / thr_b
        )
        ratios <- c(ratios, pair_ratios)
      }
    }

    if (length(ratios) > 0L) {
      persistence[i] <- min(ratios)
      mean_persistence[i] <- mean(ratios)
    }
  }

  cliques$persistence <- persistence
  cliques$mean_persistence <- mean_persistence
  cliques
}


#' Structural survival of cliques across stricter density thresholds
#'
#' Convenience wrapper that re-runs the full comparison-to-clique pipeline
#' (\code{\link{compare_neighborhoods}} -> \code{\link{summarize_comparison}}
#' -> \code{\link{comparison_to_edges}} -> \code{\link{find_cliques}}) at
#' progressively stricter thresholds. For custom threshold logic, call the
#' individual functions directly.
#'
#' No permutations are involved -- all tests are analytical
#' (hypergeometric + Storey q-values; the pi0 method is pinned to
#' \code{"storey"}, so the sweep is deterministic and does not consume
#' the global RNG).
#'
#' @param cliques Baseline output of \code{\link{find_cliques}}.
#' @param target_species Character vector of species abbreviations.
#' @param networks Named list of \code{\link{compute_network}} outputs,
#'   keyed by species abbreviation. Each element must have \code{$network}
#'   (named numeric matrix) and \code{$threshold} (scalar).
#' @param orthologs Data frame with columns \code{Species1}, \code{Species2},
#'   \code{hog} (output of \code{\link{parse_orthologs}} or
#'   \code{\link{extract_orthologs}}).
#' @param species_pairs Optional list of length-2 character vectors
#'   specifying which species pairs to compare. Defaults to all
#'   \code{combn(target_species, 2)}.
#' @param multipliers Numeric vector of threshold multipliers (each > 1).
#'   Default \code{c(1.5, 2, 3, 5, 10)}.
#' @param alternative Passed to \code{\link{summarize_comparison}} and
#'   \code{\link{comparison_to_edges}}: \code{"greater"} (default) or
#'   \code{"less"}.
#' @param alpha Significance threshold for edge classification
#'   (default 0.1).
#' @param min_species Minimum species per clique
#'   (default \code{length(target_species)}).
#' @param max_genes_per_sp Maximum genes per species per HOG (default 10).
#' @param max_missing_edges Passed to \code{\link{find_cliques}}
#'   (default 0).
#' @param edge_type Edge type filter (default \code{"conserved"}).
#' @param jaccard_threshold Minimum per-species-slot Jaccard similarity
#'   for a clique to count as "survived" (default 0.5).
#' @param n_cores Cores for \code{\link{compare_neighborhoods}}
#'   (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{survival}{Data frame with one row per (baseline clique,
#'       multiplier), including \code{multiplier = 1.0} baseline rows
#'       where all cliques trivially survive: \code{clique_idx} (1-based),
#'       \code{hog}, \code{multiplier}, \code{survived} (logical),
#'       \code{jaccard}, \code{n_species_orig}, \code{n_species_new}.}
#'     \item{sweep_cliques}{Named list of \code{find_cliques()} outputs
#'       keyed by multiplier.}
#'     \item{sweep_edges}{Named list of combined edge data frames keyed
#'       by multiplier.}
#'     \item{persistence}{Data frame with one row per baseline clique:
#'       \code{clique_idx}, \code{hog}, \code{birth} (lowest multiplier
#'       where the clique exists; 1.0 for all baseline cliques),
#'       \code{death} (first multiplier above birth where the clique is
#'       lost; \code{NA} if it survives all tested multipliers),
#'       \code{persistence} (\code{death - birth}; \code{NA} if death
#'       is \code{NA}).}
#'   }
#'
#' @examples
#' \dontrun{
#' sweep <- clique_threshold_sweep(cliques, target_species, networks,
#'   orthologs,
#'   multipliers = c(1.5, 2, 5)
#' )
#' # Survival curve
#' sapply(
#'   sort(unique(sweep$survival$multiplier)),
#'   function(m) mean(sweep$survival$survived[sweep$survival$multiplier == m])
#' )
#' }
#'
#' @export
clique_threshold_sweep <- function(
  cliques, target_species, networks, orthologs,
  species_pairs = NULL,
  multipliers = c(1.5, 2, 3, 5, 10),
  alternative = c("greater", "less"),
  alpha = 0.1,
  min_species = length(target_species),
  max_genes_per_sp = 10L,
  max_missing_edges = 0L,
  edge_type = "conserved",
  jaccard_threshold = 0.5,
  n_cores = 1L
) {
  alternative <- match.arg(alternative)

  # --- Validation ---
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques()")
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  missing_net <- setdiff(target_species, names(networks))
  if (length(missing_net) > 0) {
    stop(
      "networks missing entries for: ",
      paste(missing_net, collapse = ", ")
    )
  }
  for (sp in target_species) {
    .net_check(networks[[sp]], networks[[sp]]$threshold)
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  if (length(multipliers) == 0) {
    return(list(
      survival = data.frame(
        clique_idx = integer(0), hog = character(0),
        multiplier = numeric(0), survived = logical(0),
        jaccard = numeric(0), n_species_orig = integer(0),
        n_species_new = integer(0)
      ),
      sweep_cliques = list(),
      sweep_edges = list(),
      persistence = data.frame(
        clique_idx = integer(0), hog = character(0),
        birth = numeric(0), death = numeric(0),
        persistence = numeric(0)
      )
    ))
  }

  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(target_species, 2, simplify = FALSE)
  }

  sweep_cliques <- list()
  sweep_edges <- list()
  survival_rows <- vector("list", length(multipliers) * nrow(cliques))
  row_idx <- 0L

  for (m in sort(multipliers)) {
    m_key <- as.character(m)
    message("Threshold sweep: multiplier ", m)

    # Re-threshold all networks (shallow copy, R COW avoids matrix dup;
    # modifyList carries store_threshold etc. so the sparse guard works)
    tight_nets <- lapply(networks[target_species], function(net) {
      modifyList(net, list(threshold = net$threshold * m))
    })
    names(tight_nets) <- target_species

    # Pairwise comparisons
    pair_edges <- list()
    for (pair in species_pairs) {
      sp_a <- pair[1]
      sp_b <- pair[2]

      comparison <- tryCatch(
        compare_neighborhoods(
          tight_nets[[sp_a]], tight_nets[[sp_b]],
          orthologs, n_cores
        ),
        error = function(e) {
          warning(
            "Pair ", sp_a, "-", sp_b, " at ", m, "x failed: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      if (is.null(comparison) || nrow(comparison) == 0) next

      summary_res <- tryCatch(
        # pi0_method pinned to "storey": deterministic pre-0.2.0 q-values
        # pinned for determinism across multipliers; pass-through is a
        # P2 hand-off
        summarize_comparison(comparison, alternative, alpha,
          pi0_method = "storey"
        ),
        error = function(e) {
          warning(
            "Pair ", sp_a, "-", sp_b, " q-values at ", m,
            "x failed: ", conditionMessage(e)
          )
          NULL
        }
      )
      if (is.null(summary_res) || nrow(summary_res$results) == 0) next

      edges_df <- comparison_to_edges(
        summary_res$results, sp_a, sp_b,
        alternative, alpha
      )
      pair_edges[[length(pair_edges) + 1L]] <- edges_df
    }

    if (length(pair_edges) == 0) {
      all_edges <- data.frame(
        gene1 = character(0), gene2 = character(0),
        species1 = character(0), species2 = character(0),
        hog = character(0), q.value = numeric(0),
        effect_size = numeric(0), jaccard = numeric(0),
        power = numeric(0), type = character(0)
      )
    } else {
      all_edges <- do.call(rbind, pair_edges)
    }

    sweep_edges[[m_key]] <- all_edges

    # Find cliques at this threshold
    new_cliques <- find_cliques(all_edges, target_species,
      min_species = min_species,
      max_genes_per_sp = max_genes_per_sp,
      max_missing_edges = max_missing_edges,
      edge_type = edge_type
    )
    sweep_cliques[[m_key]] <- new_cliques

    # Match baseline cliques to new cliques
    for (i in seq_len(nrow(cliques))) {
      row_idx <- row_idx + 1L
      baseline_hog <- cliques$hog[i]
      best_jaccard <- NA_real_
      best_n_sp <- NA_integer_

      if (nrow(new_cliques) > 0) {
        candidates <- which(new_cliques$hog == baseline_hog)
        for (j in candidates) {
          jac <- jaccard_clique_match(
            cliques[i, ], new_cliques[j, ],
            target_species
          )
          if (is.na(best_jaccard) || jac > best_jaccard) {
            best_jaccard <- jac
            best_n_sp <- new_cliques$n_species[j]
          }
        }
      }

      survived <- !is.na(best_jaccard) && best_jaccard >= jaccard_threshold
      survival_rows[[row_idx]] <- data.frame(
        clique_idx = i,
        hog = baseline_hog,
        multiplier = m,
        survived = survived,
        jaccard = best_jaccard,
        n_species_orig = cliques$n_species[i],
        n_species_new = if (survived) best_n_sp else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  survival <- do.call(rbind, survival_rows[seq_len(row_idx)])
  if (is.null(survival)) {
    survival <- data.frame(
      clique_idx = integer(0), hog = character(0),
      multiplier = numeric(0), survived = logical(0),
      jaccard = numeric(0), n_species_orig = integer(0),
      n_species_new = integer(0)
    )
  }
  rownames(survival) <- NULL

  # --- Inject multiplier = 1.0 baseline rows ---
  # All baseline cliques trivially survive at their own threshold.
  if (nrow(cliques) > 0) {
    baseline_rows <- data.frame(
      clique_idx = seq_len(nrow(cliques)),
      hog = cliques$hog,
      multiplier = 1.0,
      survived = TRUE,
      jaccard = 1.0,
      n_species_orig = cliques$n_species,
      n_species_new = cliques$n_species,
      stringsAsFactors = FALSE
    )
    survival <- rbind(baseline_rows, survival)
    rownames(survival) <- NULL
  }

  # --- Compute persistence dataframe ---
  all_multipliers <- sort(unique(survival$multiplier))

  if (nrow(survival) > 0 && nrow(cliques) > 0) {
    # Safe: baseline_rows guarantees every clique_idx 1..nrow(cliques) exists
    surv_split <- split(survival, survival$clique_idx)
    persist_list <- vector("list", nrow(cliques))
    for (i in seq_len(nrow(cliques))) {
      ci_surv <- surv_split[[as.character(i)]]
      survived_at <- ci_surv$multiplier[ci_surv$survived]
      birth <- if (length(survived_at) > 0) min(survived_at) else NA_real_

      # death = first multiplier > birth where survived == FALSE
      death <- NA_real_
      if (!is.na(birth)) {
        candidates <- all_multipliers[all_multipliers > birth]
        for (cand in candidates) {
          row_match <- ci_surv[ci_surv$multiplier == cand, , drop = FALSE]
          if (nrow(row_match) > 0 && !row_match$survived[1]) {
            death <- cand
            break
          }
        }
      }

      persist_list[[i]] <- data.frame(
        clique_idx = i,
        hog = cliques$hog[i],
        birth = birth,
        death = death,
        persistence = death - birth,
        stringsAsFactors = FALSE
      )
    }
    persistence <- do.call(rbind, persist_list)
    rownames(persistence) <- NULL
  } else {
    persistence <- data.frame(
      clique_idx = integer(0), hog = character(0),
      birth = numeric(0), death = numeric(0),
      persistence = numeric(0)
    )
  }

  list(
    survival = survival,
    sweep_cliques = sweep_cliques,
    sweep_edges = sweep_edges,
    persistence = persistence
  )
}


#' Per-species-slot Jaccard similarity between two clique rows
#'
#' Compares gene assignments slot by slot across target species.
#' A slot matches if both rows have the same gene for that species.
#' @noRd
jaccard_clique_match <- function(row1, row2, target_species) {
  intersect_n <- 0L
  union_n <- 0L
  for (sp in target_species) {
    g1 <- row1[[sp]]
    g2 <- row2[[sp]]
    has1 <- !is.na(g1)
    has2 <- !is.na(g2)
    if (has1 || has2) {
      union_n <- union_n + 1L
      if (has1 && has2 && g1 == g2) intersect_n <- intersect_n + 1L
    }
  }
  if (union_n == 0L) 0 else intersect_n / union_n
}


#' Test clique robustness to edge-weight perturbation
#'
#' Adds Gaussian noise to MR scores in each species network, re-thresholds,
#' re-runs the comparison-to-cliques pipeline, and measures how often
#' baseline cliques survive. This tests robustness to measurement noise
#' without recomputing correlations from scratch.
#'
#' For sparse networks only the stored entries are perturbed: edges
#' discarded below the store threshold can never be promoted back above
#' the analysis threshold, and the RNG stream differs from the dense path
#' (one draw per stored entry instead of one per upper-triangle cell), so
#' seeded results are not comparable between the two storage classes.
#'
#' @param cliques Output from \code{\link{find_cliques}}.
#' @param target_species Character vector of target species names.
#' @param networks Named list of network objects (each with \code{$network}
#'   matrix and \code{$threshold} scalar).
#' @param orthologs Ortholog table data frame.
#' @param species_pairs Optional list of length-2 character vectors
#'   specifying which pairs to compare. Defaults to all pairwise
#'   combinations of \code{target_species}.
#' @param n_boot Number of bootstrap iterations (default 100).
#' @param noise_sd Standard deviation of Gaussian noise added to MR scores
#'   (default 0.1).
#' @param alternative Direction of test (default \code{"greater"}).
#' @param alpha Significance threshold (default 0.1).
#' @param min_species Minimum species per clique.
#' @param max_genes_per_sp Maximum genes per species per HOG.
#' @param max_missing_edges Maximum missing species-pair edges per clique
#'   (default 0). Passed to \code{\link{find_cliques}}.
#' @param edge_type Edge type filter.
#' @param jaccard_threshold Minimum Jaccard for clique matching (default 0.5).
#' @param n_cores Number of parallel cores.
#' @param seed Random seed, or \code{NULL} (default) to draw from the
#'   ambient RNG stream and leave it advanced. A seed draws from a private
#'   stream and restores the caller's on exit, the package-wide contract
#'   described under \code{\link{detect_modules}}.
#' @param cost_weights Cost weights for \code{\link{find_cliques}}
#'   (default \code{c(q = 1, effect = 0)}).
#' @param pval_combine Directional q-value combination for the edge calls
#'   in every internal rerun, passed to
#'   \code{\link{find_coexpressologs}}. Set it to whatever built the
#'   baseline \code{cliques} (default \code{"max"}): a mismatch deflates
#'   survival rates by the criterion change rather than by noise.
#' @param pi0_method pi0 estimation for every internal rerun, passed to
#'   \code{\link{find_coexpressologs}}. Default \code{"storey"} is
#'   deterministic across bootstrap iterations (\code{"randomized"} adds
#'   per-iteration Monte Carlo noise unrelated to \code{noise_sd}),
#'   matching \code{\link{clique_threshold_sweep}}. Set it to whatever
#'   built the baseline \code{cliques}.
#' @param ... Additional arguments passed to the default method.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{clique_idx}{1-based index into baseline cliques}
#'     \item{hog}{HOG identifier}
#'     \item{survival_rate}{Fraction of total bootstrap iterations where
#'       clique survived (Jaccard >= threshold)}
#'     \item{mean_jaccard}{Mean Jaccard similarity to best-matching
#'       perturbed clique, over matched iterations only (NA if no matches)}
#'     \item{n_boot}{Total number of bootstrap iterations run}
#'     \item{n_matched}{Number of iterations where the clique's HOG
#'       produced a matching clique (Jaccard > 0)}
#'   }
#'
#' @examples
#' \dontrun{
#' pert <- clique_perturbation_test(cliques, target_species, networks,
#'   orthologs,
#'   n_boot = 50, noise_sd = 0.1
#' )
#' # Cliques surviving > 80% of perturbations
#' pert[pert$survival_rate >= 0.8, ]
#' }
#'
#' @export
clique_perturbation_test <- function(cliques, ...) {
  UseMethod("clique_perturbation_test")
}

#' @rdname clique_perturbation_test
#' @export
clique_perturbation_test.default <- function(
  cliques, target_species, networks, orthologs,
  species_pairs = NULL,
  n_boot = 100L,
  noise_sd = 0.1,
  alternative = c("greater", "less"),
  alpha = 0.1,
  min_species = length(target_species),
  max_genes_per_sp = 10L,
  max_missing_edges = 0L,
  edge_type = "conserved",
  jaccard_threshold = 0.5,
  n_cores = 1L,
  seed = NULL,
  cost_weights = c(q = 1.0, effect = 0.0),
  pval_combine = c("max", "min"),
  pi0_method = c("storey", "randomized", "none"), ...
) {
  alternative <- match.arg(alternative)
  pval_combine <- match.arg(pval_combine)
  pi0_method <- match.arg(pi0_method)
  n_boot <- as.integer(n_boot)

  # --- Empty result template ---
  empty_result <- data.frame(
    clique_idx = integer(0),
    hog = character(0),
    survival_rate = numeric(0),
    mean_jaccard = numeric(0),
    n_boot = integer(0),
    n_matched = integer(0),
    stringsAsFactors = FALSE
  )

  # --- Early returns ---
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques()")
  }
  if (nrow(cliques) == 0L || n_boot == 0L) {
    return(empty_result)
  }

  # --- Validation ---
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  missing_net <- setdiff(target_species, names(networks))
  if (length(missing_net) > 0) {
    stop(
      "networks missing entries for: ",
      paste(missing_net, collapse = ", ")
    )
  }
  for (sp in target_species) {
    net <- networks[[sp]]
    if (!is.list(net) || is.null(net$network) || is.null(net$threshold)) {
      stop("each network must have 'network' and 'threshold' elements")
    }
    .net_check(net, net$threshold)
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L || noise_sd < 0) {
    stop("noise_sd must be a non-negative scalar")
  }

  # See .seed_scope() in R/rng.R for the package-wide contract.
  .seed_scope(seed)

  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(target_species, 2, simplify = FALSE)
  }

  n_cliques <- nrow(cliques)

  # Accumulators: per-clique survival count and Jaccard sum
  n_survived <- integer(n_cliques)
  n_matched <- integer(n_cliques)
  sum_jaccard <- numeric(n_cliques)

  for (b in seq_len(n_boot)) {
    # 1. Perturb each species network
    perturbed_networks <- networks
    for (sp in target_species) {
      net <- networks[[sp]]
      if (.net_is_sparse(net)) {
        # Perturb the stored entries only: edges discarded below the store
        # threshold can never be promoted back. One rnorm draw per stored
        # upper-triangle entry, so the RNG stream differs from the dense
        # path (one draw per upper-triangle cell).
        u <- Matrix::triu(net$network)
        u@x <- pmax(0, u@x + stats::rnorm(length(u@x), 0, noise_sd))
        perturbed <- methods::as(u + Matrix::t(u), "generalMatrix")
      } else {
        perturbed <- net$network
        ut <- which(upper.tri(perturbed))
        perturbed[ut] <- perturbed[ut] +
          stats::rnorm(length(ut), 0, noise_sd)
        perturbed[lower.tri(perturbed)] <- t(perturbed)[lower.tri(perturbed)]
        perturbed[perturbed < 0] <- 0
      }
      perturbed_networks[[sp]] <- modifyList(net, list(network = perturbed))
    }

    # 2. Re-run comparison pipeline (analytical = fast)
    edges_b <- tryCatch(
      find_coexpressologs(perturbed_networks, orthologs,
        species_pairs = species_pairs,
        method = "hypergeometric",
        alternative = alternative,
        alpha = alpha, n_cores = n_cores,
        pval_combine = pval_combine,
        pi0_method = pi0_method
      ),
      error = function(e) NULL
    )
    if (is.null(edges_b) || nrow(edges_b) == 0L) next

    # 3. Re-find cliques
    cliques_b <- tryCatch(
      find_cliques(edges_b, target_species,
        min_species = min_species,
        max_genes_per_sp = max_genes_per_sp,
        max_missing_edges = max_missing_edges,
        edge_type = edge_type,
        cost_weights = cost_weights
      ),
      error = function(e) NULL
    )
    if (is.null(cliques_b) || nrow(cliques_b) == 0L) next

    # 4. Match baseline cliques to perturbed cliques
    for (i in seq_len(n_cliques)) {
      best_jaccard <- -1
      candidates <- which(cliques_b$hog == cliques$hog[i])
      for (j in candidates) {
        jac <- jaccard_clique_match(
          cliques[i, ], cliques_b[j, ],
          target_species
        )
        if (jac > best_jaccard) best_jaccard <- jac
      }

      if (best_jaccard >= 0) {
        # Candidate found (sentinel is -1; Jaccard >= 0 means a HOG match)
        n_matched[i] <- n_matched[i] + 1L
        sum_jaccard[i] <- sum_jaccard[i] + best_jaccard
        if (best_jaccard >= jaccard_threshold) {
          n_survived[i] <- n_survived[i] + 1L
        }
      }
    }
  }

  if (all(n_matched == 0L) && n_boot > 0L) {
    warning(
      "no bootstrap iteration produced matching cliques; ",
      "check networks/orthologs compatibility"
    )
  }

  # Build output — mean_jaccard over matched iterations only
  data.frame(
    clique_idx = seq_len(n_cliques),
    hog = cliques$hog,
    survival_rate = n_survived / n_boot,
    mean_jaccard = ifelse(n_matched > 0L, sum_jaccard / n_matched, NA_real_),
    n_boot = rep(n_boot, n_cliques),
    n_matched = n_matched,
    stringsAsFactors = FALSE
  )
}


#' Test clique intensity against a permutation null
#'
#' Permutes the ortholog mapping and builds a null distribution of
#' clique intensity (see \code{\link{find_cliques}}: the geometric mean
#' of each clique edge's ensemble connection probability, scaled within
#' its species pair). Observed and null intensities are computed against
#' the full edge tables (\code{edges}, and every permutation's
#' \code{find_coexpressologs()} output), so pass \code{edges} unfiltered;
#' a supplied \code{edges} holding only \code{edge_type} rows warns once
#' per session (class \code{rcomplex_prefiltered_edges}).
#' By default (\code{null_model = "global"}) the null shuffles
#' \code{Species2} genes across all rows of the ortholog table,
#' destroying both the specific ortholog mapping and the within-HOG
#' gene grouping. Network topology (per-species adjacency matrices)
#' is preserved, but species structure within the ortholog table is
#' not — a gene originally in one HOG/species-pair may land in another.
#' Under \code{null_model = "within_hog"} the shuffle happens inside
#' each HOG \emph{and} species pair, so the HOG grouping survives and
#' only the gene-to-gene pairing is randomised; see \code{null_model}
#' for when that matters. Partitioning by the pair as well matters for a
#' stacked all-pairs ortholog table, where one HOG carries rows from
#' several pairs: shuffling by HOG alone would move a gene into a pair
#' whose network does not contain it, and that row is then dropped.
#' A HOG with one gene per species therefore has nothing to permute and
#' yields a null with no spread.
#' This tests whether the observed gene-to-gene correspondence produces
#' stronger co-expression conservation than random mappings.
#'
#' P-values are empirical over matched permutations only (permutations
#' where the clique's HOG produced a clique with Jaccard > 0 to the
#' baseline). Unmatched permutations are excluded, not counted as zeros.
#'
#' @param cliques Output from \code{\link{find_cliques}}.
#' @param target_species Character vector of target species names.
#' @param networks Named list of network objects (each with \code{$network}
#'   matrix and \code{$threshold} scalar).
#' @param orthologs Ortholog table data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog}.
#' @param species_pairs Optional list of length-2 character vectors
#'   specifying which pairs to compare. Defaults to all pairwise
#'   combinations of \code{target_species}.
#' @param n_perm Number of permutations (default 500).
#' @param alternative Direction of test: \code{"greater"} (default) tests
#'   whether observed intensity exceeds the null; \code{"less"} tests
#'   whether it is below the null.
#' @param alpha Significance threshold (default 0.1).
#' @param min_species Minimum species per clique.
#' @param max_genes_per_sp Maximum genes per species per HOG.
#' @param max_missing_edges Maximum missing species-pair edges per clique
#'   (default 0). Passed to \code{\link{find_cliques}}.
#' @param edge_type Edge type filter.
#' @param n_cores Number of parallel cores.
#' @param seed Random seed, or \code{NULL} (default) to draw from the
#'   ambient RNG stream and leave it advanced. A seed draws from a private
#'   stream and restores the caller's on exit, the package-wide contract
#'   described under \code{\link{detect_modules}}.
#' @param cost_weights Cost weights for \code{\link{find_cliques}}.
#' @param edges Optional edge data frame (output of
#'   \code{\link{find_coexpressologs}}). When provided, skips the
#'   baseline edge recomputation. When \code{NULL} (default), edges
#'   are computed internally.
#' @param match_clique_size For \code{null_model = "matched_edges"},
#'   draw each clique's null from edges belonging to cliques of the
#'   \emph{same size} rather than from the species pair at large.
#'   \code{TRUE} by default because pair-only pooling confounds the
#'   score with clique size: a clique exists only because all of its
#'   edges passed together, so larger cliques are built from stronger
#'   edges and beat a pool dominated by small-clique ones. An edge
#'   belonging to cliques of several sizes backs the null of each.
#'   \code{FALSE} restores the pair-only pooling.
#' @param min_pool_size Smallest pool a matched-edge draw will use;
#'   a clique whose pool is thinner is left unscored
#'   (\code{n_matched = 0}). The default of 1 is off in practice and
#'   preserves the behaviour of a pool that holds a single weight: the
#'   draw is that constant, so \code{null_sd} is 0 and \code{z} is
#'   \code{NA}, which is the right answer rather than a degenerate one.
#'   Raise it only where pools are genuinely thin -- on the
#'   eight-species Pooideae run the size-matched pools hold roughly
#'   350-575 weights, so it never binds there.
#' @param null_model How the null is built. \code{"matched_edges"} does
#'   not permute the ortholog mapping at all: it holds the mapping fixed
#'   and replaces each clique edge with one drawn from that same species
#'   pair's pool of \code{edge_type} edges, matching clique size and
#'   species-pair composition while randomising only which genes supply
#'   the weights. Nothing is re-clustered, so no clique has to be
#'   rebuilt and the null stays defined for single-copy HOGs, where the
#'   permutation models below are degenerate. It needs only \code{edges}
#'   -- \code{networks} and \code{orthologs} may be \code{NULL} -- and
#'   costs one resample per draw instead of a full
#'   \code{find_coexpressologs()} re-run. Its pools treat edges within a
#'   species pair as exchangeable, so it does not control for degree.
#'   The permutation models destroy the ortholog mapping instead.
#'   \code{"global"} (default) shuffles \code{Species2} across every HOG.
#'   That also destroys the HOGs, so a permuted run rarely contains the
#'   observed clique's HOG and there is nothing to match: on the eight-species
#'   Pooideae data every one of 204 cliques returned \code{n_matched = 0},
#'   leaving \code{z_score} and \code{p_value} \code{NA}. Raising
#'   \code{n_perm} does not help, since the match rate is near zero rather
#'   than merely small. \code{"within_hog"} permutes \code{Species2} within
#'   each HOG instead, preserving the HOG's gene multiset, so the clique can
#'   be rebuilt and the null asks whether this particular paralog combination
#'   is unusually intense given those genes.
#' @param pval_combine Directional q-value combination for the edge
#'   calls in the baseline recomputation (when \code{edges = NULL}) and
#'   in every permutation rerun, passed to
#'   \code{\link{find_coexpressologs}}. Set it to whatever built the
#'   baseline \code{cliques}/\code{edges} (default \code{"max"}): a null
#'   built under a stricter criterion than the observed edges is
#'   anti-conservative.
#' @param pi0_method pi0 estimation for the baseline recomputation and
#'   every permutation rerun, passed to
#'   \code{\link{find_coexpressologs}}. Default \code{"storey"} is
#'   deterministic across permutations (\code{"randomized"} adds
#'   per-permutation Monte Carlo noise), matching
#'   \code{\link{clique_threshold_sweep}}. Set it to whatever built the
#'   baseline \code{cliques}/\code{edges}.
#' @param ... Additional arguments passed to the default method.
#'
#' @section Comparing cliques of different sizes:
#' Under \code{null_model = "matched_edges"}, \code{z_score} and
#' \code{p_value} test whether the clique's edges are exchangeable with
#' the pool's, and the null draws its \code{n_edges} edges
#' independently, so \code{null_sd} shrinks as
#' \code{1 / sqrt(n_edges)}. A real clique's edges all belong to one
#' HOG and share its conservation level, so the observed \code{gap}
#' does not shrink: on the eight-species Pooideae run its IQR was 0.067
#' at three species and 0.080 at eight, while the IQR of
#' \code{z_score} grew from 1.6 to 5.6 and the share of |z| > 1.96 from
#' 0.12 to 0.67, in both tails. \code{z_score} therefore measures
#' evidence, and more edges are more evidence for the same effect.
#' Compare \code{z_score} and \code{p_value} within a clique size, and
#' compare \code{gap} across sizes; ranking on \code{z_score} across
#' sizes ranks largely on \code{n_edges}.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{clique_idx}{1-based index into baseline cliques}
#'     \item{hog}{HOG identifier}
#'     \item{observed_intensity}{Intensity from baseline cliques}
#'     \item{null_mean}{Mean intensity over matched permutations}
#'     \item{null_sd}{SD of intensity over matched permutations}
#'     \item{gap}{\code{observed_intensity - null_mean}: the effect on
#'       the weight scale, comparable across clique sizes (see the
#'       section above)}
#'     \item{z_score}{(observed - null_mean) / null_sd}
#'     \item{p_value}{Empirical one-sided p-value over matched
#'       permutations (direction determined by \code{alternative})}
#'     \item{n_perm}{Total number of permutations run}
#'     \item{n_matched}{Number of permutations where clique had a
#'       matching HOG with Jaccard > 0}
#'     \item{n_edges}{Number of rows of \code{edges} joining two clique
#'       members within the clique's HOG: the edges the observed
#'       intensity is taken over and, under \code{"matched_edges"}, the
#'       edges the null resamples. Equals \code{n_edges} of
#'       \code{\link{find_cliques}} when \code{edges} is the table the
#'       cliques were built from.}
#'   }
#'
#' @export
clique_intensity_test <- function(cliques, ...) {
  UseMethod("clique_intensity_test")
}

#' @rdname clique_intensity_test
#' @export
clique_intensity_test.default <- function(
  cliques, target_species, networks = NULL, orthologs = NULL,
  species_pairs = NULL,
  n_perm = 500L,
  alternative = c("greater", "less"),
  alpha = 0.1,
  min_species = length(target_species),
  max_genes_per_sp = 10L,
  max_missing_edges = 0L,
  edge_type = "conserved",
  n_cores = 1L,
  seed = NULL,
  cost_weights = c(q = 1.0, effect = 0.0),
  edges = NULL,
  pval_combine = c("max", "min"),
  pi0_method = c("storey", "randomized", "none"),
  null_model = c("global", "within_hog", "matched_edges"),
  match_clique_size = TRUE,
  min_pool_size = 1L, ...
) {
  alternative <- match.arg(alternative)
  null_model <- match.arg(null_model)
  if (!is.logical(match_clique_size) || length(match_clique_size) != 1L ||
        is.na(match_clique_size)) {
    stop("match_clique_size must be TRUE or FALSE")
  }
  # Validate before coercing: as.integer() truncates, so checking the
  # coerced value would accept 1.5 as 1 and silently change the null
  # rather than reject an argument documented as a positive integer.
  if (!is.numeric(min_pool_size) || length(min_pool_size) != 1L ||
        is.na(min_pool_size) || !is.finite(min_pool_size) ||
        min_pool_size < 1 || min_pool_size != round(min_pool_size)) {
    stop("min_pool_size must be a single positive integer")
  }
  min_pool_size <- as.integer(min_pool_size)
  pval_combine <- match.arg(pval_combine)
  pi0_method <- match.arg(pi0_method)
  n_perm <- as.integer(n_perm)

  empty_result <- data.frame(
    clique_idx = integer(0), hog = character(0),
    observed_intensity = numeric(0), null_mean = numeric(0),
    null_sd = numeric(0), gap = numeric(0), z_score = numeric(0),
    p_value = numeric(0), n_perm = integer(0),
    n_matched = integer(0), n_edges = integer(0),
    stringsAsFactors = FALSE
  )

  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques()")
  }
  if (nrow(cliques) == 0L || n_perm == 0L) {
    return(empty_result)
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  # "matched_edges" never re-runs find_coexpressologs(), so with a
  # supplied edge table it needs neither networks nor orthologs. That is
  # the point: it turns a run that had to hold every species network in
  # memory into one that reads an edge table.
  needs_networks <- is.null(edges) || null_model != "matched_edges"
  if (needs_networks) {
    if (!is.list(networks) || is.null(names(networks))) {
      stop("networks must be a named list keyed by species")
    }
    missing_net <- setdiff(target_species, names(networks))
    if (length(missing_net) > 0) {
      stop(
        "networks missing entries for: ",
        paste(missing_net, collapse = ", ")
      )
    }
    for (sp in target_species) {
      net <- networks[[sp]]
      if (!is.list(net) || is.null(net$network) || is.null(net$threshold)) {
        stop("each network must have 'network' and 'threshold' elements")
      }
      .net_check(net, net$threshold)
    }
    if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
      stop("orthologs must have columns: Species1, Species2, hog")
    }
  }

  # See .seed_scope() in R/rng.R for the package-wide contract.
  .seed_scope(seed)
  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(target_species, 2, simplify = FALSE)
  }

  # Compute baseline intensity. Only a caller-supplied table can have been
  # cut to edge_type; the one built here is every tested pair.
  if (!is.null(edges)) .warn_if_prefiltered(edges, edge_type)
  if (is.null(edges)) {
    edges <- find_coexpressologs(networks, orthologs,
      species_pairs = species_pairs,
      method = "hypergeometric",
      alternative = alternative,
      alpha = alpha, n_cores = n_cores,
      pval_combine = pval_combine,
      pi0_method = pi0_method
    )
  }
  # One Onnela fit for the observed statistic, reused by the
  # matched-edge null below. Each fit is a uniroot per species pair, and
  # the null has to score against the same weights the observation does.
  # The permutation loop still fits its own, since each permuted run has
  # its own edge table.
  ew <- .onnela_weight(edges)
  obs_stats <- compute_clique_edge_stats(
    cliques, edges,
    target_species,
    weights = ew
  )
  observed_intensity <- obs_stats$intensity
  n_cliques <- nrow(cliques)
  rows <- .clique_edge_rows(cliques, edges, target_species)
  n_edges <- lengths(rows)

  # Track null intensities and which entries are real matches
  # (not structural zeros from absent cliques)
  null_intensities <- matrix(NA_real_, nrow = n_perm, ncol = n_cliques)

  # A stacked all-pairs ortholog table puts several species pairs into
  # one HOG, so grouping the shuffle by hog alone can move a species-C
  # gene into an A-B row. compare_neighborhoods() then drops that row by
  # network membership and the mapping is lost silently (measured: 33%
  # of rows on a stacked 3-species table). Group by the species pair as
  # well. Both keys are constant across permutations, so build them once.
  shuffle_sp1 <- shuffle_sp2 <- NULL
  if (null_model == "within_hog") {
    gene_sp <- unlist(lapply(names(networks), function(s) {
      rn <- rownames(networks[[s]]$network)
      if (is.null(rn)) NULL else stats::setNames(rep(s, length(rn)), rn)
    }))
    shuffle_sp1 <- unname(gene_sp[orthologs$Species1])
    shuffle_sp2 <- unname(gene_sp[orthologs$Species2])
    shuffle_sp1[is.na(shuffle_sp1)] <- "?"
    shuffle_sp2[is.na(shuffle_sp2)] <- "?"
  }

  # Hold the ortholog mapping fixed and ask the question the statistic
  # can actually answer: is this clique's edge set unusually intense for
  # its size and species composition? Each draw replaces every clique
  # edge with one drawn from that same species pair's pool, so clique
  # size and pair composition are matched and only the genes are
  # randomised. Nothing is re-clustered, so unlike an orthology
  # permutation this null cannot collapse: no clique has to be rebuilt,
  # and it stays defined for single-copy HOGs.
  perm_iters <- if (null_model == "matched_edges") 0L else n_perm
  if (null_model == "matched_edges") {
    pool_key <- paste(
      pmin(edges$species1, edges$species2),
      pmax(edges$species1, edges$species2),
      sep = "\x01"
    )
    usable <- is.finite(ew)
    if (!is.null(edge_type) && "type" %in% names(edges)) {
      usable <- usable & edges$type %in% edge_type
    }
    # Clique size drives both halves of the statistic. The null's spread
    # shrinks as 1/sqrt(E), and a clique exists only because all of its
    # edges passed together, so larger cliques are assembled from
    # systematically stronger edges: on the eight-species Pooideae root
    # run the median edge weight runs 0.615 at three species to 0.752 at
    # eight, against a pair-only pool median of 0.638. Pooling by species
    # pair alone therefore scores a large clique against mostly
    # small-clique edges and inflates its z with size (median z ran -0.51
    # at three species to +6.32 at eight). Size is read off the row, not
    # an n_species column a caller's table need not carry.
    clique_size <- vapply(seq_len(n_cliques), function(i) {
      v <- unlist(cliques[i, target_species, drop = TRUE], use.names = FALSE)
      sum(!is.na(v))
    }, integer(1))
    if (match_clique_size) {
      # A membership, not an edge, is the unit: an edge can belong to
      # cliques of several sizes (3904 of 8212 HOGs on that run yield
      # more than one clique), and it should back the null of every
      # size class it takes part in.
      mem_idx <- unlist(rows, use.names = FALSE)
      mem_size <- rep(clique_size, lengths(rows))
      keep <- usable[mem_idx]
      pools <- split(
        ew[mem_idx][keep],
        paste(pool_key[mem_idx][keep], mem_size[keep], sep = "\x01")
      )
    } else {
      pools <- split(ew[usable], pool_key[usable])
    }
    for (i in seq_len(n_cliques)) {
      idx <- rows[[i]]
      # Same refusal as the observed statistic: an incomplete weight set
      # would silently change the denominator, so score nothing.
      if (length(idx) == 0L || anyNA(ew[idx])) next
      keys <- if (match_clique_size) {
        paste(pool_key[idx], clique_size[i], sep = "\x01")
      } else {
        pool_key[idx]
      }
      draws <- vapply(keys, function(k) {
        pool <- pools[[k]]
        # A pool this thin gives a null whose spread is a property of
        # the pool, not of the data. Refuse it rather than report an
        # inflated z; the clique surfaces as n_matched = 0.
        if (is.null(pool) || length(pool) < min_pool_size) {
          rep(NA_real_, n_perm)
        } else {
          sample(pool, n_perm, replace = TRUE)
        }
      }, numeric(n_perm))
      if (is.null(dim(draws))) draws <- matrix(draws, nrow = n_perm)
      null_intensities[, i] <- exp(rowMeans(log(draws)))
    }
  }

  for (p in seq_len(perm_iters)) {
    # How the ortholog mapping is destroyed. "global" shuffles Species2
    # across every HOG, which also destroys the HOGs themselves: a
    # permuted run then almost never contains the observed clique's HOG,
    # so there is nothing to match and n_matched stays 0 (measured: 0 of
    # 204 cliques at full scale). "within_hog" keeps each HOG's gene
    # multiset, within one species pair, and permutes only which
    # Species1 gene each Species2 gene is paired with, so the HOG
    # survives every permutation and the null asks the narrower question
    # the statistic needs: given these genes, is this paralog
    # combination unusually intense? A single-copy HOG has nothing to
    # permute, so its null has no spread and its z is NA.
    shuffled_orthologs <- orthologs
    shuffled_orthologs$Species2 <- if (null_model == "within_hog") {
      stats::ave(orthologs$Species2, orthologs$hog, shuffle_sp1,
        shuffle_sp2,
        FUN = function(g) {
          if (length(g) < 2L) g else sample(g)
        }
      )
    } else {
      sample(orthologs$Species2)
    }

    edges_p <- tryCatch(
      find_coexpressologs(networks, shuffled_orthologs,
        species_pairs = species_pairs,
        method = "hypergeometric",
        alternative = alternative,
        alpha = alpha, n_cores = n_cores,
        pval_combine = pval_combine,
        pi0_method = pi0_method
      ),
      error = function(e) NULL
    )
    if (is.null(edges_p) || nrow(edges_p) == 0L) next

    cliques_p <- tryCatch(
      find_cliques(edges_p, target_species,
        min_species = min_species,
        max_genes_per_sp = max_genes_per_sp,
        max_missing_edges = max_missing_edges,
        edge_type = edge_type,
        cost_weights = cost_weights
      ),
      error = function(e) NULL
    )
    if (is.null(cliques_p) || nrow(cliques_p) == 0L) next

    stats_p <- compute_clique_edge_stats(cliques_p, edges_p, target_species)

    for (i in seq_len(n_cliques)) {
      candidates <- which(cliques_p$hog == cliques$hog[i])
      best_jaccard <- -1
      best_intensity <- NA_real_
      for (j in candidates) {
        jac <- jaccard_clique_match(
          cliques[i, ], cliques_p[j, ],
          target_species
        )
        if (jac > best_jaccard) {
          best_jaccard <- jac
          best_intensity <- stats_p$intensity[j]
        }
      }
      if (best_jaccard >= 0 && !is.na(best_intensity)) {
        null_intensities[p, i] <- best_intensity
      }
      # Unmatched permutations (sentinel -1) stay NA
    }
  }

  # Compute statistics over matched permutations only (exclude NAs)
  n_matched <- as.integer(colSums(!is.na(null_intensities)))
  null_mean <- vapply(seq_len(n_cliques), function(i) {
    vals <- null_intensities[, i]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) NA_real_ else mean(vals)
  }, numeric(1))
  null_sd <- vapply(seq_len(n_cliques), function(i) {
    vals <- null_intensities[, i]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2L) NA_real_ else stats::sd(vals)
  }, numeric(1))

  z_score <- ifelse(!is.na(null_sd) & null_sd > 0,
    (observed_intensity - null_mean) / null_sd,
    NA_real_
  )

  # Empirical p-value over matched permutations
  upper_tail <- alternative == "greater"
  p_value <- vapply(seq_len(n_cliques), function(i) {
    vals <- null_intensities[, i]
    vals <- vals[!is.na(vals)]
    n_m <- length(vals)
    if (n_m == 0L || is.na(observed_intensity[i])) {
      return(NA_real_)
    }
    if (upper_tail) {
      (sum(vals >= observed_intensity[i]) + 1) / (n_m + 1)
    } else {
      (sum(vals <= observed_intensity[i]) + 1) / (n_m + 1)
    }
  }, numeric(1))

  data.frame(
    clique_idx = seq_len(n_cliques), hog = cliques$hog,
    observed_intensity = observed_intensity,
    null_mean = null_mean, null_sd = null_sd,
    gap = observed_intensity - null_mean,
    z_score = z_score, p_value = p_value,
    n_perm = rep(n_perm, n_cliques),
    n_matched = n_matched,
    n_edges = as.integer(n_edges),
    stringsAsFactors = FALSE
  )
}


#' Classify HOGs by clique conservation pattern
#'
#' Convenience wrapper that runs \code{\link{find_cliques}} internally
#' (once for all species, once per trait group) and applies a sequential
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
#'   \item \strong{partial}: a clique exists with \code{min_species}
#'     to \code{N-1} species (or with missing edges when
#'     \code{max_missing_edges > 0}).
#'   \item \strong{differentiated}: at least 2 trait groups each have
#'     a within-group clique, but no cross-group conserved edge exists.
#'   \item \strong{trait_specific}: exactly 1 trait group has a
#'     within-group clique.
#'   \item \strong{unclassified}: none of the above.
#' }
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
#' \code{stability_class} / \code{persistence} / \code{robust} columns
#' exist only on this side.
#'
#' \code{\link{classify_gene_cliques}} works on the \emph{gene} graph
#' built by \code{\link{gene_clique_graph}}. It asks which individual
#' gene copies are mutually conserved, so one HOG can yield several
#' overlapping cliques and the answer names paralogs rather than
#' species. It applies the taxonomy of Rodriguez et al. (2026), plus
#' \code{trait_specific}, with two explicit tolerance tiers
#' (\code{partial_significant}
#' for weak wiring, \code{partial_present} for a missing gene), and its
#' \code{lineage} split is an argument rather than the trait vector, so
#' it can be run against a clade partition the trait does not follow.
#'
#' Reach for \code{classify_gene_cliques()} when which copy sits in the
#' conserved core matters, or when the published taxonomy is what has to
#' be reported. Reach for \code{classify_cliques()} for a
#' one-row-per-HOG trait summary wired into the stability and sweep
#' machinery.
#'
#' @param edges Data frame with columns \code{gene1}, \code{gene2},
#'   \code{species1}, \code{species2}, \code{hog}, \code{q.value},
#'   \code{effect_size}, and \code{type}. Must contain ALL edges
#'   (conserved + ns + diverged), not pre-filtered, because the
#'   differentiated check needs to verify absence of cross-group
#'   conserved edges. An optional \code{power} column (from
#'   \code{\link{find_coexpressologs}}) enables the
#'   \code{"underpowered"} class; without it, or where it is \code{NA},
#'   the classification is unchanged.
#' @param target_species Character vector of all species.
#' @param species_trait Named character or factor vector mapping each
#'   species to a trait value (e.g., \code{c(SP_A = "annual",
#'   SP_B = "perennial")}).
#' @param min_species Minimum species for a partial or within-group
#'   clique (default 2).
#' @param max_genes_per_sp Passed to \code{\link{find_cliques}}
#'   (default 10).
#' @param max_missing_edges Passed to \code{\link{find_cliques}} for
#'   partial detection (default 0).
#' @param edge_type Edge types considered conserved (default
#'   \code{"conserved"}).
#' @param stability Optional output of \code{\link{clique_stability}}.
#' @param sweep Optional output of \code{\link{clique_threshold_sweep}}.
#' @param min_stability_class Minimum stability class for the
#'   \code{robust} flag (default 0).
#' @param min_persistence Minimum persistence for the \code{robust}
#'   flag (default 1.0).
#' @param min_power Detection power below which a non-conserved edge is
#'   read as uninformative rather than as evidence against conservation
#'   (default 0.8). Only used when \code{edges} has \code{power}.
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
#'     \item{trait_groups}{Comma-separated trait groups with internal
#'       cliques (NA for complete/partial/unclassified)}
#'     \item{stability_class}{From stability results (NA if not
#'       provided)}
#'     \item{persistence}{Birth/death persistence of the HOG's
#'       best clique (\code{death - birth}) from
#'       \code{\link{clique_threshold_sweep}}; falls back to highest
#'       survived multiplier for legacy sweep output (NA if not
#'       provided)}
#'     \item{robust}{Logical: passes both stability and persistence
#'       thresholds (NA if neither provided)}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- classify_cliques(edges, target_species, species_trait)
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
  edges, target_species, species_trait,
  min_species = 2L,
  max_genes_per_sp = 10L,
  max_missing_edges = 0L,
  edge_type = "conserved",
  stability = NULL,
  sweep = NULL,
  min_stability_class = 0L,
  min_persistence = 1.0,
  min_power = 0.8, ...
) {
  # --- Validation ---
  required_cols <- c(
    "gene1", "gene2", "species1", "species2", "hog",
    "q.value", "effect_size", "type"
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
  if (!is.character(species_trait) && !is.factor(species_trait)) {
    stop("species_trait must be a named character or factor vector")
  }
  if (is.null(names(species_trait))) {
    stop("species_trait must be a named vector")
  }
  missing_sp <- setdiff(target_species, names(species_trait))
  if (length(missing_sp) > 0) {
    stop(
      "species_trait missing entries for: ",
      paste(missing_sp, collapse = ", ")
    )
  }
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

  trait_char <- as.character(species_trait[target_species])
  names(trait_char) <- target_species
  trait_levels <- unique(trait_char)
  n_sp <- length(target_species)
  all_hogs <- unique(edges$hog)

  # Empty result template
  empty <- data.frame(
    hog = character(0), classification = character(0),
    n_species = integer(0), best_mean_q = numeric(0),
    trait_groups = character(0), underpowered = logical(0),
    stability_class = integer(0),
    persistence = numeric(0), robust = logical(0),
    stringsAsFactors = FALSE
  )
  if (length(all_hogs) == 0) {
    return(empty)
  }

  # --- Steps 1+2: Find all cliques (complete + partial in one pass) ---
  all_cliques <- find_cliques(edges, target_species,
    min_species = min_species,
    max_genes_per_sp = max_genes_per_sp,
    max_missing_edges = max_missing_edges,
    edge_type = edge_type
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
  within_group_hogs <- list()

  for (group in trait_levels) {
    group_sp <- names(trait_char[trait_char == group])
    if (length(group_sp) < 2L) {
      within_group_cliques[[group]] <- NULL
      within_group_hogs[[group]] <- character(0)
      next
    }
    wg <- find_cliques(edges, group_sp,
      min_species = min_species,
      max_genes_per_sp = max_genes_per_sp,
      max_missing_edges = max_missing_edges,
      edge_type = edge_type
    )
    within_group_cliques[[group]] <- wg
    within_group_hogs[[group]] <- unique(wg$hog)
  }

  # --- Step 4: Differentiated (2+ groups w/ cliques, no cross-group) ---
  remaining <- setdiff(all_hogs, c(complete_hogs, partial_hogs))

  # Identify cross-group conserved edges
  conserved <- edges[edges$type %in% edge_type, , drop = FALSE]
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

  diff_hogs <- character(0)
  diff_groups <- character(0)
  for (h in remaining) {
    groups_present <- trait_levels[vapply(trait_levels, function(g) {
      h %in% within_group_hogs[[g]]
    }, logical(1))]
    if (length(groups_present) >= 2L && !h %in% hogs_with_cross) {
      diff_hogs <- c(diff_hogs, h)
      diff_groups <- c(diff_groups, paste(groups_present, collapse = ","))
    }
  }

  # --- Step 5: Trait-specific (exactly 1 group has a clique) ---
  remaining2 <- setdiff(remaining, diff_hogs)
  ts_hogs <- character(0)
  ts_groups <- character(0)
  for (h in remaining2) {
    groups_present <- trait_levels[vapply(trait_levels, function(g) {
      h %in% within_group_hogs[[g]]
    }, logical(1))]
    if (length(groups_present) == 1L) {
      ts_hogs <- c(ts_hogs, h)
      ts_groups <- c(ts_groups, groups_present)
    }
  }

  # --- Step 5b: Underpowered specificity / divergence ---
  # Both calls rest on edges that were not conserved. One that could not
  # have been called is no evidence, so the call has to survive reading
  # it as conserved.
  up_hogs <- character(0)
  if ("power" %in% names(edges)) {
    up_hogs <- .cc_underpowered_hogs(
      edges, c(diff_hogs, ts_hogs), within_group_cliques, trait_char,
      edge_type, min_power
    )
  }

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
        best_mean_q = numeric(0)
      ))
    }
    sub <- sub[order(sub$mean_q), , drop = FALSE]
    sub <- sub[!duplicated(sub$hog), , drop = FALSE]
    data.frame(
      hog = sub$hog, n_species = sub$n_species,
      best_mean_q = sub$mean_q, stringsAsFactors = FALSE
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
      trait_groups = NA_character_, stringsAsFactors = FALSE
    )
  }

  # Partial
  if (length(partial_hogs) > 0) {
    info <- best_per_hog(all_cliques, partial_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "partial",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = NA_character_, stringsAsFactors = FALSE
    )
  }

  # Differentiated
  if (length(diff_hogs) > 0) {
    # Best within-group clique for each differentiated HOG
    wg_summary <- do.call(rbind, lapply(
      within_group_cliques[trait_levels],
      function(df) {
        if (!is.null(df) && nrow(df) > 0) {
          df[, c("hog", "n_species", "mean_q"), drop = FALSE]
        } else {
          NULL
        }
      }
    ))
    info <- best_per_hog(wg_summary, diff_hogs)
    tg <- diff_groups[match(info$hog, diff_hogs)]
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog,
      classification = "differentiated",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = tg, stringsAsFactors = FALSE
    )
  }

  # Trait-specific
  if (length(ts_hogs) > 0) {
    wg_summary2 <- do.call(rbind, lapply(
      within_group_cliques[trait_levels],
      function(df) {
        if (!is.null(df) && nrow(df) > 0) {
          df[, c("hog", "n_species", "mean_q"), drop = FALSE]
        } else {
          NULL
        }
      }
    ))
    info <- best_per_hog(wg_summary2, ts_hogs)
    tg <- ts_groups[match(info$hog, ts_hogs)]
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog,
      classification = "trait_specific",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = tg, stringsAsFactors = FALSE
    )
  }

  # Unclassified
  if (length(unclass_hogs) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      hog = unclass_hogs, classification = "unclassified",
      n_species = NA_integer_, best_mean_q = NA_real_,
      trait_groups = NA_character_, stringsAsFactors = FALSE
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
  # call survives and a caller can filter on it.
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

  # --- Sweep annotation ---
  out$persistence <- NA_real_
  if (!is.null(sweep) && "persistence" %in% names(sweep) &&
        nrow(sweep$persistence) > 0) {
    # Use formal birth/death persistence if available
    persist_df <- sweep$persistence
    # Best (max) persistence per HOG across clique indices
    best_persist <- tapply(persist_df$persistence, persist_df$hog, max,
      na.rm = FALSE
    )
    # NA-aware: if all values for a HOG are NA, tapply returns NA
    idx <- match(out$hog, names(best_persist))
    out$persistence[!is.na(idx)] <- best_persist[idx[!is.na(idx)]]
  } else if (!is.null(sweep) && "survival" %in% names(sweep)) {
    # Fallback: compute max survived multiplier from survival dataframe
    surv <- sweep$survival
    if (nrow(surv) > 0) {
      surv_ok <- surv[surv$survived, , drop = FALSE]
      if (nrow(surv_ok) > 0) {
        best_mult <- tapply(surv_ok$multiplier, surv_ok$hog, max)
        idx <- match(out$hog, names(best_mult))
        out$persistence[!is.na(idx)] <- best_mult[idx[!is.na(idx)]]
      }
    }
  }

  # --- Robust flag ---
  has_stab <- !is.null(stability$stability) && nrow(stability$stability) > 0
  has_sweep <- !is.null(sweep) &&
    (("persistence" %in% names(sweep) && nrow(sweep$persistence) > 0) ||
       ("survival" %in% names(sweep) && nrow(sweep$survival) > 0))
  if (has_stab || has_sweep) {
    stab_ok <- if (has_stab) {
      !is.na(out$stability_class) & out$stability_class >= min_stability_class
    } else {
      TRUE
    }
    sweep_ok <- if (has_sweep) {
      !is.na(out$persistence) & out$persistence >= min_persistence
    } else {
      TRUE
    }
    out$robust <- stab_ok & sweep_ok
  } else {
    out$robust <- NA
  }

  out
}


#' HOGs whose specificity or divergence call rests on underpowered edges
#'
#' A deciding edge is a non-conserved row of the HOG with `power` below
#' `min_power`, one endpoint a member of one of the HOG's within-group
#' cliques and the other in a different trait group. Such an edge could
#' not have been called, so the call cannot rule it out as conserved.
#'
#' @param edges Full edge table carrying `power`.
#' @param hogs Candidate HOGs (differentiated and trait-specific).
#' @param wg_cliques Named list of within-group [find_cliques()] tables.
#' @param trait_char Named trait of every target species.
#' @param edge_type,min_power As in [classify_cliques()].
#' @return Character vector of the HOGs to reclassify.
#' @noRd
.cc_underpowered_hogs <- function(edges, hogs, wg_cliques, trait_char,
                                  edge_type, min_power) {
  hogs <- as.character(hogs)
  e_hog <- as.character(edges$hog)
  pw <- as.numeric(edges$power)
  t1 <- unname(trait_char[as.character(edges$species1)])
  t2 <- unname(trait_char[as.character(edges$species2)])
  cand <- e_hog %in% hogs & !(edges$type %in% edge_type) &
    !is.na(pw) & pw < min_power & !is.na(t1) & !is.na(t2) & t1 != t2
  if (!any(cand)) {
    return(character(0))
  }
  sep <- "\x01"
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
