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

  all_genes <- unique(c(edges$gene1, edges$gene2))
  gene_map <- stats::setNames(seq_along(all_genes) - 1L, all_genes)

  unique_hogs <- unique(edges$hog)
  hog_map <- stats::setNames(seq_along(unique_hogs) - 1L,
                              as.character(unique_hogs))

  # Convert to 0-based integer vectors
  edge_hog <- as.integer(hog_map[as.character(edges$hog)])
  edge_g1 <- as.integer(gene_map[edges$gene1])
  edge_g2 <- as.integer(gene_map[edges$gene2])
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

  list(sp_map = sp_map, gene_map = gene_map, hog_map = hog_map,
       all_genes = all_genes, unique_hogs = unique_hogs,
       edge_hog = edge_hog, edge_g1 = edge_g1, edge_g2 = edge_g2,
       edge_sp1 = edge_sp1, edge_sp2 = edge_sp2,
       edge_qval = edge_qval, edge_effect = edge_effect,
       any_valid = any_valid)
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
#'   Optionally includes a \code{type} column for filtering.
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
#' @export
find_cliques <- function(edges, target_species,
                         min_species = length(target_species),
                         max_genes_per_sp = 10L,
                         max_missing_edges = 0L,
                         edge_type = "conserved") {
  # Validate inputs
  required_cols <- c("gene1", "gene2", "species1", "species2", "hog",
                     "q.value", "effect_size")
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

  # Empty result template
  empty_cols <- c(
    list(hog = character(0)),
    stats::setNames(lapply(target_species, \(x) character(0)), target_species),
    list(n_species = integer(0), mean_q = numeric(0), max_q = numeric(0),
         mean_effect_size = numeric(0), n_edges = integer(0),
         n_missing = integer(0))
  )
  empty_result <- as.data.frame(empty_cols)

  # Filter by edge_type if type column exists
  if ("type" %in% names(edges)) {
    edges <- edges[edges$type %in% edge_type, , drop = FALSE]
  }
  if (nrow(edges) == 0) return(empty_result)

  # Encode edges as 0-based integer vectors
  enc <- encode_clique_edges(edges, target_species)
  if (!enc$any_valid) return(empty_result)

  # Call C++
  result <- find_cliques_cpp(
    enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
    enc$edge_qval, enc$edge_effect,
    length(target_species), min_species,
    length(enc$unique_hogs), length(enc$all_genes),
    as.integer(max_genes_per_sp), as.integer(max_missing_edges)
  )

  # Map back to strings
  if (length(result$hog_idx) == 0) return(empty_result)

  # HOG names
  hog_names <- enc$unique_hogs[result$hog_idx + 1L]

  # Gene matrix: map 0-based indices back to gene names
  gene_matrix <- result$genes  # IntegerMatrix (n_cliques x n_species)
  gene_df <- as.data.frame(matrix(NA_character_, nrow = nrow(gene_matrix),
                                  ncol = ncol(gene_matrix)))
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
  out
}


#' Leave-k-out jackknife stability analysis for trait-exclusive cliques
#'
#' Tests how robust each clique's trait exclusivity is to species removal.
#' A clique is \emph{trait-exclusive} if all its species share the same
#' discrete trait value (e.g., all "annual" or all "perennial"). The analysis
#' removes k = 1, 2, \ldots, \code{max_k} species at a time, re-runs full
#' clique detection on each reduced species set, and checks whether
#' exclusivity is preserved.
#'
#' The trait system is generic: any discrete species-level trait with any
#' number of levels. For example, life habit (annual/perennial), climate
#' zone (tropical/temperate/arctic), or ploidy level (2n/4n/6n).
#'
#' @param edges Data frame (same format as \code{\link{find_cliques}}).
#' @param target_species Character vector of species that define clique
#'   membership (e.g., the 4 annuals).
#' @param species_trait Named character or factor vector mapping species to
#'   trait groups. Names must include all \code{all_species}. Example:
#'   \code{c(SP_A = "annual", SP_B = "perennial")}.
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
#'
#' @return A list with components:
#'   \describe{
#'     \item{stability}{Data frame with columns: \code{clique_idx},
#'       \code{hog}, \code{trait_value}, \code{k}, \code{n_subsets},
#'       \code{n_stable}, \code{stability_score}, \code{sole_rep}}
#'     \item{clique_disruption}{Data frame with columns: \code{species},
#'       \code{trait_value}, \code{n_cliques_disrupted} (k=1 only).
#'       One row per species in \code{all_species}.}
#'     \item{stability_class}{Named integer vector: highest k at which
#'       each exclusive clique is stable across ALL subsets (0 = unstable
#'       at k=1)}
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
#'   \item Each resulting clique is annotated by trait exclusivity
#'   \item Reduced cliques are matched to full-dataset cliques by Jaccard
#'     similarity of gene assignments (considering only non-removed species)
#'   \item For matched cliques, trait exclusivity preservation is checked
#' }
#'
#' A clique is \emph{testable} in a subset if at least 2 of its species
#' remain active. It is \emph{stable at level k} if its trait exclusivity
#' is preserved across ALL C(N, k) subsets where it is testable.
#'
#' When \code{all_species} is larger than \code{target_species}, subsets
#' may remove non-target species. With static edges, removing a non-target
#' species does not change target-species edges, so these subsets are
#' trivially stable. This becomes meaningful when edge q-values are
#' re-computed per subset (future extension).
#'
#' ## sole_rep column
#'
#' The \code{sole_rep} column in the stability data frame is \code{TRUE} if
#' this clique's trait value has only one target species representative,
#' meaning removal of that species trivially breaks exclusivity.
#'
#' ## full_cliques parameter
#'
#' When \code{full_cliques = NULL} (default), cliques are computed internally
#' via \code{\link{find_cliques}}. You can also precompute them:
#' \preformatted{
#' fc <- find_cliques(edges, target_species)
#' stab <- clique_stability(edges, target_species, species_trait,
#'                          full_cliques = fc)
#' }
#'
#' @examples
#' \dontrun{
#' # Find annual-exclusive cliques stable across all 8 Poaceae species
#' annual_sp    <- c("BDIS", "HVUL", "BMAX", "VBRO")
#' perennial_sp <- c("BSYL", "HJUB", "BMED", "FPRA")
#' all_sp       <- c(annual_sp, perennial_sp)
#' trait <- setNames(rep(c("annual", "perennial"), each = 4), all_sp)
#'
#' cliques <- find_cliques(edges, annual_sp)
#' stab    <- clique_stability(edges, annual_sp, trait,
#'                             all_species = all_sp,
#'                             full_cliques = cliques)
#' # max_k defaults to 6 (= 8 - 2), testing k = 1..6
#'
#' # Cliques surviving any single species dropout
#' k1 <- stab$stability[stab$stability$k == 1, ]
#' stable_cliques <- cliques[k1$clique_idx[k1$stability_score == 1], ]
#'
#' # Multi-level: stability_class >= 2 survives any pair of dropouts
#' deeply_stable <- which(stab$stability_class >= 2)
#' }
#'
#' @export
clique_stability <- function(edges, target_species, species_trait,
                             all_species = target_species,
                             full_cliques = NULL,
                             min_species = length(target_species),
                             max_k = length(all_species) - 2L,
                             max_genes_per_sp = 10L,
                             jaccard_threshold = 0.8,
                             edge_type = "conserved", n_cores = 1L) {
  # Validate inputs
  required_cols <- c("gene1", "gene2", "species1", "species2", "hog",
                     "q.value", "effect_size")
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
  if (!is.character(species_trait) && !is.factor(species_trait)) {
    stop("species_trait must be a named character or factor vector")
  }
  if (is.null(names(species_trait))) {
    stop("species_trait must be a named vector with species as names")
  }
  missing_sp <- setdiff(all_species, names(species_trait))
  if (length(missing_sp) > 0) {
    stop("species_trait missing entries for: ",
         paste(missing_sp, collapse = ", "))
  }
  max_k <- as.integer(max_k)
  if (max_k < 1L) {
    stop("max_k must be >= 1")
  }
  if (max_k >= length(all_species)) {
    stop("max_k must be < length(all_species)")
  }

  # Convert species_trait to 0-based integer (for ALL species)
  trait_char <- as.character(species_trait[all_species])
  trait_levels <- unique(trait_char)
  trait_int <- as.integer(match(trait_char, trait_levels) - 1L)

  # Build is_target: 1 for target species, 0 for non-target
  is_target <- as.integer(all_species %in% target_species)

  # Empty result template
  empty_stability <- data.frame(
    clique_idx = integer(0), hog = character(0), trait_value = character(0),
    k = integer(0), n_subsets = integer(0), n_stable = integer(0),
    stability_score = numeric(0), sole_rep = logical(0)
  )
  empty_disruption <- data.frame(
    species = character(0), trait_value = character(0),
    n_cliques_disrupted = integer(0)
  )
  empty_result <- list(
    stability = empty_stability,
    clique_disruption = empty_disruption,
    stability_class = integer(0),
    novel_cliques = 0L
  )

  # Filter edges by type if applicable
  if ("type" %in% names(edges)) {
    edges <- edges[edges$type %in% edge_type, , drop = FALSE]
  }
  if (nrow(edges) == 0) return(empty_result)

  # Encode edges with ALL species (full universe)
  enc <- encode_clique_edges(edges, all_species)
  if (!enc$any_valid) return(empty_result)

  # Compute full cliques if not provided
  if (is.null(full_cliques)) {
    full_cliques <- find_cliques(edges, target_species,
                                 min_species = min_species,
                                 max_genes_per_sp = max_genes_per_sp,
                                 edge_type = edge_type)
  }
  if (nrow(full_cliques) == 0) return(empty_result)

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
        fc_genes[present, j] <- as.integer(enc$gene_map[gnames[present]])
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

  # Call C++ stability function
  cpp_result <- find_cliques_stability_cpp(
    enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
    enc$edge_qval, enc$edge_effect,
    length(all_species),
    length(enc$unique_hogs), length(enc$all_genes),
    trait_int, is_target, raw_cliques,
    as.integer(max_k), as.integer(max_genes_per_sp),
    jaccard_threshold, n_cores
  )

  # Post-process: stability data frame
  stab <- cpp_result$stability
  if (nrow(stab) > 0) {
    # Convert 0-based C++ clique_idx to 1-based R indexing
    stab$clique_idx <- stab$clique_idx + 1L
    stab$hog <- enc$unique_hogs[raw_cliques$hog_idx[stab$clique_idx] + 1L]
    stab$trait_value <- trait_levels[stab$trait_value + 1L]
  } else {
    stab$hog <- character(0)
    stab$trait_value <- character(0)
  }

  # Post-process: clique_disruption (one row per all_species)
  disrupt <- cpp_result$clique_disruption
  if (nrow(disrupt) > 0) {
    disrupt$species <- all_species[disrupt$species_idx + 1L]
    disrupt$trait_value <- trait_levels[disrupt$trait_value + 1L]
    disrupt <- disrupt[, c("species", "trait_value", "n_cliques_disrupted"),
                       drop = FALSE]
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
#'   (named numeric matrix) and \code{$threshold} (scalar).
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
#' result[result$persistence > 2.0, ]  # robust to 2x threshold tightening
#' }
#'
#' @export
clique_persistence <- function(cliques, target_species, networks, edges) {
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques(); got ",
         paste(class(cliques), collapse = "/"))
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species; got ",
         paste(class(networks), collapse = "/"))
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  missing_sp <- setdiff(target_species, names(cliques))
  if (length(missing_sp) > 0) {
    stop("cliques missing columns for species: ",
         paste(missing_sp, collapse = ", "))
  }
  missing_net <- setdiff(target_species, names(networks))
  if (length(missing_net) > 0) {
    stop("networks missing entries for species: ",
         paste(missing_net, collapse = ", "))
  }
  ortho_cols <- c("gene1", "gene2", "species1", "species2")
  missing_cols <- setdiff(ortho_cols, names(edges))
  if (length(missing_cols) > 0) {
    stop("edges missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Pre-build ortholog lookup: ortho[["sp_a.sp_b"]][[gene_a]] -> c(gene_b, ...)
  ortho <- list()
  for (pair in utils::combn(target_species, 2, simplify = FALSE)) {
    sp_a <- pair[1]
    sp_b <- pair[2]
    fwd <- edges[edges$species1 == sp_a & edges$species2 == sp_b, , drop = FALSE]
    rev <- edges[edges$species1 == sp_b & edges$species2 == sp_a, , drop = FALSE]
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

        # Neighbours (excluding self)
        mr_a <- net_a[g_a, colnames(net_a) != g_a]
        neigh_a <- names(mr_a[mr_a >= thr_a])
        if (length(neigh_a) == 0L) next

        mr_b <- net_b[g_b, colnames(net_b) != g_b]
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
#' (hypergeometric + Storey q-values).
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
#'   (default 0.05).
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
#'       multiplier): \code{clique_idx} (1-based), \code{hog},
#'       \code{multiplier}, \code{survived} (logical), \code{jaccard},
#'       \code{n_species_orig}, \code{n_species_new}.}
#'     \item{sweep_cliques}{Named list of \code{find_cliques()} outputs
#'       keyed by multiplier.}
#'     \item{sweep_edges}{Named list of combined edge data frames keyed
#'       by multiplier.}
#'   }
#'
#' @examples
#' \dontrun{
#' sweep <- clique_threshold_sweep(cliques, target_species, networks,
#'                                  orthologs, multipliers = c(1.5, 2, 5))
#' # Survival curve
#' sapply(sort(unique(sweep$survival$multiplier)),
#'        function(m) mean(sweep$survival$survived[sweep$survival$multiplier == m]))
#' }
#'
#' @export
clique_threshold_sweep <- function(
    cliques, target_species, networks, orthologs,
    species_pairs = NULL,
    multipliers = c(1.5, 2, 3, 5, 10),
    alternative = c("greater", "less"),
    alpha = 0.05,
    min_species = length(target_species),
    max_genes_per_sp = 10L,
    max_missing_edges = 0L,
    edge_type = "conserved",
    jaccard_threshold = 0.5,
    n_cores = 1L) {

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
    stop("networks missing entries for: ",
         paste(missing_net, collapse = ", "))
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
        n_species_new = integer(0)),
      sweep_cliques = list(),
      sweep_edges = list()))
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

    # Re-threshold all networks (shallow copy, R COW avoids matrix dup)
    tight_nets <- lapply(networks[target_species], function(net) {
      list(network = net$network, threshold = net$threshold * m)
    })
    names(tight_nets) <- target_species

    # Pairwise comparisons
    pair_edges <- list()
    for (pair in species_pairs) {
      sp_a <- pair[1]
      sp_b <- pair[2]

      comparison <- tryCatch(
        compare_neighborhoods(tight_nets[[sp_a]], tight_nets[[sp_b]],
                              orthologs, n_cores),
        error = function(e) {
          warning("Pair ", sp_a, "-", sp_b, " at ", m, "x failed: ",
                  conditionMessage(e))
          NULL
        }
      )
      if (is.null(comparison) || nrow(comparison) == 0) next

      summary_res <- tryCatch(
        summarize_comparison(comparison, alternative, alpha),
        error = function(e) {
          warning("Pair ", sp_a, "-", sp_b, " q-values at ", m,
                  "x failed: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(summary_res) || nrow(summary_res$results) == 0) next

      edges_df <- comparison_to_edges(summary_res$results, sp_a, sp_b,
                                       alternative, alpha)
      pair_edges[[length(pair_edges) + 1L]] <- edges_df
    }

    if (length(pair_edges) == 0) {
      all_edges <- data.frame(
        gene1 = character(0), gene2 = character(0),
        species1 = character(0), species2 = character(0),
        hog = character(0), q.value = numeric(0),
        effect_size = numeric(0), type = character(0))
    } else {
      all_edges <- do.call(rbind, pair_edges)
    }

    sweep_edges[[m_key]] <- all_edges

    # Find cliques at this threshold
    new_cliques <- find_cliques(all_edges, target_species,
                                 min_species = min_species,
                                 max_genes_per_sp = max_genes_per_sp,
                                 max_missing_edges = max_missing_edges,
                                 edge_type = edge_type)
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
          jac <- jaccard_clique_match(cliques[i, ], new_cliques[j, ],
                                       target_species)
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
      n_species_new = integer(0))
  }
  rownames(survival) <- NULL

  list(
    survival = survival,
    sweep_cliques = sweep_cliques,
    sweep_edges = sweep_edges
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
#' @param edges Data frame with columns \code{gene1}, \code{gene2},
#'   \code{species1}, \code{species2}, \code{hog}, \code{q.value},
#'   \code{effect_size}, and \code{type}. Must contain ALL edges
#'   (conserved + ns + diverged), not pre-filtered, because the
#'   differentiated check needs to verify absence of cross-group
#'   conserved edges.
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
#'
#' @return A data frame with one row per HOG:
#'   \describe{
#'     \item{hog}{HOG identifier}
#'     \item{classification}{One of \code{"complete"}, \code{"partial"},
#'       \code{"differentiated"}, \code{"trait_specific"},
#'       \code{"unclassified"}}
#'     \item{n_species}{Species count in the best clique (NA for
#'       unclassified)}
#'     \item{best_mean_q}{Mean q-value of the best clique (NA for
#'       unclassified)}
#'     \item{trait_groups}{Comma-separated trait groups with internal
#'       cliques (NA for complete/partial/unclassified)}
#'     \item{stability_class}{From stability results (NA if not
#'       provided)}
#'     \item{persistence}{Highest multiplier at which the HOG's
#'       clique survived (NA if not provided)}
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
#' @export
classify_cliques <- function(
    edges, target_species, species_trait,
    min_species = 2L,
    max_genes_per_sp = 10L,
    max_missing_edges = 0L,
    edge_type = "conserved",
    stability = NULL,
    sweep = NULL,
    min_stability_class = 0L,
    min_persistence = 1.0) {

  # --- Validation ---
  required_cols <- c("gene1", "gene2", "species1", "species2", "hog",
                     "q.value", "effect_size", "type")
  missing_cols <- setdiff(required_cols, names(edges))
  if (length(missing_cols) > 0) {
    stop("edges missing required columns: ",
         paste(missing_cols, collapse = ", "))
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
    stop("species_trait missing entries for: ",
         paste(missing_sp, collapse = ", "))
  }
  min_species <- as.integer(min_species)
  if (min_species < 2L) stop("min_species must be >= 2")

  if (!is.null(stability)) {
    if (!is.list(stability) || is.null(stability$stability))
      stop("stability must be output of clique_stability()")
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
    trait_groups = character(0), stability_class = integer(0),
    persistence = numeric(0), robust = logical(0),
    stringsAsFactors = FALSE
  )
  if (length(all_hogs) == 0) return(empty)

  # --- Steps 1+2: Find all cliques (complete + partial in one pass) ---
  all_cliques <- find_cliques(edges, target_species,
                               min_species = min_species,
                               max_genes_per_sp = max_genes_per_sp,
                               max_missing_edges = max_missing_edges,
                               edge_type = edge_type)

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
      if (length(traits_in) >= 2L) { any_cross <- TRUE; break }
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
                        edge_type = edge_type)
    within_group_cliques[[group]] <- wg
    within_group_hogs[[group]] <- unique(wg$hog)
  }

  # --- Step 4: Differentiated (2+ groups with cliques, no cross-group conserved) ---
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

  # --- Step 6: Unclassified ---
  classified <- c(complete_hogs, partial_hogs, diff_hogs, ts_hogs)
  unclass_hogs <- setdiff(all_hogs, classified)

  # --- Build output ---
  # Helper: best clique per HOG from a cliques df
  best_per_hog <- function(cliques_df, hogs) {
    sub <- cliques_df[cliques_df$hog %in% hogs, , drop = FALSE]
    if (nrow(sub) == 0) {
      return(data.frame(hog = character(0), n_species = integer(0),
                        best_mean_q = numeric(0)))
    }
    sub <- sub[order(sub$mean_q), , drop = FALSE]
    sub <- sub[!duplicated(sub$hog), , drop = FALSE]
    data.frame(hog = sub$hog, n_species = sub$n_species,
               best_mean_q = sub$mean_q, stringsAsFactors = FALSE)
  }

  rows <- list()

  # Complete
  if (length(complete_hogs) > 0) {
    info <- best_per_hog(all_cliques[is_complete, , drop = FALSE], complete_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "complete",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = NA_character_, stringsAsFactors = FALSE)
  }

  # Partial
  if (length(partial_hogs) > 0) {
    info <- best_per_hog(all_cliques, partial_hogs)
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "partial",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = NA_character_, stringsAsFactors = FALSE)
  }

  # Differentiated
  if (length(diff_hogs) > 0) {
    # Best within-group clique for each differentiated HOG
    wg_summary <- do.call(rbind, lapply(within_group_cliques[trait_levels],
      function(df) if (!is.null(df) && nrow(df) > 0)
        df[, c("hog", "n_species", "mean_q"), drop = FALSE] else NULL))
    info <- best_per_hog(wg_summary, diff_hogs)
    tg <- diff_groups[match(info$hog, diff_hogs)]
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "differentiated",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = tg, stringsAsFactors = FALSE)
  }

  # Trait-specific
  if (length(ts_hogs) > 0) {
    wg_summary2 <- do.call(rbind, lapply(within_group_cliques[trait_levels],
      function(df) if (!is.null(df) && nrow(df) > 0)
        df[, c("hog", "n_species", "mean_q"), drop = FALSE] else NULL))
    info <- best_per_hog(wg_summary2, ts_hogs)
    tg <- ts_groups[match(info$hog, ts_hogs)]
    rows[[length(rows) + 1L]] <- data.frame(
      hog = info$hog, classification = "trait_specific",
      n_species = info$n_species, best_mean_q = info$best_mean_q,
      trait_groups = tg, stringsAsFactors = FALSE)
  }

  # Unclassified
  if (length(unclass_hogs) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      hog = unclass_hogs, classification = "unclassified",
      n_species = NA_integer_, best_mean_q = NA_real_,
      trait_groups = NA_character_, stringsAsFactors = FALSE)
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) return(empty)
  rownames(out) <- NULL

  # --- Stability annotation ---
  out$stability_class <- NA_integer_
  if (!is.null(stability$stability) && nrow(stability$stability) > 0) {
    stab_df <- stability$stability
    sc_vec <- stability$stability_class
    if ("hog" %in% names(stab_df) && length(sc_vec) > 0) {
      # stability_class is per-exclusive-clique (same order as unique
      # clique_idx values in the stability DF). Map to HOG via the DF.
      excl_cliques <- unique(stab_df$clique_idx)
      excl_hogs <- stab_df$hog[match(excl_cliques, stab_df$clique_idx)]
      # Best (max) stability_class per HOG: a HOG's classification is
      # determined by its best clique, so we take the most stable one.
      # Use max (optimistic) not min (conservative) — consistent with
      # best_per_hog() which picks the lowest-mean-q clique per HOG.
      best_sc <- tapply(sc_vec, excl_hogs, max)
      idx <- match(out$hog, names(best_sc))
      out$stability_class[!is.na(idx)] <- as.integer(best_sc[idx[!is.na(idx)]])
    }
  }

  # --- Sweep annotation ---
  out$persistence <- NA_real_
  if (!is.null(sweep) && "survival" %in% names(sweep)) {
    surv <- sweep$survival
    if (nrow(surv) > 0) {
      # For each HOG, find the highest multiplier where it survived
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
  has_sweep <- !is.null(sweep) && "survival" %in% names(sweep) &&
    nrow(sweep$survival) > 0
  if (has_stab || has_sweep) {
    stab_ok <- if (has_stab) {
      !is.na(out$stability_class) & out$stability_class >= min_stability_class
    } else TRUE
    sweep_ok <- if (has_sweep) {
      !is.na(out$persistence) & out$persistence >= min_persistence
    } else TRUE
    out$robust <- stab_ok & sweep_ok
  } else {
    out$robust <- NA
  }

  out
}
