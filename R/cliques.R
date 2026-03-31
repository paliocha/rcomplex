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
#' (minimising mean q-value across all edges).
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
#' @param edge_type If \code{edges} has a \code{type} column, keep only
#'   edges with \code{type \%in\% edge_type} (default \code{"conserved"}).
#'
#' @return A data frame with one row per clique:
#'   \describe{
#'     \item{hog}{Ortholog group identifier}
#'     \item{<species>}{One column per target species (gene ID or NA)}
#'     \item{n_species}{Number of species in the clique}
#'     \item{mean_q}{Mean q-value across clique edges}
#'     \item{max_q}{Maximum q-value across clique edges}
#'     \item{mean_effect_size}{Mean effect size across clique edges}
#'     \item{n_edges}{Number of edges (= C(n_species, 2))}
#'   }
#'
#' @export
find_cliques <- function(edges, target_species,
                         min_species = length(target_species),
                         max_genes_per_sp = 10L,
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
         mean_effect_size = numeric(0), n_edges = integer(0))
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
    as.integer(max_genes_per_sp)
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
#' stable_cliques <- cliques[k1$clique_idx[k1$stability_score == 1] + 1L, ]
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
    stab$hog <- enc$unique_hogs[raw_cliques$hog_idx[stab$clique_idx + 1L] + 1L]
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


#' Identify hub genes that recur across cliques
#'
#' For each gene appearing in the output of \code{\link{find_cliques}},
#' counts how many cliques it participates in (clique degree). A HOG can
#' produce multiple cliques with different species subsets; \code{n_cliques}
#' counts total clique appearances, while \code{n_hogs} counts distinct
#' HOGs. When stability results are provided, each clique appearance is
#' weighted by its phylogenetic stability score.
#'
#' @param cliques Output of \code{\link{find_cliques}}.
#' @param target_species Character vector of species (must match column
#'   names in \code{cliques}).
#' @param species_trait Optional named character or factor vector mapping
#'   species to trait groups (same format as \code{\link{clique_stability}}).
#'   When provided, cliques are annotated by trait exclusivity and per-trait
#'   counts are included in the output.
#' @param stability Optional output of \code{\link{clique_stability}}.
#'   When provided, adds stability-weighted columns: \code{mean_stability}
#'   (mean k=1 stability score across the gene's exclusive cliques) and
#'   \code{n_stable} (count of exclusive cliques with stability_score = 1.0
#'   at k=1). Requires \code{species_trait} to also be provided.
#' @param min_cliques Minimum number of clique appearances to include a
#'   gene in the output (default 2).
#'
#' @return A data frame sorted by \code{n_stable} (if stability provided),
#'   \code{n_exclusive} (if trait provided), or \code{n_cliques}, with:
#'   \describe{
#'     \item{gene}{Gene identifier}
#'     \item{species}{Species the gene belongs to}
#'     \item{n_cliques}{Number of cliques this gene participates in}
#'     \item{n_hogs}{Number of distinct HOGs (secondary; equals
#'       \code{n_cliques} when each HOG produces one clique)}
#'     \item{n_exclusive}{Number of trait-exclusive clique appearances
#'       (only present when \code{species_trait} is provided)}
#'     \item{n_<trait>}{One column per trait level counting exclusive
#'       clique appearances (only present when \code{species_trait} is
#'       provided)}
#'     \item{mean_stability}{Mean stability score (k=1) across exclusive
#'       cliques (only present when \code{stability} is provided)}
#'     \item{n_stable}{Count of exclusive cliques with perfect stability
#'       at k=1 (only present when \code{stability} is provided)}
#'   }
#'
#' @export
clique_hubs <- function(cliques, target_species,
                        species_trait = NULL, stability = NULL,
                        min_cliques = 2L) {
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques()")
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  missing_sp <- setdiff(target_species, names(cliques))
  if (length(missing_sp) > 0) {
    stop("cliques missing columns for species: ",
         paste(missing_sp, collapse = ", "))
  }
  if (!is.null(stability)) {
    if (is.null(species_trait))
      stop("species_trait is required when stability is provided")
    if (!is.list(stability) || is.null(stability$stability))
      stop("stability must be output of clique_stability()")
  }

  # Unpivot: one row per (clique row, gene, species)
  # Use row index (not hog) as join key — a HOG can produce multiple cliques
  # with different species compositions and therefore different trait annotations.
  long <- do.call(rbind, lapply(target_species, function(sp) {
    genes <- cliques[[sp]]
    keep <- !is.na(genes)
    if (!any(keep)) return(NULL)
    data.frame(row_idx = which(keep), hog = cliques$hog[keep],
               gene = genes[keep], species = sp)
  }))

  if (is.null(long) || nrow(long) == 0) {
    cols <- list(gene = character(0), species = character(0),
                 n_cliques = integer(0), n_hogs = integer(0))
    if (!is.null(species_trait)) {
      cols$n_exclusive <- integer(0)
    }
    return(as.data.frame(cols))
  }

  # Annotate cliques with trait exclusivity if trait provided
  if (!is.null(species_trait)) {
    trait_char <- as.character(species_trait)
    trait_levels <- unique(trait_char)

    # Per clique row: exclusive trait value or NA (mixed)
    clique_trait <- vapply(seq_len(nrow(cliques)), function(i) {
      spp <- target_species[!is.na(cliques[i, target_species])]
      traits <- trait_char[match(spp, names(species_trait))]
      if (length(unique(traits)) == 1L) traits[1] else NA_character_
    }, character(1))

    long$trait_value <- clique_trait[long$row_idx]
    long$exclusive <- !is.na(long$trait_value)
  }

  # Aggregate per gene: clique count (primary) and HOG count (secondary)
  gene_sp <- unique(long[, c("gene", "species")])
  gene_cliques <- tapply(long$row_idx, long$gene, length)
  gene_hogs <- tapply(long$hog, long$gene, function(h) length(unique(h)))
  gene_sp$n_cliques <- as.integer(gene_cliques[gene_sp$gene])
  gene_sp$n_hogs <- as.integer(gene_hogs[gene_sp$gene])

  if (!is.null(species_trait)) {
    excl <- long[long$exclusive, , drop = FALSE]

    gene_excl <- tapply(excl$row_idx, excl$gene, length)
    gene_sp$n_exclusive <- as.integer(gene_excl[gene_sp$gene])
    gene_sp$n_exclusive[is.na(gene_sp$n_exclusive)] <- 0L

    # Per-trait columns (count clique appearances, not unique HOGs)
    for (tl in trait_levels) {
      sub <- excl[excl$trait_value == tl, , drop = FALSE]
      counts <- tapply(sub$row_idx, sub$gene, length)
      col_name <- paste0("n_", tl)
      gene_sp[[col_name]] <- as.integer(counts[gene_sp$gene])
      gene_sp[[col_name]][is.na(gene_sp[[col_name]])] <- 0L
    }
  }

  # Stability-weighted columns
  if (!is.null(stability) && !is.null(species_trait)) {
    stab_k1 <- stability$stability[stability$stability$k == 1,
                                   c("clique_idx", "stability_score"),
                                   drop = FALSE]
    if (nrow(stab_k1) > 0) {
      # clique_idx is 0-based; row_idx in long is 1-based
      stab_k1$row_idx <- stab_k1$clique_idx + 1L
      long_stab <- merge(long[long$exclusive, , drop = FALSE],
                         stab_k1[, c("row_idx", "stability_score")],
                         by = "row_idx", all.x = TRUE)
      long_stab$stability_score[is.na(long_stab$stability_score)] <- 0

      gene_mean_stab <- tapply(long_stab$stability_score, long_stab$gene, mean)
      gene_n_stable <- tapply(long_stab$stability_score, long_stab$gene,
                              function(s) sum(s >= 1.0))
      gene_sp$mean_stability <- as.numeric(gene_mean_stab[gene_sp$gene])
      gene_sp$mean_stability[is.na(gene_sp$mean_stability)] <- 0
      gene_sp$n_stable <- as.integer(gene_n_stable[gene_sp$gene])
      gene_sp$n_stable[is.na(gene_sp$n_stable)] <- 0L
    } else {
      gene_sp$mean_stability <- 0
      gene_sp$n_stable <- 0L
    }
  }

  # Filter and sort
  gene_sp <- gene_sp[gene_sp$n_cliques >= min_cliques, , drop = FALSE]

  if (!is.null(stability) && "n_stable" %in% names(gene_sp)) {
    sort_ord <- order(-gene_sp$n_stable, -gene_sp$mean_stability,
                      -gene_sp$n_cliques)
  } else if (!is.null(species_trait)) {
    sort_ord <- order(-gene_sp$n_exclusive, -gene_sp$n_cliques)
  } else {
    sort_ord <- order(-gene_sp$n_cliques)
  }
  gene_sp <- gene_sp[sort_ord, , drop = FALSE]

  rownames(gene_sp) <- NULL
  gene_sp
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
#' @export
clique_persistence <- function(cliques, target_species, networks, edges) {
  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data frame from find_cliques()")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
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
