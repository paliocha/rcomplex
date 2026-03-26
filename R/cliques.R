#' Find co-expression cliques from a multi-species edge list
#'
#' For each Hierarchical Ortholog Group (HOG), builds an undirected graph from
#' qualifying edges and uses [igraph::cliques()] to find fully connected
#' subgraphs containing exactly one gene per target species.
#'
#' @param edges Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (first gene in pair)}
#'     \item{gene2}{Gene identifier (second gene in pair)}
#'     \item{species1}{Species abbreviation for gene1}
#'     \item{species2}{Species abbreviation for gene2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{type}{Edge classification (e.g. "conserved", "diverged", "n.s.")}
#'     \item{effect_size}{Numeric effect size for the edge}
#'   }
#' @param target_species Character vector of species abbreviations that must
#'   each contribute exactly one gene to each clique (e.g.
#'   `c("BDIS", "HVUL", "VBRO")`).
#' @param edge_type Character string specifying which edge type(s) qualify for
#'   clique membership (default `"conserved"`). Only edges with
#'   `type %in% edge_type` are used.
#'
#' @return A data frame with one row per clique, containing columns:
#'   \describe{
#'     \item{hog}{Ortholog group identifier}
#'     \item{<species>}{One column per target species with the gene ID}
#'     \item{min_effect_size}{Minimum effect size across all edges
#'       in the clique}
#'     \item{mean_effect_size}{Mean effect size across all edges in the clique}
#'   }
#'   Returns a zero-row data frame with correct columns if no cliques are found.
#'
#' @export
find_coexpression_cliques <- function(edges, target_species,
                                      edge_type = "conserved") {
  # Validate inputs
  required_cols <- c("gene1", "gene2", "species1", "species2", "hog",
                     "type", "effect_size")
  missing <- setdiff(required_cols, names(edges))
  if (length(missing) > 0) {
    stop("edges missing required columns: ", paste(missing, collapse = ", "))
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }

  # Empty result template
  empty_cols <- c(
    list(hog = character(0)),
    stats::setNames(lapply(target_species, \(x) character(0)), target_species),
    list(min_effect_size = numeric(0), mean_effect_size = numeric(0))
  )
  empty_result <- as.data.frame(empty_cols, stringsAsFactors = FALSE)

  if (nrow(edges) == 0) return(empty_result)

  # Filter to qualifying edges
  edges <- edges[edges$type %in% edge_type, , drop = FALSE]
  if (nrow(edges) == 0) return(empty_result)

  # Gene-to-species lookup (from both edge columns)
  gene_species <- c(
    stats::setNames(edges$species1, edges$gene1),
    stats::setNames(edges$species2, edges$gene2)
  )
  gene_species <- gene_species[!duplicated(names(gene_species))]

  k <- length(target_species)

  # Per-HOG: build graph, find k-cliques with one gene per species
  result <- split(edges, edges$hog) |>
    lapply(\(hog_df) {
      g <- igraph::graph_from_data_frame(
        hog_df[, c("gene1", "gene2")], directed = FALSE
      )
      igraph::E(g)$effect_size <- hog_df$effect_size
      igraph::V(g)$species <- gene_species[igraph::V(g)$name]

      # Restrict to target-species vertices
      target_v <- which(igraph::V(g)$species %in% target_species)
      if (length(target_v) < k) return(NULL)
      g <- igraph::induced_subgraph(g, target_v)
      g <- igraph::simplify(g, edge.attr.comb = list(effect_size = "first"))

      # Find k-cliques, keep those spanning all target species
      valid <- Filter(
        \(cl) setequal(igraph::V(g)$species[cl], target_species),
        igraph::cliques(g, min = k, max = k)
      )
      if (length(valid) == 0) return(NULL)

      # Extract gene IDs and effect sizes per clique
      valid |>
        lapply(\(cl) {
          vnames <- igraph::V(g)$name[cl]
          vspp <- igraph::V(g)$species[cl]
          genes <- stats::setNames(vnames, vspp)
          effs <- igraph::E(g)[cl %--% cl]$effect_size
          c(list(hog = hog_df$hog[1]),
            as.list(genes[target_species]),
            list(min_effect_size = min(effs), mean_effect_size = mean(effs)))
        }) |>
        dplyr::bind_rows()
    }) |>
    dplyr::bind_rows()

  if (nrow(result) == 0) return(empty_result)
  as.data.frame(result)
}


#' Encode clique edge data as 0-based integer vectors for C++
#'
#' Shared helper for [find_cliques()] and [clique_stability()]. Builds
#' string-to-int maps, converts edge columns to 0-based integer vectors,
#' and filters out edges where either species is not in `target_species`.
#'
#' @param edges Data frame with columns gene1, gene2, species1, species2,
#'   hog, fdr, effect_size (already type-filtered).
#' @param target_species Character vector of species abbreviations.
#' @return A list with components: sp_map, gene_map, hog_map, all_genes,
#'   unique_hogs, edge_hog, edge_g1, edge_g2, edge_sp1, edge_sp2,
#'   edge_fdr, edge_effect, and a logical `any_valid` flag.
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
  edge_fdr <- as.numeric(edges$fdr[valid])
  edge_effect <- as.numeric(edges$effect_size[valid])

  list(sp_map = sp_map, gene_map = gene_map, hog_map = hog_map,
       all_genes = all_genes, unique_hogs = unique_hogs,
       edge_hog = edge_hog, edge_g1 = edge_g1, edge_g2 = edge_g2,
       edge_sp1 = edge_sp1, edge_sp2 = edge_sp2,
       edge_fdr = edge_fdr, edge_effect = edge_effect,
       any_valid = any_valid)
}


#' Find co-expression cliques using C++ two-level decomposition
#'
#' For each Hierarchical Ortholog Group (HOG), uses Bron-Kerbosch with
#' pivoting on the species-level adjacency graph to find maximal species
#' cliques, then backtracking to assign the best gene per species
#' (minimising mean FDR across all edges).
#'
#' @param edges Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (first gene in pair)}
#'     \item{gene2}{Gene identifier (second gene in pair)}
#'     \item{species1}{Species for gene1}
#'     \item{species2}{Species for gene2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{fdr}{FDR (or q-value) for the edge}
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
#'     \item{mean_fdr}{Mean FDR across clique edges}
#'     \item{max_fdr}{Maximum FDR across clique edges}
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
                     "fdr", "effect_size")
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
    list(n_species = integer(0), mean_fdr = numeric(0), max_fdr = numeric(0),
         mean_effect_size = numeric(0), n_edges = integer(0))
  )
  empty_result <- as.data.frame(empty_cols, stringsAsFactors = FALSE)

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
    enc$edge_fdr, enc$edge_effect,
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
                                  ncol = ncol(gene_matrix)),
                           stringsAsFactors = FALSE)
  names(gene_df) <- target_species
  for (j in seq_len(ncol(gene_matrix))) {
    idx <- gene_matrix[, j]
    present <- !is.na(idx)
    gene_df[present, j] <- enc$all_genes[idx[present] + 1L]
  }

  # Build result data frame
  out <- data.frame(hog = hog_names, stringsAsFactors = FALSE)
  out <- cbind(out, gene_df)
  out$n_species <- as.integer(result$n_species)
  out$mean_fdr <- result$mean_fdr
  out$max_fdr <- result$max_fdr
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
#' @param target_species Character vector of species abbreviations.
#' @param species_trait Named character or factor vector mapping species to
#'   trait groups. Names must include all \code{target_species}. Example:
#'   \code{c(SP_A = "annual", SP_B = "perennial")} or
#'   \code{c(SP_A = "tropical", SP_B = "temperate", SP_C = "arctic")}.
#' @param full_cliques Output of \code{\link{find_cliques}}, or \code{NULL} to
#'   compute internally (default).
#' @param min_species Minimum species per clique
#'   (default: \code{length(target_species)}).
#' @param max_k Maximum number of species to leave out (default 3).
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
#'       \code{trait_value}, \code{n_cliques_disrupted} (k=1 only)}
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
#' For each combination of k species removed (k = 1, ..., max_k):
#' \enumerate{
#'   \item All edges involving the removed species are dropped
#'   \item Full clique detection re-runs on the reduced edge set
#'   \item Each resulting clique is annotated by trait exclusivity
#'   \item Reduced cliques are matched to full-dataset cliques by Jaccard
#'     similarity of gene assignments (considering only non-removed species)
#'   \item For matched cliques, trait exclusivity preservation is checked
#' }
#'
#' A clique is \emph{stable at level k} if its trait exclusivity is preserved
#' across ALL C(N, k) species-removal subsets where it remains testable.
#'
#' ## sole_rep column
#'
#' The \code{sole_rep} column in the stability data frame is \code{TRUE} if
#' this clique's trait value has only one species representative in the full
#' target set, meaning removal of that species trivially breaks exclusivity.
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
#' @export
clique_stability <- function(edges, target_species, species_trait,
                             full_cliques = NULL,
                             min_species = length(target_species),
                             max_k = 3L, max_genes_per_sp = 10L,
                             jaccard_threshold = 0.8,
                             edge_type = "conserved", n_cores = 1L) {
  # Validate inputs
  required_cols <- c("gene1", "gene2", "species1", "species2", "hog",
                     "fdr", "effect_size")
  missing <- setdiff(required_cols, names(edges))
  if (length(missing) > 0) {
    stop("edges missing required columns: ", paste(missing, collapse = ", "))
  }
  if (length(target_species) < 2) {
    stop("target_species must have at least 2 species")
  }
  if (!is.character(species_trait) && !is.factor(species_trait)) {
    stop("species_trait must be a named character or factor vector")
  }
  if (is.null(names(species_trait))) {
    stop("species_trait must be a named vector with species as names")
  }
  missing_sp <- setdiff(target_species, names(species_trait))
  if (length(missing_sp) > 0) {
    stop("species_trait missing entries for: ",
         paste(missing_sp, collapse = ", "))
  }
  max_k <- as.integer(max_k)
  n_cores <- as.integer(n_cores)
  if (max_k < 1L) {
    stop("max_k must be >= 1")
  }
  if (max_k >= length(target_species)) {
    stop("max_k must be < length(target_species)")
  }

  # Convert species_trait to 0-based integer
  trait_char <- as.character(species_trait[target_species])
  trait_levels <- unique(trait_char)
  trait_int <- as.integer(match(trait_char, trait_levels) - 1L)

  # Empty result template
  empty_stability <- data.frame(
    clique_idx = integer(0), hog = character(0), trait_value = character(0),
    k = integer(0), n_subsets = integer(0), n_stable = integer(0),
    stability_score = numeric(0), sole_rep = logical(0),
    stringsAsFactors = FALSE
  )
  empty_disruption <- data.frame(
    species = character(0), trait_value = character(0),
    n_cliques_disrupted = integer(0), stringsAsFactors = FALSE
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

  # Encode edges as 0-based integer vectors
  enc <- encode_clique_edges(edges, target_species)
  if (!enc$any_valid) return(empty_result)

  # Compute or re-encode full_cliques
  if (is.null(full_cliques)) {
    raw_cliques <- find_cliques_cpp(
      enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
      enc$edge_fdr, enc$edge_effect,
      length(target_species), as.integer(min_species),
      length(enc$unique_hogs), length(enc$all_genes),
      as.integer(max_genes_per_sp)
    )
  } else {
    # Re-encode from string-based find_cliques() output to integer lists
    fc_hog_idx <- as.integer(enc$hog_map[as.character(full_cliques$hog)])
    n_fc <- nrow(full_cliques)
    fc_genes <- matrix(NA_integer_, nrow = n_fc, ncol = length(target_species))
    for (j in seq_along(target_species)) {
      sp <- target_species[j]
      if (sp %in% names(full_cliques)) {
        gnames <- full_cliques[[sp]]
        mapped <- enc$gene_map[gnames[!is.na(gnames)]]
        fc_genes[!is.na(gnames), j] <- as.integer(mapped)
      }
    }
    raw_cliques <- list(
      hog_idx = fc_hog_idx,
      genes = fc_genes,
      n_species = as.integer(full_cliques$n_species),
      mean_fdr = full_cliques$mean_fdr,
      max_fdr = full_cliques$max_fdr,
      mean_effect_size = full_cliques$mean_effect_size,
      n_edges = as.integer(full_cliques$n_edges)
    )
  }

  if (length(raw_cliques$hog_idx) == 0) return(empty_result)

  # Call C++ stability function
  cpp_result <- find_cliques_stability_cpp(
    enc$edge_hog, enc$edge_g1, enc$edge_g2, enc$edge_sp1, enc$edge_sp2,
    enc$edge_fdr, enc$edge_effect,
    length(target_species), as.integer(min_species),
    length(enc$unique_hogs), length(enc$all_genes),
    trait_int, raw_cliques,
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

  # Post-process: clique_disruption data frame (k=1 only)
  disrupt <- cpp_result$clique_disruption
  if (nrow(disrupt) > 0) {
    disrupt$species <- target_species[disrupt$species_idx + 1L]
    disrupt$trait_value <- trait_levels[disrupt$trait_value + 1L]
    disrupt <- disrupt[, c("species", "trait_value", "n_cliques_disrupted"),
                       drop = FALSE]
  } else {
    disrupt <- empty_disruption
  }

  # Post-process: stability_class — highest k at which each exclusive
  # clique is stable across ALL subsets (0 = unstable at k=1)
  sc <- cpp_result$stability_class

  list(
    stability = stab,
    clique_disruption = disrupt,
    stability_class = sc,
    novel_cliques = cpp_result$novel_cliques
  )
}
