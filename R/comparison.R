#' Compare co-expression neighborhoods across species
#'
#' For each ortholog pair, tests whether co-expression neighborhoods are
#' conserved between species using hypergeometric tests in both directions.
#'
#' @section Ortholog preparation for de novo transcriptomes:
#' When working with de novo transcriptome assemblies, it is **critical** to
#' filter ortholog tables to expressed genes **before** applying any
#' maximum-paralogs filter. De novo assemblies can produce 10--50x more gene
#' models than are actually expressed (e.g. 205K models but only 14K expressed).
#' Unexpressed fragments inflate paralog counts per ortholog group, causing the
#' paralogs filter to discard groups that would pass if only expressed members
#' were counted. Failing to do this can reduce the usable gene set by 40--70%
#' and eliminate all statistical power from the hypergeometric test.
#'
#' Recommended workflow:
#' 1. Build expression matrices for both species
#' 2. Filter ortholog table to genes present in expression data
#' 3. Apply maximum-paralogs filter on the reduced ortholog table
#' 4. Compute networks and run comparison
#'
#' @param net1 Network object for species 1 (output of [compute_network()]).
#' @param net2 Network object for species 2 (output of [compute_network()]).
#' @param orthologs Data frame with columns `Species1`, `Species2`, and
#'   `hog` (output of [parse_orthologs()]).
#' @param n_cores Number of threads for parallel computation (default 1).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{Species1}{Gene identifier for species 1}
#'     \item{Species2}{Gene identifier for species 2}
#'     \item{hog}{Ortholog group identifier}
#'     \item{Species1.neigh}{Number of neighbors of Species1 gene in net1}
#'     \item{Species1.ortho.neigh}{Number of ortholog-mapped neighbors
#'       from net2}
#'     \item{Species1.neigh.overlap}{Intersection size}
#'     \item{Species1.p.val.con}{Upper-tail hypergeometric p-value for
#'       conservation (direction 1)}
#'     \item{Species1.p.val.div}{Lower-tail hypergeometric p-value for
#'       divergence (direction 1)}
#'     \item{Species1.effect.size}{Fold enrichment (direction 1). Values > 1
#'       indicate conservation, < 1 indicate divergence.}
#'     \item{Species2.neigh}{Number of neighbors of Species2 gene in net2}
#'     \item{Species2.ortho.neigh}{Number of ortholog-mapped neighbors
#'       from net1}
#'     \item{Species2.neigh.overlap}{Intersection size}
#'     \item{Species2.p.val.con}{Upper-tail hypergeometric p-value for
#'       conservation (direction 2)}
#'     \item{Species2.p.val.div}{Lower-tail hypergeometric p-value for
#'       divergence (direction 2)}
#'     \item{Species2.effect.size}{Fold enrichment (direction 2). Values > 1
#'       indicate conservation, < 1 indicate divergence.}
#'   }
#'
#' @examples
#' \dontrun{
#' comparison <- compare_neighborhoods(net_A, net_B, orthologs)
#' head(comparison[, c("Species1", "Species2", "hog",
#'                      "Species1.effect.size")])
#' }
#'
#' @export
compare_neighborhoods <- function(net1, net2, orthologs, n_cores = 1L) {
  # Validate inputs
  if (!is.list(net1) || is.null(net1$network)) {
    stop("net1 must be a network object from compute_network()")
  }
  if (!is.list(net2) || is.null(net2$network)) {
    stop("net2 must be a network object from compute_network()")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }

  net1_mat <- net1$network
  net2_mat <- net2$network
  thr1 <- net1$threshold
  thr2 <- net2$threshold

  net1_genes <- rownames(net1_mat)
  net2_genes <- rownames(net2_mat)

  # Filter orthologs to genes present in both networks and deduplicate
  orthologs <- orthologs[orthologs$Species1 %in% net1_genes &
                           orthologs$Species2 %in% net2_genes, ,
                         drop = FALSE]
  orthologs <- unique(orthologs[, c("Species1", "Species2", "hog"),
                                drop = FALSE])

  if (nrow(orthologs) == 0) {
    stop("No orthologs found in both networks")
  }

  # Build gene name -> 0-based index maps
  idx1 <- stats::setNames(seq_along(net1_genes) - 1L, net1_genes)
  idx2 <- stats::setNames(seq_along(net2_genes) - 1L, net2_genes)

  sp1_idx <- as.integer(idx1[orthologs$Species1])
  sp2_idx <- as.integer(idx2[orthologs$Species2])

  # Call C++
  result <- compare_neighborhoods_cpp(
    net1 = net1_mat,
    net2 = net2_mat,
    thr1 = thr1,
    thr2 = thr2,
    pair_sp1_idx = sp1_idx,
    pair_sp2_idx = sp2_idx,
    ortho_sp1_idx = sp1_idx,
    ortho_sp2_idx = sp2_idx,
    n_cores = n_cores
  )

  # Combine with ortholog info
  cbind(
    orthologs[, c("Species1", "Species2", "hog"), drop = FALSE],
    result
  )
}


#' Build clique edges from pairwise comparison results
#'
#' Converts output from \code{\link{summarize_comparison}} into the edge
#' format expected by \code{\link{find_cliques}} and
#' \code{\link{clique_persistence}}. Renames columns, injects species
#' identity, and computes per-pair effect sizes and classification.
#'
#' @param comparison Data frame from \code{\link{summarize_comparison}}
#'   (the \code{$results} element). Must contain columns \code{Species1},
#'   \code{Species2}, \code{hog}, plus q-value and effect-size columns
#'   from both directions.
#' @param sp1 Species abbreviation for \code{Species1} genes (e.g.,
#'   \code{"BDIS"}).
#' @param sp2 Species abbreviation for \code{Species2} genes.
#' @param alternative Which test direction to use for q-values and
#'   classification: \code{"greater"} (conservation, default) or
#'   \code{"less"} (divergence).
#' @param alpha Significance threshold for the \code{type} column
#'   (default 0.05).
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{gene1}{Gene identifier (Species1)}
#'     \item{gene2}{Gene identifier (Species2)}
#'     \item{species1}{Species abbreviation for gene1 (\code{sp1})}
#'     \item{species2}{Species abbreviation for gene2 (\code{sp2})}
#'     \item{hog}{Ortholog group identifier}
#'     \item{q.value}{Minimum of the two directional q-values}
#'     \item{effect_size}{Geometric mean of directional effect sizes}
#'     \item{type}{\code{"conserved"} or \code{"diverged"} if
#'       \code{q.value < alpha}; \code{"ns"} otherwise}
#'   }
#'
#' @examples
#' \dontrun{
#' summary_AB <- summarize_comparison(comparison_AB)
#' edges_AB <- comparison_to_edges(summary_AB$results, "SP_A", "SP_B")
#'
#' # Combine multiple species pairs for find_cliques()
#' edges <- rbind(edges_AB, edges_AC, edges_BC)
#' }
#'
#' @export
comparison_to_edges <- function(comparison, sp1, sp2,
                                alternative = c("greater", "less"),
                                alpha = 0.05) {
  alternative <- match.arg(alternative)

  suffix <- if (alternative == "greater") "con" else "div"
  q1_col <- paste0("Species1.q.val.", suffix)
  q2_col <- paste0("Species2.q.val.", suffix)

  required <- c("Species1", "Species2", "hog",
                 "Species1.effect.size", "Species2.effect.size",
                 q1_col, q2_col)
  missing_cols <- setdiff(required, names(comparison))
  if (length(missing_cols) > 0) {
    stop("comparison missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Did you pass summarize_comparison()$results?")
  }

  q_min <- pmin(comparison[[q1_col]], comparison[[q2_col]], na.rm = TRUE)
  q_min[is.infinite(q_min)] <- NA_real_
  eff_geo <- sqrt(comparison$Species1.effect.size *
                  comparison$Species2.effect.size)

  type_label <- if (alternative == "greater") "conserved" else "diverged"
  type <- ifelse(q_min < alpha, type_label, "ns")

  data.frame(
    gene1 = comparison$Species1,
    gene2 = comparison$Species2,
    species1 = sp1,
    species2 = sp2,
    hog = comparison$hog,
    q.value = q_min,
    effect_size = eff_geo,
    type = type
  )
}


#' Run pairwise comparisons across all species pairs
#'
#' Convenience function that runs the full comparison pipeline
#' (\code{\link{compare_neighborhoods}} -> \code{\link{summarize_comparison}}
#' -> \code{\link{comparison_to_edges}}) for each species pair and combines
#' the results into a single edge data frame ready for
#' \code{\link{find_cliques}}.
#'
#' @param networks Named list of \code{\link{compute_network}} outputs,
#'   keyed by species abbreviation.
#' @param orthologs Data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog} (from \code{\link{parse_orthologs}} or
#'   \code{\link{extract_orthologs}}).
#' @param species_pairs Optional list of length-2 character vectors
#'   specifying which pairs to compare. Defaults to all
#'   \code{combn(names(networks), 2)}.
#' @param alternative Passed to \code{\link{summarize_comparison}}:
#'   \code{"greater"} (conservation, default) or \code{"less"} (divergence).
#' @param alpha Significance threshold (default 0.05).
#' @param n_cores Cores for \code{\link{compare_neighborhoods}}
#'   (default 1).
#'
#' @return Data frame with columns \code{gene1}, \code{gene2},
#'   \code{species1}, \code{species2}, \code{hog}, \code{q.value},
#'   \code{effect_size}, \code{type}. Ready for \code{\link{find_cliques}}
#'   or \code{\link{classify_cliques}}.
#'
#' @examples
#' \dontrun{
#' edges <- run_pairwise_comparisons(
#'   networks = list(SP_A = net_a, SP_B = net_b, SP_C = net_c),
#'   orthologs = ortho
#' )
#' cliques <- find_cliques(edges, c("SP_A", "SP_B", "SP_C"))
#' }
#'
#' @export
run_pairwise_comparisons <- function(
    networks, orthologs,
    species_pairs = NULL,
    alternative = c("greater", "less"),
    alpha = 0.05,
    n_cores = 1L) {

  alternative <- match.arg(alternative)

  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  if (length(networks) < 2) {
    stop("networks must contain at least 2 species")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }

  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(names(networks), 2, simplify = FALSE)
  }

  n_pairs <- length(species_pairs)
  pair_edges <- vector("list", n_pairs)
  idx <- 0L

  for (pair in species_pairs) {
    sp_a <- pair[1]
    sp_b <- pair[2]

    if (!sp_a %in% names(networks)) {
      stop("species '", sp_a, "' not found in networks")
    }
    if (!sp_b %in% names(networks)) {
      stop("species '", sp_b, "' not found in networks")
    }

    comparison <- tryCatch(
      compare_neighborhoods(networks[[sp_a]], networks[[sp_b]],
                            orthologs, n_cores),
      error = function(e) {
        warning("Pair ", sp_a, "-", sp_b, " failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(comparison) || nrow(comparison) == 0) next

    summary_res <- tryCatch(
      summarize_comparison(comparison, alternative, alpha),
      error = function(e) {
        warning("Pair ", sp_a, "-", sp_b, " q-value computation failed: ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(summary_res) || nrow(summary_res$results) == 0) next

    edges_df <- comparison_to_edges(summary_res$results, sp_a, sp_b,
                                     alternative, alpha)
    idx <- idx + 1L
    pair_edges[[idx]] <- edges_df
  }

  if (idx == 0L) {
    return(data.frame(
      gene1 = character(0), gene2 = character(0),
      species1 = character(0), species2 = character(0),
      hog = character(0), q.value = numeric(0),
      effect_size = numeric(0), type = character(0)))
  }

  do.call(rbind, pair_edges[seq_len(idx)])
}


#' Query co-expression partners of a candidate HOG across species
#'
#' For a given HOG, finds which other HOGs co-express with it in each
#' species network, then aggregates across species. Useful after
#' \code{\link{identify_module_hubs}} and \code{\link{classify_hub_conservation}}
#' to explore the co-expression neighborhood of a hub gene.
#'
#' @param candidate_hog Character string: the HOG ID to query (e.g.,
#'   \code{"HOG42"}).
#' @param networks Named list of \code{\link{compute_network}} outputs,
#'   keyed by species abbreviation.
#' @param orthologs Data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog} (from \code{\link{parse_orthologs}} or
#'   \code{\link{extract_orthologs}}).
#' @param species Character vector of species to query (default: all in
#'   \code{networks}).
#' @param species_trait Optional named character vector mapping species to
#'   trait labels (e.g., \code{c(SP_A = "annual", SP_B = "perennial")}).
#'   Enables the \code{coexpressed_traits} output column.
#' @param min_species Minimum number of species where co-expression must
#'   occur for a partner to be reported (default 2).
#' @param edges Optional stacked edge data frame from
#'   \code{\link{run_pairwise_comparisons}}. When provided, attaches
#'   conservation statistics for each partner HOG.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{partner_hog}{HOG ID of the co-expression partner}
#'     \item{n_species}{Number of species where co-expression occurs}
#'     \item{coexpressed_species}{Comma-separated species names}
#'     \item{coexpressed_traits}{Comma-separated unique trait values of
#'       those species (sorted); \code{NA} if \code{species_trait} is
#'       \code{NULL}}
#'     \item{mean_weight}{Mean co-expression weight across species}
#'   }
#'   When \code{edges} is provided, three additional columns are appended:
#'   \describe{
#'     \item{partner_conserved}{Logical: does the partner HOG have any
#'       \code{"conserved"} edge?}
#'     \item{partner_mean_effect}{Mean effect size from edges for the
#'       partner HOG}
#'     \item{partner_min_q}{Minimum q-value from edges for the partner HOG}
#'   }
#'   Rows are sorted by \code{n_species} (descending), then
#'   \code{mean_weight} (descending).
#'
#' @examples
#' \dontrun{
#' trait <- c(SP_A = "annual", SP_B = "annual",
#'            SP_C = "perennial", SP_D = "perennial")
#' partners <- get_coexpressed_hogs("HOG42", networks, orthologs,
#'                                   species_trait = trait,
#'                                   min_species = 2L, edges = edges)
#' # Annual-only partners
#' partners[partners$coexpressed_traits == "annual", ]
#' }
#'
#' @export
get_coexpressed_hogs <- function(candidate_hog, networks, orthologs,
                                  species = names(networks),
                                  species_trait = NULL,
                                  min_species = 2L,
                                  edges = NULL) {
  # --- Input validation ---
  if (!is.character(candidate_hog) || length(candidate_hog) != 1L) {
    stop("candidate_hog must be a single character string")
  }
  if (!is.list(networks) || is.null(names(networks))) {
    stop("networks must be a named list keyed by species")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  species <- as.character(species)
  bad_sp <- setdiff(species, names(networks))
  if (length(bad_sp) > 0L) {
    stop("species not found in networks: ", paste(bad_sp, collapse = ", "))
  }
  min_species <- as.integer(min_species)
  if (!is.null(species_trait)) {
    missing_trait <- setdiff(species, names(species_trait))
    if (length(missing_trait) > 0L) {
      stop("species_trait missing entries for: ",
           paste(missing_trait, collapse = ", "))
    }
  }

  # --- Build gene -> HOG lookup (stack both Species columns) ---
  g1 <- orthologs[, c("Species1", "hog"), drop = FALSE]
  g2 <- orthologs[, c("Species2", "hog"), drop = FALSE]
  names(g1) <- c("gene", "hog")
  names(g2) <- c("gene", "hog")
  gene_hog_df <- unique(rbind(g1, g2))
  gene_hog_df <- gene_hog_df[!duplicated(gene_hog_df$gene), , drop = FALSE]
  gene_to_hog <- stats::setNames(gene_hog_df$hog, gene_hog_df$gene)

  # --- Build HOG -> genes-per-species lookup ---
  # For each species, find genes in the network that belong to each HOG
  hog_genes_by_sp <- list()
  for (sp in species) {
    net_genes <- rownames(networks[[sp]]$network)
    mapped <- gene_to_hog[intersect(names(gene_to_hog), net_genes)]
    hog_genes_by_sp[[sp]] <- split(names(mapped), mapped)
  }

  # --- Per-species neighbor extraction ---
  # Collect list of data frames: one row per (partner_hog, species) with weight
  records <- vector("list", length(species))
  names(records) <- species

  for (sp in species) {
    candidate_genes <- hog_genes_by_sp[[sp]][[candidate_hog]]
    if (is.null(candidate_genes) || length(candidate_genes) == 0L) next

    net_mat <- networks[[sp]]$network
    thr <- networks[[sp]]$threshold

    # Union of neighbors across all candidate genes (multi-copy)
    all_neighbors <- character(0)
    for (cg in candidate_genes) {
      row_vals <- net_mat[cg, ]
      nbrs <- names(row_vals[row_vals >= thr])
      # Exclude self
      nbrs <- setdiff(nbrs, candidate_genes)
      all_neighbors <- union(all_neighbors, nbrs)
    }
    if (length(all_neighbors) == 0L) next

    # Map neighbors to HOGs
    nbr_hogs <- gene_to_hog[all_neighbors]
    nbr_hogs <- nbr_hogs[!is.na(nbr_hogs)]
    if (length(nbr_hogs) == 0L) next

    # Per-partner-HOG mean weight
    partner_hogs <- unique(nbr_hogs)
    weights <- vapply(partner_hogs, function(ph) {
      ph_genes <- names(nbr_hogs[nbr_hogs == ph])
      mean(vapply(candidate_genes, function(cg) {
        mean(net_mat[cg, ph_genes])
      }, numeric(1)))
    }, numeric(1))

    records[[sp]] <- data.frame(
      partner_hog = partner_hogs,
      species = sp,
      weight = weights,
      stringsAsFactors = FALSE
    )
  }

  # Stack all records
  all_records <- do.call(rbind, records[!vapply(records, is.null, logical(1))])
  if (is.null(all_records) || nrow(all_records) == 0L) {
    out <- data.frame(
      partner_hog = character(0), n_species = integer(0),
      coexpressed_species = character(0), coexpressed_traits = character(0),
      mean_weight = numeric(0), stringsAsFactors = FALSE
    )
    if (!is.null(edges)) {
      out$partner_conserved <- logical(0)
      out$partner_mean_effect <- numeric(0)
      out$partner_min_q <- numeric(0)
    }
    return(out)
  }

  # --- Aggregate per partner HOG ---
  partner_split <- split(all_records, all_records$partner_hog)

  agg <- do.call(rbind, lapply(names(partner_split), function(ph) {
    df <- partner_split[[ph]]
    sp_list <- sort(unique(df$species))
    traits_val <- if (!is.null(species_trait)) {
      paste(sort(unique(species_trait[sp_list])), collapse = ",")
    } else {
      NA_character_
    }
    data.frame(
      partner_hog = ph,
      n_species = length(sp_list),
      coexpressed_species = paste(sp_list, collapse = ","),
      coexpressed_traits = traits_val,
      mean_weight = mean(df$weight),
      stringsAsFactors = FALSE
    )
  }))

  # Filter by min_species
  agg <- agg[agg$n_species >= min_species, , drop = FALSE]
  if (nrow(agg) == 0L) {
    out <- agg
    if (!is.null(edges)) {
      out$partner_conserved <- logical(0)
      out$partner_mean_effect <- numeric(0)
      out$partner_min_q <- numeric(0)
    }
    return(out)
  }

  # --- Optional edges join ---
  if (!is.null(edges)) {
    if (!all(c("hog", "q.value", "effect_size", "type") %in% names(edges))) {
      stop("edges must have columns: hog, q.value, effect_size, type")
    }
    # Summarize edges per HOG
    edge_by_hog <- split(edges, edges$hog)
    partner_conserved <- vapply(agg$partner_hog, function(ph) {
      e <- edge_by_hog[[ph]]
      if (is.null(e)) return(FALSE)
      any(e$type == "conserved", na.rm = TRUE)
    }, logical(1))
    partner_mean_effect <- vapply(agg$partner_hog, function(ph) {
      e <- edge_by_hog[[ph]]
      if (is.null(e)) return(NA_real_)
      mean(e$effect_size, na.rm = TRUE)
    }, numeric(1))
    partner_min_q <- vapply(agg$partner_hog, function(ph) {
      e <- edge_by_hog[[ph]]
      if (is.null(e)) return(NA_real_)
      min(e$q.value, na.rm = TRUE)
    }, numeric(1))

    agg$partner_conserved <- partner_conserved
    agg$partner_mean_effect <- partner_mean_effect
    agg$partner_min_q <- partner_min_q
  }

  # Sort by n_species desc, mean_weight desc
  agg <- agg[order(-agg$n_species, -agg$mean_weight), , drop = FALSE]
  rownames(agg) <- NULL
  agg
}
