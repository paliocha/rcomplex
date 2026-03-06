#' Detect co-expression modules in a network
#'
#' Applies community detection to a thresholded co-expression network using
#' either the Leiden algorithm or label propagation.
#'
#' @param net Network object from [compute_network()].
#' @param method Community detection method: `"leiden"` (default) or
#'   `"label_prop"`.
#' @param resolution Resolution parameter for Leiden algorithm (default 1.0).
#'   Higher values produce more, smaller modules.
#' @param objective_function Leiden objective: `"CPM"` (default) or
#'   `"modularity"`.
#' @param n_iterations Number of Leiden iterations (default 2).
#' @param seed Random seed for reproducibility (default `NULL`).
#'
#' @return A list with components:
#'   \describe{
#'     \item{modules}{Named integer vector of module assignments (gene -> module ID)}
#'     \item{module_genes}{Named list: module ID -> character vector of gene names}
#'     \item{n_modules}{Number of modules detected}
#'     \item{modularity}{Modularity score of the partition}
#'     \item{graph}{The igraph graph object used for community detection}
#'     \item{method}{Method used}
#'     \item{params}{List of parameters used}
#'   }
#'
#' @export
detect_modules <- function(net,
                           method = c("leiden", "label_prop"),
                           resolution = 1.0,
                           objective_function = c("CPM", "modularity"),
                           n_iterations = 2L,
                           seed = NULL) {
  method <- match.arg(method)
  objective_function <- match.arg(objective_function)
  n_iterations <- as.integer(n_iterations)

  if (!is.list(net) || is.null(net$network)) {
    stop("net must be a network object from compute_network()")
  }

  mat <- net$network
  thr <- net$threshold
  genes <- rownames(mat)

  # Extract edges above threshold (upper triangle only, vectorized)
  edges <- which(upper.tri(mat) & mat >= thr, arr.ind = TRUE)

  if (nrow(edges) == 0L) {
    stop("No edges above threshold; cannot detect modules")
  }

  g <- igraph::graph_from_data_frame(
    data.frame(from = genes[edges[, 1L]], to = genes[edges[, 2L]],
               stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = data.frame(name = genes)
  )
  igraph::E(g)$weight <- mat[edges]

  if (!is.null(seed)) set.seed(seed)

  if (method == "leiden") {
    comm <- igraph::cluster_leiden(
      g,
      resolution = resolution,
      objective_function = objective_function,
      n_iterations = n_iterations
    )
  } else {
    comm <- igraph::cluster_label_prop(g, weights = igraph::E(g)$weight)
  }

  membership <- igraph::membership(comm)
  names(membership) <- igraph::V(g)$name
  module_genes <- split(names(membership), membership)

  list(
    modules = membership,
    module_genes = module_genes,
    n_modules = length(module_genes),
    modularity = igraph::modularity(g, membership),
    graph = g,
    method = method,
    params = list(
      resolution = resolution,
      objective_function = objective_function,
      n_iterations = n_iterations,
      seed = seed
    )
  )
}


#' Compare co-expression modules across species
#'
#' Tests all pairs of modules between two species for significant overlap of
#' ortholog-mapped genes, using either hypergeometric or Jaccard + permutation
#' tests.
#'
#' @section Hypergeometric method:
#' For each (module_i, module_j) pair, maps module_i genes to species 2 via
#' orthologs, then tests overlap with module_j using [stats::phyper()].
#' Q-values are computed via [qvalue::qvalue()].
#'
#' @section Jaccard + permutation method:
#' Computes observed Jaccard index between ortholog-mapped module_i genes and
#' module_j genes (restricted to ortholog-mappable genes). The null permutes
#' the ortholog mapping (which sp1 genes map to which sp2 genes) using
#' Besag-Clifford adaptive stopping. Q-values use the Liang (2016) discrete
#' method via [DiscreteQvalue::DQ()].
#'
#' @param modules1 Module detection result for species 1
#'   (output of [detect_modules()]).
#' @param modules2 Module detection result for species 2
#'   (output of [detect_modules()]).
#' @param orthologs Data frame with columns `Species1`, `Species2`, and
#'   `OrthoGroup` (output of [parse_orthologs()]).
#' @param method Test method: `"hypergeometric"` (default) or `"jaccard"`.
#' @param min_exceedances Besag-Clifford stopping parameter for Jaccard method
#'   (default 50).
#' @param max_permutations Maximum permutations for Jaccard method
#'   (default 10000).
#' @param n_cores Number of threads for Jaccard permutation (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{pairs}{Data frame with one row per module pair, including
#'       `module_sp1`, `module_sp2`, `size_sp1`, `size_sp2`, `overlap`,
#'       `jaccard`, `p.value`, `q.value`}
#'     \item{best_matches}{Data frame with best match per module
#'       (from both species directions)}
#'   }
#'
#' @export
compare_modules <- function(modules1, modules2, orthologs,
                            method = c("hypergeometric", "jaccard"),
                            min_exceedances = 50L,
                            max_permutations = 10000L,
                            n_cores = 1L) {
  method <- match.arg(method)
  min_exceedances <- as.integer(min_exceedances)
  max_permutations <- as.integer(max_permutations)
  n_cores <- as.integer(n_cores)

  if (!is.list(modules1) || is.null(modules1$module_genes)) {
    stop("modules1 must be output from detect_modules()")
  }
  if (!is.list(modules2) || is.null(modules2$module_genes)) {
    stop("modules2 must be output from detect_modules()")
  }
  if (!all(c("Species1", "Species2") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2")
  }

  mg1 <- modules1$module_genes
  mg2 <- modules2$module_genes

  # Filter orthologs to genes present in both module sets
  ortho <- orthologs[orthologs$Species1 %in% names(modules1$modules) &
                       orthologs$Species2 %in% names(modules2$modules), ,
                     drop = FALSE]
  ortho <- unique(ortho[, c("Species1", "Species2"), drop = FALSE])

  if (nrow(ortho) == 0L) {
    stop("No orthologs found in both module sets")
  }

  sp2_mappable <- unique(ortho$Species2)
  sp1_to_sp2 <- split(ortho$Species2, ortho$Species1)

  if (method == "hypergeometric") {
    pairs <- compare_modules_hypergeometric(
      mg1, mg2, sp1_to_sp2, sp2_mappable
    )
  } else {
    pairs <- compare_modules_jaccard(
      mg1, mg2, ortho, sp2_mappable, sp1_to_sp2,
      min_exceedances, max_permutations, n_cores
    )
  }

  best <- compute_best_matches(pairs)
  list(pairs = pairs, best_matches = best)
}


#' Hypergeometric module comparison (internal)
#' @noRd
compare_modules_hypergeometric <- function(mg1, mg2, sp1_to_sp2, sp2_mappable) {
  mod_names1 <- names(mg1)
  mod_names2 <- names(mg2)
  N <- length(sp2_mappable)
  n_pairs <- length(mod_names1) * length(mod_names2)

  # Pre-allocate output vectors
  out_mod1 <- character(n_pairs)
  out_mod2 <- character(n_pairs)
  out_size1 <- integer(n_pairs)
  out_size2 <- integer(n_pairs)
  out_overlap <- integer(n_pairs)
  out_jaccard <- numeric(n_pairs)
  out_pval <- numeric(n_pairs)
  idx <- 0L

  for (i in seq_along(mod_names1)) {
    genes_i <- mg1[[i]]
    mapped_i <- unique(unlist(sp1_to_sp2[genes_i], use.names = FALSE))
    if (is.null(mapped_i)) mapped_i <- character(0)
    k <- length(mapped_i)

    for (j in seq_along(mod_names2)) {
      idx <- idx + 1L
      genes_j_mappable <- intersect(mg2[[j]], sp2_mappable)
      m <- length(genes_j_mappable)
      ol <- length(intersect(mapped_i, genes_j_mappable))
      un <- k + m - ol

      out_mod1[idx] <- mod_names1[i]
      out_mod2[idx] <- mod_names2[j]
      out_size1[idx] <- length(genes_i)
      out_size2[idx] <- length(mg2[[j]])
      out_overlap[idx] <- ol
      out_jaccard[idx] <- if (un > 0L) ol / un else 0
      out_pval[idx] <- if (ol > 0L && k > 0L && m > 0L) {
        stats::phyper(ol - 1L, m, N - m, k, lower.tail = FALSE)
      } else {
        1
      }
    }
  }

  pairs <- data.frame(
    module_sp1 = out_mod1, module_sp2 = out_mod2,
    size_sp1 = out_size1, size_sp2 = out_size2,
    overlap = out_overlap, jaccard = out_jaccard,
    p.value = out_pval, stringsAsFactors = FALSE
  )

  pairs$q.value <- if (n_pairs < 2L) pairs$p.value else {
    compute_qvalues(pairs$p.value)
  }
  pairs
}


#' Jaccard + permutation module comparison (internal)
#' @noRd
compare_modules_jaccard <- function(mg1, mg2, ortho, sp2_mappable,
                                    sp1_to_sp2,
                                    min_exceedances, max_permutations,
                                    n_cores) {
  mod_names1 <- names(mg1)
  mod_names2 <- names(mg2)

  # Integer-index mapping for sp2 genes (mappable only)
  sp2_idx_map <- stats::setNames(seq_along(sp2_mappable) - 1L, sp2_mappable)
  n_sp2 <- length(sp2_mappable)

  # Unique sp1 genes and their indices
  sp1_unique <- unique(ortho$Species1)
  sp1_idx_map <- stats::setNames(seq_along(sp1_unique) - 1L, sp1_unique)

  # Per ortholog row: sp1 unique gene index and sp2 gene index
  ortho_sp1_int <- as.integer(sp1_idx_map[ortho$Species1])
  ortho_sp2_int <- as.integer(sp2_idx_map[ortho$Species2])

  # Per sp2 module: mappable gene indices
  mod_sp2_sets <- lapply(mg2, function(genes) {
    as.integer(sp2_idx_map[intersect(genes, sp2_mappable)])
  })

  # Per sp1 module: mapped sp2 indices (for observed Jaccard)
  mod_sp1_mapped <- lapply(mg1, function(genes) {
    mapped <- unique(unlist(sp1_to_sp2[genes], use.names = FALSE))
    if (is.null(mapped)) return(integer(0))
    as.integer(sp2_idx_map[mapped])
  })

  # Per sp1 module: unique sp1 gene indices (for C++ permutation)
  mod1_sp1_genes <- lapply(mg1, function(genes) {
    as.integer(sp1_idx_map[intersect(genes, sp1_unique)])
  })

  # Compute observed Jaccard for all pairs
  n_pairs <- length(mod_names1) * length(mod_names2)
  mod_i_idx <- integer(n_pairs)
  mod_j_idx <- integer(n_pairs)
  obs_jaccard <- numeric(n_pairs)
  obs_overlap <- integer(n_pairs)
  size_sp1 <- integer(n_pairs)
  size_sp2 <- integer(n_pairs)

  idx <- 0L
  for (i in seq_along(mod_names1)) {
    a_set <- mod_sp1_mapped[[i]]
    for (j in seq_along(mod_names2)) {
      idx <- idx + 1L
      mod_i_idx[idx] <- i
      mod_j_idx[idx] <- j
      b_set <- mod_sp2_sets[[j]]
      ol <- length(intersect(a_set, b_set))
      un <- length(a_set) + length(b_set) - ol
      obs_overlap[idx] <- ol
      obs_jaccard[idx] <- if (un > 0L) ol / un else 0
      size_sp1[idx] <- length(mg1[[i]])
      size_sp2[idx] <- length(mg2[[j]])
    }
  }

  perm_result <- module_jaccard_permutation_cpp(
    ortho_sp1_gene = ortho_sp1_int,
    ortho_sp2_gene = ortho_sp2_int,
    n_sp1_unique = length(sp1_unique),
    n_sp2_universe = n_sp2,
    mod1_sp1_genes = mod1_sp1_genes,
    mod_sp2_sets = mod_sp2_sets,
    mod_i_idx = mod_i_idx - 1L,
    mod_j_idx = mod_j_idx - 1L,
    obs_jaccard = obs_jaccard,
    min_exceedances = min_exceedances,
    max_permutations = max_permutations,
    n_cores = n_cores
  )

  pairs <- data.frame(
    module_sp1 = mod_names1[mod_i_idx],
    module_sp2 = mod_names2[mod_j_idx],
    size_sp1 = size_sp1, size_sp2 = size_sp2,
    overlap = obs_overlap, jaccard = obs_jaccard,
    p.value = perm_result$p_value,
    stringsAsFactors = FALSE
  )

  # Discrete q-values (Liang 2016) for Besag-Clifford p-values
  early_support <- (min_exceedances + 1) /
    (seq.int(min_exceedances, max_permutations) + 1)
  maxp_support <- seq_len(min_exceedances) / (max_permutations + 1)
  bc_support <- sort(unique(c(early_support, maxp_support)))
  bc_support <- bc_support[bc_support <= 0.5]

  pairs$q.value <- if (n_pairs < 2L) pairs$p.value else {
    DiscreteQvalue::DQ(
      pairs$p.value, ss = bc_support, method = "Liang"
    )$q.values
  }
  pairs
}


#' Extract best match per module in one direction (internal)
#' @noRd
best_match_direction <- function(pairs, from_col, to_col, species_label) {
  do.call(rbind, lapply(
    split(pairs, pairs[[from_col]]),
    function(df) {
      best <- df[which.min(df$p.value), , drop = FALSE]
      data.frame(
        module = best[[from_col]],
        species = species_label,
        best_match = best[[to_col]],
        overlap = best$overlap,
        jaccard = best$jaccard,
        p.value = best$p.value,
        q.value = best$q.value,
        n_significant = sum(df$q.value < 0.05),
        stringsAsFactors = FALSE
      )
    }
  ))
}


#' Compute best matches for each module (internal)
#' @noRd
compute_best_matches <- function(pairs) {
  if (nrow(pairs) == 0L) {
    return(data.frame(
      module = character(0), species = character(0),
      best_match = character(0), overlap = integer(0),
      jaccard = numeric(0), p.value = numeric(0),
      q.value = numeric(0), n_significant = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  result <- rbind(
    best_match_direction(pairs, "module_sp1", "module_sp2", "sp1"),
    best_match_direction(pairs, "module_sp2", "module_sp1", "sp2")
  )
  rownames(result) <- NULL
  result
}


#' Classify modules by conservation status
#'
#' Classifies each module as conserved, partially conserved, or
#' species-specific based on cross-species module comparison results.
#'
#' @section Classification criteria:
#' \describe{
#'   \item{Conserved}{Best-match q-value < `alpha` AND Jaccard >= `jaccard_threshold`}
#'   \item{Partially conserved}{Best-match q-value < `alpha` AND
#'     Jaccard < `jaccard_threshold` (module split or merged across species)}
#'   \item{Species-specific}{No significant match (best-match q-value >= `alpha`)}
#' }
#'
#' @param module_comparison Output of [compare_modules()].
#' @param alpha Significance threshold for q-values (default 0.05).
#' @param jaccard_threshold Minimum Jaccard index to call a module conserved
#'   (default 0.1).
#'
#' @return A data frame with one row per module (from both species), including:
#'   \describe{
#'     \item{module}{Module identifier}
#'     \item{species}{Which species this module belongs to (`"sp1"` or `"sp2"`)}
#'     \item{classification}{One of `"conserved"`, `"partially_conserved"`,
#'       or `"species_specific"`}
#'     \item{best_match}{Best matching module in the other species}
#'     \item{best_jaccard}{Jaccard index with best match}
#'     \item{best_q}{Q-value for best match}
#'     \item{n_significant}{Number of significant matches in the other species}
#'   }
#'
#' @export
classify_modules <- function(module_comparison,
                             alpha = 0.05,
                             jaccard_threshold = 0.1) {
  if (!is.list(module_comparison) ||
        !all(c("pairs", "best_matches") %in% names(module_comparison))) {
    stop("module_comparison must be output from compare_modules()")
  }

  best <- module_comparison$best_matches

  if (nrow(best) == 0L) {
    return(data.frame(
      module = character(0), species = character(0),
      classification = character(0), best_match = character(0),
      best_jaccard = numeric(0), best_q = numeric(0),
      n_significant = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  significant <- best$q.value < alpha
  classification <- ifelse(
    !significant, "species_specific",
    ifelse(best$jaccard >= jaccard_threshold,
           "conserved", "partially_conserved")
  )

  data.frame(
    module = best$module,
    species = best$species,
    classification = classification,
    best_match = best$best_match,
    best_jaccard = best$jaccard,
    best_q = best$q.value,
    n_significant = best$n_significant,
    stringsAsFactors = FALSE
  )
}
