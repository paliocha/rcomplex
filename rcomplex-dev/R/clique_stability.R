#' Leave-k-out jackknife structural stability for cliques
#'
#' Tests how structurally robust each clique is to species removal.
#' The analysis removes k = 1, 2, \ldots, \code{max_k} species at a time,
#' re-runs full clique detection on each reduced species set, and checks
#' whether matching gene assignments are preserved (Jaccard similarity).
#' ALL cliques are tested, regardless of trait composition. Trait
#' annotations are added post-hoc if \code{clades} is provided.
#'
#' @param edges Data frame (same format as \code{\link{find_cliques}}).
#' @param target_species Character vector of species that define clique
#'   membership.
#' @param clades Optional named list of species vectors, one per clade.
#'   Clades may nest but must not cross. The trait groups are the
#'   top-level clades; a species in no clade is a group of its own. If
#'   provided, the output gains trait annotations (\code{traits},
#'   \code{sole_rep}). If \code{NULL} (default), trait columns are
#'   \code{NA}.
#' @param all_species Character vector of ALL species in the analysis
#'   universe (default: \code{target_species}). Leave-k-out subsets are
#'   drawn from \code{all_species}. Must be a superset of
#'   \code{target_species}. When larger than \code{target_species},
#'   removing a non-target species tests whether the clique signal is
#'   robust to changes in the broader phylogenetic context.
#' @param full_cliques Output of \code{\link{find_cliques}}, or \code{NULL} to
#'   compute internally (default).
#' @param max_k Maximum number of species to leave out
#'   (default: \code{length(all_species) - 2}, leaving at least 2 species).
#' @param n_cores Number of OpenMP threads (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{stability}{Data frame with columns: \code{clique_idx} (1-based),
#'       \code{hog}, \code{k}, \code{n_subsets}, \code{n_stable},
#'       \code{stability_score}, \code{species_present}, \code{traits},
#'       \code{sole_rep}. One row per (clique, k) pair.}
#'     \item{clique_disruption}{Data frame with columns: \code{species},
#'       \code{n_cliques_disrupted} (k=1 only), and optionally
#'       \code{trait_value} if \code{clades} is provided.
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
#' added post-hoc from \code{clades} if provided.
#'
#' ## sole_rep column
#'
#' The \code{sole_rep} column is \code{TRUE} if this clique has a single
#' trait value and that trait value has only one target species
#' representative. Only populated when \code{clades} is provided.
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
#'   clades = split(names(trait), trait),
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
  clades = NULL,
  all_species = target_species,
  full_cliques = NULL,
  max_k = length(all_species) - 2L,
  n_cores = 1L, ...
) {
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
  if (!all(target_species %in% all_species)) {
    stop("target_species must be a subset of all_species")
  }
  trait_of <- NULL
  if (!is.null(clades)) {
    trait_of <- .clade_groups(.check_clades(clades, all_species), all_species)
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
    edges <- edges[edges$type %in% "conserved", , drop = FALSE]
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
    full_cliques <- find_cliques(edges_all, target_species)
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
    enc$edge_hog, enc$edge_g1, enc$edge_g2,
    enc$edge_species1, enc$edge_species2,
    enc$edge_qval, enc$edge_effect,
    length(all_species),
    length(enc$unique_hogs), length(enc$all_genes),
    is_target, raw_cliques,
    as.integer(max_k), 10L,
    0.8, n_cores,
    1, 0
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

  if (!is.null(trait_of)) {
    trait_char <- unname(trait_of[sp_cols])
    trait_counts <- table(trait_of[target_species])

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
    if (!is.null(trait_of)) {
      disrupt$trait_value <- unname(trait_of[disrupt$species])
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
#' @param orthologs Data frame with columns \code{gene1}, \code{gene2},
#'   \code{hog}.
#' @param multipliers Numeric vector of threshold multipliers (each > 1).
#'   Default \code{c(1.5, 2, 3, 5, 10)}.
#' @param min_species Minimum species per clique
#'   (default \code{length(target_species)}).
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
  multipliers = c(1.5, 2, 3, 5, 10),
  min_species = length(target_species),
  n_cores = 1L
) {
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
  if (!all(c("gene1", "gene2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: gene1, gene2, hog")
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

  species_pairs <- utils::combn(target_species, 2, simplify = FALSE)

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
        summarize_comparison(comparison, "greater", 0.1,
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
        "greater", 0.1
      )
      pair_edges[[length(pair_edges) + 1L]] <- edges_df
    }

    if (length(pair_edges) == 0) {
      all_edges <- data.frame(
        gene1 = character(0), gene2 = character(0),
        species1 = character(0), species2 = character(0),
        hog = character(0), q_value = numeric(0),
        effect_size = numeric(0), jaccard = numeric(0),
        power = numeric(0), type = character(0)
      )
    } else {
      all_edges <- do.call(rbind, pair_edges)
    }

    sweep_edges[[m_key]] <- all_edges

    # Find cliques at this threshold
    new_cliques <- find_cliques(all_edges, target_species,
      min_species = min_species
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

      survived <- !is.na(best_jaccard) && best_jaccard >= 0.5
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
