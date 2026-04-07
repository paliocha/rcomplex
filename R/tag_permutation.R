#' Permutation test for trait-specific module recurrence
#'
#' Tests whether HOGs recur in trait-specific modules across independent
#' species pairs more than expected by chance. The null shuffles trait
#' labels across species (preserving marginal frequencies), re-tags
#' species-specific modules, and counts HOG recurrence under the
#' permuted labelling.
#'
#' @section Null model:
#' Each permutation randomly reassigns trait labels to species,
#' preserving the number of species per trait value. Module structure
#' and gene content are held fixed --- only the trait-to-species
#' mapping changes. A pair contributes HOGs for \code{target_group}
#' only when exactly one of its two species carries that trait under
#' the (permuted) labelling. If both or neither species carry the
#' target trait, the pair contributes nothing. This ensures the test
#' is sensitive to the trait-species association, not to module size
#' or gene overlap alone.
#'
#' The test generalises to any number of trait values with arbitrary
#' frequencies.
#'
#' @param classification Data frame from
#'   \code{\link{compare_modules_paired}()$classification}. Must
#'   contain columns \code{pair_name}, \code{module}, \code{species}
#'   (\code{"sp1"} or \code{"sp2"}), and \code{classification}.
#' @param modules Named list of \code{\link{detect_modules}} outputs.
#'   Names must include all species referenced by \code{pairs}.
#' @param orthologs Data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog}.
#' @param pairs Data frame with columns \code{sp1}, \code{sp2}, and
#'   \code{pair_name}, matching the contrasts in \code{classification}.
#' @param group Named character vector mapping species identifiers to
#'   trait values (e.g., \code{c(BDIS = "annual", BSYL = "perennial")}).
#'   All species in \code{pairs} must have entries.
#' @param target_group Character string: the trait value to test
#'   recurrence for (e.g., \code{"annual"}).
#' @param n_perm Number of permutations (default 1000).
#' @param min_recurrence Minimum number of pairs in which a HOG must
#'   appear to be counted as recurring (default 2).
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Integer: number of HOGs recurring in
#'       \code{>= min_recurrence} pairs for \code{target_group}.}
#'     \item{null_distribution}{Integer vector of length
#'       \code{n_perm}: recurrence counts under the null.}
#'     \item{p_value}{One-sided p-value (conservative):
#'       \code{(sum(null >= observed) + 1) / (n_perm + 1)}.}
#'     \item{recurrence_table}{Data frame with columns \code{hog} and
#'       \code{n_pairs}: observed per-HOG recurrence counts (only
#'       HOGs appearing in at least 1 pair).}
#'     \item{target_group}{Echo of the input.}
#'     \item{min_recurrence}{Echo of the input.}
#'     \item{n_perm}{Echo of the input.}
#'   }
#'
#' @examples
#' \dontrun{
#' mod_results <- compare_modules_paired(modules, orthologs,
#'   pairs = data.frame(sp1 = c("BDIS", "HVUL"),
#'                      sp2 = c("BSYL", "HJUB"),
#'                      pair_name = c("Brachypodium", "Hordeum")),
#'   group = c(BDIS = "annual", BSYL = "perennial",
#'             HVUL = "annual", HJUB = "perennial"))
#'
#' result <- tag_permutation(
#'   mod_results$classification, modules, orthologs,
#'   pairs = data.frame(sp1 = c("BDIS", "HVUL"),
#'                      sp2 = c("BSYL", "HJUB"),
#'                      pair_name = c("Brachypodium", "Hordeum")),
#'   group = c(BDIS = "annual", BSYL = "perennial",
#'             HVUL = "annual", HJUB = "perennial"),
#'   target_group = "annual"
#' )
#' result$p_value
#' }
#'
#' @export
tag_permutation <- function(classification, modules, orthologs, pairs,
                            group, target_group,
                            n_perm = 1000L,
                            min_recurrence = 2L) {
  # --- Validation ---
  req_cls <- c("pair_name", "module", "species", "classification")
  missing_cls <- setdiff(req_cls, names(classification))
  if (length(missing_cls) > 0L) {
    stop("classification missing columns: ",
         paste(missing_cls, collapse = ", "))
  }
  if (!is.list(modules) || is.null(names(modules))) {
    stop("modules must be a named list keyed by species")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  if (!all(c("sp1", "sp2", "pair_name") %in% names(pairs))) {
    stop("pairs must have columns: sp1, sp2, pair_name")
  }
  if (!is.character(group) || is.null(names(group))) {
    stop("group must be a named character vector")
  }
  all_sp <- unique(c(pairs$sp1, pairs$sp2))
  missing_grp <- setdiff(all_sp, names(group))
  if (length(missing_grp) > 0L) {
    stop("group missing entries for: ",
         paste(missing_grp, collapse = ", "))
  }
  missing_mod <- setdiff(all_sp, names(modules))
  if (length(missing_mod) > 0L) {
    stop("modules missing entries for: ",
         paste(missing_mod, collapse = ", "))
  }
  if (!is.character(target_group) || length(target_group) != 1L) {
    stop("target_group must be a single character string")
  }
  if (!target_group %in% group) {
    stop("target_group '", target_group,
         "' not found in group values: ",
         paste(unique(group), collapse = ", "))
  }
  n_perm <- as.integer(n_perm)
  min_recurrence <- as.integer(min_recurrence)

  # --- Build gene -> HOG lookup ---
  g1 <- orthologs[, c("Species1", "hog"), drop = FALSE]
  g2 <- orthologs[, c("Species2", "hog"), drop = FALSE]
  names(g1) <- c("gene", "hog")
  names(g2) <- c("gene", "hog")
  gene_hog_df <- unique(rbind(g1, g2))
  multi_hog <- duplicated(gene_hog_df$gene)
  if (any(multi_hog)) {
    n_multi <- length(unique(gene_hog_df$gene[multi_hog]))
    warning(n_multi, " gene(s) map to multiple HOGs; ",
            "keeping first occurrence for each gene")
  }
  gene_hog_df <- gene_hog_df[!duplicated(gene_hog_df$gene), , drop = FALSE]
  gene_to_hog <- stats::setNames(gene_hog_df$hog, gene_hog_df$gene)

  # --- Pre-compute HOG sets per (pair, side) ---
  ss <- classification[classification$classification == "species_specific", ,
                       drop = FALSE]

  n_pairs <- nrow(pairs)
  hog_pool <- vector("list", n_pairs)
  names(hog_pool) <- pairs$pair_name

  for (p in seq_len(n_pairs)) {
    pn <- pairs$pair_name[p]
    s1 <- pairs$sp1[p]
    s2 <- pairs$sp2[p]

    ss_pair <- ss[ss$pair_name == pn, , drop = FALSE]

    # sp1 side
    mods_sp1 <- ss_pair$module[ss_pair$species == "sp1"]
    genes_sp1 <- if (length(mods_sp1) > 0L) {
      unlist(modules[[s1]]$module_genes[as.character(mods_sp1)],
             use.names = FALSE)
    } else {
      character(0)
    }
    hogs_sp1 <- unique(stats::na.omit(gene_to_hog[genes_sp1]))

    # sp2 side
    mods_sp2 <- ss_pair$module[ss_pair$species == "sp2"]
    genes_sp2 <- if (length(mods_sp2) > 0L) {
      unlist(modules[[s2]]$module_genes[as.character(mods_sp2)],
             use.names = FALSE)
    } else {
      character(0)
    }
    hogs_sp2 <- unique(stats::na.omit(gene_to_hog[genes_sp2]))

    hog_pool[[pn]] <- list(sp1 = hogs_sp1, sp2 = hogs_sp2,
                           sp1_name = s1, sp2_name = s2)
  }

  # --- Pair-level HOG selection (shared by counter and table builder) ---
  get_pair_hogs <- function(grp) {
    pair_hogs <- vector("list", n_pairs)
    for (p in seq_len(n_pairs)) {
      pool <- hog_pool[[p]]
      t1 <- grp[pool$sp1_name]
      t2 <- grp[pool$sp2_name]

      # Contribute only when exactly one side has the target group
      if (t1 == target_group && t2 != target_group) {
        pair_hogs[[p]] <- pool$sp1
      } else if (t2 == target_group && t1 != target_group) {
        pair_hogs[[p]] <- pool$sp2
      }
      # Both or neither: skip
    }
    pair_hogs
  }

  # --- Recurrence counter ---
  count_recurring <- function(grp) {
    all_hogs <- unlist(get_pair_hogs(grp))
    if (length(all_hogs) == 0L) return(0L)
    sum(tabulate(match(all_hogs, unique(all_hogs))) >= min_recurrence)
  }

  # --- Observed ---
  obs <- count_recurring(group)

  # --- Build observed recurrence table ---
  obs_all <- unlist(get_pair_hogs(group))
  if (length(obs_all) > 0L) {
    recurrence_table <- as.data.frame(table(obs_all),
                                       stringsAsFactors = FALSE)
    names(recurrence_table) <- c("hog", "n_pairs")
    recurrence_table$n_pairs <- as.integer(recurrence_table$n_pairs)
    recurrence_table <- recurrence_table[order(-recurrence_table$n_pairs), ,
                                          drop = FALSE]
    rownames(recurrence_table) <- NULL
  } else {
    recurrence_table <- data.frame(hog = character(0),
                                    n_pairs = integer(0))
  }

  # --- Permutation loop ---
  group_names <- names(group)
  group_vals <- unname(group)
  null_dist <- integer(n_perm)

  for (i in seq_len(n_perm)) {
    perm_group <- stats::setNames(sample(group_vals), group_names)
    null_dist[i] <- count_recurring(perm_group)
  }

  # --- P-value (conservative one-sided) ---
  p_value <- (sum(null_dist >= obs) + 1L) / (n_perm + 1L)

  list(
    observed = obs,
    null_distribution = null_dist,
    p_value = p_value,
    recurrence_table = recurrence_table,
    target_group = target_group,
    min_recurrence = min_recurrence,
    n_perm = n_perm
  )
}
