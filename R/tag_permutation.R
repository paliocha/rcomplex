#' Permutation test for trait-specific module recurrence
#'
#' Tests whether HOGs recur in trait-specific modules across independent
#' species pairs more than expected by chance. The null swaps the two
#' trait labels within each pair, re-tags diverged modules -- those whose
#' topology is not preserved in the partner network -- and counts HOG
#' recurrence under the relabelled design.
#'
#' @section Null model:
#' The design is a randomised block: the pair is the block, the two
#' species in it are the two levels of the trait, and the pairs are the
#' unit of replication. The null hypothesis is that which member of a
#' pair carries \code{target_group} is unrelated to which of its modules
#' diverge. Under that hypothesis the trait label is exchangeable
#' \emph{within} a pair and not across pairs, so each draw independently
#' either swaps a pair's two labels or leaves them alone. Module
#' structure and gene content are held fixed.
#'
#' A pair contributes HOGs for \code{target_group} only when exactly one
#' of its two species carries that trait. Every swap preserves each
#' pair's pair of labels, so a pair that contributes under the observed
#' labelling contributes under every draw --- the null stays inside the
#' design that was actually run.
#'
#' Only pairs with exactly one side in \code{target_group} can be
#' swapped. A swap within a pair whose species share a label is the
#' identity, and with three or more trait values a pair whose labels
#' differ but where neither is the target contributes nothing under
#' either orientation. With \code{k} swappable pairs the null therefore
#' has exactly \code{2^k} distinct labellings. When
#' \code{2^k <= n_perm} the null is
#' \strong{enumerated exactly} and \code{n_perm} is ignored; otherwise
#' \code{n_perm} independent swap vectors are drawn. Either way the
#' observed labelling is one of the points, so the smallest attainable
#' p-value is \code{p_min} (\code{2^-k} when enumerated). \strong{A
#' design with fewer than 5 swappable pairs cannot reach p < 0.05 no
#' matter how strong the signal}, because \code{2^-4 = 0.0625}; the
#' function warns when \code{p_min > 0.05}.
#'
#' Because a swap of \emph{every} pair maps the statistic for one trait
#' value onto the statistic for the other, running the test for two
#' complementary trait values reads two entries of the same null
#' distribution. Those are not independent tests and must not be
#' corrected as if they were.
#'
#' The test generalises to any number of trait values with arbitrary
#' frequencies. It requires a disjoint pairing: each species may appear
#' in at most one pair, otherwise a within-pair swap would change another
#' pair's labels.
#'
#' An earlier version permuted trait labels across all species without
#' conditioning on the pairing. That null mixed the observed design with
#' designs having fewer contributing pairs (for four pairs, 77\% of its
#' support), which inflated its variance and made it anti-conservative.
#'
#' @param classification Data frame from
#'   \code{\link{preservation_paired}()$classification}. Must contain
#'   columns \code{pair_name}, \code{module}, \code{reference},
#'   \code{test} and \code{classification}. \code{reference} names the
#'   species whose modules the row describes and \code{test} the partner
#'   network they were tested in; \code{classification} is one of
#'   \code{"conserved"}, \code{"moderate"}, \code{"diverged"} or
#'   \code{"untested"}, and only \code{"diverged"} rows contribute.
#'   Preservation only reports modules with at least
#'   \code{min_module_size} mapped genes, so the HOG pool is smaller than
#'   the retired gene-overlap engine's.
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
#' @param n_perm Number of swap vectors to draw when the null is too
#'   large to enumerate (default 1000). Ignored, with a message, when
#'   \code{2^k <= n_perm} for \code{k} swappable pairs, because the null
#'   is then enumerated exactly and sampling could only add noise.
#' @param min_recurrence Minimum number of pairs in which a HOG must
#'   appear to be counted as recurring (default 2). \strong{This does not
#'   scale with the number of pairs and should be raised as pairs are
#'   added.} A HOG lands in at least 2 of \code{k} sides by chance with
#'   probability \code{1 - (1 - p)^k - k p (1 - p)^(k - 1)}, which is
#'   about 0.03 at \code{k = 4} but 0.19 at \code{k = 10}: the statistic
#'   saturates on chance recurrence, and simulation shows power becoming
#'   non-monotone in \code{k} and collapsing by \code{k = 10}. Scaling it
#'   as \code{max(2, round(k / 2))} restores monotone power. The default
#'   is left at 2 for continuity, not because it is right at every
#'   \code{k}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Integer: number of HOGs recurring in
#'       \code{>= min_recurrence} pairs for \code{target_group}.}
#'     \item{null_distribution}{Integer vector of recurrence counts under
#'       the null: all \code{2^k} labellings when \code{exact} is
#'       \code{TRUE}, otherwise \code{n_perm} sampled ones.}
#'     \item{p_value}{One-sided p-value. When \code{exact} is
#'       \code{TRUE} this is \code{mean(null >= observed)} over the
#'       complete null, which already includes the observed labelling.
#'       When sampled it is \code{(sum(null >= observed) + 1) /
#'       (n_perm + 1)}.}
#'     \item{p_min}{Smallest p-value this design can produce
#'       (\code{2^-k} when enumerated, \code{1 / (n_perm + 1)} when
#'       sampled). A \code{p_value} above \code{alpha} is uninformative
#'       when \code{p_min} is also above it.}
#'     \item{exact}{Logical: was the null enumerated?}
#'     \item{n_swappable}{Number of pairs with exactly one side in
#'       \code{target_group} --- the \code{k} above. For a binary trait
#'       this is every pair whose two labels differ.}
#'     \item{pair_sizes}{Data frame of per-pair diverged-HOG set sizes,
#'       with the target and partner sides named. The within-pair swap is
#'       exchangeable only if the target side is not systematically the
#'       larger one; if it is, the statistic reads set size rather than
#'       recurrence. A warning is emitted when the target side is larger
#'       in every swappable pair (and there are at least three).}
#'     \item{recurrence_table}{Data frame with columns \code{hog} and
#'       \code{n_pairs}: observed per-HOG recurrence counts (only
#'       HOGs appearing in at least 1 pair).}
#'     \item{target_group}{Echo of the input.}
#'     \item{min_recurrence}{Echo of the input.}
#'     \item{n_perm}{Number of null values actually used.}
#'   }
#'
#' @examples
#' \dontrun{
#' mod_results <- preservation_paired(modules, networks, orthologs,
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
  req_cls <- c("pair_name", "module", "reference", "test",
               "classification")
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
  # The within-pair swap null needs a disjoint pairing: swapping one
  # pair's labels must not change another pair's. A species appearing
  # twice makes the swaps dependent and the 2^k support wrong.
  # A self-pair is also a duplicate, so diagnose it first.
  if (any(pairs$sp1 == pairs$sp2)) {
    stop("pairs must have two distinct species per row")
  }
  sp_all <- c(pairs$sp1, pairs$sp2)
  dup_sp <- unique(sp_all[duplicated(sp_all)])
  if (length(dup_sp) > 0L) {
    stop("pairs must be disjoint (each species in at most one pair); ",
         "repeated: ", paste(dup_sp, collapse = ", "))
  }
  all_sp <- unique(sp_all)
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
  # "conserved" belongs to both vocabularies, so testing for known levels
  # would pass any old table carrying one. Test for the old-only levels.
  retired_cls <- c("species_specific", "partially_conserved")
  if (any(classification$classification %in% retired_cls)) {
    stop("classification uses the retired gene-overlap vocabulary (",
         paste(retired_cls, collapse = "/"), "); expected the output of ",
         "preservation_paired()")
  }
  # Both checks are needed. Testing only for retired levels lets a foreign
  # vocabulary through (relabelled rows, "Diverged", another tool's output),
  # and testing only for known levels lets an old table through on its
  # "conserved" rows, which belong to both vocabularies. Either way the filter
  # below finds no diverged rows and returns observed = 0 with no error.
  known_cls <- c("conserved", "moderate", "diverged", "untested")
  if (!all(classification$classification %in% known_cls)) {
    stop("classification$classification holds levels outside ",
         paste(known_cls, collapse = "/"), "; expected the output of ",
         "preservation_paired()")
  }
  miss_sp <- setdiff(unique(c(pairs$sp1, pairs$sp2)),
                     unique(c(classification$reference,
                              classification$test)))
  if (length(miss_sp) > 0L) {
    stop("classification has no rows for: ",
         paste(miss_sp, collapse = ", "))
  }
  # preservation_paired() defaults pair_name to "sp1.sp2" while this function
  # requires the caller to supply it, so a mismatch is easy to produce and
  # would otherwise yield observed = 0 with no error.
  miss_pn <- setdiff(pairs$pair_name, unique(classification$pair_name))
  if (length(miss_pn) > 0L) {
    stop("classification has no rows for pair_name: ",
         paste(miss_pn, collapse = ", "))
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
  ss <- classification[classification$classification == "diverged", ,
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
    mods_sp1 <- ss_pair$module[ss_pair$reference == s1 &
                                 ss_pair$test == s2]
    genes_sp1 <- if (length(mods_sp1) > 0L) {
      unlist(modules[[s1]]$module_genes[as.character(mods_sp1)],
             use.names = FALSE)
    } else {
      character(0)
    }
    hogs_sp1 <- unique(stats::na.omit(gene_to_hog[genes_sp1]))

    # sp2 side
    mods_sp2 <- ss_pair$module[ss_pair$reference == s2 &
                                 ss_pair$test == s1]
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

  # --- Null: swap the two trait labels within a pair, or not ---
  # Only pairs with exactly one side in target_group can change the
  # statistic under a swap. A same-label pair swaps to itself, and with
  # three or more trait values a pair whose labels differ but where
  # neither is the target contributes nothing under either orientation.
  # Counting either kind would duplicate every labelling and so report a
  # p_min below the floor the design can actually reach -- which is the
  # guard this function exists to provide. For a binary trait this is
  # exactly the set of pairs whose two labels differ.
  swappable <- which(xor(group[pairs$sp1] == target_group,
                         group[pairs$sp2] == target_group))
  k <- length(swappable)

  apply_swaps <- function(flip) {
    grp <- group
    idx <- swappable[flip]
    if (length(idx) > 0L) {
      a <- pairs$sp1[idx]
      b <- pairs$sp2[idx]
      grp[a] <- group[b]
      grp[b] <- group[a]
    }
    grp
  }

  # 2^k fits in an integer for k <= 30; beyond that it is far past any
  # n_perm a caller would pass, so the sampled branch takes over.
  n_labellings <- if (k <= 30L) 2^k else Inf
  exact <- n_labellings <= n_perm
  if (exact) {
    n_draw <- as.integer(n_labellings)
    if (n_perm != n_draw) {
      message("null enumerated exactly over ", n_draw,
              " labellings (2^", k, "); n_perm = ", n_perm, " ignored")
    }
    null_dist <- integer(n_draw)
    for (i in seq_len(n_draw)) {
      # bit i-1 of the counter selects which swappable pairs flip
      bits <- as.logical(bitwAnd(i - 1L, bitwShiftL(1L, seq_len(k) - 1L)))
      null_dist[i] <- count_recurring(apply_swaps(bits))
    }
    # The enumeration already contains the observed labelling, so no +1.
    p_value <- sum(null_dist >= obs) / n_draw
    p_min <- 1 / n_draw
  } else {
    n_draw <- n_perm
    null_dist <- integer(n_draw)
    for (i in seq_len(n_draw)) {
      flip <- stats::runif(k) < 0.5
      null_dist[i] <- count_recurring(apply_swaps(flip))
    }
    p_value <- (sum(null_dist >= obs) + 1L) / (n_draw + 1L)
    p_min <- 1 / (n_draw + 1L)
  }

  if (p_min > 0.05) {
    warning("only ", k, " swappable pair(s): the smallest attainable ",
            "p-value is ", signif(p_min, 3), ", so p < 0.05 is ",
            "unreachable for any signal. At least 5 swappable pairs ",
            "are needed.")
  }

  # The within-pair swap is only exchangeable if, under the null, the two
  # sides of a pair are interchangeable in what the statistic reads. A
  # target side that is systematically the larger one breaks that, and
  # the test then measures set size. Simulation puts the rejection rate
  # at 0.74 for a 13% systematic size excess with no recurrence signal at
  # all, so the caller needs to be able to see it. Per-pair sizes vary
  # for ordinary reasons; what matters is whether the target side is
  # consistently the larger across pairs.
  pair_sizes <- data.frame(
    pair_name = pairs$pair_name,
    sp1 = pairs$sp1,
    sp2 = pairs$sp2,
    n_hogs_sp1 = vapply(hog_pool, function(z) length(z$sp1), integer(1)),
    n_hogs_sp2 = vapply(hog_pool, function(z) length(z$sp2), integer(1)),
    swappable = seq_len(n_pairs) %in% swappable,
    stringsAsFactors = FALSE
  )
  target_is_sp1 <- group[pairs$sp1] == target_group
  pair_sizes$n_hogs_target <- ifelse(target_is_sp1,
                                     pair_sizes$n_hogs_sp1,
                                     pair_sizes$n_hogs_sp2)
  pair_sizes$n_hogs_partner <- ifelse(target_is_sp1,
                                      pair_sizes$n_hogs_sp2,
                                      pair_sizes$n_hogs_sp1)
  pair_sizes$n_hogs_target[!pair_sizes$swappable] <- NA_integer_
  pair_sizes$n_hogs_partner[!pair_sizes$swappable] <- NA_integer_
  rownames(pair_sizes) <- NULL
  n_larger <- sum(pair_sizes$n_hogs_target > pair_sizes$n_hogs_partner,
                  na.rm = TRUE)
  if (k >= 3L && n_larger == k) {
    warning("the ", target_group, " side holds the larger HOG set in all ",
            k, " swappable pairs, so the within-pair swap is not ",
            "exchangeable and this test may be reading set size rather ",
            "than recurrence; see $pair_sizes")
  }

  list(
    observed = obs,
    null_distribution = null_dist,
    p_value = p_value,
    p_min = p_min,
    exact = exact,
    n_swappable = k,
    pair_sizes = pair_sizes,
    recurrence_table = recurrence_table,
    target_group = target_group,
    min_recurrence = min_recurrence,
    n_perm = n_draw
  )
}
