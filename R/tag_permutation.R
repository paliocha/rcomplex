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
#' A component counts toward \code{k} only when relabelling it changes
#' the statistic. A contrast whose two species share a label is pinned;
#' with three or more trait values a contrast whose labels differ but
#' where neither is the target contributes nothing under either
#' orientation; and a component whose labellings all select the same HOG
#' sets --- most often one with no diverged module above
#' \code{min_module_size} on either side --- adds no resolution.
#' Counting any of these would multiply the label space without
#' separating anything, and \code{p_min} would then claim a floor the
#' design cannot reach. The label space is enumerated exactly while it
#' holds at most \code{2^enum_max} labellings, and \code{n_perm} is then
#' ignored; beyond that \code{n_perm} labellings are drawn, one
#' admissible labelling per component. Either way the
#' observed labelling is one of the points, so the smallest attainable
#' p-value is \code{p_min} (\code{2^-k} when enumerated). \strong{A
#' design with fewer than 5 independent contrast groups cannot reach
#' p < 0.05 no matter how strong the signal}, because
#' \code{2^-4 = 0.0625}. The function warns whenever
#' \code{max(p_min, p_attainable)} exceeds 0.05, naming whichever of the
#' two --- too small a label space, or ties within it --- actually
#' binds.
#'
#' Because a swap of \emph{every} pair maps the statistic for one trait
#' value onto the statistic for the other, running the test for two
#' complementary trait values reads two entries of the same null
#' distribution. Those are not independent tests and must not be
#' corrected as if they were.
#'
#' Contrasts sharing a species are \emph{coupled}: relabelling one
#' changes the other, so they cannot be swapped independently. The unit of
#' independence is therefore the connected component of the graph whose
#' nodes are species and whose edges are contrasts, and the null is the
#' product over components of each component's admissible labellings. A
#' component's labelling is fixed by the label given to any one of its
#' species, since each contrast then forces its partner, so enumerating
#' the alphabet for one species and propagating finds them all. A disjoint
#' pairing gives one component per contrast with two labellings each, which
#' is where \code{2^k} comes from --- the special case, not the
#' assumption. The test generalises to any number of trait values with
#' arbitrary frequencies, to unbalanced designs, and to species appearing
#' in several contrasts.
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
#' @param n_perm Number of labellings to draw when the null is too large
#'   to enumerate (default 1000). Only used when the label space exceeds
#'   \code{2^enum_max}, which needs 21 independent contrast groups at the
#'   default --- 42 species if they are disjoint --- so in practice the
#'   null is always enumerated and this argument is ignored, with a
#'   message when it was supplied explicitly. Enumeration is decided on
#'   cost rather than on \code{n_perm} so that raising \code{n_perm} for
#'   precision cannot demote the null from exact to sampled. See
#'   \code{enum_max} to lower that ceiling.
#' @param seed Integer seed for the sampled branch, or \code{NULL} (default)
#'   to draw from the ambient RNG stream and leave it advanced. A seed draws
#'   from a private stream and restores the caller's on exit, the
#'   package-wide contract described under \code{\link{detect_modules}}. An
#'   enumerated null draws nothing, so a seed changes nothing there.
#' @param enum_max Enumerate while the label space holds at most
#'   \code{2^enum_max} labellings (default 20). That is
#'   \code{2^enum_max} serial evaluations of the statistic in the worst
#'   case --- about a million at the default,
#'   which on Pooideae-sized HOG pools is tens of minutes --- so lower it
#'   to fall back to \code{n_perm} sampled draws when that cost is not
#'   worth an exact p-value. Values above 30 are refused.
#' @param statistic Either \code{"count"} (default) or \code{"excess"}.
#'   Inference is exact under both --- the permutation null is recomputed
#'   on whichever is chosen --- so this is a power choice, not a validity
#'   one. \code{observed} is the recurrence count either way and
#'   \code{statistic_observed} is the value actually tested.
#'
#'   \code{"count"} is the raw number of recurring HOGs. It scales with
#'   how many HOGs the selected sides happen to hold, so a contrast whose
#'   two sides differ greatly in size dominates the null, which then ranks
#'   labellings largely by set size. Measured on the eight-species
#'   Pooideae set, 98\% of the variance of the count null is explained by
#'   the total size of the selected sides.
#'
#'   \code{"excess"} subtracts the count expected from independent sides
#'   of exactly those sizes, as a Poisson-binomial upper tail over
#'   \code{universe}. \strong{Whether this helps depends on the regime,
#'   and it can make matters worse.} It works when the sides are small
#'   relative to the universe --- on the Pooideae set it cuts the variance
#'   explained by set size from 0.98 to 0.33 and moves the observed
#'   labelling from 6th to 3rd of 16. It over-corrects when the sides
#'   between them account for most of the universe, because the
#'   independence model then predicts far more overlap than disjoint sides
#'   can produce: on a simulated design whose sides nearly partition the
#'   universe, the size dependence rises rather than falls (measured 0.97
#'   to 0.999 as the sides were made more disjoint). Check
#'   \code{$pair_sizes} against \code{universe} before trusting it, and
#'   treat a large negative \code{statistic_observed} as a sign that the
#'   universe is too small.
#' @param universe Reference universe for \code{statistic = "excess"}:
#'   a count, or a vector of HOG identifiers whose distinct values are
#'   counted. Defaults to every HOG in \code{orthologs}. \strong{This is
#'   not identifiable from the data} --- it is the set of HOGs that could
#'   have appeared in a diverged module, which the data cannot report ---
#'   and the correction is sensitive to it, so it is an argument rather
#'   than something inferred. A universe at or below the union of the
#'   selected sides forces the expectation above what disjoint sides can
#'   deliver and drives the statistic systematically negative. Ignored for
#'   \code{statistic = "count"}.
#' @param min_recurrence Minimum number of pairs in which a HOG must
#'   appear to be counted as recurring. \code{NULL} (default) uses half
#'   the contributing contrasts, at least 2 --- which is 2 for a
#'   four-contrast design, so small designs behave as they did when the
#'   default was the constant 2. A constant does not describe a design of
#'   arbitrary size. Writing \code{p} for the chance that one HOG falls in a
#'   diverged module on one side of one pair, a HOG reaches at least 2 of
#'   \code{k} target sides by chance alone with probability
#'   \code{1 - (1 - p)^k - k p (1 - p)^(k - 1)}. At the \code{p = 0.08}
#'   of the eight-species Pooideae set that is about 0.03 at
#'   \code{k = 4} and 0.19 at \code{k = 10}: the statistic
#'   saturates on chance recurrence, and simulation shows power becoming
#'   non-monotone in \code{k} and collapsing by \code{k = 10}. Scaling it
#'   as \code{max(2, round(k / 2))} restores monotone power, and is what
#'   \code{NULL} does. Exceeding the number of contributing pairs is an
#'   error:
#'   no HOG could recur in that many, so the statistic would be 0 under
#'   every labelling. At \code{min_recurrence = 1} the statistic
#'   degenerates to the size of the union of the chosen sides, a pure
#'   set-size measure.
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Integer: number of HOGs recurring in
#'       \code{>= min_recurrence} pairs for \code{target_group}.}
#'     \item{null_distribution}{Integer vector of recurrence counts under
#'       the null: all \code{2^k} labellings when \code{exact} is
#'       \code{TRUE}, otherwise \code{n_perm} sampled ones.}
#'     \item{p_value}{One-sided p-value, comparing
#'       \code{statistic_observed} --- not \code{observed} --- against
#'       \code{null_distribution}, which is a double vector. When
#'       \code{exact} is \code{TRUE} this is
#'       \code{mean(null >= statistic_observed)} over the complete null,
#'       which already includes the observed labelling; when sampled it
#'       is \code{(sum(null >= statistic_observed) + 1) / (n_perm + 1)}.
#'       The comparison carries a tolerance, since under
#'       \code{statistic = "excess"} mathematically tied labellings need
#'       not be bitwise equal.}
#'     \item{p_min}{Smallest p-value this design can produce
#'       (\code{2^-k} when enumerated, \code{1 / (n_perm + 1)} when
#'       sampled). A \code{p_value} above \code{alpha} is uninformative
#'       when \code{p_min} is also above it.}
#'     \item{p_attainable}{The floor after ties: labellings sharing a
#'       statistic cannot be separated, so nothing scores below the tail
#'       of the maximum. Equals \code{p_min} when the maximum is unique,
#'       and whenever the null was sampled.}
#'     \item{statistic}{Echo of the input.}
#'     \item{statistic_observed}{The observed value of the chosen
#'       statistic --- what \code{null_distribution} is compared against.
#'       Equals \code{observed} when \code{statistic = "count"}.}
#'     \item{n_labellings}{Size of the admissible label space, the
#'       product over components, which is \code{2^k}: a component with
#'       at least one contrast admits at most two labellings, since the
#'       seed species can only carry one of that contrast's two labels.}
#'     \item{n_contributing}{Number of pairs feeding the statistic ---
#'       those with exactly one side in \code{target_group}. Exceeds
#'       \code{n_swappable} both when contrasts are coupled into one
#'       component and when a component's labellings all select the same
#'       HOG sets.}
#'     \item{exact}{Logical: was the null enumerated?}
#'     \item{n_swappable}{The \code{k} above: the number of
#'       \emph{components} some relabelling moves, not a count of pairs.
#'       For a disjoint design the two coincide.}
#'     \item{pair_sizes}{Data frame with one row per contrast and
#'       columns \code{pair_name}, \code{sp1}, \code{sp2},
#'       \code{n_hogs_sp1}, \code{n_hogs_sp2}, \code{contributes},
#'       \code{block}, \code{swappable}, \code{n_hogs_target} and
#'       \code{n_hogs_partner}. The last two are \code{NA} on rows that
#'       do not contribute, since those have no target side.
#'       \code{block} gives the
#'       component each contrast belongs to and \code{swappable} is a
#'       property of that component, so it is \code{TRUE} for every
#'       contrast in a component some relabelling moves. The relabelling
#'       is exchangeable only if the target side is not systematically
#'       the larger one; if it is, the statistic reads set size rather
#'       than recurrence.}
#'     \item{size_asymmetry_p}{One-sided sign-test p-value for the target
#'       side being the larger one across the non-tied contrasts in
#'       components some relabelling moves.
#'       A warning is emitted when it is at or below 0.10 --- an advisory
#'       threshold, since like \code{p_min} this cannot reach 0.05 below
#'       five pairs. \code{NA} when
#'       every swappable pair is tied. This is a diagnostic on the null's
#'       exchangeability assumption, not part of the test.}
#'     \item{recurrence_table}{Data frame with columns \code{hog} and
#'       \code{n_pairs}: observed per-HOG recurrence counts (only
#'       HOGs appearing in at least 1 pair).}
#'     \item{target_group}{Echo of the input.}
#'     \item{min_recurrence}{The threshold actually used --- the
#'       resolved value, not an echo, when the input was \code{NULL}.}
#'     \item{n_perm}{Number of null values actually used.}
#'   }
#'
#' @examples
#' \dontrun{
#' mod_results <- preservation_paired(modules, networks, orthologs,
#'   pairs = data.frame(
#'     sp1 = c("BDIS", "HVUL"),
#'     sp2 = c("BSYL", "HJUB"),
#'     pair_name = c("Brachypodium", "Hordeum")
#'   ),
#'   group = c(
#'     BDIS = "annual", BSYL = "perennial",
#'     HVUL = "annual", HJUB = "perennial"
#'   )
#' )
#'
#' result <- tag_permutation(
#'   mod_results$classification, modules, orthologs,
#'   pairs = data.frame(
#'     sp1 = c("BDIS", "HVUL"),
#'     sp2 = c("BSYL", "HJUB"),
#'     pair_name = c("Brachypodium", "Hordeum")
#'   ),
#'   group = c(
#'     BDIS = "annual", BSYL = "perennial",
#'     HVUL = "annual", HJUB = "perennial"
#'   ),
#'   target_group = "annual"
#' )
#' result$p_value
#' }
#'
#' @export
tag_permutation <- function(classification, modules, orthologs, pairs,
                            group, target_group,
                            n_perm = 1000L,
                            min_recurrence = NULL,
                            statistic = c("count", "excess"),
                            universe = NULL,
                            enum_max = 20L,
                            seed = NULL) {
  # Only the sampled branch draws, but the scope is opened unconditionally:
  # enum_max decides which branch runs, so a caller cannot tell from the call
  # site whether a seed matters. See .seed_scope() in R/rng.R.
  .seed_scope(seed)

  # Capture before n_perm is reassigned: missing() reports FALSE once an
  # argument has been written to, so this cannot be asked for later.
  n_perm_supplied <- !missing(n_perm)
  min_recurrence_in <- min_recurrence
  statistic <- match.arg(statistic)
  # --- Validation ---
  req_cls <- c(
    "pair_name", "module", "reference", "test",
    "classification"
  )
  missing_cls <- setdiff(req_cls, names(classification))
  if (length(missing_cls) > 0L) {
    stop(
      "classification missing columns: ",
      paste(missing_cls, collapse = ", ")
    )
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
  if (any(pairs$sp1 == pairs$sp2)) {
    stop("pairs must have two distinct species per row")
  }
  # hog_pool is named by pair_name but indexed positionally, so a
  # duplicate silently leaves later slots NULL and dies later with an
  # error naming neither pairs nor pair_name.
  if (anyDuplicated(pairs$pair_name)) {
    stop(
      "pairs$pair_name must be unique; repeated: ",
      paste(unique(pairs$pair_name[duplicated(pairs$pair_name)]),
        collapse = ", "
      )
    )
  }
  all_sp <- unique(c(pairs$sp1, pairs$sp2))
  missing_grp <- setdiff(all_sp, names(group))
  if (length(missing_grp) > 0L) {
    stop(
      "group missing entries for: ",
      paste(missing_grp, collapse = ", ")
    )
  }
  if (anyNA(group[all_sp])) {
    stop(
      "group has NA trait values for: ",
      paste(all_sp[is.na(group[all_sp])], collapse = ", ")
    )
  }
  missing_mod <- setdiff(all_sp, names(modules))
  if (length(missing_mod) > 0L) {
    stop(
      "modules missing entries for: ",
      paste(missing_mod, collapse = ", ")
    )
  }
  if (!is.character(target_group) || length(target_group) != 1L) {
    stop("target_group must be a single character string")
  }
  if (!target_group %in% group) {
    stop(
      "target_group '", target_group,
      "' not found in group values: ",
      paste(unique(group), collapse = ", ")
    )
  }
  # "conserved" belongs to both vocabularies, so testing for known levels
  # would pass any old table carrying one. Test for the old-only levels.
  retired_cls <- c("species_specific", "partially_conserved")
  if (any(classification$classification %in% retired_cls)) {
    stop(
      "classification uses the retired gene-overlap vocabulary (",
      paste(retired_cls, collapse = "/"), "); expected the output of ",
      "preservation_paired()"
    )
  }
  # Both checks are needed. Testing only for retired levels lets a foreign
  # vocabulary through (relabelled rows, "Diverged", another tool's output),
  # and testing only for known levels lets an old table through on its
  # "conserved" rows, which belong to both vocabularies. Either way the filter
  # below finds no diverged rows and returns observed = 0 with no error.
  known_cls <- c("conserved", "moderate", "diverged", "untested")
  if (!all(classification$classification %in% known_cls)) {
    stop(
      "classification$classification holds levels outside ",
      paste(known_cls, collapse = "/"), "; expected the output of ",
      "preservation_paired()"
    )
  }
  miss_sp <- setdiff(
    unique(c(pairs$sp1, pairs$sp2)),
    unique(c(
      classification$reference,
      classification$test
    ))
  )
  if (length(miss_sp) > 0L) {
    stop(
      "classification has no rows for: ",
      paste(miss_sp, collapse = ", ")
    )
  }
  # preservation_paired() defaults pair_name to "sp1.sp2" while this function
  # requires the caller to supply it, so a mismatch is easy to produce and
  # would otherwise yield observed = 0 with no error.
  miss_pn <- setdiff(pairs$pair_name, unique(classification$pair_name))
  if (length(miss_pn) > 0L) {
    stop(
      "classification has no rows for pair_name: ",
      paste(miss_pn, collapse = ", ")
    )
  }
  if (!is.numeric(n_perm) || length(n_perm) != 1L || is.na(n_perm) ||
    n_perm < 1) {
    stop("n_perm must be a single positive number")
  }
  # as.integer() overflows to NA above 2^31, which used to surface as
  # "missing value where TRUE/FALSE needed" from an unrelated `if`.
  n_perm <- as.integer(min(n_perm, .Machine$integer.max))
  if (!is.null(min_recurrence)) {
    if (!is.numeric(min_recurrence) || length(min_recurrence) != 1L ||
      is.na(min_recurrence) || min_recurrence < 1) {
      stop("min_recurrence must be a single positive number or NULL")
    }
    min_recurrence <- as.integer(min_recurrence)
  }
  if (!is.numeric(enum_max) || length(enum_max) != 1L ||
    is.na(enum_max) || enum_max < 0 || enum_max > 30 ||
    enum_max != round(enum_max)) {
    stop("enum_max must be a single whole number between 0 and 30")
  }
  enum_max <- as.integer(enum_max)

  # --- Build gene -> HOG lookup ---
  g1 <- orthologs[, c("Species1", "hog"), drop = FALSE]
  g2 <- orthologs[, c("Species2", "hog"), drop = FALSE]
  names(g1) <- c("gene", "hog")
  names(g2) <- c("gene", "hog")
  gene_hog_df <- unique(rbind(g1, g2))
  multi_hog <- duplicated(gene_hog_df$gene)
  if (any(multi_hog)) {
    n_multi <- length(unique(gene_hog_df$gene[multi_hog]))
    warning(
      n_multi, " gene(s) map to multiple HOGs; ",
      "keeping first occurrence for each gene"
    )
  }
  gene_hog_df <- gene_hog_df[!duplicated(gene_hog_df$gene), , drop = FALSE]
  gene_to_hog <- stats::setNames(gene_hog_df$hog, gene_hog_df$gene)

  # --- Pre-compute HOG sets per (pair, side) ---
  ss <- classification[classification$classification == "diverged", ,
    drop = FALSE
  ]

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
        use.names = FALSE
      )
    } else {
      character(0)
    }
    hogs_sp1 <- unique(stats::na.omit(gene_to_hog[genes_sp1]))

    # sp2 side
    mods_sp2 <- ss_pair$module[ss_pair$reference == s2 &
      ss_pair$test == s1]
    genes_sp2 <- if (length(mods_sp2) > 0L) {
      unlist(modules[[s2]]$module_genes[as.character(mods_sp2)],
        use.names = FALSE
      )
    } else {
      character(0)
    }
    hogs_sp2 <- unique(stats::na.omit(gene_to_hog[genes_sp2]))

    hog_pool[[pn]] <- list(
      sp1 = hogs_sp1, sp2 = hogs_sp2,
      sp1_name = s1, sp2_name = s2
    )
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

  # Which contrasts feed the statistic, and how deep a HOG must recur.
  contributes <- unname(xor(
    group[pairs$sp1] == target_group,
    group[pairs$sp2] == target_group
  ))
  n_contributing <- sum(contributes)
  # A fixed min_recurrence does not describe a design of arbitrary size.
  # The chance a HOG reaches two of k sides on its own grows steeply with
  # k -- about 0.03 at k = 4 but 0.25 at k = 12 for a per-side rate of
  # 0.08 -- so the statistic saturates on chance recurrence and power
  # stops rising with k. Half the contributing contrasts holds that
  # chance near constant. At four contrasts this is 2, the previous
  # fixed default, so small designs are unaffected.
  # Scale on the contrasts that can actually supply a HOG, not on those
  # that merely have a target side. A contrast whose target side holds no
  # diverged module contributes the empty set, so counting it raises the
  # threshold above what the data can reach: the observed statistic then
  # collapses to 0 while the null still scores, and the guard below
  # cannot see it because the threshold is under n_contributing.
  n_supplying <- sum(lengths(get_pair_hogs(group)) > 0L)
  if (is.null(min_recurrence)) {
    # ceiling, not round: round() is half-to-even, which makes "half"
    # flat across k = 7, 8, 9 and again across 11, 12, 13.
    min_recurrence <- max(2L, as.integer(ceiling(n_supplying / 2)))
  }
  if (min_recurrence > n_contributing) {
    stop(
      "min_recurrence (", min_recurrence, ") exceeds the number of ",
      "pairs contributing to the statistic (", n_contributing,
      "): no HOG can recur in that many pairs, so the statistic is 0 ",
      "under every labelling and the p-value is 1 by construction"
    )
  }
  # The message belongs after the guard, or a design that is about to
  # error first announces a threshold it will never use.
  if (is.null(min_recurrence_in) && n_supplying > 0L) {
    message(
      "min_recurrence = ", min_recurrence, " (half of ",
      n_supplying, " contrast(s) able to supply a HOG); pass it ",
      "explicitly to override"
    )
  }

  # --- Recurrence counter ---
  count_recurring <- function(grp) {
    all_hogs <- unlist(get_pair_hogs(grp))
    if (length(all_hogs) == 0L) {
      return(0L)
    }
    sum(tabulate(match(all_hogs, unique(all_hogs))) >= min_recurrence)
  }

  # The reference universe for the "excess" statistic. It is NOT
  # identifiable from the data and the correction is sensitive to it, so
  # it is an argument rather than an inference: see @param universe.
  n_universe <- if (is.null(universe)) {
    length(unique(stats::na.omit(gene_to_hog)))
  } else if (is.numeric(universe) && length(universe) == 1L) {
    if (is.na(universe) || !is.finite(universe) || universe < 1 ||
      universe != round(universe) ||
      universe > .Machine$integer.max) {
      stop(
        "universe must be a whole number of at least 1 and at most ",
        .Machine$integer.max, ", or a vector of HOG identifiers"
      )
    }
    as.integer(universe)
  } else {
    if (length(universe) == 0L || all(is.na(universe))) {
      stop(
        "universe must be a whole number or a non-empty vector of ",
        "HOG identifiers"
      )
    }
    length(unique(stats::na.omit(universe)))
  }
  if (statistic == "excess" && n_universe <= 0L) {
    stop("statistic = \"excess\" needs a positive universe size")
  }

  # --- The statistic ---
  # "count" is the raw recurrence count, which scales with how many HOGs
  # the chosen sides happen to hold: a contrast whose two sides differ
  # greatly in size then dominates, and the relabelling null measures set
  # size rather than shared identity. "excess" subtracts the count
  # expected from independent sides of exactly those sizes, which removes
  # that term without assuming anything about the trait.
  statistic_of <- function(grp) {
    sel <- get_pair_hogs(grp)
    all_hogs <- unlist(sel)
    obs_n <- if (length(all_hogs) == 0L) {
      0L
    } else {
      sum(tabulate(match(all_hogs, unique(all_hogs))) >= min_recurrence)
    }
    if (statistic == "count") {
      return(as.numeric(obs_n))
    }
    sizes <- vapply(sel, length, integer(1))
    sizes <- sizes[sizes > 0L]
    obs_n - .tp_expected(sizes, n_universe, min_recurrence)
  }

  # --- Observed ---
  obs <- statistic_of(group)
  observed_count <- count_recurring(group)

  # --- Build observed recurrence table ---
  obs_all <- unlist(get_pair_hogs(group))
  if (length(obs_all) > 0L) {
    recurrence_table <- as.data.frame(table(obs_all),
      stringsAsFactors = FALSE
    )
    names(recurrence_table) <- c("hog", "n_pairs")
    recurrence_table$n_pairs <- as.integer(recurrence_table$n_pairs)
    recurrence_table <- recurrence_table[order(-recurrence_table$n_pairs), ,
      drop = FALSE
    ]
    rownames(recurrence_table) <- NULL
  } else {
    recurrence_table <- data.frame(
      hog = character(0),
      n_pairs = integer(0)
    )
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
  # Contrasts sharing a species are coupled, so the unit of independence
  # is the connected component of the species-by-contrast graph, not the
  # contrast. A disjoint pairing gives one component per contrast with
  # two labellings each, recovering the 2^k space as a special case.
  blocks <- .tp_blocks(pairs, group)
  # A component whose labellings all select the same HOG sets adds no
  # resolution: keeping them would multiply the label space without
  # separating anything, and p_min would then claim a floor the design
  # cannot reach. Reduce such a component to one labelling.
  informative <- vapply(seq_along(blocks$labellings), function(ci) {
    labs <- blocks$labellings[[ci]]
    if (length(labs) < 2L) {
      return(FALSE)
    }
    sigs <- vapply(labs, function(l) {
      grp <- group
      grp[names(l)] <- l
      paste(vapply(get_pair_hogs(grp), function(h) {
        paste(sort(h), collapse = "\r")
      }, character(1)), collapse = "\n")
    }, character(1))
    length(unique(sigs)) > 1L
  }, logical(1))
  for (ci in which(!informative)) {
    blocks$labellings[[ci]] <- blocks$labellings[[ci]][1L]
  }
  blocks$n <- vapply(blocks$labellings, length, integer(1))

  n_labellings <- prod(as.numeric(blocks$n))
  # k is reported as the number of components that actually move the
  # statistic; for a disjoint binary design this is the familiar count of
  # swappable pairs and n_labellings is 2^k.
  k <- sum(blocks$n > 1L)

  apply_labelling <- function(pick) {
    grp <- group
    for (ci in seq_along(pick)) {
      l <- blocks$labellings[[ci]][[pick[ci]]]
      grp[names(l)] <- l
    }
    grp
  }

  # Enumeration is a cost question, not a Monte Carlo budget question.
  # The ceiling is on the number of labellings, so a design with many
  # coupled contrasts is judged by the size of its actual label space
  # rather than by a pair count that no longer describes it.
  enum_cap <- 2^enum_max
  exact <- n_labellings <= enum_cap
  # Mixed-radix counter over the components: digit ci ranges over that
  # component's admissible labellings, so the product space is walked
  # without materialising it.
  radix <- blocks$n
  pick_at <- function(i) {
    rem <- i
    pick <- integer(length(radix))
    for (ci in seq_along(radix)) {
      pick[ci] <- (rem %% radix[ci]) + 1L
      rem <- rem %/% radix[ci]
    }
    pick
  }

  if (exact) {
    n_draw <- as.integer(n_labellings)
    if (n_perm_supplied) {
      message(
        "null enumerated exactly over ", n_draw,
        " labellings; n_perm = ", n_perm, " ignored"
      )
    }
    null_dist <- numeric(n_draw)
    for (i in seq_len(n_draw)) {
      null_dist[i] <- statistic_of(apply_labelling(pick_at(i - 1L)))
    }
    # The enumeration already contains the observed labelling, so no +1.
    # "excess" makes the null doubles whose expectation term is a
    # convolution taken in contrast order, so two labellings that are
    # mathematically tied need not be bitwise equal. An exact >= would
    # drop such points from the tail (anti-conservative) and undercount
    # tied maxima, weakening the p_attainable guard.
    tol <- 1e-9 * max(1, abs(obs), max(abs(null_dist)))
    p_value <- sum(null_dist >= obs - tol) / n_draw
    p_min <- 1 / n_draw
    # The floor this data actually reaches. Ties at the maximum raise it
    # above 1 / n_labellings: labellings sharing a statistic cannot be
    # separated, so nothing scores below their shared tail.
    p_attainable <- sum(null_dist >= max(null_dist) - tol) / n_draw
  } else {
    n_draw <- n_perm
    null_dist <- numeric(n_draw)
    for (i in seq_len(n_draw)) {
      pick <- vapply(radix, function(r) sample.int(r, 1L), integer(1))
      null_dist[i] <- statistic_of(apply_labelling(pick))
    }
    tol <- 1e-9 * max(1, abs(obs), max(abs(null_dist)))
    p_value <- (sum(null_dist >= obs - tol) + 1L) / (n_draw + 1L)
    p_min <- 1 / (n_draw + 1L)
    # Ties bind here too. Setting this to p_min unconditionally meant
    # lowering enum_max to cap runtime silently switched the whole tie
    # diagnostic off, even with most draws sharing the maximum.
    p_attainable <- (sum(null_dist >= max(null_dist) - tol) + 1L) /
      (n_draw + 1L)
  }

  p_floor <- max(p_min, p_attainable)
  if (p_floor > 0.05) {
    # Contributing contrasts sitting in a component that no relabelling
    # moves. n_contributing - k would absorb the coupling as well, and
    # would then claim pairs are degenerate when they are merely joined.
    pinned <- sum(contributes & !informative[blocks$membership[pairs$sp1]])
    extra <- if (pinned > 0L) {
      paste0(
        " (", pinned, " contributing contrast(s) sit in a group no ",
        "relabelling moves)"
      )
    } else {
      ""
    }
    # Two things can put the floor above 0.05, and they call for
    # different remedies, so name whichever actually binds -- and both
    # when both do. Blaming the label space in the sampled branch would
    # prescribe more contrasts when the knob is n_perm; blaming ties when
    # p_min alone already exceeds 0.05 would deny that more contrasts
    # help, which they do.
    too_few <- p_min > 0.05
    tied <- p_attainable > p_min && p_attainable > 0.05
    space <- if (exact) {
      paste0(
        k, " independent contrast group(s)", extra, " over ",
        n_contributing, " contributing contrast(s) give ",
        n_draw, " labellings"
      )
    } else {
      paste0("a sampled null of ", n_draw, " draws")
    }
    remedy <- character(0)
    if (too_few) {
      remedy <- c(remedy, if (exact) {
        "at least 5 independent contrast groups are needed"
      } else {
        paste0("n_perm must be at least 20 (it is ", n_perm, ")")
      })
    }
    if (tied) {
      shared <- round(p_attainable * n_draw)
      remedy <- c(
        remedy,
        paste0(
          "the maximum is shared by ", shared, " of ",
          n_draw, " draws, so separating them matters ",
          "more than adding contrasts"
        )
      )
    }
    warning(
      space, ", so the smallest attainable p-value is ",
      signif(p_floor, 3), " and p < 0.05 is unreachable for any ",
      "signal: ", paste(remedy, collapse = "; "),
      ". See $null_distribution."
    )
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
    contributes = contributes,
    block = unname(blocks$membership[pairs$sp1]),
    swappable = blocks$n[blocks$membership[pairs$sp1]] > 1L,
    stringsAsFactors = FALSE
  )
  target_is_sp1 <- group[pairs$sp1] == target_group
  pair_sizes$n_hogs_target <- ifelse(target_is_sp1,
    pair_sizes$n_hogs_sp1,
    pair_sizes$n_hogs_sp2
  )
  pair_sizes$n_hogs_partner <- ifelse(target_is_sp1,
    pair_sizes$n_hogs_sp2,
    pair_sizes$n_hogs_sp1
  )
  # Keyed on `contributes`, not `swappable`: a pair can feed the
  # statistic while its swap is the identity, and it still has a
  # well-defined target side.
  pair_sizes$n_hogs_target[!pair_sizes$contributes] <- NA_integer_
  pair_sizes$n_hogs_partner[!pair_sizes$contributes] <- NA_integer_
  rownames(pair_sizes) <- NULL
  # Score the asymmetry rather than demanding a clean sweep. A sweep rule
  # has sensitivity P(target larger)^k, which *falls* as pairs are added
  # -- quietest at the k >= 5 where p_min finally allows a significant
  # result, and a single tie disarms it outright. A sign test over the
  # non-tied swappable pairs is the same question asked properly.
  # Both counts must range over the same rows. Once contrasts can be
  # coupled, a pinned component can hold a contributing contrast with
  # unequal sides: that row would feed n_larger but not n_nontied, and
  # binom.test() then stops with x > n, killing a complete analysis for
  # the sake of an advisory diagnostic.
  larger <- pair_sizes$n_hogs_target > pair_sizes$n_hogs_partner
  tied <- pair_sizes$n_hogs_target == pair_sizes$n_hogs_partner
  countable <- !tied & pair_sizes$swappable & !is.na(tied)
  n_larger <- sum(larger & countable, na.rm = TRUE)
  n_nontied <- sum(countable, na.rm = TRUE)
  size_p <- if (n_nontied > 0L) {
    stats::binom.test(n_larger, n_nontied, 0.5,
      alternative = "greater"
    )$p.value
  } else {
    NA_real_
  }
  # 0.10, not 0.05: this is an advisory diagnostic on the null's
  # assumption, not a test, and like p_min it cannot reach 0.05 below
  # five pairs -- a clean sweep is p = 0.125 at k = 3 and 0.0625 at
  # k = 4. Read $size_asymmetry_p directly at small k.
  if (!is.na(size_p) && size_p <= 0.10) {
    warning(
      "the ", target_group, " side holds the larger HOG set in ",
      n_larger, " of ", n_nontied, " swappable pairs with a size ",
      "difference (sign test p = ", signif(size_p, 3), "), so the ",
      "within-pair swap may not be exchangeable and this test may ",
      "be reading set size rather than recurrence; see $pair_sizes"
    )
  }

  list(
    observed = observed_count,
    statistic = statistic,
    statistic_observed = obs,
    n_labellings = n_labellings,
    null_distribution = null_dist,
    p_value = p_value,
    p_min = p_min,
    p_attainable = p_attainable,
    exact = exact,
    n_swappable = k,
    n_contributing = n_contributing,
    pair_sizes = pair_sizes,
    size_asymmetry_p = size_p,
    recurrence_table = recurrence_table,
    target_group = target_group,
    min_recurrence = min_recurrence,
    n_perm = n_draw
  )
}
