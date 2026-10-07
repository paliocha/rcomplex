#' Identify hub genes within co-expression modules
#'
#' Computes within-module centrality for each gene and flags the top-ranked
#' genes as hubs.  Optionally maps genes to ortholog groups (HOGs) for
#' downstream conservation analysis with [classify_hub_conservation()].
#'
#' @section Centrality measures:
#' Centrality is computed on the **within-module subgraph** (edges between
#' genes in the same module only):
#' \describe{
#'   \item{degree}{Weighted degree (`igraph::strength`): sum of edge weights
#'     to other genes in the same module.}
#'   \item{betweenness}{Shortest-path betweenness using inverse edge weights
#'     as distances.  Identifies genes that bridge sub-clusters within a
#'     module.}
#'   \item{eigenvector}{Eigenvector centrality (`igraph::eigen_centrality`):
#'     high for genes connected to other high-centrality genes.}
#' }
#'
#' @section Tie-breaking cascade:
#' When genes share the same primary centrality score, hub selection uses a
#' biologically informed cascade (all available tiers evaluated):
#' \enumerate{
#'   \item Primary centrality (user-selected measure)
#'   \item Global weighted degree across the full network
#'   \item Alternative within-module centrality (betweenness if primary is
#'     degree; degree otherwise)
#'   \item Mean within-module edge weight (strength / degree)
#'   \item Per-gene conservation effect size (requires `comparison`)
#'   \item Per-HOG minimum q-value (requires `comparison`; lower = better)
#' }
#'
#' @param modules Output of [detect_modules()].
#' @param net Output of [compute_network()].
#' @param orthologs Optional data frame from [parse_orthologs()] with columns
#'   `gene1`, `gene2`, `hog`.  The function auto-detects which column
#'   matches the gene names in `modules`.  If `NULL`, the `hog` column in the
#'   result is all `NA`.
#' @param comparison Optional data frame: the `$results` element from
#'   [summarize_comparison()].  When provided, enables conservation-informed
#'   tie-breaking (tiers 5--6).  Must contain columns `gene1`, `gene2`,
#'   `hog`, `species1.effect_size`, `species2.effect_size`, plus at least one
#'   pair of q-value columns (`species1.q_value_con`/`species2.q_value_con` or
#'   the `.div` variants).
#' @param centrality Centrality measure: `"degree"` (default), `"betweenness"`,
#'   or `"eigenvector"`.
#' @param top_n Integer: flag the top N genes per module as hubs.  If `NULL`
#'   (default), uses `top_fraction` instead.
#' @param top_fraction Numeric in (0, 1): fraction of genes per module to flag
#'   as hubs (default 0.1).  Ignored when `top_n` is non-NULL.
#' @param min_module_size Integer: modules with fewer genes get
#'   `is_hub = FALSE` for all genes (default 3).
#'
#' @return A data frame with one row per gene, ordered by module then rank:
#'   \describe{
#'     \item{gene}{Gene identifier}
#'     \item{module}{Module ID (integer)}
#'     \item{degree}{Within-module weighted degree (`igraph::strength`)}
#'     \item{betweenness}{Within-module betweenness centrality}
#'     \item{eigenvector}{Within-module eigenvector centrality}
#'     \item{mean_edge_weight}{Mean weight of edges to other module members}
#'     \item{global_degree}{Weighted degree in the full (thresholded) network}
#'     \item{rank}{Rank within module by primary centrality
#'       (1 = highest; ties use `"min"`)}
#'     \item{is_hub}{`TRUE` if the gene is in the top slice after the
#'       6-tier tie-breaking cascade}
#'     \item{hog}{HOG identifier (`NA` if `orthologs` not provided or gene
#'       not in the ortholog table)}
#'   }
#'
#' @examples
#' \dontrun{
#' hubs <- identify_module_hubs(modules, net, orthologs,
#'   comparison = summary$results
#' )
#' hubs[hubs$is_hub, ]
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
identify_module_hubs <- function(modules, ...) {
  UseMethod("identify_module_hubs")
}

#' @rdname identify_module_hubs
#' @export
identify_module_hubs.default <- function(modules, net, orthologs = NULL,
                                         comparison = NULL,
                                         centrality = c(
                                           "degree", "betweenness",
                                           "eigenvector"
                                         ),
                                         top_n = NULL,
                                         top_fraction = 0.1,
                                         min_module_size = 3L, ...) {
  centrality <- match.arg(centrality)

  if (!is.list(modules) || is.null(modules$module_genes) ||
        is.null(modules$graph) || is.null(modules$modules)) {
    stop("modules must be output from detect_modules()")
  }
  if (!is.list(net) || is.null(net$network)) {
    stop("net must be output from compute_network()")
  }
  .net_check(net, net$threshold)
  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (top_n < 1L) stop("top_n must be >= 1")
  } else {
    if (top_fraction <= 0 || top_fraction >= 1) {
      stop("top_fraction must be in (0, 1)")
    }
  }
  min_module_size <- as.integer(min_module_size)

  g <- modules$graph

  # Pre-compute global weighted degree (tie-breaker tier 2)
  global_str <- igraph::strength(g)

  # Pre-compute conservation lookups if comparison provided (tiers 5-6)
  gene_conserv <- NULL
  hog_min_q <- NULL
  if (!is.null(comparison)) {
    if (!all(c(
      "gene1", "gene2", "hog",
      "species1.effect_size", "species2.effect_size"
    ) %in%
      names(comparison))) {
      stop("comparison must be $results from summarize_comparison()")
    }
    # Auto-detect which column has our genes
    all_genes <- names(modules$modules)
    in_sp1 <- sum(all_genes %in% comparison$gene1)
    in_sp2 <- sum(all_genes %in% comparison$gene2)
    comp_col <- if (in_sp1 >= in_sp2) "gene1" else "gene2"

    # Per-row geometric mean of effect sizes
    geo_eff <- sqrt(comparison$species1.effect_size *
                      comparison$species2.effect_size)

    # Per-gene mean conservation effect (higher = more conserved)
    comp_genes <- comparison[[comp_col]]
    gene_conserv <- vapply(
      split(geo_eff, comp_genes), mean, numeric(1),
      na.rm = TRUE
    )

    # Per-HOG minimum q-value (lower = more conserved)
    q1_col <- if ("species1.q_value_con" %in% names(comparison)) {
      "species1.q_value_con"
    } else if ("species1.q_value_div" %in% names(comparison)) {
      "species1.q_value_div"
    } else {
      NULL
    }
    q2_col <- if ("species2.q_value_con" %in% names(comparison)) {
      "species2.q_value_con"
    } else if ("species2.q_value_div" %in% names(comparison)) {
      "species2.q_value_div"
    } else {
      NULL
    }
    if (!is.null(q1_col) && !is.null(q2_col)) {
      pair_q <- pmin(comparison[[q1_col]], comparison[[q2_col]], na.rm = TRUE)
      hog_min_q <- vapply(
        split(pair_q, comparison$hog), min, numeric(1),
        na.rm = TRUE
      )
    }
  }

  # HOG mapping (needed for tier-6 tie-breaking and output)
  hog_lookup <- NULL
  if (!is.null(orthologs)) {
    if (!all(c("gene1", "gene2", "hog") %in% names(orthologs))) {
      stop("orthologs must have columns: gene1, gene2, hog")
    }
    all_genes <- names(modules$modules)
    in_sp1 <- sum(all_genes %in% orthologs$gene1)
    in_sp2 <- sum(all_genes %in% orthologs$gene2)
    gene_col <- if (in_sp1 >= in_sp2) "gene1" else "gene2"

    gene_hog <- unique(orthologs[, c(gene_col, "hog"), drop = FALSE])
    gene_hog <- gene_hog[!duplicated(gene_hog[[gene_col]]), , drop = FALSE]
    hog_lookup <- stats::setNames(
      as.character(gene_hog$hog), gene_hog[[gene_col]]
    )
  }

  rows <- vector("list", length(modules$module_genes))

  for (i in seq_along(modules$module_genes)) {
    mod_id <- names(modules$module_genes)[i]
    genes <- modules$module_genes[[i]]
    n_genes <- length(genes)

    if (n_genes < min_module_size) {
      rows[[i]] <- data.frame(
        gene = genes, module = as.integer(mod_id),
        degree = NA_real_, betweenness = NA_real_, eigenvector = NA_real_,
        mean_edge_weight = NA_real_,
        global_degree = global_str[genes],
        rank = NA_integer_, is_hub = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    sub <- igraph::induced_subgraph(g, genes)
    w <- igraph::E(sub)$weight
    inv_w <- if (!is.null(w)) 1 / w else NULL

    # Compute all three centrality measures once
    sub_str <- igraph::strength(sub)
    sub_btw <- igraph::betweenness(sub, weights = inv_w)
    sub_eig <- tryCatch(
      igraph::eigen_centrality(sub, weights = w)$vector,
      error = function(e) {
        warning(
          "eigen_centrality failed for module ", mod_id, ": ",
          conditionMessage(e), "; using zero fallback"
        )
        stats::setNames(rep(0, length(genes)), genes)
      }
    )

    # Primary centrality for ranking/tie-breaking (tier 1)
    cent_vals <- switch(centrality,
      degree = sub_str,
      betweenness = sub_btw,
      eigenvector = sub_eig
    )
    # Mean within-module edge weight (tier 4): strength / degree
    sub_deg <- igraph::degree(sub)
    mean_ew <- ifelse(sub_deg > 0, sub_str / sub_deg, 0)

    rnk <- rank(-cent_vals, ties.method = "min")

    rows[[i]] <- data.frame(
      gene = names(cent_vals), module = as.integer(mod_id),
      degree = as.numeric(sub_str),
      betweenness = as.numeric(sub_btw),
      eigenvector = as.numeric(sub_eig),
      mean_edge_weight = as.numeric(mean_ew),
      global_degree = as.numeric(global_str[genes]),
      rank = as.integer(rnk),
      is_hub = FALSE, # filled below
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  # HOG mapping (vectorized, once for all genes)
  result$hog <- NA_character_
  if (!is.null(hog_lookup)) {
    matched <- match(result$gene, names(hog_lookup))
    result$hog[!is.na(matched)] <- hog_lookup[matched[!is.na(matched)]]
  }

  # Conservation lookups (vectorized, once for all genes — tiers 5-6)
  result$conserv_eff <- 0
  result$hog_q <- 1
  if (!is.null(gene_conserv)) {
    matched <- match(result$gene, names(gene_conserv))
    result$conserv_eff[!is.na(matched)] <-
      gene_conserv[matched[!is.na(matched)]]
  }
  if (!is.null(hog_min_q) && !is.null(hog_lookup)) {
    matched <- match(result$hog, names(hog_min_q))
    result$hog_q[!is.na(matched)] <- hog_min_q[matched[!is.na(matched)]]
  }

  # Primary and alternative centrality column names for tie-breaking
  primary_col <- centrality # "degree", "betweenness", or "eigenvector"
  alt_col <- if (centrality == "degree") "betweenness" else "degree"

  # Hub selection per module: tie-breaking cascade across all 6 tiers
  # Tiers 1-5 descending (higher = better), tier 6 ascending (lower = better)
  mod_ids <- unique(result$module[!is.na(result$degree)])
  for (m in mod_ids) {
    idx <- which(result$module == m & !is.na(result$degree))
    n_mod <- length(idx)
    hub_cutoff <- if (!is.null(top_n)) {
      min(top_n, n_mod)
    } else {
      max(1L, ceiling(top_fraction * n_mod))
    }
    ord <- order(
      -result[[primary_col]][idx], -result$global_degree[idx],
      -result[[alt_col]][idx], -result$mean_edge_weight[idx],
      -result$conserv_eff[idx], result$hog_q[idx]
    )
    result$is_hub[idx[ord[seq_len(hub_cutoff)]]] <- TRUE
  }

  # Drop internal tie-breaking columns (conservation lookups)
  result$conserv_eff <- NULL
  result$hog_q <- NULL

  attr(result, "primary_centrality") <- centrality
  result
}


#' Classify hub gene conservation across species and traits
#'
#' Given per-species hub identification results (from
#' [identify_module_hubs()]), maps hub genes to HOGs and classifies each HOG
#' by its hub conservation pattern relative to a discrete trait (e.g.
#' annual / perennial).
#'
#' @section Classification waterfall:
#' For each HOG that appears in at least one species:
#' \describe{
#'   \item{conserved_hub}{Hub in multiple trait groups **and** the hub modules
#'     correspond across traits (checked via `module_comparisons`).}
#'   \item{rewired_hub}{Hub in multiple trait groups but in
#'     **non-corresponding** modules -- the gene kept its centrality but
#'     changed regulatory context.}
#'   \item{multi_trait_hub}{Hub in multiple trait groups; module correspondence
#'     unknown (`module_comparisons` not provided).}
#'   \item{\emph{trait}_specific_hub}{Hub in exactly one trait group (e.g.
#'     `"annual_specific_hub"`).}
#'   \item{sporadic_hub}{Hub in some species but does not reach
#'     `min_trait_fraction` in any trait group.}
#'   \item{non_hub}{Present in modules but not a hub in any species.}
#' }
#'
#' @param hub_results Named list keyed by species name.  Each element is the
#'   data frame output of [identify_module_hubs()] (with `orthologs`
#'   provided so the `hog` column is populated).
#' @param species_trait Named character or factor vector mapping species to
#'   trait groups, e.g. `c(SP_A = "annual", SP_B = "annual",
#'   SP_C = "perennial", SP_D = "perennial")`.
#' @param module_comparisons Optional named list of
#'   [module_correspondence()] outputs keyed by alphabetically sorted species
#'   pair (e.g. `"SP_A.SP_C"`). Required for the conserved_hub vs rewired_hub
#'   distinction. Two things the caller must now satisfy themselves: each
#'   element must be built with the alphabetically first species as
#'   `modules_ref`, and [module_correspondence()] needs a map from
#'   [resolve_ortholog_map()], so the networks are required for the gene
#'   universes. Pass `sp_ref` / `sp_test` to [module_correspondence()] and
#'   that orientation is checked here instead of taken on trust. A module
#'   pair absent from the table counts as not corresponding.
#' @param alpha Significance threshold for module correspondence (default 0.1).
#' @param jaccard_threshold Jaccard threshold for module correspondence
#'   (default 0.1). [module_correspondence()] computes this over the
#'   one-to-one paralog-resolved projection, whereas the retired gene-overlap
#'   engine used the paralog-expanded mappable set, so values now run
#'   systematically higher and the unchanged default is slightly more
#'   permissive.
#' @param min_trait_fraction Minimum fraction of species (within a trait group)
#'   where the HOG must be a hub for the group to count (default 0.5). The
#'   denominator is the size of the trait group -- every species of that
#'   group in `hub_results` -- not just the ones carrying the HOG, so a HOG
#'   confined to one of four annuals cannot reach 0.5 there. Lower the
#'   threshold to admit accessory HOGs.
#' @param correspondence_threshold Fraction of cross-trait hub pairs that must
#'   have corresponding modules for the HOG to be classified as
#'   `conserved_hub` rather than `rewired_hub` (default 0.5).
#'
#' @return A data frame with one row per HOG:
#'   \describe{
#'     \item{hog}{HOG identifier}
#'     \item{classification}{Conservation category (see Classification
#'       waterfall)}
#'     \item{n_species_hub}{Number of species where the HOG is a hub}
#'     \item{n_species_present}{Number of species where the HOG has genes.
#'       Reported for context; it is not the `min_trait_fraction`
#'       denominator.}
#'     \item{hub_trait_groups}{Comma-separated trait groups where it qualifies
#'       as hub (`NA` for non_hub)}
#'     \item{n_corresponding}{Cross-trait hub pairs with corresponding modules
#'       (`NA` without `module_comparisons`)}
#'     \item{n_cross_pairs}{Total cross-trait hub pairs checked (`NA` without
#'       `module_comparisons`)}
#'     \item{max_centrality}{Highest centrality score across species}
#'     \item{best_hub_species}{Species with highest centrality}
#'   }
#'
#' @examples
#' \dontrun{
#' hub_list <- list(
#'   SP_A = identify_module_hubs(mods_A, net_A, ortho_A),
#'   SP_B = identify_module_hubs(mods_B, net_B, ortho_B)
#' )
#' trait <- c(SP_A = "annual", SP_B = "perennial")
#' classify_hub_conservation(hub_list, trait)
#'
#' # With module correspondence, for the conserved_hub / rewired_hub split.
#' # The list key must be the alphabetically sorted species pair.
#' map <- resolve_ortholog_map(
#'   ortho_AB, rownames(net_A$network), rownames(net_B$network)
#' )
#' corr <- list(SP_A.SP_B = module_correspondence(
#'   mods_A, mods_B, map,
#'   sp_ref = "SP_A", sp_test = "SP_B"
#' ))
#' classify_hub_conservation(hub_list, trait, module_comparisons = corr)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
classify_hub_conservation <- function(hub_results, ...) {
  UseMethod("classify_hub_conservation")
}

#' @rdname classify_hub_conservation
#' @export
classify_hub_conservation.default <- function(hub_results, species_trait,
                                              module_comparisons = NULL,
                                              alpha = 0.1,
                                              jaccard_threshold = 0.1,
                                              min_trait_fraction = 0.5,
                                              correspondence_threshold = 0.5,
                                              ...) {
  # --- Validation ---
  if (!is.list(hub_results) || is.null(names(hub_results))) {
    stop("hub_results must be a named list keyed by species")
  }
  if (!is.character(species_trait) && !is.factor(species_trait)) {
    stop("species_trait must be a named character or factor vector")
  }
  if (is.null(names(species_trait))) {
    stop("species_trait must be a named vector")
  }
  missing_sp <- setdiff(names(hub_results), names(species_trait))
  if (length(missing_sp) > 0) {
    stop(
      "species_trait missing entries for: ",
      paste(missing_sp, collapse = ", ")
    )
  }
  req_cols <- c("gene", "module", "is_hub", "hog", "degree")
  for (sp in names(hub_results)) {
    if (!is.data.frame(hub_results[[sp]]) ||
          !all(req_cols %in% names(hub_results[[sp]]))) {
      stop(
        "hub_results[['", sp,
        "']] must be output from identify_module_hubs() with orthologs"
      )
    }
  }

  trait_char <- as.character(species_trait[names(hub_results)])
  names(trait_char) <- names(hub_results)
  species_by_trait <- split(names(trait_char), trait_char)

  # Determine which centrality column to use for max_centrality / hub_module.
  # Reads the attribute set by identify_module_hubs(); falls back to "degree".
  primary_col <- unique(vapply(hub_results, function(hr) {
    pc <- attr(hr, "primary_centrality")
    if (is.null(pc)) "degree" else pc
  }, character(1)))
  if (length(primary_col) != 1L) primary_col <- "degree"

  # --- Build HOG-level summary: stack all results, aggregate per (hog, sp) ---
  tagged <- lapply(names(hub_results), function(sp) {
    hr <- hub_results[[sp]]
    hr <- hr[!is.na(hr$hog), , drop = FALSE]
    if (nrow(hr) == 0L) {
      return(NULL)
    }
    hr$species <- sp
    hr
  })
  stacked <- do.call(rbind, tagged)

  # Empty result template
  empty <- data.frame(
    hog = character(0), classification = character(0),
    n_species_hub = integer(0), n_species_present = integer(0),
    hub_trait_groups = character(0),
    n_corresponding = integer(0), n_cross_pairs = integer(0),
    max_centrality = numeric(0), best_hub_species = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(stacked) || nrow(stacked) == 0L) {
    return(empty)
  }

  # A species with no HOG-mapped gene still sits in the min_trait_fraction
  # denominator, where it counts as non-hub for every HOG and depresses the
  # whole trait group -- silently, if the caller forgot its ortholog table.
  no_hog <- setdiff(names(hub_results), unique(stacked$species))
  if (length(no_hog) > 0L) {
    warning(
      "no HOG-mapped genes in hub_results for: ",
      paste(no_hog, collapse = ", "),
      "; they count as non-hub in their trait group"
    )
  }

  # One row per (hog, species): is_hub (OR), hub_module, max_centrality
  hog_df <- do.call(rbind, lapply(
    split(stacked, paste(stacked$hog, stacked$species, sep = "\x01")),
    function(df) {
      any_hub <- any(df$is_hub)
      hub_mod <- if (any_hub) {
        hub_rows <- df[df$is_hub, , drop = FALSE]
        hub_rows$module[which.max(hub_rows[[primary_col]])]
      } else {
        NA_integer_
      }
      cent <- df[[primary_col]][!is.na(df[[primary_col]])]
      data.frame(
        hog = df$hog[1], species = df$species[1],
        is_hub = any_hub, hub_module = hub_mod,
        max_centrality = if (length(cent) == 0L) NA_real_ else max(cent),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(hog_df) <- NULL

  # --- Pre-compute (hog x trait) hub fraction matrix ---
  # The denominator is the size of the trait group, not the number of its
  # species that carry the HOG. Dividing by the species present scores a HOG
  # seen in one annual and a hub there as 1.0 -- the same as a hub in all
  # four annuals -- so accessory HOGs are not comparable with core ones and
  # sporadic_hub is unreachable for a HOG present in a single species.
  hog_df$trait <- trait_char[hog_df$species]
  hub_n <- tapply(hog_df$is_hub, list(hog_df$hog, hog_df$trait), sum)
  hub_n[is.na(hub_n)] <- 0
  group_n <- lengths(species_by_trait)[colnames(hub_n)]
  hub_frac <- sweep(hub_n, 2L, group_n, "/")
  is_hub_group <- hub_frac >= min_trait_fraction # logical matrix

  # --- Pre-compute per-HOG aggregates ---
  hog_n_present <- tapply(hog_df$species, hog_df$hog, length)
  hog_n_hub <- tapply(hog_df$is_hub, hog_df$hog, sum)
  hog_max_cent <- tapply(hog_df$max_centrality, hog_df$hog, function(x) {
    cx <- x[!is.na(x)]
    if (length(cx) == 0L) NA_real_ else max(cx)
  })
  hog_best_sp <- tapply(
    seq_len(nrow(hog_df)), hog_df$hog,
    function(idx) {
      sub <- hog_df[idx, , drop = FALSE]
      cx <- sub$max_centrality
      if (all(is.na(cx))) sub$species[1] else sub$species[which.max(cx)]
    }
  )

  # --- Pre-build module correspondence lookup per species pair ---
  # Validate first: an element without a $pairs data frame (what
  # preservation_paired()$raw gives you) would otherwise leave is_match as
  # logical(0) and report every HOG as NA, indistinguishable from having
  # supplied no comparison at all.
  if (!is.null(module_comparisons)) {
    # Keys first. An unnamed list makes the loops below iterate over NULL and
    # a wrongly-ordered key never matches the sorted lookup, and both leave
    # every HOG at NA -- indistinguishable from supplying no comparison, which
    # is the failure this guard exists to prevent.
    nm <- names(module_comparisons)
    if (is.null(nm) || !all(nzchar(nm))) {
      stop(
        "module_comparisons must be a named list keyed by ",
        "alphabetically sorted species pair (e.g. \"SP_A.SP_C\")"
      )
    }
    known_sp <- names(species_trait)
    valid_keys <- if (length(known_sp) >= 2L) {
      apply(utils::combn(sort(known_sp), 2L), 2L, paste, collapse = ".")
    } else {
      character(0)
    }
    bad_keys <- setdiff(nm, valid_keys)
    if (length(bad_keys) > 0L) {
      stop(
        "module_comparisons keys must be alphabetically sorted species ",
        "pairs drawn from species_trait; unusable: ",
        paste(bad_keys, collapse = ", ")
      )
    }
    # Orientation, when the producer recorded it. A transposed call --
    # module_correspondence(mods_B, mods_A, ...) filed under "A.B" -- passes
    # the name and shape checks and then matches lookups with module_sp1 and
    # module_sp2 swapped, giving wrong verdicts rather than a detectable NA.
    for (k in nm) {
      ref <- module_comparisons[[k]]$sp_ref
      if (is.null(ref)) next
      tst <- module_comparisons[[k]]$sp_test
      # Rebuild the key from the recorded labels rather than splitting it.
      # Splitting on "." mangles species names that contain one, and
      # comparing the pair as a whole also catches a sp_test naming a third
      # species, which a first-element check would pass.
      rebuilt <- if (is.null(tst)) NULL else paste(c(ref, tst), collapse = ".")
      ok <- if (is.null(rebuilt)) {
        startsWith(k, paste0(ref, "."))
      } else {
        identical(k, rebuilt)
      }
      if (!ok) {
        stop(
          "module_comparisons[[\"", k, "\"]] was built with sp_ref = \"",
          ref, "\"", if (!is.null(tst)) paste0(", sp_test = \"", tst, "\""),
          "; module_sp1 must belong to the first species of the key, so ",
          "the arguments or the key are wrong"
        )
      }
    }

    req_corr <- c("module_sp1", "module_sp2", "jaccard", "q_value")
    for (k in nm) {
      pk <- module_comparisons[[k]]$pairs
      if (!is.data.frame(pk) || !all(req_corr %in% names(pk))) {
        stop(
          "module_comparisons[[\"", k, "\"]] must be a ",
          "module_correspondence() result: a list with a `pairs` data ",
          "frame carrying ", paste(req_corr, collapse = ", ")
        )
      }
    }
  }

  corresp_lookup <- list() # keyed by "SP_A.SP_C", values = named logical
  if (!is.null(module_comparisons)) {
    for (pair_key in names(module_comparisons)) {
      pairs <- module_comparisons[[pair_key]]$pairs
      is_match <- pairs$q_value < alpha & pairs$jaccard >= jaccard_threshold
      keys <- paste(pairs$module_sp1, pairs$module_sp2, sep = "\x01")
      corresp_lookup[[pair_key]] <- stats::setNames(is_match, keys)
    }
  }

  # O(1) module correspondence check
  check_correspondence <- function(sp_a, mod_a, sp_b, mod_b) {
    if (length(corresp_lookup) == 0L) {
      return(NA)
    }
    pair_key <- paste(sort(c(sp_a, sp_b)), collapse = ".")
    lkp <- corresp_lookup[[pair_key]]
    if (is.null(lkp)) {
      return(NA)
    }
    sorted <- sort(c(sp_a, sp_b))
    mod_key <- if (sp_a == sorted[1]) {
      paste(mod_a, mod_b, sep = "\x01")
    } else {
      paste(mod_b, mod_a, sep = "\x01")
    }
    val <- lkp[mod_key]
    if (is.na(val)) FALSE else val
  }

  # --- Classify each HOG ---
  hog_groups <- split(hog_df, hog_df$hog)
  all_hogs <- names(hog_groups)

  out_rows <- lapply(all_hogs, function(h) {
    h_df <- hog_groups[[h]]
    n_present <- hog_n_present[[h]]
    n_hub <- hog_n_hub[[h]]
    max_cent <- hog_max_cent[[h]]
    best_sp <- hog_best_sp[[h]]

    # Trait-group hub status from pre-computed matrix (fix #6)
    hub_group_names <- colnames(is_hub_group)[is_hub_group[h, ]]

    n_corresponding <- NA_integer_
    n_cross_pairs <- NA_integer_

    if (n_hub == 0L) {
      classification <- "non_hub"
    } else if (length(hub_group_names) >= 2L) {
      # Hub in multiple trait groups -- check module correspondence
      hub_sp_by_group <- lapply(hub_group_names, function(g) {
        h_df$species[h_df$species %in% species_by_trait[[g]] & h_df$is_hub]
      })
      group_indices <- seq_along(hub_group_names)
      pair_mat <- if (length(group_indices) == 2L) {
        matrix(group_indices, nrow = 2)
      } else {
        utils::combn(group_indices, 2)
      }
      cross_pairs <- do.call(rbind, lapply(
        seq_len(ncol(pair_mat)), function(k) {
          expand.grid(
            sp_a = hub_sp_by_group[[pair_mat[1, k]]],
            sp_b = hub_sp_by_group[[pair_mat[2, k]]],
            stringsAsFactors = FALSE
          )
        }
      ))

      n_cross_pairs <- nrow(cross_pairs)

      corresp <- vapply(seq_len(n_cross_pairs), function(j) {
        mod_a <- h_df$hub_module[h_df$species == cross_pairs$sp_a[j]]
        mod_b <- h_df$hub_module[h_df$species == cross_pairs$sp_b[j]]
        check_correspondence(
          cross_pairs$sp_a[j], mod_a,
          cross_pairs$sp_b[j], mod_b
        )
      }, logical(1))

      if (all(is.na(corresp))) {
        classification <- "multi_trait_hub"
        n_corresponding <- NA_integer_
      } else {
        n_corresponding <- sum(corresp, na.rm = TRUE)
        n_available <- sum(!is.na(corresp))
        classification <- if (n_corresponding / n_available >=
                                correspondence_threshold) {
          "conserved_hub"
        } else {
          "rewired_hub"
        }
      }
    } else if (length(hub_group_names) == 1L) {
      classification <- paste0(hub_group_names, "_specific_hub")
    } else {
      classification <- "sporadic_hub"
    }

    data.frame(
      hog = h,
      classification = classification,
      n_species_hub = as.integer(n_hub),
      n_species_present = as.integer(n_present),
      hub_trait_groups = if (length(hub_group_names) > 0L) {
        paste(sort(hub_group_names), collapse = ",")
      } else {
        NA_character_
      },
      n_corresponding = n_corresponding,
      n_cross_pairs = n_cross_pairs,
      max_centrality = max_cent,
      best_hub_species = best_sp,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out_rows)
  rownames(result) <- NULL
  result
}
