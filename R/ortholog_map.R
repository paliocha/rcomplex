# Paralog-resolved ortholog maps for cross-species module projection.
#
# A multi-copy HOG offers several candidate gene pairs between two species.
# Projecting module labels through all of them lets co-expressed paralogs
# inflate within-module connectivity, so the copies are resolved to a single
# counterpart wherever the pipeline already carries evidence for one:
# clique membership first (globally consistent across every species at once),
# then a mutual-best coexpressolog pair, then nothing.
#
# The resolution layers may only choose WHICH copy carries a label; they must
# never change which genes are mappable at all.  Coexpressologs are defined by
# conserved co-expression neighbourhoods and module preservation measures
# conserved local topology, so filtering the mappable set on coexpressolog
# evidence would select the tested genes on the statistic being tested.
# `retain_uncovered()` below is what enforces that invariant.


#' Resolve paralog copies in an ortholog map
#'
#' Reduces the candidate gene pairs of multi-copy ortholog groups (HOGs) to a
#' single counterpart per gene, using evidence the pipeline has already
#' produced.  Intended for [module_preservation()], which projects module
#' labels from one species onto another and is sensitive to the paralog
#' expansion that a raw ortholog table implies.
#'
#' @section Resolution waterfall:
#' \describe{
#'   \item{cliques}{A clique fixes one gene per species simultaneously, so the
#'     copy choices it implies cannot contradict each other across a
#'     multi-species run the way independent pairwise matchings can.  A HOG
#'     appearing in several cliques resolves by `n_species` descending, then
#'     `mean_q` ascending.}
#'   \item{coexpressologs}{For genes no clique reached, the mutual-best
#'     significant pair: `gene1`'s highest-ranked partner must also rank
#'     `gene1` highest.}
#'   \item{unresolved}{Everything else keeps all its candidate pairs, left for
#'     the consumer to resolve (see the `source` column).}
#' }
#'
#' @section Preserved gene set:
#' Resolution never removes a species-2 gene from the map.  A candidate pair
#' survives when it was resolved, or when its species-2 gene was not claimed by
#' any resolved pair, so the set of mappable species-2 genes is identical to
#' the set implied by `orthologs` alone.  Only which species-1 copy points at
#' each gene changes.
#'
#' @section Ranking column:
#' `find_coexpressologs(method = "permutation")` computes q-values at the HOG
#' level and broadcasts them to every gene pair of the HOG, so `q.value` is
#' constant within a HOG and cannot discriminate paralogs.  `effect_size` and
#' `jaccard` stay pair-level under both methods and are the usable ranks;
#' `rank_by = "q.value"` errors when the supplied table has no within-HOG
#' q-value variation.
#'
#' @param orthologs Data frame with columns `Species1`, `Species2`, `hog`
#'   (output of [parse_orthologs()]).  `Species1` / `Species2` hold gene
#'   identifiers; species membership is resolved against `genes1` / `genes2`.
#' @param genes1,genes2 Character vectors giving the gene universes of the two
#'   species, e.g. `rownames(net$network)`.
#' @param sp1,sp2 Species labels.  Required only when `edges` or `cliques` is
#'   supplied, to select clique columns and filter edges.
#' @param edges Optional coexpressolog edge table from
#'   [find_coexpressologs()], with columns `gene1`, `gene2`, `species1`,
#'   `species2`, `hog` and the column named by `rank_by`.
#' @param cliques Optional clique table from [find_cliques()]: `hog`, one
#'   column per species holding a gene identifier or `NA`, and `n_species`.
#' @param rank_by Column ranking coexpressolog partners: `"effect_size"`
#'   (default), `"jaccard"`, or `"q.value"` (lower is better).
#' @param alpha Significance threshold applied to `edges$q.value` when the
#'   table has no `type` column (default 0.05).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{gene1}{Species-1 gene identifier}
#'     \item{gene2}{Species-2 gene identifier}
#'     \item{hog}{Ortholog group identifier}
#'     \item{source}{`"clique"`, `"coexpressolog"`, or `"unresolved"`}
#'   }
#'
#' @examples
#' \dontrun{
#' map <- resolve_ortholog_map(
#'   ortho, rownames(net_a$network), rownames(net_b$network),
#'   sp1 = "SP_A", sp2 = "SP_B", edges = rcx$edges, cliques = rcx$cliques
#' )
#' table(map$source)
#' }
#'
#' @export
resolve_ortholog_map <- function(orthologs, genes1, genes2,
                                 sp1 = NULL, sp2 = NULL,
                                 edges = NULL, cliques = NULL,
                                 rank_by = c(
                                   "effect_size", "jaccard",
                                   "q.value"
                                 ),
                                 alpha = 0.05) {
  rank_by <- match.arg(rank_by)

  if (!is.data.frame(orthologs)) {
    stop("orthologs must be a data.frame")
  }
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2, hog")
  }
  if (!is.character(genes1) || !is.character(genes2)) {
    stop("genes1 and genes2 must be character vectors")
  }
  if ((!is.null(edges) || !is.null(cliques)) &&
    (is.null(sp1) || is.null(sp2))) {
    stop("sp1 and sp2 are required when edges or cliques is supplied")
  }

  # Candidate pairs: the same orientation rule the rest of the package uses
  # (compare_neighborhoods() filters Species1 against net1, Species2 against
  # net2), so a table whose orientation is flipped contributes nothing.
  cand <- orthologs[orthologs$Species1 %in% genes1 &
    orthologs$Species2 %in% genes2, , drop = FALSE]
  cand <- unique(cand[, c("Species1", "Species2", "hog"), drop = FALSE])
  names(cand) <- c("gene1", "gene2", "hog")
  cand$hog <- as.character(cand$hog)

  if (nrow(cand) == 0L) {
    stop("No orthologs found in both gene universes")
  }

  cand_key <- paste(cand$gene1, cand$gene2, sep = "\x01")

  resolved <- .map_clique_layer(cliques, sp1, sp2, cand, cand_key)

  coexpr <- .map_coexpressolog_layer(
    edges, sp1, sp2, cand, cand_key, resolved$gene1, rank_by, alpha
  )
  resolved <- rbind(resolved, coexpr)

  # Preserved-gene-set invariant: keep every candidate pair whose species-2
  # gene no resolved pair claims, so resolution only redirects copies.
  uncovered <- cand[!cand$gene2 %in% resolved$gene2, , drop = FALSE]
  if (nrow(uncovered) > 0L) {
    uncovered$source <- "unresolved"
  }

  out <- rbind(resolved, uncovered)
  out <- out[order(out$hog, out$gene1, out$gene2), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Clique layer of the resolution waterfall (internal)
#'
#' One gene per species per clique, so the copy choice is consistent across
#' every species at once.  Only pairs that are genuine ortholog candidates are
#' kept: a clique edge whose gene pair is absent from `cand` (different HOG
#' assignment, or a gene outside the network) is dropped.
#'
#' @noRd
.map_clique_layer <- function(cliques, sp1, sp2, cand, cand_key) {
  empty <- cand[0, , drop = FALSE]
  empty$source <- character(0)
  if (is.null(cliques)) {
    return(empty)
  }

  if (!is.data.frame(cliques) || !"hog" %in% names(cliques)) {
    stop("cliques must be a data.frame with a 'hog' column")
  }
  missing_cols <- setdiff(c(sp1, sp2), names(cliques))
  if (length(missing_cols) > 0L) {
    stop(
      "cliques has no column for species: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  cl <- cliques[!is.na(cliques[[sp1]]) & !is.na(cliques[[sp2]]), ,
    drop = FALSE
  ]
  if (nrow(cl) == 0L) {
    return(empty)
  }

  # Best clique per HOG: most species, then lowest mean q-value.
  n_species <- if ("n_species" %in% names(cl)) cl$n_species else 0L
  mean_q <- if ("mean_q" %in% names(cl)) cl$mean_q else 0
  cl <- cl[order(-n_species, mean_q, cl[[sp1]], cl[[sp2]]), , drop = FALSE]
  cl <- cl[!duplicated(as.character(cl$hog)), , drop = FALSE]

  hits <- data.frame(
    gene1 = as.character(cl[[sp1]]),
    gene2 = as.character(cl[[sp2]]),
    hog = as.character(cl$hog),
    stringsAsFactors = FALSE
  )
  hits <- hits[paste(hits$gene1, hits$gene2, sep = "\x01") %in% cand_key, ,
    drop = FALSE
  ]
  if (nrow(hits) == 0L) {
    return(empty)
  }

  hits$source <- "clique"
  hits
}


#' Coexpressolog layer of the resolution waterfall (internal)
#'
#' Mutual-best significant pairs among the genes no clique reached.  Mutual
#' rather than one-sided best so the choice does not depend on which species is
#' treated as the reference.
#'
#' @noRd
.map_coexpressolog_layer <- function(edges, sp1, sp2, cand, cand_key,
                                     done_genes, rank_by, alpha) {
  empty <- cand[0, , drop = FALSE]
  empty$source <- character(0)
  if (is.null(edges)) {
    return(empty)
  }

  if (!is.data.frame(edges)) {
    stop("edges must be a data.frame")
  }
  req <- c("gene1", "gene2", "species1", "species2", rank_by)
  missing_cols <- setdiff(req, names(edges))
  if (length(missing_cols) > 0L) {
    stop("edges missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Orient every edge so gene1 is the species-1 side.
  fwd <- edges$species1 == sp1 & edges$species2 == sp2
  rev <- edges$species1 == sp2 & edges$species2 == sp1
  e <- rbind(
    edges[fwd, , drop = FALSE],
    .swap_edge_genes(edges[rev, , drop = FALSE])
  )
  if (nrow(e) == 0L) {
    return(empty)
  }

  # Significance: prefer the explicit type label, fall back to q.value.
  if ("type" %in% names(e)) {
    e <- e[e$type == "conserved", , drop = FALSE]
  } else if ("q.value" %in% names(e)) {
    e <- e[!is.na(e$q.value) & e$q.value < alpha, , drop = FALSE]
  }

  # Only candidate pairs, and only genes the clique layer left open.
  e <- e[!e$gene1 %in% done_genes, , drop = FALSE]
  e <- e[paste(e$gene1, e$gene2, sep = "\x01") %in% cand_key, , drop = FALSE]
  e <- e[!is.na(e[[rank_by]]), , drop = FALSE]
  if (nrow(e) == 0L) {
    return(empty)
  }

  if (rank_by == "q.value") .check_pairwise_qvalues(e)

  # Higher is better for effect_size and jaccard, lower for q.value.
  score <- if (rank_by == "q.value") -e[[rank_by]] else e[[rank_by]]
  e <- e[order(-score, e$gene1, e$gene2), , drop = FALSE]

  best1 <- paste(e$gene1, e$gene2, sep = "\x01")[!duplicated(e$gene1)]
  best2 <- paste(e$gene1, e$gene2, sep = "\x01")[!duplicated(e$gene2)]
  mutual <- intersect(best1, best2)
  if (length(mutual) == 0L) {
    return(empty)
  }

  keep <- cand[cand_key %in% mutual, , drop = FALSE]
  if (nrow(keep) == 0L) {
    return(empty)
  }
  keep$source <- "coexpressolog"
  keep
}


#' Swap the two gene/species columns of an edge table (internal)
#' @noRd
.swap_edge_genes <- function(e) {
  if (nrow(e) == 0L) {
    return(e)
  }
  g1 <- e$gene1
  s1 <- e$species1
  e$gene1 <- e$gene2
  e$species1 <- e$species2
  e$gene2 <- g1
  e$species2 <- s1
  e
}


#' Refuse q-value ranking when q-values cannot discriminate paralogs (internal)
#'
#' `find_coexpressologs(method = "permutation")` assigns one q-value per HOG
#' and broadcasts it to every gene pair in that HOG, so ranking copies by
#' q-value would pick an arbitrary one.
#'
#' @noRd
.check_pairwise_qvalues <- function(e) {
  multi <- split(e$q.value, e$hog)
  multi <- multi[vapply(multi, length, integer(1)) > 1L]
  if (length(multi) == 0L) {
    return(invisible(NULL))
  }

  varies <- vapply(multi, function(q) {
    q <- q[!is.na(q)]
    length(q) > 1L && diff(range(q)) > 0
  }, logical(1))

  if (!any(varies)) {
    stop(
      "rank_by = \"q.value\" cannot resolve paralogs: q-values are ",
      "constant within every multi-copy HOG. This is expected from ",
      "find_coexpressologs(method = \"permutation\"), which computes ",
      "q-values at the HOG level. Use rank_by = \"effect_size\" or ",
      "\"jaccard\"."
    )
  }
  invisible(NULL)
}
