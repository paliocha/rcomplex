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
      valid <- igraph::cliques(g, min = k, max = k) |>
        Filter(\(cl) setequal(igraph::V(g)$species[cl], target_species), x = _)
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
