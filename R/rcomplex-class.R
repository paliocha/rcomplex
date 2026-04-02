# S3 container class for rcomplex pipeline
#
# A validated named list that threads networks, orthologs, species,
# traits, and results between pipeline functions.


#' Create an rcomplex analysis object
#'
#' Constructs a container that carries species metadata, networks,
#' orthologs, and analysis results through the rcomplex pipeline.
#' Pipeline functions gain S3 methods that read inputs from the object
#' and store results back, eliminating repetitive argument passing.
#'
#' @param species Character vector of species abbreviations.
#' @param traits Named character vector mapping species to trait labels
#'   (e.g., \code{c(BDIS = "annual", BSYL = "perennial")}).
#' @param networks Named list of \code{\link{compute_network}} outputs,
#'   keyed by species.
#' @param orthologs Data frame with columns \code{Species1},
#'   \code{Species2}, \code{hog}.
#' @param species_pairs Optional list of length-2 character vectors
#'   specifying species pairs for gene-level comparison. Defaults to
#'   all pairwise combinations.
#' @param phylo_pairs Optional data frame with columns \code{sp1},
#'   \code{sp2} (and optionally \code{pair_name}) for module-level
#'   comparison.
#'
#' @return An S3 object of class \code{"rcomplex"}.
#'
#' @examples
#' \dontrun{
#' rcx <- rcomplex(
#'   species = c("SP_A", "SP_B", "SP_C"),
#'   traits = c(SP_A = "annual", SP_B = "perennial", SP_C = "annual"),
#'   networks = list(SP_A = net_a, SP_B = net_b, SP_C = net_c),
#'   orthologs = ortho
#' )
#' rcx <- find_coexpressologs(rcx, method = "permutation")
#' rcx <- find_cliques(rcx)
#' rcx$edges     # access results directly
#' rcx$cliques
#' }
#'
#' @export
rcomplex <- function(species, traits, networks, orthologs,
                     species_pairs = NULL, phylo_pairs = NULL) {

  # --- Validate species ---
  if (!is.character(species) || length(species) < 2L)
    stop("species must be a character vector with at least 2 elements")
  if (anyDuplicated(species))
    stop("species must not contain duplicates")

  # --- Validate traits ---
  if (!is.character(traits) && !is.factor(traits))
    stop("traits must be a named character or factor vector")
  if (is.null(names(traits)))
    stop("traits must be named")
  missing_sp <- setdiff(species, names(traits))
  if (length(missing_sp) > 0L)
    stop("traits missing species: ", paste(missing_sp, collapse = ", "))

  # --- Validate networks ---
  if (!is.list(networks) || is.null(names(networks)))
    stop("networks must be a named list")
  if (!setequal(names(networks), species))
    stop("network names must match species")
  for (sp in species) {
    net <- networks[[sp]]
    if (!is.list(net) || is.null(net$network) || is.null(net$threshold))
      stop("networks[['", sp, "']] must have 'network' and 'threshold'")
  }

  # --- Validate orthologs ---
  if (!is.data.frame(orthologs))
    stop("orthologs must be a data.frame")
  if (!all(c("Species1", "Species2", "hog") %in% names(orthologs)))
    stop("orthologs must have columns: Species1, Species2, hog")

  # --- Validate / default species_pairs ---
  if (is.null(species_pairs)) {
    species_pairs <- utils::combn(species, 2, simplify = FALSE)
  } else {
    if (!is.list(species_pairs))
      stop("species_pairs must be a list of length-2 character vectors")
    for (i in seq_along(species_pairs)) {
      p <- species_pairs[[i]]
      if (!is.character(p) || length(p) != 2L)
        stop("each species_pair must be a length-2 character vector")
      bad <- setdiff(p, species)
      if (length(bad) > 0L)
        stop("species_pairs references unknown species: ",
             paste(bad, collapse = ", "))
    }
  }

  # --- Validate phylo_pairs ---
  if (!is.null(phylo_pairs)) {
    if (!is.data.frame(phylo_pairs))
      stop("phylo_pairs must be a data.frame")
    if (!all(c("sp1", "sp2") %in% names(phylo_pairs)))
      stop("phylo_pairs must have columns: sp1, sp2")
    bad <- setdiff(c(phylo_pairs$sp1, phylo_pairs$sp2), species)
    if (length(bad) > 0L)
      stop("phylo_pairs references unknown species: ",
           paste(unique(bad), collapse = ", "))
    if (!"pair_name" %in% names(phylo_pairs))
      phylo_pairs$pair_name <- paste(phylo_pairs$sp1, phylo_pairs$sp2,
                                      sep = ".")
  }

  structure(
    list(
      species       = species,
      traits        = traits[species],
      networks      = networks[species],  # enforce species order
      orthologs     = orthologs,
      species_pairs = species_pairs,
      phylo_pairs   = phylo_pairs,
      # Result slots — populated by pipeline steps
      edges              = NULL,
      sweep              = NULL,
      modules            = NULL,
      module_comparisons = NULL,
      hubs               = NULL,
      hub_classification = NULL,
      cliques            = NULL,
      stability          = NULL,
      classification     = NULL
    ),
    class = "rcomplex"
  )
}


#' @export
print.rcomplex <- function(x, ...) {
  trait_levels <- sort(unique(x$traits[x$species]))
  n_sp <- length(x$species)

  cat("rcomplex:", n_sp, "species,",
      length(trait_levels), "traits (",
      paste(trait_levels, collapse = ", "), ")\n")
  cat("  Networks:       ", n_sp, "/", n_sp, "species\n", sep = "")
  cat("  Orthologs:      ", nrow(x$orthologs), " pairs\n", sep = "")
  cat("  Contrasts:      ", length(x$species_pairs), " species pairs",
      if (!is.null(x$phylo_pairs))
        paste0(", ", nrow(x$phylo_pairs), " phylo pairs")
      else "",
      "\n", sep = "")

  slot_info <- function(val, label, detail_fn) {
    if (is.null(val)) {
      cat("  ", label, " [pending]\n", sep = "")
    } else {
      cat("  ", label, " ", detail_fn(val), "\n", sep = "")
    }
  }

  slot_info(x$edges, "Edges:         ", function(e)
    paste0(nrow(e), " (", sum(e$type == "conserved", na.rm = TRUE),
           " conserved)"))
  slot_info(x$modules, "Modules:       ", function(m)
    paste0(length(m), " species"))
  slot_info(x$module_comparisons, "Mod. comp.:    ", function(mc)
    paste0(nrow(mc$summary), " pairs"))
  slot_info(x$hubs, "Hubs:          ", function(h)
    paste0(length(h), " species"))
  slot_info(x$hub_classification, "Hub class.:    ", function(hc)
    paste0(nrow(hc), " HOGs"))
  slot_info(x$cliques, "Cliques:       ", function(cl)
    paste0(nrow(cl)))
  slot_info(x$stability, "Stability:     ", function(s) {
    n <- length(s$stability_class)
    mx <- if (n > 0L) max(s$stability_class) else 0L
    paste0(n, " cliques, max class=", mx)
  })
  slot_info(x$classification, "Classification:", function(cl)
    paste0(nrow(cl), " HOGs"))
  slot_info(x$sweep, "Sweep:         ", function(sw)
    paste0(nrow(sw), " multipliers"))

  invisible(x)
}


#' @export
summary.rcomplex <- function(object, ...) {
  x <- object
  out <- list(
    n_species   = length(x$species),
    trait_levels = sort(unique(x$traits[x$species])),
    n_orthologs = nrow(x$orthologs),
    n_pairs     = length(x$species_pairs)
  )
  if (!is.null(x$edges)) {
    out$n_edges     <- nrow(x$edges)
    out$n_conserved <- sum(x$edges$type == "conserved", na.rm = TRUE)
  }
  if (!is.null(x$cliques))
    out$n_cliques <- nrow(x$cliques)
  if (!is.null(x$classification))
    out$classification_table <- table(x$classification$classification)
  if (!is.null(x$modules))
    out$n_modules <- vapply(x$modules, function(m) m$n_modules, integer(1))

  class(out) <- "summary.rcomplex"
  out
}


#' @export
print.summary.rcomplex <- function(x, ...) {
  cat("rcomplex summary:", x$n_species, "species,",
      x$n_orthologs, "ortholog pairs\n")
  if (!is.null(x$n_edges))
    cat("  Edges:", x$n_edges, "(", x$n_conserved, "conserved )\n")
  if (!is.null(x$n_cliques))
    cat("  Cliques:", x$n_cliques, "\n")
  if (!is.null(x$classification_table)) {
    cat("  Classification:\n")
    print(x$classification_table)
  }
  if (!is.null(x$n_modules)) {
    cat("  Modules per species:",
        paste(names(x$n_modules), x$n_modules, sep = "=", collapse = ", "),
        "\n")
  }
  invisible(x)
}


# ---- S3 methods for pipeline functions ----

#' @export
find_coexpressologs.rcomplex <- function(networks, ...) {
  x <- networks
  x$edges <- find_coexpressologs.default(
    x$networks, x$orthologs,
    species_pairs = x$species_pairs,
    ...
  )
  x
}


#' @export
density_sweep.rcomplex <- function(networks, ...) {
  x <- networks
  x$sweep <- density_sweep.default(x$networks, x$orthologs, ...)
  x
}


#' @export
detect_modules.rcomplex <- function(net, ...) {
  x <- net
  x$modules <- stats::setNames(
    lapply(x$species, function(sp)
      detect_modules.default(x$networks[[sp]], ...)),
    x$species
  )
  x
}


#' @export
compare_modules_paired.rcomplex <- function(modules, ..., group = NULL) {
  x <- modules
  if (is.null(x$modules))
    stop("run detect_modules() first")
  if (is.null(x$phylo_pairs))
    stop("phylo_pairs not set; pass to rcomplex() constructor")
  x$module_comparisons <- compare_modules_paired.default(
    x$modules, x$orthologs,
    pairs = x$phylo_pairs,
    group = if (is.null(group)) x$traits else group,
    ...
  )
  x
}


#' @export
identify_module_hubs.rcomplex <- function(modules, ...) {
  x <- modules
  if (is.null(x$modules))
    stop("run detect_modules() first")
  x$hubs <- stats::setNames(
    lapply(x$species, function(sp)
      identify_module_hubs.default(
        x$modules[[sp]], x$networks[[sp]],
        orthologs = x$orthologs, ...)),
    x$species
  )
  x
}


#' @export
classify_hub_conservation.rcomplex <- function(hub_results, ...) {
  x <- hub_results
  if (is.null(x$hubs))
    stop("run identify_module_hubs() first")
  mod_raw <- if (!is.null(x$module_comparisons)) {
    x$module_comparisons$raw
  } else {
    NULL
  }
  x$hub_classification <- classify_hub_conservation.default(
    x$hubs, x$traits,
    module_comparisons = mod_raw,
    ...
  )
  x
}


#' @export
find_cliques.rcomplex <- function(edges, ...) {
  x <- edges
  if (is.null(x$edges))
    stop("run find_coexpressologs() first")
  x$cliques <- find_cliques.default(x$edges, x$species, ...)
  x
}


#' @export
clique_stability.rcomplex <- function(edges, ...) {
  x <- edges
  if (is.null(x$edges))
    stop("run find_coexpressologs() first")
  x$stability <- clique_stability.default(
    x$edges, x$species,
    species_trait = x$traits,
    full_cliques = x$cliques,
    ...
  )
  x
}


#' @export
classify_cliques.rcomplex <- function(edges, ...) {
  x <- edges
  if (is.null(x$edges))
    stop("run find_coexpressologs() first")
  x$classification <- classify_cliques.default(
    x$edges, x$species, x$traits,
    stability = x$stability,
    sweep = x$sweep,
    ...
  )
  x
}
