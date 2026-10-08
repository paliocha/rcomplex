#' Clades from a species tree
#'
#' Turns a species tree into the clade list that the classifiers take.
#' Each internal node gives one clade: the tip labels below it. The root
#' holds every species, so it carries no contrast and is dropped.
#'
#' @param phy A tree of class `phylo` from the ape package.
#' @param min_size Smallest clade to keep. Use `1L` to keep the tips.
#' @return A named list of character vectors, one per clade. A clade
#'   takes its node label. An empty or repeated label becomes
#'   `"node<k>"`, with `k` the node number. A tip takes its own label.
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- ape::read.tree(text = "((A,B)AB,(C,D)CD);")
#'   clades_from_tree(phy)
#' }
#' @export
clades_from_tree <- function(phy, min_size = 2L) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Install the ape package to read a tree.")
  }
  if (!inherits(phy, "phylo")) stop("phy must be an ape phylo object.")
  pp <- ape::prop.part(phy)
  tips <- attr(pp, "labels")
  lab <- phy$node.label
  if (is.null(lab)) lab <- character(length(pp))
  bad <- !nzchar(lab) | duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[bad] <- paste0("node", length(tips) + seq_along(pp))[bad]
  cl <- c(as.list(tips), lapply(pp, function(i) tips[sort(i)]))
  names(cl) <- c(tips, lab)
  n <- lengths(cl)
  cl[n >= min_size & n < length(tips)]
}

#' Check a clade list and restrict it to the analysed species
#'
#' Clades must be pairwise nested or disjoint. After restriction to
#' `species`, empty clades and repeated species sets are dropped; the
#' first name of a repeated set is kept. A species in no clade forms its
#' own clade, named by itself, so no clade may take that name. A message
#' names those species.
#' @noRd
.check_clades <- function(clades, species) {
  nm <- names(clades)
  if (!is.list(clades) || length(clades) == 0L || is.null(nm) ||
        !all(nzchar(nm)) || anyDuplicated(nm) > 0L) {
    stop("clades must be a named list with unique names.")
  }
  ok <- vapply(clades, function(v) {
    is.character(v) && length(v) > 0L && !anyNA(v) && !anyDuplicated(v)
  }, logical(1))
  if (!all(ok)) stop("clade ", nm[!ok][1L], " must hold unique species.")
  all_sp <- unique(unlist(clades))
  inc <- matrix(unlist(lapply(clades, `%in%`, x = all_sp)), length(all_sp))
  common <- crossprod(inc)
  size <- diag(common)
  bad <- common > 0 & common < outer(size, size, pmin) & upper.tri(common)
  if (any(bad)) {
    ij <- which(bad, arr.ind = TRUE)[1L, ]
    stop(
      "clades ", nm[ij[1L]], " and ", nm[ij[2L]],
      " overlap, but neither holds the other."
    )
  }
  clades <- Filter(length, lapply(clades, intersect, species))
  clades <- clades[!duplicated(lapply(clades, sort))]
  alone <- setdiff(species, unlist(clades))
  if (length(alone) > 0L) {
    message("Species ", toString(alone), " are in no clade; ",
            "each forms its own clade.")
  }
  clash <- intersect(names(clades), alone)
  if (length(clash) > 0L) {
    stop("clade names equal species outside every clade: ", clash[1L])
  }
  clades
}

#' Name of the smallest clade that holds every species in `sp`
#' @noRd
.clade_home <- function(clades, sp) {
  hit <- vapply(clades, function(v) all(sp %in% v), logical(1))
  if (!any(hit)) return(NA_character_)
  names(clades)[hit][which.min(lengths(clades)[hit])]
}

#' Top-level group of each species: its largest clade, or its own name
#' @noRd
.clade_groups <- function(clades, species) {
  out <- stats::setNames(species, species)
  for (k in names(clades)[order(lengths(clades))]) {
    out[intersect(clades[[k]], species)] <- k
  }
  out
}
