#' Build a module assignment from any gene partition
#'
#' Turns gene sets from any source -- anchored regulons, curated pathways,
#' clusters from another tool -- into the module assignment that
#' [module_preservation()], [module_correspondence()] and
#' [preservation_paired()] take, so they can be tested without
#' [detect_modules()].
#'
#' The preservation tests give every gene one module label, so the sets must
#' not overlap. Split overlapping sets (for example regulons that share genes)
#' into disjoint batches and test each batch separately.
#'
#' Module labels matter for reproducibility: [module_preservation()] draws
#' its permutation null module by module in the C-locale (radix) order of
#' the labels as strings (so "10" before "2"), so under a
#' fixed seed the permutation p- and q-values are reproducible for a given
#' labelling but change when the same partition is relabelled. The observed
#' statistics do not depend on the labels.
#'
#' @param x One of: a named vector mapping gene to module label (`NA` for an
#'   unassigned gene); a named list mapping module label to a character
#'   vector of genes; or an existing module assignment (a [detect_modules()]
#'   or `as_modules()` result), which is returned unchanged.
#'
#' @return A list with `modules` (named character vector, gene -> module
#'   label, unassigned genes omitted), `module_genes` (named list, module
#'   label -> genes), `n_modules`, `method = "external"` and `params`.
#'   `module_genes` follows the list order for list input, numeric order for
#'   numeric labels, and order of first appearance for other labels; the
#'   order changes only listings such as `coverage`, not statistics.
#' @export
#' @examples
#' mods <- as_modules(list(
#'   photosynthesis = c("g1", "g2", "g3"),
#'   ribosome = c("g4", "g5", "g6")
#' ))
#' mods$modules
as_modules <- function(x) {
  if (is.data.frame(x)) {
    stop(
      "pass a named vector, e.g. setNames(df$module, df$gene), ",
      "not a data frame"
    )
  }
  if (is.list(x) && !is.null(x$modules) && is.list(x$module_genes)) {
    return(x)
  }
  if (is.list(x)) {
    if (is.null(names(x)) || any(!nzchar(names(x))) || anyNA(names(x))) {
      stop("a list of gene sets must be named by module label")
    }
    if (anyDuplicated(names(x))) {
      stop("module label used twice: ", names(x)[anyDuplicated(names(x))])
    }
    if (!all(vapply(x, is.character, logical(1)))) {
      stop("every gene set must be a character vector of gene names")
    }
    x <- lapply(x, unique)
    genes <- unlist(x, use.names = FALSE)
    if (anyNA(genes) || any(!nzchar(genes))) {
      stop("gene names must not be NA or empty")
    }
    if (anyDuplicated(genes)) {
      dup <- unique(genes[duplicated(genes)])
      where <- vapply(utils::head(dup, 5L), function(g) {
        holders <- names(x)[vapply(x, `%in%`, x = g, logical(1))]
        paste0(g, " (", paste(holders, collapse = ", "), ")")
      }, character(1))
      stop(
        "gene sets overlap; the preservation tests need a partition. ",
        length(dup), " shared gene(s), e.g. ", paste(where, collapse = "; ")
      )
    }
    modules <- stats::setNames(rep(names(x), lengths(x)), genes)
    levels <- names(x)
  } else {
    if (is.null(names(x)) || any(!nzchar(names(x))) || anyNA(names(x))) {
      stop("a module vector must be named by gene")
    }
    if (anyDuplicated(names(x))) {
      stop("gene assigned twice: ", names(x)[anyDuplicated(names(x))])
    }
    # numeric labels sort numerically, as detect_modules() orders them;
    # other labels keep their order of first appearance
    lab <- x[!is.na(x)]
    lab <- if (is.numeric(lab)) sort(unique(lab)) else unique(lab)
    # unique after conversion: doubles that print alike share a label
    levels <- unique(as.character(lab))
    modules <- stats::setNames(as.character(x), names(x))
    modules <- modules[!is.na(modules)]
    if (any(!nzchar(modules))) {
      stop("empty module label for gene ", names(modules)[!nzchar(modules)][1])
    }
  }
  levels <- levels[levels %in% modules]
  module_genes <- split(names(modules), factor(modules, levels = levels))
  list(
    modules = modules,
    module_genes = module_genes,
    n_modules = length(module_genes),
    method = "external",
    params = list()
  )
}
