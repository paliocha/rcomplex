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
#' @param x One of: a named vector mapping gene to module label (`NA` for an
#'   unassigned gene); a named list mapping module label to a character
#'   vector of genes; or an existing module assignment (a [detect_modules()]
#'   or `as_modules()` result), which is returned unchanged.
#' @param min_size Modules with fewer genes than this become unassigned and
#'   are dropped from `module_genes`.
#'
#' @return A list with `modules` (named character vector, gene -> module
#'   label, unassigned genes omitted), `module_genes` (named list, module
#'   label -> genes), `n_modules`, `method = "external"` and `params`.
#' @export
#' @examples
#' mods <- as_modules(list(
#'   photosynthesis = c("g1", "g2", "g3"),
#'   ribosome = c("g4", "g5", "g6")
#' ))
#' mods$modules
as_modules <- function(x, min_size = 1L) {
  ok_size <- is.numeric(min_size) && length(min_size) == 1L &&
    !is.na(min_size) && min_size >= 1
  if (!ok_size) {
    stop("min_size must be a single number >= 1")
  }
  if (is.list(x) && !is.null(x$modules) && is.list(x$module_genes)) {
    if (min_size != 1) {
      stop("min_size applies to a vector or list of gene sets, ",
           "not to an existing module object")
    }
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
    levels <- as.character(lab)
    modules <- stats::setNames(as.character(x), names(x))
    modules <- modules[!is.na(modules)]
    if (any(!nzchar(modules))) {
      stop("empty module label for gene ", names(modules)[!nzchar(modules)][1])
    }
  }
  size <- table(factor(modules, levels = levels))
  levels <- levels[size >= min_size]
  modules <- modules[modules %in% levels]
  module_genes <- split(names(modules), factor(modules, levels = levels))
  list(
    modules = modules,
    module_genes = module_genes,
    n_modules = length(module_genes),
    method = "external",
    params = list(min_size = min_size)
  )
}
