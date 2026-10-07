#' rcomplex: Comparative Co-Expression Network Analysis Across Species
#'
#' Compares gene co-expression networks across species by mapping orthologous
#' genes, building co-expression networks independently per species, then
#' testing conservation at gene, module, and clique levels.
#'
#' @docType package
#' @name rcomplex-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib rcomplex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods setGeneric setMethod is new
#' @importFrom rlang .data .env %||%
#' @importFrom stats setNames
#' @importFrom utils read.delim modifyList
## usethis namespace: end
NULL
