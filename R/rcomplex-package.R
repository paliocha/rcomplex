#' rcomplex: Comparative Co-Expression Network Analysis Across Species
#'
#' Compares gene co-expression networks across plant species by mapping
#' orthologous genes, building co-expression networks independently per species,
#' then testing whether network neighborhoods are significantly preserved using
#' hypergeometric tests with FDR correction and effect sizes.
#'
#' @section Main functions:
#' - [parse_orthologs()]: Parse PLAZA ortholog group files
#' - [compute_network()]: Build co-expression network from expression matrix
#' - [compare_neighborhoods()]: Test conservation of network neighborhoods
#' - [summarize_comparison()]: FDR correction and summary statistics
#' - [permutation_hog_test()]: Permutation-based HOG-level conservation test
#'
#' @section Typical workflow:
#' 1. Parse ortholog groups with [parse_orthologs()]
#' 2. Build networks for each species with [compute_network()]
#' 3. Compare neighborhoods with [compare_neighborhoods()]
#' 4. Summarize results with [summarize_comparison()] or test HOG-level
#'    conservation with [permutation_hog_test()]
#'
#' @docType package
#' @name rcomplex-package
#' @aliases rcomplex
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib rcomplex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom igraph %--%
#' @importFrom rlang .data .env
#' @importFrom stats phyper p.adjust pchisq setNames
#' @importFrom utils read.delim
## usethis namespace: end
NULL
