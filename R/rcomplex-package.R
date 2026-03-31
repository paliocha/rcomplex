#' rcomplex: Comparative Co-Expression Network Analysis Across Species
#'
#' Compares gene co-expression networks across species by mapping orthologous
#' genes, building co-expression networks independently per species, then
#' testing conservation at gene, module, and clique levels.
#'
#' @section Main functions:
#' \describe{
#'   \item{[parse_orthologs()]}{Parse ortholog group files}
#'   \item{[reduce_orthogroups()]}{Merge correlated paralogs within HOGs}
#'   \item{[extract_orthologs()]}{Derive ortholog pairs from SummarizedExperiment objects}
#'   \item{[compute_network()]}{Build co-expression network (matrix or SummarizedExperiment)}
#'   \item{[compare_neighborhoods()]}{Pair-level hypergeometric tests}
#'   \item{[summarize_comparison()]}{Q-value correction and summary}
#'   \item{[comparison_to_edges()]}{Convert comparison results to clique edge format}
#'   \item{[run_pairwise_comparisons()]}{Batch pairwise comparison pipeline}
#'   \item{[permutation_hog_test()]}{HOG-level permutation test}
#'   \item{[detect_modules()]}{Community detection with consensus}
#'   \item{[compare_modules()]}{Cross-species module comparison}
#'   \item{[classify_modules()]}{Module conservation classification}
#'   \item{[find_cliques()]}{C++ clique detection (Bron-Kerbosch)}
#'   \item{[clique_stability()]}{Leave-k-out jackknife stability}
#'   \item{[clique_persistence()]}{Co-expressolog persistence scores}
#'   \item{[clique_threshold_sweep()]}{Threshold sweep (convenience wrapper)}
#'   \item{[classify_cliques()]}{HOG classification (convenience wrapper)}
#' }
#'
#' @docType package
#' @name rcomplex-package
#' @aliases rcomplex
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib rcomplex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods setGeneric setMethod is new
#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom utils read.delim
## usethis namespace: end
NULL
