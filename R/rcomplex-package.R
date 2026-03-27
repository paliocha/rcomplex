#' rcomplex: Comparative Co-Expression Network Analysis Across Species
#'
#' Compares gene co-expression networks across species by mapping orthologous
#' genes, building co-expression networks independently per species, then
#' testing conservation at gene, module, and clique levels.
#'
#' @section Main functions:
#' \subsection{Preprocessing}{
#' - [parse_orthologs()]: Parse ortholog group files (tab-delimited)
#' - [reduce_orthogroups()]: Merge correlated paralogs within HOGs (Ward.D2)
#' - [compute_network()]: Build co-expression network from expression matrix
#' }
#' \subsection{Gene / HOG-level analysis}{
#' - [compare_neighborhoods()]: Test conservation of network neighborhoods
#' - [summarize_comparison()]: Q-value correction and summary statistics
#' - [permutation_hog_test()]: Permutation-based HOG-level conservation test
#' }
#' \subsection{Module-level analysis}{
#' - [detect_modules()]: Community detection (Leiden / Infomap / SBM;
#'   multi-resolution consensus)
#' - [compare_modules()]: Cross-species module comparison
#' - [classify_modules()]: Classify modules as conserved/species-specific
#' }
#' \subsection{Clique-level analysis}{
#' - [find_cliques()]: C++ clique detection via Bron-Kerbosch decomposition
#' - [clique_stability()]: Leave-k-out jackknife stability for trait-exclusive
#'   cliques
#' - [clique_hubs()]: Rank genes by recurrence across trait-exclusive cliques
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
#' @importFrom igraph %--%
#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom utils read.delim
## usethis namespace: end
NULL
