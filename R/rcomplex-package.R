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
#'   \item{[extract_orthologs()]}{Derive ortholog pairs from
#'     SummarizedExperiment objects}
#'   \item{[compute_network()]}{Build co-expression network (matrix or
#'     SummarizedExperiment)}
#'   \item{[compare_neighborhoods()]}{Pair-level hypergeometric tests}
#'   \item{[summarize_comparison()]}{Q-value correction and summary}
#'   \item{[comparison_to_edges()]}{Convert comparison results to clique
#'     edge format}
#'   \item{[compare_specificity()]}{Pair-level specificity score (opt-in
#'     `method = "rank"`)}
#'   \item{[null_network()]}{Shuffled-partner null network}
#'   \item{[summarize_specificity()]}{Empirical calibration and q-values
#'     for the specificity score}
#'   \item{[run_pairwise_comparisons()]}{Batch pairwise comparison pipeline}
#'   \item{[get_coexpressed_hogs()]}{Query co-expression partners of a
#'     candidate HOG}
#'   \item{[permutation_hog_test()]}{HOG-level permutation test}
#'   \item{[detect_modules()]}{Community detection with consensus}
#'   \item{[resolve_ortholog_map()]}{Paralog-resolved ortholog map}
#'   \item{[module_preservation()]}{Cross-species module topology
#'     preservation}
#'   \item{[classify_preservation()]}{Module preservation classification}
#'   \item{[module_correspondence()]}{Cross-species module matching}
#'   \item{[preservation_paired()]}{Module preservation over species pairs}
#'   \item{[identify_module_hubs()]}{Within-module hub gene identification}
#'   \item{[classify_hub_conservation()]}{Hub conservation across traits}
#'   \item{[find_cliques()]}{C++ clique detection (Bron-Kerbosch)}
#'   \item{[clique_stability()]}{Leave-k-out jackknife stability}
#'   \item{[clique_persistence()]}{Co-expressolog persistence scores}
#'   \item{[clique_threshold_sweep()]}{Threshold sweep (convenience wrapper)}
#'   \item{[classify_cliques()]}{HOG classification on the species graph
#'     (convenience wrapper): one row per HOG, trait-aware, and the input
#'     the stability and sweep rankings consume}
#'   \item{[gene_clique_graph()]}{Maximal cliques of the per-orthogroup
#'     gene graph, one row per clique member}
#'   \item{[classify_gene_cliques()]}{Six-tier gene-clique taxonomy
#'     (Rodriguez et al. 2026): copy-level rather than HOG-level, so use
#'     it when which paralog sits in the conserved core matters, or when
#'     that published taxonomy has to be reported. Neither it nor
#'     [classify_cliques()] replaces the other}
#'   \item{[characterize_hubs()]}{Regulatory potential metrics for hub genes}
#'   \item{[tag_permutation()]}{Permutation test for trait-specific module
#'     recurrence}
#' }
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
#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom utils read.delim modifyList
## usethis namespace: end
NULL
