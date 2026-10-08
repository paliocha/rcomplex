#' rcomplex: Comparative Co-Expression Network Analysis Across Species
#'
#' Compares gene co-expression networks across species. Maps orthologous
#' genes, builds one network per species, and tests conservation. The tests
#' work at the gene, hog, module and clique levels.
#'
#' @details
#' Start with [rcomplex()], which runs the whole analysis. The building
#' blocks below run the same steps one at a time.
#'
#' **Run**
#' * [rcomplex()]: the driver, from expression and orthologs to cliques.
#' * [write_rcomplex()]: write the result as tab-separated files.
#'
#' **Inputs**
#' * [read_orthologs()]: read an ortholog group file.
#' * [as_network()]: import a network built outside rcomplex.
#' * [as_modules()]: import modules from any gene partition.
#' * [clades_from_tree()]: derive clades from a species tree.
#'
#' **Harmonise**
#' * [reduce_orthogroups()]: merge correlated paralogs.
#' * [prepare_orthologs()]: pair orthologs for every species pair.
#' * [resolve_ortholog_map()]: pick one paralog copy per gene.
#' * [split_layers()]: split expression into wiring and deployment parts.
#'
#' **Networks**
#' * [compute_network()]: build a mutual-rank co-expression network.
#'
#' **Co-expressologs**
#' * [find_coexpressologs()]: test every ortholog pair across species pairs.
#' * [density_sweep()]: repeat the test at several density thresholds.
#' * [get_coexpressed_hogs()]: query the partners of a candidate hog.
#'
#' **Cliques**
#' * [find_cliques()]: cliques of the species graph.
#' * [gene_clique_graph()]: every maximal clique of the gene graph.
#' * [classify_cliques()], [classify_gene_cliques()]: conservation tiers.
#' * [clique_stability()], [clique_threshold_sweep()]: robustness checks.
#'
#' **Modules**
#' * [detect_modules()]: find modules with Leiden.
#' * [module_preservation()], [preservation_paired()],
#'   [classify_preservation()]: test whether modules are preserved.
#' * [module_correspondence()]: match modules across species.
#' * [identify_module_hubs()], [classify_hub_conservation()]: hub genes.
#'
#' **Traits**
#' * [preservation_matrix_test()]: relabelling test on all species pairs.
#'
#' **Nulls and diagnostics**
#' * [coexpressolog_null()]: degree-preserving edge-swap null.
#' * [null_network()]: shuffled-partner null network.
#' * [pvalue_resolution()]: how many distinct p-values a test can give.
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
