#!/usr/bin/env Rscript

# Degree-preserving test of regulatory HOG recruitment around rewiring
# hotspot anchors. This is an exploratory conditional null: hotspot anchors
# remain fixed while their within-species network edges are rewired.

required_packages <- c(
  "AnnotationDbi", "devtools", "GO.db", "igraph", "Matrix",
  "SummarizedExperiment"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}
if (
  !file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION")[1L, "Package"]), "rcomplex")
) {
  stop("Run this script from the rcomplex repository root")
}

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  if (length(hit) != 1L) {
    stop("Supply at most one ", prefix, " argument")
  }
  sub(paste0("^", prefix), "", hit)
}

annotation_file <- arg_value("--annotation-file=")
if (is.null(annotation_file) || !file.exists(annotation_file)) {
  stop("Supply --annotation-file= with the B. distachyon annotation file")
}
n_perm <- as.integer(arg_value("--n-perm=", "499"))
n_cores <- as.integer(arg_value("--n-cores=", "1"))
swap_factor <- as.integer(arg_value("--swap-factor=", "10"))
seed <- as.integer(arg_value("--seed=", "93481"))
if (
  anyNA(c(n_perm, n_cores, swap_factor, seed)) ||
    n_perm < 1L || n_cores < 1L || swap_factor < 1L
) {
  stop("n-perm, n-cores, swap-factor, and seed must be positive integers")
}

cache_file <- file.path(
  "analysis", "cache",
  "pooideae-probe-v1-density-030-alpha-050.rds"
)
screen_dir <- file.path("analysis", "cache", "rewiring-screen")
output_dir <- file.path("analysis", "output")
input_files <- c(
  cache_file,
  file.path(screen_dir, "selected-anchors.csv"),
  file.path(output_dir, "rewiring-hot-hogs.csv")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}

read_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
write_result <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = ""
  )
}
collapse_text <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  paste(x, collapse = "; ")
}

cache <- readRDS(cache_file)
se_list <- readRDS(file.path("inst", "extdata", "pooideae_vignette.rds"))
selected <- read_csv(file.path(screen_dir, "selected-anchors.csv"))
hotspots <- read_csv(file.path(output_dir, "rewiring-hot-hogs.csv"))
hotspots <- hotspots[hotspots$topology_supported, , drop = FALSE]
if (
  nrow(hotspots) == 0L ||
    anyDuplicated(hotspots[c("anchor_hog", "direction")])
) {
  stop("Expected unique topology-supported hotspot HOGs")
}

species <- names(cache$networks)
if (!identical(species, names(se_list))) {
  stop("Network and expression species are not aligned")
}
species_metadata <- do.call(rbind, lapply(species, function(sp) {
  metadata <- as.data.frame(SummarizedExperiment::colData(se_list[[sp]]))
  if (
    length(unique(metadata$life_cycle)) != 1L ||
      length(unique(metadata$pair)) != 1L
  ) {
    stop("Species metadata are not constant for ", sp)
  }
  data.frame(
    species = sp,
    life_cycle = as.character(unique(metadata$life_cycle)),
    phylogenetic_pair = as.character(unique(metadata$pair)),
    stringsAsFactors = FALSE
  )
}))
life_cycle <- stats::setNames(
  species_metadata$life_cycle,
  species_metadata$species
)
phylogenetic_pair <- stats::setNames(
  species_metadata$phylogenetic_pair,
  species_metadata$species
)

build_gene_hog_map <- function(se, reduction) {
  annotation <- data.frame(
    original = rownames(se),
    hog = as.character(SummarizedExperiment::rowData(se)$hog),
    stringsAsFactors = FALSE
  )
  mapped <- merge(
    reduction$gene_map,
    annotation,
    by = "original",
    all.x = TRUE,
    sort = FALSE
  )
  mapped <- unique(mapped[c("representative", "hog")])
  mapped <- mapped[!is.na(mapped$hog) & nzchar(mapped$hog), , drop = FALSE]
  conflicts <- aggregate(
    hog ~ representative,
    mapped,
    function(x) length(unique(x))
  )
  if (any(conflicts$hog != 1L)) {
    stop("Reduced genes map to multiple HOGs")
  }
  mapped <- mapped[!duplicated(mapped$representative), , drop = FALSE]
  stats::setNames(mapped$hog, mapped$representative)
}

gene_hog <- Map(
  build_gene_hog_map,
  se_list,
  cache$reductions
)

annotation_all <- utils::read.delim(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation_columns <- c(
  "locusName", "GO", "Pfam", "Panther", "best_arabi_gene",
  "best_arabi_defline", "best_rice_gene", "best_rice_defline"
)
if (!all(annotation_columns %in% names(annotation_all))) {
  stop(
    "Annotation file must contain: ",
    paste(annotation_columns, collapse = ", ")
  )
}

bdis_map <- data.frame(
  bdis_gene = rownames(se_list$BDIS),
  locus = sub("[.]v3[.]2$", "", rownames(se_list$BDIS)),
  hog = as.character(SummarizedExperiment::rowData(se_list$BDIS)$hog),
  stringsAsFactors = FALSE
)
annotation <- merge(
  bdis_map,
  annotation_all[annotation_columns],
  by.x = "locus",
  by.y = "locusName",
  all.x = TRUE,
  sort = FALSE
)

direct_go_rows <- lapply(seq_len(nrow(annotation)), function(index) {
  value <- annotation$GO[index]
  if (is.na(value) || !nzchar(value)) {
    return(data.frame())
  }
  terms <- strsplit(value, "\\s+")[[1L]]
  terms <- unique(terms[grepl("^GO:[0-9]+$", terms)])
  if (length(terms) == 0L) {
    return(data.frame())
  }
  data.frame(
    hog = annotation$hog[index],
    go_id = terms,
    stringsAsFactors = FALSE
  )
})
direct_go <- unique(do.call(rbind, direct_go_rows))
go_info <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(direct_go$go_id),
  columns = "ONTOLOGY",
  keytype = "GOID"
)
go_info <- go_info[!is.na(go_info$ONTOLOGY), , drop = FALSE]
ontology_by_go <- stats::setNames(go_info$ONTOLOGY, go_info$GOID)
ancestor_maps <- list(
  BP = as.list(GO.db::GOBPANCESTOR),
  MF = as.list(GO.db::GOMFANCESTOR),
  CC = as.list(GO.db::GOCCANCESTOR)
)
expanded_go_rows <- lapply(seq_len(nrow(direct_go)), function(index) {
  go_id <- direct_go$go_id[index]
  ontology <- unname(ontology_by_go[go_id])
  if (is.na(ontology)) {
    return(data.frame())
  }
  ancestors <- ancestor_maps[[ontology]][[go_id]]
  terms <- unique(c(go_id, ancestors))
  terms <- terms[!is.na(terms) & terms != "all"]
  data.frame(
    hog = direct_go$hog[index],
    go_id = terms,
    stringsAsFactors = FALSE
  )
})
hog_go <- unique(do.call(rbind, expanded_go_rows))
go_by_hog <- split(hog_go$go_id, hog_go$hog)

text_columns <- c(
  "Pfam", "Panther", "best_arabi_defline", "best_rice_defline"
)
annotation_text <- aggregate(
  annotation[text_columns],
  list(hog = annotation$hog),
  collapse_text
)
annotation_text$combined <- apply(
  annotation_text[text_columns],
  1L,
  paste,
  collapse = "; "
)
annotation_text$combined <- tolower(annotation_text$combined)

all_hogs <- sort(unique(annotation$hog))
has_go <- function(hog, go_id) {
  terms <- go_by_hog[[hog]]
  !is.null(terms) && go_id %in% terms
}
go_flag <- function(go_id) {
  vapply(all_hogs, has_go, logical(1), go_id = go_id)
}
text_by_hog <- stats::setNames(
  annotation_text$combined,
  annotation_text$hog
)
annotation_strings <- unname(text_by_hog[all_hogs])
annotation_strings[is.na(annotation_strings)] <- ""

transcription_text <- grepl(
  paste(
    "transcription factor|transcriptional regulator|",
    "transcription regulator|transcriptional coactivator|",
    "transcriptional corepressor",
    sep = ""
  ),
  annotation_strings
)
chromatin_text <- grepl(
  paste(
    "chromatin|histone (acetyltransferase|deacetylase|methyltransferase)|",
    "nucleosome assembly|bromodomain|phd[- ]finger|snf2|swi/snf",
    sep = ""
  ),
  annotation_strings
)
rna_text <- grepl(
  "splicing factor|spliceosom|pre-mrna processing|mrna processing",
  annotation_strings
)

regulator_classes <- data.frame(
  hog = all_hogs,
  transcription_regulator = go_flag("GO:0140110") | transcription_text,
  dna_binding_tf = go_flag("GO:0003700"),
  chromatin_regulator = go_flag("GO:0006325") | chromatin_text,
  protein_kinase = go_flag("GO:0004672"),
  ubiquitin_transferase = go_flag("GO:0004842"),
  rna_processing_factor = rna_text,
  stringsAsFactors = FALSE
)
base_classes <- setdiff(names(regulator_classes), "hog")
regulator_classes$gene_expression_regulator <-
  regulator_classes$transcription_regulator |
  regulator_classes$dna_binding_tf |
  regulator_classes$chromatin_regulator
regulator_classes$any_regulator <-
  regulator_classes$gene_expression_regulator |
  regulator_classes$protein_kinase |
  regulator_classes$ubiquitin_transferase
class_names <- c(
  "any_regulator", "gene_expression_regulator", base_classes
)
class_matrix <- as.matrix(regulator_classes[class_names])
rownames(class_matrix) <- regulator_classes$hog
annotated_hogs <- regulator_classes$hog[
  nzchar(annotation_strings) |
    regulator_classes$hog %in% names(go_by_hog)
]
regulator_hogs <- regulator_classes$hog[
  regulator_classes$any_regulator
]

graphs <- lapply(cache$networks, function(net) {
  adjacency <- net$network
  adjacency@x <- as.numeric(adjacency@x >= net$threshold)
  adjacency <- Matrix::drop0(adjacency)
  igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    diag = FALSE
  )
})

candidate_rows <- merge(
  hotspots[c("anchor_hog", "direction")],
  selected,
  by.x = "anchor_hog",
  by.y = "hog",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(candidate_rows$anchor_id)) {
  stop("Some hotspot HOGs are absent from selected anchors")
}

tasks <- do.call(rbind, lapply(seq_len(nrow(candidate_rows)), function(i) {
  do.call(rbind, lapply(species, function(sp) {
    data.frame(
      direction = candidate_rows$direction[i],
      anchor_id = candidate_rows$anchor_id[i],
      anchor_hog = candidate_rows$anchor_hog[i],
      species = sp,
      anchor_gene = as.character(candidate_rows[[sp]][i]),
      life_cycle = unname(life_cycle[sp]),
      phylogenetic_pair = unname(phylogenetic_pair[sp]),
      target = unname(life_cycle[sp]) == candidate_rows$direction[i],
      stringsAsFactors = FALSE
    )
  }))
}))
if (anyNA(tasks$anchor_gene)) {
  stop("Hotspot anchors must have one selected gene in every species")
}

task_gene_index <- vapply(seq_len(nrow(tasks)), function(i) {
  match(tasks$anchor_gene[i], igraph::V(graphs[[tasks$species[i]]])$name)
}, integer(1))
if (anyNA(task_gene_index)) {
  stop("Some selected anchor genes are absent from their species graph")
}

neighbor_hogs <- function(graph, task_index) {
  sp <- tasks$species[task_index]
  gene_indices <- as.integer(
    igraph::neighbors(graph, task_gene_index[task_index])
  )
  genes <- igraph::V(graph)$name[gene_indices]
  hogs <- unique(unname(gene_hog[[sp]][genes]))
  hogs <- hogs[!is.na(hogs) & hogs != tasks$anchor_hog[task_index]]
  hogs
}

profile_graphs <- function(graph_list, retain_hogs = FALSE) {
  proportions <- matrix(
    NA_real_,
    nrow = nrow(tasks),
    ncol = length(class_names),
    dimnames = list(NULL, class_names)
  )
  regulator_occurrence <- matrix(
    0L,
    nrow = nrow(tasks),
    ncol = length(regulator_hogs),
    dimnames = list(NULL, regulator_hogs)
  )
  all_neighbors <- if (retain_hogs) vector("list", nrow(tasks)) else NULL

  for (i in seq_len(nrow(tasks))) {
    hogs <- neighbor_hogs(graph_list[[tasks$species[i]]], i)
    if (retain_hogs) {
      all_neighbors[[i]] <- hogs
    }
    annotated <- intersect(hogs, annotated_hogs)
    if (length(annotated) > 0L) {
      proportions[i, ] <- colMeans(class_matrix[annotated, , drop = FALSE])
    }
    regulators <- intersect(hogs, regulator_hogs)
    if (length(regulators) > 0L) {
      regulator_occurrence[i, regulators] <- 1L
    }
  }

  summary_values <- numeric()
  occurrence_values <- numeric()
  for (direction in c("annual", "perennial")) {
    in_direction <- tasks$direction == direction
    target <- in_direction & tasks$target
    other <- in_direction & !tasks$target
    target_mean <- colMeans(proportions[target, , drop = FALSE], na.rm = TRUE)
    other_mean <- colMeans(proportions[other, , drop = FALSE], na.rm = TRUE)
    summary_values <- c(
      summary_values,
      stats::setNames(
        target_mean,
        paste(direction, class_names, "target", sep = "|")
      ),
      stats::setNames(
        target_mean - other_mean,
        paste(direction, class_names, "contrast", sep = "|")
      )
    )
    target_counts <- colSums(
      regulator_occurrence[target, , drop = FALSE]
    )
    other_counts <- colSums(
      regulator_occurrence[other, , drop = FALSE]
    )
    occurrence_values <- c(
      occurrence_values,
      stats::setNames(
        target_counts,
        paste(direction, regulator_hogs, "target", sep = "|")
      ),
      stats::setNames(
        target_counts - other_counts,
        paste(direction, regulator_hogs, "contrast", sep = "|")
      )
    )
  }
  list(
    summary = summary_values,
    occurrence = occurrence_values,
    proportions = proportions,
    regulator_occurrence = regulator_occurrence,
    neighbors = all_neighbors
  )
}

observed <- profile_graphs(graphs, retain_hogs = TRUE)

rewire_graph <- function(graph) {
  igraph::rewire(
    graph,
    igraph::keeping_degseq(
      loops = FALSE,
      niter = swap_factor * igraph::ecount(graph)
    )
  )
}
one_permutation <- function(index) {
  set.seed(rcomplex:::.task_seed(seed, 17L, index))
  rewired <- lapply(graphs, rewire_graph)
  result <- profile_graphs(rewired)
  c(result$summary, result$occurrence)
}

use_mc <- .Platform$OS.type == "unix" && n_cores > 1L
if (use_mc) {
  null_list <- parallel::mclapply(
    seq_len(n_perm),
    one_permutation,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
} else {
  null_list <- lapply(seq_len(n_perm), one_permutation)
}
null <- do.call(rbind, null_list)
if (
  nrow(null) != n_perm ||
    anyNA(null) ||
    any(!is.finite(null))
) {
  stop("Edge-swap null produced invalid statistics")
}

summary_rows <- list()
for (direction in c("annual", "perennial")) {
  in_direction <- tasks$direction == direction
  target <- in_direction & tasks$target
  other <- in_direction & !tasks$target
  for (class_name in class_names) {
    target_name <- paste(direction, class_name, "target", sep = "|")
    contrast_name <- paste(direction, class_name, "contrast", sep = "|")
    observed_target <- unname(observed$summary[target_name])
    observed_contrast <- unname(observed$summary[contrast_name])
    null_target <- null[, target_name]
    null_contrast <- null[, contrast_name]

    pair_values <- vapply(
      unique(tasks$phylogenetic_pair),
      function(pair) {
        pair_rows <- in_direction & tasks$phylogenetic_pair == pair
        target_value <- mean(
          observed$proportions[
            pair_rows & tasks$target,
            class_name
          ],
          na.rm = TRUE
        )
        other_value <- mean(
          observed$proportions[
            pair_rows & !tasks$target,
            class_name
          ],
          na.rm = TRUE
        )
        target_value - other_value
      },
      numeric(1)
    )

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      direction = direction,
      regulator_class = class_name,
      class_hogs = sum(class_matrix[, class_name]),
      n_candidate_anchors = length(unique(tasks$anchor_hog[in_direction])),
      observed_target_fraction = observed_target,
      observed_other_fraction = observed_target - observed_contrast,
      observed_target_minus_other = observed_contrast,
      n_pairs_target_higher = sum(pair_values > 0),
      null_target_mean = mean(null_target),
      target_fold_over_null = if (mean(null_target) == 0) {
        NA_real_
      } else {
        observed_target / mean(null_target)
      },
      p_edge_swap_target = (
        1 + sum(null_target >= observed_target)
      ) / (n_perm + 1),
      null_contrast_mean = mean(null_contrast),
      p_edge_swap_contrast = (
        1 + sum(null_contrast >= observed_contrast)
      ) / (n_perm + 1),
      stringsAsFactors = FALSE
    )
  }
}
regulator_enrichment <- do.call(rbind, summary_rows)
regulator_enrichment$q_edge_swap_target <- stats::p.adjust(
  regulator_enrichment$p_edge_swap_target,
  method = "BH"
)
regulator_enrichment$q_edge_swap_contrast <- stats::p.adjust(
  regulator_enrichment$p_edge_swap_contrast,
  method = "BH"
)
regulator_enrichment <- regulator_enrichment[
  order(
    regulator_enrichment$direction,
    regulator_enrichment$p_edge_swap_contrast,
    regulator_enrichment$p_edge_swap_target,
    -regulator_enrichment$target_fold_over_null
  ),
  ,
  drop = FALSE
]

hog_annotation <- aggregate(
  annotation[c(
    "bdis_gene", "best_arabi_gene", "best_arabi_defline",
    "best_rice_gene", "best_rice_defline"
  )],
  list(hog = annotation$hog),
  collapse_text
)

cooccurrence_rows <- list()
for (direction in c("annual", "perennial")) {
  in_direction <- tasks$direction == direction
  target <- in_direction & tasks$target
  other <- in_direction & !tasks$target
  n_target_tasks <- sum(target)
  n_other_tasks <- sum(other)
  for (hog in regulator_hogs) {
    target_name <- paste(direction, hog, "target", sep = "|")
    contrast_name <- paste(direction, hog, "contrast", sep = "|")
    target_count <- unname(observed$occurrence[target_name])
    other_count <- target_count -
      unname(observed$occurrence[contrast_name])
    target_by_anchor <- tapply(
      observed$regulator_occurrence[target, hog],
      tasks$anchor_hog[target],
      sum
    )
    pair_values <- vapply(
      unique(tasks$phylogenetic_pair),
      function(pair) {
        pair_rows <- in_direction & tasks$phylogenetic_pair == pair
        sum(
          observed$regulator_occurrence[
            pair_rows & tasks$target,
            hog
          ]
        ) - sum(
          observed$regulator_occurrence[
            pair_rows & !tasks$target,
            hog
          ]
        )
      },
      numeric(1)
    )
    null_target <- null[, target_name]
    null_contrast <- null[, contrast_name]
    classes <- class_names[class_matrix[hog, class_names]]
    cooccurrence_rows[[length(cooccurrence_rows) + 1L]] <- data.frame(
      direction = direction,
      hog = hog,
      regulator_classes = paste(classes, collapse = ";"),
      target_occurrences = target_count,
      other_occurrences = other_count,
      target_prevalence = target_count / n_target_tasks,
      other_prevalence = other_count / n_other_tasks,
      target_minus_other = target_count / n_target_tasks -
        other_count / n_other_tasks,
      n_target_anchors = sum(target_by_anchor > 0),
      n_target_anchors_3_of_4 = sum(target_by_anchor >= 3L),
      n_pairs_target_higher = sum(pair_values > 0),
      null_target_mean = mean(null_target),
      fold_over_edge_swap = if (mean(null_target) == 0) {
        NA_real_
      } else {
        target_count / mean(null_target)
      },
      p_edge_swap_target = (
        1 + sum(null_target >= target_count)
      ) / (n_perm + 1),
      p_edge_swap_contrast = (
        1 + sum(null_contrast >= target_count - other_count)
      ) / (n_perm + 1),
      stringsAsFactors = FALSE
    )
  }
}
regulator_cooccurrence <- do.call(rbind, cooccurrence_rows)
regulator_cooccurrence$q_edge_swap_target <- ave(
  regulator_cooccurrence$p_edge_swap_target,
  regulator_cooccurrence$direction,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)
regulator_cooccurrence$q_edge_swap_contrast <- ave(
  regulator_cooccurrence$p_edge_swap_contrast,
  regulator_cooccurrence$direction,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)
regulator_cooccurrence <- merge(
  regulator_cooccurrence,
  hog_annotation,
  by = "hog",
  all.x = TRUE,
  sort = FALSE
)
regulator_cooccurrence <- regulator_cooccurrence[
  order(
    regulator_cooccurrence$direction,
    regulator_cooccurrence$q_edge_swap_contrast,
    regulator_cooccurrence$p_edge_swap_contrast,
    -regulator_cooccurrence$n_target_anchors_3_of_4,
    -regulator_cooccurrence$target_minus_other,
    regulator_cooccurrence$hog
  ),
  ,
  drop = FALSE
]

parameters <- data.frame(
  n_perm = n_perm,
  n_cores = n_cores,
  swap_factor = swap_factor,
  seed = seed,
  n_hotspot_anchors = nrow(candidate_rows),
  n_annotated_hogs = length(annotated_hogs),
  n_regulator_hogs = length(regulator_hogs),
  null = "degree-preserving edge swap within each species",
  selection = paste(
    "conditional on fixed hotspots selected from observed deployment",
    "and topology"
  ),
  stringsAsFactors = FALSE
)

write_result(
  regulator_enrichment,
  "rewiring-regulator-enrichment.csv"
)
write_result(
  regulator_cooccurrence,
  "rewiring-regulator-cooccurrence.csv"
)
write_result(
  parameters,
  "rewiring-regulator-null-parameters.csv"
)

cat(
  "Edge-swap permutations:", n_perm,
  "| swap factor:", swap_factor,
  "| seed:", seed, "\n"
)
cat(
  "Hotspot anchors:", nrow(candidate_rows),
  "| annotated HOGs:", length(annotated_hogs),
  "| regulator HOGs:", length(regulator_hogs), "\n"
)
cat("\nRegulator-class enrichment\n")
print(regulator_enrichment, row.names = FALSE)
cat("\nTop regulator co-occurrences by direction\n")
for (direction in c("annual", "perennial")) {
  rows <- regulator_cooccurrence[
    regulator_cooccurrence$direction == direction &
      regulator_cooccurrence$target_occurrences > 0L,
    ,
    drop = FALSE
  ]
  cat("\n", direction, "\n", sep = "")
  print(
    utils::head(
      rows[c(
        "hog", "regulator_classes", "target_occurrences",
        "other_occurrences", "n_target_anchors_3_of_4",
        "n_pairs_target_higher", "fold_over_edge_swap",
        "p_edge_swap_target", "p_edge_swap_contrast",
        "best_arabi_gene", "best_arabi_defline"
      )],
      15L
    ),
    row.names = FALSE
  )
}
cat(
  "\nGuardrail: edge-swap p-values condition on hotspot anchors selected",
  " from these observed networks; they are exploratory, not an",
  " independent confirmation.\n"
)
