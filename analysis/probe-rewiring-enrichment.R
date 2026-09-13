#!/usr/bin/env Rscript

# Experimental preranked GO enrichment of annual/perennial module rewiring.
# Positive scores are annual-leaning; negative scores are perennial-leaning.

required_packages <- c("AnnotationDbi", "fgsea", "GO.db")
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

args <- commandArgs(trailingOnly = TRUE)
annotation_arg <- grep(
  "^--annotation-file=",
  args,
  value = TRUE
)
if (length(annotation_arg) != 1L) {
  stop(
    "Supply one --annotation-file= path to the B. distachyon ",
    "v3.2 annotation_info.txt file"
  )
}
annotation_file <- sub("^--annotation-file=", "", annotation_arg)
if (!file.exists(annotation_file)) {
  stop("Annotation file does not exist: ", annotation_file)
}

screen_dir <- file.path("analysis", "cache", "rewiring-screen")
output_dir <- file.path("analysis", "output")
required_inputs <- c(
  "selected-anchors.csv",
  "life-history-topology-exact.csv",
  "life-history-topology-sister-contrasts.csv",
  "life-history-deployment-sister-contrasts.csv"
)
missing_inputs <- required_inputs[
  !file.exists(file.path(screen_dir, required_inputs))
]
if (length(missing_inputs) > 0L) {
  stop(
    "Run probe-clique-module-deployment.R --rewiring-screen first; ",
    "missing: ",
    paste(missing_inputs, collapse = ", ")
  )
}

read_screen <- function(filename) {
  utils::read.csv(
    file.path(screen_dir, filename),
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

selected <- read_screen("selected-anchors.csv")
topology_exact <- read_screen("life-history-topology-exact.csv")
topology_contrasts <- read_screen(
  "life-history-topology-sister-contrasts.csv"
)
deployment_contrasts <- read_screen(
  "life-history-deployment-sister-contrasts.csv"
)

topology_features <- c(
  "anchor_module_mean_strength_ratio",
  "internal_edge_density",
  "recurrent_coverage"
)
topology_contrasts <- topology_contrasts[
  topology_contrasts$feature %in% topology_features,
  ,
  drop = FALSE
]
topology_contrasts$ranking <- topology_contrasts$feature
topology_contrasts$tissue <- NA_character_

deployment_contrasts$ranking <- paste(
  deployment_contrasts$tissue,
  deployment_contrasts$feature,
  sep = "_"
)

rank_contrasts <- rbind(
  topology_contrasts[c(
    "anchor_hog", "ranking", "tissue", "feature",
    "phylogenetic_pair", "annual_minus_perennial"
  )],
  deployment_contrasts[c(
    "anchor_hog", "ranking", "tissue", "feature",
    "phylogenetic_pair", "annual_minus_perennial"
  )]
)

rank_groups <- split(
  seq_len(nrow(rank_contrasts)),
  rank_contrasts$ranking
)
ranking_rows <- lapply(rank_groups, function(index) {
  rows <- rank_contrasts[index, , drop = FALSE]
  values <- split(rows$annual_minus_perennial, rows$anchor_hog)
  data.frame(
    ranking = rows$ranking[1L],
    tissue = rows$tissue[1L],
    feature = rows$feature[1L],
    anchor_hog = names(values),
    score = vapply(values, stats::median, numeric(1)),
    n_annual_higher = vapply(
      values,
      function(x) sum(x > 0),
      integer(1)
    ),
    n_perennial_higher = vapply(
      values,
      function(x) sum(x < 0),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
})
rankings <- do.call(rbind, ranking_rows)
rownames(rankings) <- NULL
rankings$direction <- ifelse(
  rankings$score > 0,
  "annual",
  ifelse(rankings$score < 0, "perennial", "neutral")
)
rankings$rank_within_profile <- ave(
  -rankings$score,
  rankings$ranking,
  FUN = function(x) rank(x, ties.method = "first")
)

ranking_sizes <- table(rankings$ranking)
stopifnot(
  length(ranking_sizes) == 7L,
  all(ranking_sizes == nrow(selected)),
  !anyDuplicated(rankings[c("ranking", "anchor_hog")])
)

annotation_all <- utils::read.delim(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_annotation_columns <- c("locusName", "GO")
if (!all(required_annotation_columns %in% names(annotation_all))) {
  stop(
    "Annotation file must contain columns: ",
    paste(required_annotation_columns, collapse = ", ")
  )
}

selected$locus <- sub("[.]v3[.]2$", "", selected$BDIS)
annotation <- annotation_all[
  annotation_all$locusName %in% selected$locus &
    !is.na(annotation_all$GO) &
    nzchar(annotation_all$GO),
  c("locusName", "GO"),
  drop = FALSE
]

direct_go_rows <- lapply(seq_len(nrow(annotation)), function(index) {
  terms <- strsplit(annotation$GO[index], "\\s+")[[1L]]
  terms <- terms[grepl("^GO:[0-9]+$", terms)]
  if (length(terms) == 0L) {
    return(data.frame())
  }
  data.frame(
    locus = annotation$locusName[index],
    go_id = terms,
    stringsAsFactors = FALSE
  )
})
direct_go <- unique(do.call(rbind, direct_go_rows))
direct_go <- merge(
  selected[c("hog", "locus")],
  direct_go,
  by = "locus",
  all = FALSE,
  sort = FALSE
)

go_info <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(direct_go$go_id),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)
go_info <- go_info[
  !is.na(go_info$ONTOLOGY),
  ,
  drop = FALSE
]
direct_go <- direct_go[
  direct_go$go_id %in% go_info$GOID,
  ,
  drop = FALSE
]

ancestor_maps <- list(
  BP = as.list(GO.db::GOBPANCESTOR),
  MF = as.list(GO.db::GOMFANCESTOR),
  CC = as.list(GO.db::GOCCANCESTOR)
)
ontology_by_go <- stats::setNames(go_info$ONTOLOGY, go_info$GOID)
expanded_go_rows <- lapply(seq_len(nrow(direct_go)), function(index) {
  go_id <- direct_go$go_id[index]
  ontology <- unname(ontology_by_go[go_id])
  ancestors <- ancestor_maps[[ontology]][[go_id]]
  terms <- unique(c(go_id, ancestors))
  terms <- terms[!is.na(terms) & terms != "all"]
  data.frame(
    anchor_hog = direct_go$hog[index],
    go_id = terms,
    stringsAsFactors = FALSE
  )
})
hog_go <- unique(do.call(rbind, expanded_go_rows))

all_go_info <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = unique(hog_go$go_id),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)
all_go_info <- unique(all_go_info[
  !is.na(all_go_info$TERM) &
    !is.na(all_go_info$ONTOLOGY),
  c("GOID", "TERM", "ONTOLOGY")
])
hog_go <- hog_go[hog_go$go_id %in% all_go_info$GOID, , drop = FALSE]

pathways <- split(hog_go$anchor_hog, hog_go$go_id)
pathways <- lapply(pathways, unique)
pathway_sizes <- lengths(pathways)
pathways <- pathways[pathway_sizes >= 10L & pathway_sizes <= 100L]
if (length(pathways) == 0L) {
  stop("No GO pathways contain 10-100 ranked HOGs")
}

pair_order <- unique(rank_contrasts$phylogenetic_pair)
add_pair_consistency <- function(result, contrast_rows) {
  result$n_directional_pairs <- integer(nrow(result))
  result$pair_consistency <- numeric(nrow(result))
  result$pair_median_effects <- character(nrow(result))

  for (index in seq_len(nrow(result))) {
    members <- pathways[[result$go_id[index]]]
    rows <- contrast_rows[
      contrast_rows$anchor_hog %in% members,
      ,
      drop = FALSE
    ]
    pair_effects <- vapply(pair_order, function(pair_name) {
      values <- rows$annual_minus_perennial[
        rows$phylogenetic_pair == pair_name
      ]
      if (length(values) == 0L) {
        return(NA_real_)
      }
      stats::median(values)
    }, numeric(1))
    expected_sign <- sign(result$NES[index])
    result$n_directional_pairs[index] <- sum(
      sign(pair_effects) == expected_sign,
      na.rm = TRUE
    )
    result$pair_consistency[index] <-
      result$n_directional_pairs[index] / sum(is.finite(pair_effects))
    result$pair_median_effects[index] <- paste(
      paste0(
        names(pair_effects),
        "=",
        formatC(pair_effects, digits = 4L, format = "fg")
      ),
      collapse = ";"
    )
  }
  result
}

ranking_names <- names(rank_groups)
enrichment_rows <- lapply(seq_along(ranking_names), function(ranking_index) {
  ranking_name <- ranking_names[ranking_index]
  ranking <- rankings[
    rankings$ranking == ranking_name,
    ,
    drop = FALSE
  ]
  ranking <- ranking[
    order(-ranking$score, ranking$anchor_hog),
    ,
    drop = FALSE
  ]
  stats <- stats::setNames(ranking$score, ranking$anchor_hog)
  set.seed(1000L + ranking_index)
  result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = 10L,
    maxSize = 100L,
    eps = 0
  )
  significant <- result[result$padj < 0.05, ]
  main_pathways <- character()
  if (nrow(significant) > 0L) {
    set.seed(2000L + ranking_index)
    collapsed <- fgsea::collapsePathways(
      fgseaRes = significant,
      pathways = pathways,
      stats = stats,
      pval.threshold = 0.05
    )
    main_pathways <- collapsed$mainPathways
  }
  result$main_pathway <- result$pathway %in% main_pathways
  result <- as.data.frame(result)
  if (nrow(result) == 0L) {
    return(data.frame())
  }
  names(result)[names(result) == "pathway"] <- "go_id"
  result$leadingEdge <- vapply(
    result$leadingEdge,
    paste,
    collapse = ";",
    character(1)
  )
  result <- merge(
    result,
    all_go_info,
    by.x = "go_id",
    by.y = "GOID",
    all.x = TRUE,
    sort = FALSE
  )
  result$ranking <- ranking_name
  result$direction <- ifelse(result$NES > 0, "annual", "perennial")
  result <- add_pair_consistency(
    result,
    rank_contrasts[rank_contrasts$ranking == ranking_name, ]
  )
  result
})
enrichment <- do.call(rbind, enrichment_rows)
rownames(enrichment) <- NULL
enrichment$padj_global <- stats::p.adjust(enrichment$pval, method = "BH")
names(enrichment)[names(enrichment) == "padj"] <-
  "padj_within_ranking"
enrichment <- enrichment[
  order(
    enrichment$padj_within_ranking,
    enrichment$pval,
    -abs(enrichment$NES),
    enrichment$ranking,
    enrichment$go_id
  ),
  ,
  drop = FALSE
]

annual <- enrichment[enrichment$direction == "annual", , drop = FALSE]
perennial <- enrichment[
  enrichment$direction == "perennial",
  ,
  drop = FALSE
]
annual_main <- annual[annual$main_pathway, , drop = FALSE]
perennial_main <- perennial[
  perennial$main_pathway,
  ,
  drop = FALSE
]

main_deployment <- enrichment[
  enrichment$main_pathway &
    enrichment$padj_within_ranking < 0.05 &
    grepl("^(leaf|root)_", enrichment$ranking),
  ,
  drop = FALSE
]
leading_edge_rows <- lapply(seq_len(nrow(main_deployment)), function(index) {
  hogs <- strsplit(
    main_deployment$leadingEdge[index],
    ";",
    fixed = TRUE
  )[[1L]]
  data.frame(
    anchor_hog = hogs,
    ranking = main_deployment$ranking[index],
    go_id = main_deployment$go_id[index],
    term = main_deployment$TERM[index],
    pathway_direction = main_deployment$direction[index],
    pathway_fdr = main_deployment$padj_within_ranking[index],
    stringsAsFactors = FALSE
  )
})
leading_edges <- do.call(rbind, leading_edge_rows)
leading_edges <- merge(
  leading_edges,
  rankings[c(
    "ranking", "anchor_hog", "score",
    "n_annual_higher", "n_perennial_higher"
  )],
  by = c("ranking", "anchor_hog"),
  all.x = TRUE,
  sort = FALSE
)

topology_exact <- topology_exact[
  order(
    topology_exact$p_value,
    -topology_exact$signed_statistic,
    topology_exact$anchor_hog
  ),
  ,
  drop = FALSE
]
best_exact <- topology_exact[
  !duplicated(topology_exact$anchor_hog),
  ,
  drop = FALSE
]
topology_rankings <- rankings[
  rankings$ranking %in% topology_features,
  ,
  drop = FALSE
]
topology_rankings$pair_agreement <- pmax(
  topology_rankings$n_annual_higher,
  topology_rankings$n_perennial_higher
)
topology_rankings <- topology_rankings[
  order(
    -topology_rankings$pair_agreement,
    -abs(topology_rankings$score),
    topology_rankings$anchor_hog,
    topology_rankings$ranking
  ),
  ,
  drop = FALSE
]
best_topology <- topology_rankings[
  !duplicated(topology_rankings$anchor_hog),
  ,
  drop = FALSE
]

hotspot_groups <- split(
  seq_len(nrow(leading_edges)),
  paste(
    leading_edges$pathway_direction,
    leading_edges$anchor_hog,
    sep = "\x1f"
  )
)
hotspot_rows <- lapply(hotspot_groups, function(index) {
  rows <- leading_edges[index, , drop = FALSE]
  data.frame(
    direction = rows$pathway_direction[1L],
    anchor_hog = rows$anchor_hog[1L],
    n_main_pathways = length(unique(paste(rows$ranking, rows$go_id))),
    n_deployment_rankings = length(unique(rows$ranking)),
    best_pathway_fdr = min(rows$pathway_fdr),
    max_deployment_pair_agreement = max(pmax(
      rows$n_annual_higher,
      rows$n_perennial_higher
    )),
    pathway_terms = paste(sort(unique(rows$term)), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
hotspots <- do.call(rbind, hotspot_rows)
rownames(hotspots) <- NULL
hotspots <- merge(
  hotspots,
  best_exact[c(
    "anchor_hog", "metric", "signed_statistic", "p_value"
  )],
  by = "anchor_hog",
  all.x = TRUE,
  sort = FALSE
)
names(hotspots)[names(hotspots) == "metric"] <- "topology_metric"
names(hotspots)[names(hotspots) == "signed_statistic"] <-
  "topology_group_separation"
names(hotspots)[names(hotspots) == "p_value"] <- "topology_exact_p"
hotspots <- merge(
  hotspots,
  best_topology[c(
    "anchor_hog", "ranking", "score", "pair_agreement", "direction"
  )],
  by = "anchor_hog",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_topology")
)
names(hotspots)[names(hotspots) == "ranking"] <- "topology_feature"
names(hotspots)[names(hotspots) == "score"] <- "topology_median"
names(hotspots)[names(hotspots) == "pair_agreement"] <-
  "topology_pair_agreement"
names(hotspots)[names(hotspots) == "direction_topology"] <-
  "topology_direction"

hotspots$topology_direction <- ifelse(
  hotspots$topology_median > 0,
  "annual",
  ifelse(hotspots$topology_median < 0, "perennial", "neutral")
)
hotspots$topology_supported <-
  hotspots$topology_exact_p <= 4 / 35 &
  hotspots$topology_pair_agreement >= 3L &
  hotspots$direction == hotspots$topology_direction

species_columns <- c(
  "BDIS", "BSYL", "HVUL", "HJUB", "BMAX", "BMED", "VBRO", "FPRA"
)
annotation_columns <- c(
  "locusName", "best_arabi_gene", "best_arabi_defline",
  "best_rice_gene", "best_rice_defline"
)
annotation_summary <- unique(annotation_all[
  annotation_all$locusName %in% selected$locus,
  annotation_columns,
  drop = FALSE
])
annotation_summary <- annotation_summary[
  !duplicated(annotation_summary$locusName),
  ,
  drop = FALSE
]
hotspots <- merge(
  hotspots,
  selected[c(
    "hog", species_columns, "locus", "mean_effect_size",
    "max_q", "n_multicopy_species"
  )],
  by.x = "anchor_hog",
  by.y = "hog",
  all.x = TRUE,
  sort = FALSE
)
hotspots <- merge(
  hotspots,
  annotation_summary,
  by.x = "locus",
  by.y = "locusName",
  all.x = TRUE,
  sort = FALSE
)
hotspots <- hotspots[
  order(
    hotspots$direction,
    -hotspots$topology_supported,
    -hotspots$n_deployment_rankings,
    -hotspots$n_main_pathways,
    hotspots$best_pathway_fdr,
    -abs(hotspots$topology_median),
    hotspots$anchor_hog
  ),
  ,
  drop = FALSE
]
hotspots$hotspot_rank_within_direction <- ave(
  seq_len(nrow(hotspots)),
  hotspots$direction,
  FUN = seq_along
)

write_result(rankings, "rewiring-hog-rankings.csv")
write_result(enrichment, "rewiring-gsea-all.csv")
write_result(annual, "rewiring-gsea-annual.csv")
write_result(perennial, "rewiring-gsea-perennial.csv")
write_result(annual_main, "rewiring-gsea-annual-main.csv")
write_result(perennial_main, "rewiring-gsea-perennial-main.csv")
write_result(hotspots, "rewiring-hot-hogs.csv")

cat("Rewiring HOG enrichment\n")
cat("========================\n")
cat("Eligible HOG universe:", nrow(selected), "\n")
cat(
  "HOGs with valid direct BDIS GO annotation:",
  length(unique(direct_go$hog)),
  "\n"
)
cat("Expanded GO pathways of size 10-100:", length(pathways), "\n")
cat("Signed rankings:", length(rank_groups), "\n")
cat(
  "Within-ranking FDR < 0.05:",
  sum(enrichment$padj_within_ranking < 0.05),
  "\n"
)
cat(
  "Global FDR < 0.05:",
  sum(enrichment$padj_global < 0.05),
  "\n"
)
cat(
  "Directionally matched pathway-leading HOGs with topology support:",
  sum(hotspots$topology_supported),
  "\n"
)

show_results <- function(x, direction_label) {
  cat("\n", direction_label, "-leaning pathways\n", sep = "")
  columns <- c(
    "ranking", "go_id", "TERM", "ONTOLOGY", "NES", "pval",
    "padj_within_ranking", "size", "n_directional_pairs"
  )
  print(utils::head(x[columns], 15L), row.names = FALSE)
}

show_results(annual_main, "Annual")
show_results(perennial_main, "Perennial")

show_hotspots <- function(direction_label) {
  rows <- hotspots[
    hotspots$direction == tolower(direction_label) &
      hotspots$topology_supported,
    ,
    drop = FALSE
  ]
  cat("\n", direction_label, " cross-layer HOG hotspots\n", sep = "")
  print(
    utils::head(
      rows[c(
        "anchor_hog", "BDIS", "best_arabi_defline",
        "n_main_pathways", "n_deployment_rankings",
        "topology_feature", "topology_pair_agreement"
      )],
      10L
    ),
    row.names = FALSE
  )
}

show_hotspots("Annual")
show_hotspots("Perennial")

cat(
  "\nGuardrails:\n",
  "- Positive NES is annual-leaning; negative NES is perennial-leaning.\n",
  "- Both directions come from one signed ranking and share its null.\n",
  "- Rankings use median sister-pair contrasts, not discrete p-values.\n",
  "- Each HOG appears once; GO mapping uses its BDIS representative.\n",
  "- GO ancestors are included, and tested sets contain 10-100 HOGs.\n",
  "- Hotspots require a significant deployment leading edge plus",
  " same-direction topology in at least 3/4 sister pairs and an exact",
  " topology rank no worse than 8/70.\n",
  "- fgsea is a pathway-discovery layer. Pair consistency is descriptive",
  " and does not create independent replication.\n",
  sep = ""
)
