#!/usr/bin/env Rscript

# Experimental feasibility probe: clique-resolved local module families and
# their topology versus expression-context deployment in the bundled Pooideae
# data. This is intentionally an analysis script, not a public package API.

script_start <- proc.time()[["elapsed"]]

required_packages <- c(
  "devtools", "SummarizedExperiment", "Matrix", "igraph"
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
rebuild <- "--rebuild" %in% args

config <- list(
  seed = 1L,
  density = 0.03,
  store_density = 0.06,
  alpha = 0.05,
  min_clique_species = 4L,
  max_anchors = 12L,
  target_multicopy_anchors = 6L,
  recurrent_fraction = 0.50,
  min_recurrent_species = 3L,
  min_informative_edges = 6L
)

dir.create("analysis/cache", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/output", recursive = TRUE, showWarnings = FALSE)

cache_file <- file.path(
  "analysis/cache",
  sprintf(
    "pooideae-probe-v1-density-%03d-alpha-%03d.rds",
    round(1000 * config$density),
    round(1000 * config$alpha)
  )
)

timing_state <- new.env(parent = emptyenv())
timing_state$timings <- data.frame(stage = character(), seconds = numeric())
append_timing <- function(stage, seconds) {
  timing_state$timings <- rbind(
    timing_state$timings,
    data.frame(stage = stage, seconds = round(seconds, 3))
  )
}

run_stage <- function(label, value) {
  start <- proc.time()[["elapsed"]]
  result <- force(value)
  elapsed <- proc.time()[["elapsed"]] - start
  append_timing(label, elapsed)
  result
}

write_table <- function(x, filename) {
  utils::write.csv(
    x,
    file.path("analysis/output", filename),
    row.names = FALSE,
    na = ""
  )
}

safe_rbind <- function(x) {
  x <- x[lengths(x) > 0L]
  if (length(x) == 0L) {
    return(data.frame())
  }
  do.call(rbind, x)
}

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
  stopifnot(all(conflicts$hog == 1L))

  mapped <- mapped[!duplicated(mapped$representative), , drop = FALSE]
  stats::setNames(mapped$hog, mapped$representative)
}

network_neighbors <- function(net, anchor_gene) {
  matrix <- net$network
  gene_index <- match(anchor_gene, colnames(matrix))
  if (is.na(gene_index)) {
    stop("Anchor gene is absent from its network: ", anchor_gene)
  }

  if (inherits(matrix, "dgCMatrix")) {
    first <- matrix@p[gene_index] + 1L
    last <- matrix@p[gene_index + 1L]
    if (first > last) {
      return(setNames(numeric(), character()))
    }
    positions <- seq.int(first, last)
    values <- matrix@x[positions]
    genes <- rownames(matrix)[matrix@i[positions] + 1L]
  } else {
    values <- matrix[, gene_index]
    genes <- rownames(matrix)
  }

  keep <- genes != anchor_gene & values >= net$threshold
  stats::setNames(values[keep], genes[keep])
}

local_hog_neighborhood <- function(net, anchor_gene, gene_hog) {
  neighbors <- network_neighbors(net, anchor_gene)
  if (length(neighbors) == 0L) {
    return(data.frame(
      gene = character(),
      hog = character(),
      strength = numeric()
    ))
  }

  hog <- unname(gene_hog[names(neighbors)])
  result <- data.frame(
    gene = names(neighbors),
    hog = hog,
    strength = unname(neighbors),
    stringsAsFactors = FALSE
  )
  result <- result[!is.na(result$hog), , drop = FALSE]
  result[order(result$hog, -result$strength, result$gene), , drop = FALSE]
}

selected_copy_edge_profile <- function(net, selected_genes) {
  if (length(selected_genes) < 2L) {
    return(setNames(numeric(), character()))
  }
  hogs <- names(selected_genes)
  pairs <- utils::combn(hogs, 2L)
  keys <- paste(pairs[1L, ], pairs[2L, ], sep = "\x1f")
  matrix <- as.matrix(
    net$network[selected_genes, selected_genes, drop = FALSE]
  )
  values <- matrix[cbind(
    match(pairs[1L, ], hogs),
    match(pairs[2L, ], hogs)
  )]
  values[values < net$threshold] <- 0
  stats::setNames(as.numeric(values), keys)
}

set_jaccard <- function(x, y) {
  union_size <- length(union(x, y))
  if (union_size == 0L) {
    return(NA_real_)
  }
  length(intersect(x, y)) / union_size
}

robust_scale_rows <- function(matrix) {
  centers <- apply(matrix, 1L, stats::median, na.rm = TRUE)
  scales <- apply(
    matrix,
    1L,
    stats::mad,
    constant = 1.4826,
    na.rm = TRUE
  )
  bad <- !is.finite(scales) | scales <= .Machine$double.eps
  if (any(bad)) {
    scales[bad] <- apply(
      matrix[bad, , drop = FALSE],
      1L,
      stats::sd,
      na.rm = TRUE
    )
  }
  scales[!is.finite(scales) | scales <= .Machine$double.eps] <- 1
  sweep(sweep(matrix, 1L, centers, "-"), 1L, scales, "/")
}

context_plan <- function(se_list) {
  required <- c("tissue", "time.point")
  metadata <- lapply(names(se_list), function(species) {
    col_data <- as.data.frame(SummarizedExperiment::colData(se_list[[species]]))
    if (!all(required %in% names(col_data))) {
      return(NULL)
    }
    data.frame(
      species = species,
      sample = rownames(col_data),
      tissue = as.character(col_data$tissue),
      time_point = as.character(col_data$time.point),
      day = if ("day" %in% names(col_data)) {
        as.numeric(col_data$day)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  if (any(vapply(metadata, is.null, logical(1)))) {
    return(list(
      ok = FALSE,
      reason = "tissue/time.point metadata are not available for every species"
    ))
  }

  metadata <- do.call(rbind, metadata)
  metadata$context <- paste(metadata$tissue, metadata$time_point, sep = "::")
  coverage <- aggregate(
    species ~ context,
    unique(metadata[c("species", "context")]),
    length
  )
  common <- coverage$context[coverage$species == length(se_list)]
  if (length(common) < 2L) {
    return(list(
      ok = FALSE,
      reason = "fewer than two tissue/time-point contexts match all species"
    ))
  }

  metadata <- metadata[metadata$context %in% common, , drop = FALSE]
  replicate_counts <- aggregate(
    sample ~ species + tissue + time_point + context,
    metadata,
    length
  )
  names(replicate_counts)[names(replicate_counts) == "sample"] <-
    "n_replicates"

  day_by_time <- aggregate(day ~ time_point, metadata, function(x) {
    if (any(!is.finite(x)) || diff(range(x)) > 1e-8) {
      return(NA_real_)
    }
    x[1L]
  })
  ordered_axis <- all(is.finite(day_by_time$day))
  if (ordered_axis) {
    time_levels <- day_by_time$time_point[order(day_by_time$day)]
  } else {
    time_levels <- sort(unique(metadata$time_point))
  }

  list(
    ok = TRUE,
    reason = if (ordered_axis) {
      paste(
        "Exact tissue x time-point categories match all species;",
        "day orders the shared grid, so no registration or warping is used."
      )
    } else {
      paste(
        "Tissue x time-point categories match all species, but no common",
        "numeric ordered axis is available; contexts are categorical."
      )
    },
    metadata = metadata,
    replicate_counts = replicate_counts,
    common_contexts = sort(common),
    ordered_axis = ordered_axis,
    time_levels = time_levels,
    day_by_time = day_by_time
  )
}

se_list <- run_stage(
  "load_data",
  readRDS("inst/extdata/pooideae_vignette.rds")
)
species <- names(se_list)
stopifnot(length(species) == 8L)

n_cores <- max(
  1L,
  min(4L, parallel::detectCores(logical = FALSE))
)

if (!rebuild && file.exists(cache_file)) {
  core <- run_stage("load_cache", readRDS(cache_file))
} else {
  reductions <- run_stage("reduce_orthogroups", {
    lapply(se_list, function(se) {
      expression <- SummarizedExperiment::assay(se, 1L)
      annotation <- data.frame(
        gene = rownames(se),
        hog = as.character(SummarizedExperiment::rowData(se)$hog),
        stringsAsFactors = FALSE
      )
      annotation <- annotation[!is.na(annotation$hog), , drop = FALSE]
      reduce_orthogroups(
        expression,
        annotation,
        gene_col = "gene",
        cor_threshold = 0.7
      )
    })
  })

  networks <- run_stage("compute_networks", {
    lapply(reductions, function(reduction) {
      compute_network(
        reduction$expr_matrix,
        norm_method = "MR",
        density = config$density,
        store_density = config$store_density,
        n_cores = n_cores
      )
    })
  })

  orthologs <- run_stage(
    "prepare_orthologs",
    prepare_orthologs(se_list, reductions)
  )

  edges <- run_stage("analytical_coexpressologs", {
    find_coexpressologs(
      networks,
      orthologs,
      method = "analytical",
      alpha = config$alpha,
      n_cores = n_cores,
      pi0_method = "randomized",
      pval_combine = "max",
      seed = config$seed
    )
  })

  core <- list(
    reductions = reductions,
    networks = networks,
    orthologs = orthologs,
    edges = edges
  )
  saveRDS(core, cache_file, compress = FALSE)
}

reductions <- core$reductions
networks <- core$networks
orthologs <- core$orthologs
edges <- core$edges
names(reductions) <- names(networks) <- species

gene_hog <- Map(build_gene_hog_map, se_list, reductions)
hog_genes <- lapply(species, function(sp) {
  mapping <- gene_hog[[sp]]
  mapping <- mapping[names(mapping) %in% rownames(networks[[sp]]$network)]
  split(names(mapping), unname(mapping))
})
names(gene_hog) <- names(hog_genes) <- species

candidate_cliques <- run_stage(
  "find_cliques",
  find_cliques(
    edges,
    target_species = species,
    min_species = config$min_clique_species
  )
)
if (nrow(candidate_cliques) == 0L) {
  stop(
    "No cliques span at least ",
    config$min_clique_species,
    " species at alpha = ",
    config$alpha
  )
}

candidate_cliques$gene_signature <- apply(
  candidate_cliques[species],
  1L,
  function(x) paste(ifelse(is.na(x), "", x), collapse = "|")
)
candidate_cliques$n_multicopy_species <- vapply(
  seq_len(nrow(candidate_cliques)),
  function(i) {
    represented <- species[!is.na(candidate_cliques[i, species])]
    sum(vapply(represented, function(sp) {
      length(hog_genes[[sp]][[candidate_cliques$hog[i]]]) > 1L
    }, logical(1)))
  },
  integer(1)
)
candidate_cliques$max_reduced_copies <- vapply(
  seq_len(nrow(candidate_cliques)),
  function(i) {
    represented <- species[!is.na(candidate_cliques[i, species])]
    max(vapply(represented, function(sp) {
      length(hog_genes[[sp]][[candidate_cliques$hog[i]]])
    }, integer(1)))
  },
  integer(1)
)

candidate_order <- order(
  -candidate_cliques$n_species,
  candidate_cliques$max_q,
  -candidate_cliques$mean_effect_size,
  -candidate_cliques$intensity,
  candidate_cliques$hog,
  candidate_cliques$gene_signature,
  na.last = TRUE
)
ranked_candidates <- candidate_cliques[candidate_order, , drop = FALSE]
ranked_candidates <- ranked_candidates[
  !duplicated(ranked_candidates$hog),
  ,
  drop = FALSE
]

multicopy_rows <- which(ranked_candidates$n_multicopy_species > 0L)
take_multicopy <- head(
  multicopy_rows,
  config$target_multicopy_anchors
)
remaining <- setdiff(seq_len(nrow(ranked_candidates)), take_multicopy)
take <- c(
  take_multicopy,
  head(remaining, config$max_anchors - length(take_multicopy))
)
selected <- ranked_candidates[take, , drop = FALSE]
selected <- selected[order(match(take, seq_len(nrow(ranked_candidates)))), ]
selected$anchor_id <- sprintf("anchor_%02d", seq_len(nrow(selected)))
selected$selection_stratum <- ifelse(
  selected$n_multicopy_species > 0L,
  "multi-copy-resolved",
  "strongest-remainder"
)
rownames(selected) <- NULL

local_tables <- list()
module_hog_rows <- list()
membership_rows <- list()
module_rows <- list()
module_start <- proc.time()[["elapsed"]]

for (anchor_index in seq_len(nrow(selected))) {
  anchor <- selected[anchor_index, , drop = FALSE]
  represented <- species[!is.na(anchor[1L, species])]
  anchor_hog <- as.character(anchor$hog)
  local_for_anchor <- list()

  for (sp in represented) {
    anchor_gene <- as.character(anchor[[sp]])
    stopifnot(
      anchor_gene %in% rownames(networks[[sp]]$network),
      identical(unname(gene_hog[[sp]][anchor_gene]), anchor_hog)
    )
    local <- local_hog_neighborhood(
      networks[[sp]],
      anchor_gene,
      gene_hog[[sp]]
    )
    local$anchor_id <- anchor$anchor_id
    local$species <- sp
    local$anchor_gene <- anchor_gene
    local_for_anchor[[sp]] <- local
    local_tables[[paste(anchor$anchor_id, sp)]] <- local
  }

  observed <- unique(safe_rbind(lapply(names(local_for_anchor), function(sp) {
    local <- local_for_anchor[[sp]]
    if (nrow(local) == 0L) {
      return(data.frame())
    }
    data.frame(species = sp, hog = local$hog)
  })))
  support <- if (nrow(observed) == 0L) {
    data.frame(hog = character(), n_species_neighbor = integer())
  } else {
    aggregate(species ~ hog, observed, length)
  }
  names(support)[names(support) == "species"] <- "n_species_neighbor"
  support <- support[support$hog != anchor_hog, , drop = FALSE]

  recurrence_cutoff <- max(
    config$min_recurrent_species,
    ceiling(config$recurrent_fraction * length(represented))
  )
  core_hogs <- sort(
    support$hog[support$n_species_neighbor >= recurrence_cutoff]
  )

  if (nrow(support) > 0L) {
    support$n_species_testable <- vapply(support$hog, function(hog) {
      sum(vapply(represented, function(sp) {
        length(hog_genes[[sp]][[hog]]) > 0L
      }, logical(1)))
    }, integer(1))
    support$fraction_represented <- support$n_species_neighbor /
      length(represented)
    support$is_recurrent <- support$hog %in% core_hogs
    support$role <- "neighbor"
    support$anchor_id <- anchor$anchor_id
    module_hog_rows[[anchor$anchor_id]] <- support
  }
  module_hog_rows[[paste0(anchor$anchor_id, "_anchor")]] <- data.frame(
    hog = anchor_hog,
    n_species_neighbor = length(represented),
    n_species_testable = length(represented),
    fraction_represented = 1,
    is_recurrent = TRUE,
    role = "anchor",
    anchor_id = anchor$anchor_id
  )

  for (sp in represented) {
    local <- local_for_anchor[[sp]]
    local_hogs <- unique(local$hog)
    for (hog in core_hogs) {
      genes <- hog_genes[[sp]][[hog]]
      candidates <- local[local$hog == hog, , drop = FALSE]
      selected_gene <- if (nrow(candidates) == 0L) {
        NA_character_
      } else {
        candidates$gene[1L]
      }
      anchor_edge_strength <- if (nrow(candidates) == 0L) {
        NA_real_
      } else {
        candidates$strength[1L]
      }
      status <- if (length(genes) == 0L) {
        "not_in_expression_panel"
      } else if (hog %in% local_hogs) {
        "anchor_neighbor"
      } else {
        "present_not_neighbor"
      }
      membership_rows[[length(membership_rows) + 1L]] <- data.frame(
        anchor_id = anchor$anchor_id,
        anchor_hog = anchor_hog,
        species = sp,
        anchor_gene = as.character(anchor[[sp]]),
        module_hog = hog,
        status = status,
        n_reduced_genes = length(genes),
        n_anchor_neighbor_copies = nrow(candidates),
        selected_gene = selected_gene,
        anchor_edge_strength = anchor_edge_strength
      )
    }
  }

  module_rows[[anchor$anchor_id]] <- data.frame(
    anchor_id = anchor$anchor_id,
    anchor_hog = anchor_hog,
    n_species = length(represented),
    recurrence_cutoff = recurrence_cutoff,
    n_candidate_neighbor_hogs = nrow(support),
    n_recurrent_hogs = length(core_hogs),
    median_species_support = if (length(core_hogs) == 0L) {
      NA_real_
    } else {
      stats::median(
        support$n_species_neighbor[support$hog %in% core_hogs]
      )
    }
  )
}

local_neighborhoods <- safe_rbind(local_tables)
module_hogs <- safe_rbind(module_hog_rows)
module_membership <- safe_rbind(membership_rows)
modules <- safe_rbind(module_rows)
append_timing(
  "construct_recurrent_modules",
  proc.time()[["elapsed"]] - module_start
)

stopifnot(
  nrow(modules) == nrow(selected),
  all(modules$n_recurrent_hogs >= 0L),
  all(
    module_hogs$n_species_neighbor <=
      selected$n_species[match(module_hogs$anchor_id, selected$anchor_id)]
  ),
  all(
    is.na(module_membership$selected_gene) ==
      (module_membership$status != "anchor_neighbor")
  ),
  all(
    is.na(module_membership$anchor_edge_strength) ==
      (module_membership$status != "anchor_neighbor")
  )
)
for (key in names(local_tables)) {
  local <- local_tables[[key]]
  expected <- local[!duplicated(local$hog), c("hog", "gene", "strength")]
  parts <- strsplit(key, " ", fixed = TRUE)[[1L]]
  actual <- module_membership[
    module_membership$anchor_id == parts[1L] &
      module_membership$species == parts[2L] &
      module_membership$status == "anchor_neighbor",
    c("module_hog", "selected_gene", "anchor_edge_strength")
  ]
  expected <- expected[expected$hog %in% actual$module_hog, , drop = FALSE]
  expected <- expected[order(expected$hog), , drop = FALSE]
  actual <- actual[order(actual$module_hog), , drop = FALSE]
  stopifnot(
    identical(expected$hog, actual$module_hog),
    identical(expected$gene, actual$selected_gene),
    isTRUE(all.equal(
      expected$strength,
      actual$anchor_edge_strength,
      tolerance = 0
    ))
  )
}

topology_rows <- list()
topology_edge_profiles <- list()
topology_neighbor_sets <- list()
topology_start <- proc.time()[["elapsed"]]

for (anchor_index in seq_len(nrow(selected))) {
  anchor <- selected[anchor_index, , drop = FALSE]
  represented <- species[!is.na(anchor[1L, species])]
  core_hogs <- module_hogs$hog[
    module_hogs$anchor_id == anchor$anchor_id &
      module_hogs$role == "neighbor" &
      module_hogs$is_recurrent
  ]

  for (sp in represented) {
    anchor_gene <- as.character(anchor[[sp]])
    local <- local_tables[[paste(anchor$anchor_id, sp)]]
    neighbor_hogs <- unique(local$hog)
    membership <- module_membership[
      module_membership$anchor_id == anchor$anchor_id &
        module_membership$species == sp,
      ,
      drop = FALSE
    ]
    testable_hogs <- membership$module_hog[
      membership$status != "not_in_expression_panel"
    ]
    resolved <- membership[membership$status == "anchor_neighbor", ]
    selected_genes <- stats::setNames(
      resolved$selected_gene,
      resolved$module_hog
    )
    anchor_weights <- stats::setNames(
      resolved$anchor_edge_strength,
      resolved$module_hog
    )
    edge_weights <- selected_copy_edge_profile(
      networks[[sp]],
      selected_genes
    )
    present_edges <- edge_weights[is.finite(edge_weights) & edge_weights > 0]

    if (length(selected_genes) == 0L) {
      n_components <- NA_integer_
      largest_component_fraction <- NA_real_
    } else if (length(selected_genes) == 1L) {
      n_components <- 1L
      largest_component_fraction <- 1
    } else {
      positive_keys <- names(present_edges)
      edge_frame <- if (length(positive_keys) == 0L) {
        data.frame(from = character(), to = character())
      } else {
        parts <- strsplit(positive_keys, "\x1f", fixed = TRUE)
        data.frame(
          from = vapply(parts, `[`, character(1), 1L),
          to = vapply(parts, `[`, character(1), 2L)
        )
      }
      graph <- igraph::graph_from_data_frame(
        edge_frame,
        directed = FALSE,
        vertices = data.frame(name = names(selected_genes))
      )
      components <- igraph::components(graph)
      n_components <- components$no
      largest_component_fraction <- max(components$csize) /
        length(selected_genes)
    }

    possible_edges <- choose(length(selected_genes), 2L)
    topology_rows[[length(topology_rows) + 1L]] <- data.frame(
      anchor_id = anchor$anchor_id,
      anchor_hog = anchor$hog,
      species = sp,
      anchor_gene = anchor_gene,
      local_neighbor_genes = nrow(local),
      local_neighbor_hogs = length(neighbor_hogs),
      recurrent_hogs = length(core_hogs),
      testable_recurrent_hogs = length(testable_hogs),
      resolved_recurrent_hogs = length(selected_genes),
      unresolved_present_hogs = sum(
        membership$status == "present_not_neighbor"
      ),
      recurrent_coverage = if (length(core_hogs) == 0L) {
        NA_real_
      } else {
        length(selected_genes) / length(core_hogs)
      },
      testable_recurrent_coverage = if (length(testable_hogs) == 0L) {
        NA_real_
      } else {
        length(selected_genes) / length(testable_hogs)
      },
      anchor_module_mean_strength_ratio = if (length(testable_hogs) == 0L) {
        NA_real_
      } else {
        sum(anchor_weights / networks[[sp]]$threshold) /
          length(testable_hogs)
      },
      anchor_module_present_strength_ratio = if (
        sum(anchor_weights > 0) == 0L
      ) {
        NA_real_
      } else {
        mean(
          anchor_weights[anchor_weights > 0] /
            networks[[sp]]$threshold
        )
      },
      internal_possible_edges = possible_edges,
      internal_present_edges = length(present_edges),
      internal_edge_density = if (possible_edges == 0L) {
        NA_real_
      } else {
        length(present_edges) / possible_edges
      },
      internal_mean_strength_ratio = if (length(edge_weights) == 0L) {
        NA_real_
      } else {
        mean(edge_weights / networks[[sp]]$threshold)
      },
      internal_present_strength_ratio = if (length(present_edges) == 0L) {
        NA_real_
      } else {
        mean(present_edges / networks[[sp]]$threshold)
      },
      n_components = n_components,
      largest_component_fraction = largest_component_fraction
    )
    key <- paste(anchor$anchor_id, sp)
    topology_edge_profiles[[key]] <- edge_weights
    topology_neighbor_sets[[key]] <- neighbor_hogs
  }
}

topology <- safe_rbind(topology_rows)
stopifnot(
  !anyDuplicated(topology[c("anchor_id", "species")]),
  all(
    is.na(topology$recurrent_coverage) |
      topology$recurrent_coverage >= 0 &
        topology$recurrent_coverage <= 1
  )
)

topology_pair_rows <- list()
for (anchor_index in seq_len(nrow(selected))) {
  anchor <- selected[anchor_index, , drop = FALSE]
  represented <- species[!is.na(anchor[1L, species])]
  pairs <- utils::combn(represented, 2L, simplify = FALSE)

  for (pair in pairs) {
    sp1 <- pair[1L]
    sp2 <- pair[2L]
    key1 <- paste(anchor$anchor_id, sp1)
    key2 <- paste(anchor$anchor_id, sp2)
    weights1 <- topology_edge_profiles[[key1]]
    weights2 <- topology_edge_profiles[[key2]]
    matched_keys <- intersect(names(weights1), names(weights2))
    values1 <- weights1[matched_keys]
    values2 <- weights2[matched_keys]
    informative <- values1 > 0 | values2 > 0
    shared <- values1 > 0 & values2 > 0
    n_informative <- sum(informative)
    n_shared <- sum(shared)
    edge_union <- sum(values1 > 0 | values2 > 0)

    scaled1 <- log(values1[shared] / networks[[sp1]]$threshold)
    scaled2 <- log(values2[shared] / networks[[sp2]]$threshold)
    enough <- n_shared >= config$min_informative_edges &&
      stats::sd(scaled1) > 0 &&
      stats::sd(scaled2) > 0
    topology_pair_rows[[length(topology_pair_rows) + 1L]] <- data.frame(
      anchor_id = anchor$anchor_id,
      anchor_hog = anchor$hog,
      species1 = sp1,
      species2 = sp2,
      neighborhood_jaccard = set_jaccard(
        topology_neighbor_sets[[key1]],
        topology_neighbor_sets[[key2]]
      ),
      n_matched_module_edges = length(matched_keys),
      n_informative_module_edges = n_informative,
      n_shared_present_edges = n_shared,
      module_edge_jaccard = if (edge_union == 0L) {
        NA_real_
      } else {
        n_shared / edge_union
      },
      matched_edge_correlation = if (enough) {
        stats::cor(scaled1, scaled2)
      } else {
        NA_real_
      },
      matched_edge_mean_abs_log_threshold_distance = if (n_shared == 0L) {
        NA_real_
      } else {
        mean(abs(scaled1 - scaled2))
      },
      enough_edges_for_correlation = enough
    )
  }
}
topology_pairs <- safe_rbind(topology_pair_rows)
append_timing(
  "topology_profiles",
  proc.time()[["elapsed"]] - topology_start
)

context <- context_plan(se_list)
deployment <- data.frame()
deployment_pairs <- data.frame()
background_summary <- data.frame()
background_kernel <- data.frame()
background_pcs <- data.frame()
background_eigenvalues <- data.frame()
context_diagnostic <- data.frame(
  context_matching_available = context$ok,
  ordered_axis_available = isTRUE(context$ordered_axis),
  deployment_primary = "tissue_specific",
  combined_deployment = "within_tissue_centered",
  recurrent_hog_copy_rule = paste(
    "strongest thresholded anchor edge;",
    "stable gene-ID tie break; unresolved if no anchor edge"
  ),
  anchor_selection_note = paste(
    "12 strongest candidates, including six multi-copy-resolved anchors;",
    "all selected anchors span all eight species"
  ),
  decision = context$reason
)

species_metadata <- safe_rbind(lapply(species, function(sp) {
  col_data <- as.data.frame(SummarizedExperiment::colData(se_list[[sp]]))
  stopifnot(
    length(unique(col_data$life_cycle)) == 1L,
    length(unique(col_data$pair)) == 1L
  )
  data.frame(
    species = sp,
    life_cycle = as.character(unique(col_data$life_cycle)),
    phylogenetic_pair = as.character(unique(col_data$pair))
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
stopifnot(
  all(table(life_cycle)[c("annual", "perennial")] == 4L),
  all(table(phylogenetic_pair) == 2L),
  all(vapply(
    split(species_metadata, species_metadata$phylogenetic_pair),
    \(x) setequal(x$life_cycle, c("annual", "perennial")),
    logical(1)
  ))
)

if (context$ok) {
  deployment_start <- proc.time()[["elapsed"]]
  scaled_expression <- run_stage(
    "robust_expression_scaling",
    lapply(reductions, function(x) robust_scale_rows(x$expr_matrix))
  )
  names(scaled_expression) <- species

  deployment_rows <- list()
  for (anchor_index in seq_len(nrow(selected))) {
    anchor <- selected[anchor_index, , drop = FALSE]
    represented <- species[!is.na(anchor[1L, species])]
    core_hogs <- module_hogs$hog[
      module_hogs$anchor_id == anchor$anchor_id &
        module_hogs$role == "neighbor" &
        module_hogs$is_recurrent
    ]
    if (length(core_hogs) == 0L) {
      next
    }

    for (sp in represented) {
      expression <- scaled_expression[[sp]]
      membership <- module_membership[
        module_membership$anchor_id == anchor$anchor_id &
          module_membership$species == sp,
        ,
        drop = FALSE
      ]
      resolved <- membership[membership$status == "anchor_neighbor", ]
      if (nrow(resolved) == 0L) {
        next
      }
      stopifnot(all(resolved$selected_gene %in% rownames(expression)))
      module_score <- colMeans(
        expression[resolved$selected_gene, , drop = FALSE]
      )
      anchor_gene <- as.character(anchor[[sp]])
      anchor_score <- as.numeric(expression[anchor_gene, ])

      metadata <- context$metadata[
        context$metadata$species == sp,
        ,
        drop = FALSE
      ]
      metadata <- metadata[
        match(colnames(expression), metadata$sample),
        ,
        drop = FALSE
      ]
      stopifnot(identical(metadata$sample, colnames(expression)))
      sample_scores <- data.frame(
        metadata,
        module_score = module_score,
        anchor_score = anchor_score
      )

      groups <- split(seq_len(nrow(sample_scores)), sample_scores$context)
      for (context_name in intersect(
        context$common_contexts,
        names(groups)
      )) {
        index <- groups[[context_name]]
        module_values <- sample_scores$module_score[index]
        anchor_values <- sample_scores$anchor_score[index]
        deployment_rows[[length(deployment_rows) + 1L]] <- data.frame(
          anchor_id = anchor$anchor_id,
          anchor_hog = anchor$hog,
          species = sp,
          tissue = sample_scores$tissue[index[1L]],
          time_point = sample_scores$time_point[index[1L]],
          day = sample_scores$day[index[1L]],
          context = context_name,
          n_replicates = length(index),
          n_core_hogs_present = sum(
            membership$status != "not_in_expression_panel"
          ),
          n_core_hogs_scored = nrow(resolved),
          n_core_hogs_unresolved_present = sum(
            membership$status == "present_not_neighbor"
          ),
          core_hog_coverage = nrow(resolved) / length(core_hogs),
          module_score_mean = mean(module_values),
          module_score_sd = stats::sd(module_values),
          module_score_se = stats::sd(module_values) / sqrt(length(index)),
          anchor_score_mean = mean(anchor_values),
          anchor_score_sd = stats::sd(anchor_values)
        )
      }
    }
  }
  deployment <- safe_rbind(deployment_rows)

  deployment_pair_rows <- list()
  for (anchor_index in seq_len(nrow(selected))) {
    anchor <- selected[anchor_index, , drop = FALSE]
    anchor_deployment <- deployment[
      deployment$anchor_id == anchor$anchor_id,
      ,
      drop = FALSE
    ]
    represented <- unique(anchor_deployment$species)
    if (length(represented) < 2L) {
      next
    }

    pairs <- utils::combn(represented, 2L, simplify = FALSE)
    for (pair in pairs) {
      first <- anchor_deployment[
        anchor_deployment$species == pair[1L],
        ,
        drop = FALSE
      ]
      second <- anchor_deployment[
        anchor_deployment$species == pair[2L],
        ,
        drop = FALSE
      ]
      matched_all <- merge(
        first[c("context", "tissue", "module_score_mean")],
        second[c("context", "tissue", "module_score_mean")],
        by = c("context", "tissue"),
        suffixes = c("_1", "_2")
      )

      for (scope in sort(unique(anchor_deployment$tissue))) {
        matched <- matched_all[matched_all$tissue == scope, , drop = FALSE]
        enough <- nrow(matched) >= 4L &&
          stats::sd(matched$module_score_mean_1) > 0 &&
          stats::sd(matched$module_score_mean_2) > 0
        deployment_pair_rows[[
          length(deployment_pair_rows) + 1L
        ]] <- data.frame(
          anchor_id = anchor$anchor_id,
          anchor_hog = anchor$hog,
          species1 = pair[1L],
          species2 = pair[2L],
          scope = scope,
          comparison_type = "tissue_specific",
          n_matched_contexts = nrow(matched),
          deployment_correlation = if (enough) {
            stats::cor(
              matched$module_score_mean_1,
              matched$module_score_mean_2
            )
          } else {
            NA_real_
          },
          deployment_mean_abs_distance = if (nrow(matched) == 0L) {
            NA_real_
          } else {
            mean(abs(
              matched$module_score_mean_1 -
                matched$module_score_mean_2
            ))
          }
        )
      }

      matched <- matched_all
      matched$module_score_mean_1 <- ave(
        matched$module_score_mean_1,
        matched$tissue,
        FUN = function(x) x - mean(x)
      )
      matched$module_score_mean_2 <- ave(
        matched$module_score_mean_2,
        matched$tissue,
        FUN = function(x) x - mean(x)
      )
      enough <- nrow(matched) >= 4L &&
        stats::sd(matched$module_score_mean_1) > 0 &&
        stats::sd(matched$module_score_mean_2) > 0
      deployment_pair_rows[[
        length(deployment_pair_rows) + 1L
      ]] <- data.frame(
        anchor_id = anchor$anchor_id,
        anchor_hog = anchor$hog,
        species1 = pair[1L],
        species2 = pair[2L],
        scope = "within_tissue_centered",
        comparison_type = "combined_sensitivity",
        n_matched_contexts = nrow(matched),
        deployment_correlation = if (enough) {
          stats::cor(
            matched$module_score_mean_1,
            matched$module_score_mean_2
          )
        } else {
          NA_real_
        },
        deployment_mean_abs_distance = if (nrow(matched) == 0L) {
          NA_real_
        } else {
          mean(abs(
            matched$module_score_mean_1 -
              matched$module_score_mean_2
          ))
        }
      )
    }
  }
  deployment_pairs <- safe_rbind(deployment_pair_rows)
  deployment_coverage <- unique(deployment[c(
    "anchor_id", "species", "n_core_hogs_scored"
  )])
  topology_coverage <- topology[c(
    "anchor_id", "species", "resolved_recurrent_hogs"
  )]
  coverage_check <- merge(
    deployment_coverage,
    topology_coverage,
    by = c("anchor_id", "species"),
    all = TRUE
  )
  stopifnot(
    !"all" %in% deployment_pairs$scope,
    all(
      coverage_check$n_core_hogs_scored ==
        coverage_check$resolved_recurrent_hogs
    )
  )
  append_timing(
    "deployment_profiles",
    proc.time()[["elapsed"]] - deployment_start
  )

  # Replicates are averaged within exact matched contexts. Each common HOG
  # profile is standardized over time separately within each tissue, then
  # the flattened HOG-context features are standardized across species.
  background_start <- proc.time()[["elapsed"]]
  common_background_hogs <- sort(Reduce(
    intersect,
    lapply(hog_genes, names)
  ))
  context_order <- unique(
    context$metadata[c("tissue", "time_point", "day", "context")]
  )
  context_order <- context_order[
    order(
      context_order$tissue,
      match(context_order$time_point, context$time_levels)
    ),
    ,
    drop = FALSE
  ]
  context_names <- context_order$context

  background_profiles <- array(
    NA_real_,
    dim = c(
      length(species),
      length(common_background_hogs),
      length(context_names)
    ),
    dimnames = list(species, common_background_hogs, context_names)
  )

  for (sp in species) {
    expression <- reductions[[sp]]$expr_matrix
    sample_hog <- vapply(common_background_hogs, function(hog) {
      genes <- intersect(hog_genes[[sp]][[hog]], rownames(expression))
      colMeans(expression[genes, , drop = FALSE])
    }, numeric(ncol(expression)))

    metadata <- context$metadata[
      context$metadata$species == sp,
      ,
      drop = FALSE
    ]
    metadata <- metadata[
      match(colnames(expression), metadata$sample),
      ,
      drop = FALSE
    ]
    stopifnot(identical(metadata$sample, colnames(expression)))
    group <- factor(metadata$context, levels = context_names)
    context_hog <- rowsum(sample_hog, group, reorder = FALSE)
    context_hog <- sweep(
      context_hog,
      1L,
      as.numeric(table(group)),
      "/"
    )
    hog_context <- t(context_hog)

    standardized <- matrix(
      NA_real_,
      nrow = nrow(hog_context),
      ncol = ncol(hog_context),
      dimnames = dimnames(hog_context)
    )
    valid <- rep(TRUE, nrow(hog_context))
    for (tissue in unique(context_order$tissue)) {
      columns <- which(context_order$tissue == tissue)
      profile_mean <- rowMeans(hog_context[, columns, drop = FALSE])
      profile_sd <- apply(
        hog_context[, columns, drop = FALSE],
        1L,
        stats::sd
      )
      tissue_valid <- is.finite(profile_sd) &
        profile_sd > .Machine$double.eps
      standardized[tissue_valid, columns] <- sweep(
        sweep(
          hog_context[tissue_valid, columns, drop = FALSE],
          1L,
          profile_mean[tissue_valid],
          "-"
        ),
        1L,
        profile_sd[tissue_valid],
        "/"
      )
      valid <- valid & tissue_valid
    }
    hog_context <- standardized
    hog_context[!valid, ] <- NA_real_
    background_profiles[sp, , ] <- hog_context
  }

  valid_background_hogs <- common_background_hogs[
    apply(is.finite(background_profiles), 2L, all)
  ]
  stopifnot(length(valid_background_hogs) > 2L)

  background_summary_rows <- list()
  background_kernel_rows <- list()
  background_pc_rows <- list()
  background_eigen_rows <- list()

  for (anchor_index in seq_len(nrow(selected))) {
    anchor <- selected[anchor_index, , drop = FALSE]
    recurrent_hogs <- module_hogs$hog[
      module_hogs$anchor_id == anchor$anchor_id &
        module_hogs$role == "neighbor" &
        module_hogs$is_recurrent
    ]
    excluded_hogs <- unique(c(as.character(anchor$hog), recurrent_hogs))
    retained_hogs <- setdiff(valid_background_hogs, excluded_hogs)
    profile_slice <- background_profiles[
      species,
      retained_hogs,
      context_names,
      drop = FALSE
    ]

    z <- matrix(
      NA_real_,
      nrow = length(species),
      ncol = length(retained_hogs) * length(context_names),
      dimnames = list(species, NULL)
    )
    for (i in seq_along(species)) {
      z[i, ] <- as.vector(t(profile_slice[i, , , drop = TRUE]))
    }
    feature_sd <- apply(z, 2L, stats::sd)
    variable_features <- is.finite(feature_sd) &
      feature_sd > .Machine$double.eps
    z <- scale(
      z[, variable_features, drop = FALSE],
      center = TRUE,
      scale = TRUE
    )
    kernel <- tcrossprod(z) / ncol(z)
    dimnames(kernel) <- list(species, species)

    eig <- eigen(kernel, symmetric = TRUE)
    eig$values[eig$values < 0 & abs(eig$values) < 1e-10] <- 0
    positive <- eig$values > max(eig$values) * 1e-8
    kernel_rank <- sum(positive)
    stopifnot(
      ncol(z) > 0L,
      max(abs(kernel - t(kernel))) < 1e-10,
      min(eig$values) > -1e-8,
      kernel_rank <= length(species) - 1L
    )
    positive_sum <- sum(eig$values[eig$values > 0])
    pc_scores <- sweep(
      eig$vectors[, seq_len(min(2L, ncol(eig$vectors))), drop = FALSE],
      2L,
      sqrt(pmax(eig$values[seq_len(min(2L, length(eig$values)))], 0)),
      "*"
    )

    pair_index <- utils::combn(species, 2L)
    similarity <- kernel[cbind(
      match(pair_index[1L, ], species),
      match(pair_index[2L, ], species)
    )]
    same_trait <- life_cycle[pair_index[1L, ]] ==
      life_cycle[pair_index[2L, ]]
    same_phylogenetic_pair <- phylogenetic_pair[pair_index[1L, ]] ==
      phylogenetic_pair[pair_index[2L, ]]

    background_summary_rows[[anchor$anchor_id]] <- data.frame(
      anchor_id = anchor$anchor_id,
      anchor_hog = anchor$hog,
      profile_standardization = "within_tissue",
      n_matched_contexts = length(context_names),
      n_common_profiled_hogs = length(valid_background_hogs),
      n_excluded_module_family_hogs = sum(
        valid_background_hogs %in% excluded_hogs
      ),
      n_background_hogs = length(retained_hogs),
      n_standardized_features = ncol(z),
      kernel_rank = kernel_rank,
      eigenvalue_1 = eig$values[1L],
      eigenvalue_2 = eig$values[2L],
      variance_fraction_pc1 = eig$values[1L] / positive_sum,
      variance_fraction_pc2 = eig$values[2L] / positive_sum,
      similarity_trait_correlation = stats::cor(
        similarity,
        as.numeric(same_trait)
      ),
      similarity_phylogenetic_pair_correlation = stats::cor(
        similarity,
        as.numeric(same_phylogenetic_pair)
      ),
      mean_similarity_same_trait = mean(similarity[same_trait]),
      mean_similarity_different_trait = mean(similarity[!same_trait]),
      mean_similarity_same_phylogenetic_pair = mean(
        similarity[same_phylogenetic_pair]
      ),
      mean_similarity_other_pair = mean(
        similarity[!same_phylogenetic_pair]
      )
    )

    upper <- which(upper.tri(kernel, diag = TRUE), arr.ind = TRUE)
    background_kernel_rows[[anchor$anchor_id]] <- data.frame(
      anchor_id = anchor$anchor_id,
      species1 = species[upper[, 1L]],
      species2 = species[upper[, 2L]],
      similarity = kernel[upper],
      same_trait = life_cycle[species[upper[, 1L]]] ==
        life_cycle[species[upper[, 2L]]],
      same_phylogenetic_pair =
        phylogenetic_pair[species[upper[, 1L]]] ==
        phylogenetic_pair[species[upper[, 2L]]]
    )
    background_pc_rows[[anchor$anchor_id]] <- data.frame(
      anchor_id = anchor$anchor_id,
      species = species,
      pc1 = pc_scores[, 1L],
      pc2 = pc_scores[, 2L],
      life_cycle = unname(life_cycle[species]),
      phylogenetic_pair = unname(phylogenetic_pair[species])
    )
    background_eigen_rows[[anchor$anchor_id]] <- data.frame(
      anchor_id = anchor$anchor_id,
      component = seq_along(eig$values),
      eigenvalue = eig$values,
      variance_fraction = if (positive_sum > 0) {
        pmax(eig$values, 0) / positive_sum
      } else {
        NA_real_
      }
    )
  }

  background_summary <- safe_rbind(background_summary_rows)
  background_kernel <- safe_rbind(background_kernel_rows)
  background_pcs <- safe_rbind(background_pc_rows)
  background_eigenvalues <- safe_rbind(background_eigen_rows)
  append_timing(
    "background_expression_diagnostic",
    proc.time()[["elapsed"]] - background_start
  )
}

life_history_start <- proc.time()[["elapsed"]]
balanced_annual_sets <- utils::combn(species, 4L, simplify = FALSE)
stopifnot(length(balanced_annual_sets) == 70L)

topology_metric_info <- data.frame(
  metric = c(
    "neighborhood_jaccard",
    "module_edge_jaccard",
    "matched_edge_correlation",
    "matched_edge_mean_abs_log_threshold_distance"
  ),
  metric_type = c("similarity", "similarity", "similarity", "distance"),
  stringsAsFactors = FALSE
)

topology_exact_rows <- list()
topology_phylo_rows <- list()
for (anchor_id in selected$anchor_id) {
  anchor_pairs <- topology_pairs[
    topology_pairs$anchor_id == anchor_id,
    ,
    drop = FALSE
  ]
  stopifnot(nrow(anchor_pairs) == choose(length(species), 2L))

  for (metric_index in seq_len(nrow(topology_metric_info))) {
    metric <- topology_metric_info$metric[metric_index]
    metric_type <- topology_metric_info$metric_type[metric_index]
    values <- anchor_pairs[[metric]]
    stopifnot(all(is.finite(values)))

    statistic <- function(annual_species) {
      concordant <- (anchor_pairs$species1 %in% annual_species) ==
        (anchor_pairs$species2 %in% annual_species)
      if (metric_type == "distance") {
        mean(values[!concordant]) - mean(values[concordant])
      } else {
        mean(values[concordant]) - mean(values[!concordant])
      }
    }

    observed <- statistic(names(life_cycle)[life_cycle == "annual"])
    null <- vapply(balanced_annual_sets, statistic, numeric(1))
    tolerance <- 1e-9 * max(1, abs(observed), max(abs(null)))
    n_ge <- sum(null >= observed - tolerance)
    n_tied_observed <- sum(abs(null - observed) <= tolerance)
    n_tied_max <- sum(null >= max(null) - tolerance)
    observed_concordant <- life_cycle[anchor_pairs$species1] ==
      life_cycle[anchor_pairs$species2]

    topology_exact_rows[[length(topology_exact_rows) + 1L]] <- data.frame(
      anchor_id = anchor_id,
      anchor_hog = unique(anchor_pairs$anchor_hog),
      metric = metric,
      metric_type = metric_type,
      statistic_orientation = if (metric_type == "distance") {
        "discordant_minus_concordant"
      } else {
        "concordant_minus_discordant"
      },
      concordant_mean = mean(values[observed_concordant]),
      discordant_mean = mean(values[!observed_concordant]),
      signed_statistic = observed,
      p_value = n_ge / length(null),
      p_min = 1 / length(null),
      p_attainable = n_tied_max / length(null),
      rank_from_top = n_ge,
      n_tied_observed = n_tied_observed,
      n_tied_max = n_tied_max,
      n_distinct_statistics = length(unique(signif(null, 12L))),
      n_labellings = length(null),
      global_flip_multiplicity = 2L,
      blocked_label_space = 16L,
      blocked_p_min = 1 / 16,
      blocked_p_attainable_at_best = 2 / 16,
      blocked_p_below_0.05_attainable = FALSE
    )
    stopifnot(n_tied_observed >= 2L, n_tied_max >= 2L)

    for (pair_name in unique(species_metadata$phylogenetic_pair)) {
      pair_metadata <- species_metadata[
        species_metadata$phylogenetic_pair == pair_name,
        ,
        drop = FALSE
      ]
      annual_species <- pair_metadata$species[
        pair_metadata$life_cycle == "annual"
      ]
      perennial_species <- pair_metadata$species[
        pair_metadata$life_cycle == "perennial"
      ]
      row <- (
        anchor_pairs$species1 == annual_species &
          anchor_pairs$species2 == perennial_species
      ) | (
        anchor_pairs$species1 == perennial_species &
          anchor_pairs$species2 == annual_species
      )
      stopifnot(sum(row) == 1L)
      value <- values[row]
      topology_phylo_rows[[length(topology_phylo_rows) + 1L]] <- data.frame(
        anchor_id = anchor_id,
        anchor_hog = unique(anchor_pairs$anchor_hog),
        metric = metric,
        metric_type = metric_type,
        phylogenetic_pair = pair_name,
        annual_species = annual_species,
        perennial_species = perennial_species,
        metric_value = value,
        direction_normalized_value = if (metric_type == "distance") {
          -value
        } else {
          value
        },
        blocked_label_space = 16L,
        blocked_p_attainable_at_best = 2 / 16
      )
    }
  }
}
topology_trait_exact <- safe_rbind(topology_exact_rows)
topology_trait_exact$q_bh <- stats::p.adjust(
  topology_trait_exact$p_value,
  method = "BH"
)
topology_phylo_pairs <- safe_rbind(topology_phylo_rows)
within_pair_median <- aggregate(
  metric_value ~ anchor_id + metric,
  topology_phylo_pairs,
  stats::median
)
names(within_pair_median)[3L] <- "within_phylogenetic_pair_median"
topology_trait_exact <- merge(
  topology_trait_exact,
  within_pair_median,
  by = c("anchor_id", "metric"),
  all.x = TRUE,
  sort = FALSE
)
topology_trait_exact <- topology_trait_exact[order(
  match(topology_trait_exact$anchor_id, selected$anchor_id),
  match(topology_trait_exact$metric, topology_metric_info$metric)
), ]

summarize_signed_contrasts <- function(data, group_columns) {
  key <- do.call(paste, c(data[group_columns], sep = "\x1e"))
  rows <- lapply(split(seq_len(nrow(data)), key), function(index) {
    values <- data$annual_minus_perennial[index]
    n_positive <- sum(values > 0)
    n_negative <- sum(values < 0)
    n_equal <- sum(values == 0)
    majority_direction <- if (n_positive > n_negative) {
      "annual_higher"
    } else if (n_negative > n_positive) {
      "perennial_higher"
    } else {
      "mixed_or_tied"
    }
    result <- data[index[1L], group_columns, drop = FALSE]
    result$median_annual_minus_perennial <- stats::median(values)
    result$n_annual_higher <- n_positive
    result$n_perennial_higher <- n_negative
    result$n_equal <- n_equal
    result$majority_direction <- majority_direction
    result$n_agreeing_with_majority <- max(n_positive, n_negative)
    result
  })
  result <- safe_rbind(rows)
  rownames(result) <- NULL
  result
}

topology_feature_info <- data.frame(
  feature = c(
    "recurrent_coverage",
    "anchor_module_mean_strength_ratio",
    "internal_edge_density",
    "largest_component_fraction"
  ),
  interpretation = c(
    "fraction of recurrent HOGs resolved as anchor neighbors",
    "mean threshold-normalized anchor strength over present recurrent HOGs",
    "edge density among anchor-resolved recurrent-HOG copies",
    "fraction of resolved HOGs in the largest connected component"
  ),
  stringsAsFactors = FALSE
)

topology_sister_rows <- list()
for (anchor_id in selected$anchor_id) {
  anchor_topology <- topology[topology$anchor_id == anchor_id, ]
  for (feature_index in seq_len(nrow(topology_feature_info))) {
    feature <- topology_feature_info$feature[feature_index]
    for (pair_name in unique(species_metadata$phylogenetic_pair)) {
      pair_metadata <- species_metadata[
        species_metadata$phylogenetic_pair == pair_name,
        ,
        drop = FALSE
      ]
      annual_species <- pair_metadata$species[
        pair_metadata$life_cycle == "annual"
      ]
      perennial_species <- pair_metadata$species[
        pair_metadata$life_cycle == "perennial"
      ]
      annual_value <- anchor_topology[
        anchor_topology$species == annual_species,
        feature
      ]
      perennial_value <- anchor_topology[
        anchor_topology$species == perennial_species,
        feature
      ]
      stopifnot(
        length(annual_value) == 1L,
        length(perennial_value) == 1L
      )
      topology_sister_rows[[length(topology_sister_rows) + 1L]] <-
        data.frame(
          anchor_id = anchor_id,
          anchor_hog = unique(anchor_topology$anchor_hog),
          feature = feature,
          interpretation = topology_feature_info$interpretation[
            feature_index
          ],
          phylogenetic_pair = pair_name,
          annual_species = annual_species,
          perennial_species = perennial_species,
          annual_value = annual_value,
          perennial_value = perennial_value,
          annual_minus_perennial = annual_value - perennial_value
        )
    }
  }
}
topology_sister_contrasts <- safe_rbind(topology_sister_rows)
topology_sister_summary <- summarize_signed_contrasts(
  topology_sister_contrasts,
  c("anchor_id", "anchor_hog", "feature", "interpretation")
)
stopifnot(
  all(topology_trait_exact$n_labellings == 70L),
  all(topology_trait_exact$n_tied_observed >= 2L),
  all(topology_trait_exact$p_attainable >= 2 / 70),
  all(table(
    topology_sister_contrasts$anchor_id,
    topology_sister_contrasts$feature
  ) == 4L)
)

deployment_species <- data.frame()
deployment_sister_contrasts <- data.frame()
deployment_sister_summary <- data.frame()
if (context$ok) {
  deployment_species_rows <- list()
  deployment_groups <- split(
    seq_len(nrow(deployment)),
    paste(deployment$anchor_id, deployment$species, deployment$tissue)
  )
  for (index in deployment_groups) {
    rows <- deployment[index, , drop = FALSE]
    rows <- rows[order(rows$day), , drop = FALSE]
    stopifnot(
      nrow(rows) == length(context$time_levels),
      length(unique(rows$day)) == nrow(rows)
    )
    deployment_species_rows[[length(deployment_species_rows) + 1L]] <-
      data.frame(
        anchor_id = rows$anchor_id[1L],
        anchor_hog = rows$anchor_hog[1L],
        species = rows$species[1L],
        tissue = rows$tissue[1L],
        n_time_points = nrow(rows),
        mean_activity = mean(rows$module_score_mean),
        early_to_late = rows$module_score_mean[nrow(rows)] -
          rows$module_score_mean[1L],
        early_day = rows$day[1L],
        late_day = rows$day[nrow(rows)]
      )
  }
  deployment_species <- safe_rbind(deployment_species_rows)

  deployment_feature_info <- data.frame(
    feature = c("mean_activity", "early_to_late"),
    interpretation = c(
      "equal-time-point mean module activity within tissue",
      "last minus first matched time-point module activity"
    ),
    stringsAsFactors = FALSE
  )
  deployment_sister_rows <- list()
  for (anchor_id in selected$anchor_id) {
    anchor_deployment <- deployment_species[
      deployment_species$anchor_id == anchor_id,
      ,
      drop = FALSE
    ]
    for (tissue in sort(unique(anchor_deployment$tissue))) {
      tissue_deployment <- anchor_deployment[
        anchor_deployment$tissue == tissue,
        ,
        drop = FALSE
      ]
      for (feature_index in seq_len(nrow(deployment_feature_info))) {
        feature <- deployment_feature_info$feature[feature_index]
        for (pair_name in unique(species_metadata$phylogenetic_pair)) {
          pair_metadata <- species_metadata[
            species_metadata$phylogenetic_pair == pair_name,
            ,
            drop = FALSE
          ]
          annual_species <- pair_metadata$species[
            pair_metadata$life_cycle == "annual"
          ]
          perennial_species <- pair_metadata$species[
            pair_metadata$life_cycle == "perennial"
          ]
          annual_value <- tissue_deployment[
            tissue_deployment$species == annual_species,
            feature
          ]
          perennial_value <- tissue_deployment[
            tissue_deployment$species == perennial_species,
            feature
          ]
          stopifnot(
            length(annual_value) == 1L,
            length(perennial_value) == 1L
          )
          deployment_sister_rows[[
            length(deployment_sister_rows) + 1L
          ]] <- data.frame(
            anchor_id = anchor_id,
            anchor_hog = unique(tissue_deployment$anchor_hog),
            tissue = tissue,
            feature = feature,
            interpretation = deployment_feature_info$interpretation[
              feature_index
            ],
            phylogenetic_pair = pair_name,
            annual_species = annual_species,
            perennial_species = perennial_species,
            annual_value = annual_value,
            perennial_value = perennial_value,
            annual_minus_perennial = annual_value - perennial_value
          )
        }
      }
    }
  }
  deployment_sister_contrasts <- safe_rbind(deployment_sister_rows)
  deployment_sister_summary <- summarize_signed_contrasts(
    deployment_sister_contrasts,
    c("anchor_id", "anchor_hog", "tissue", "feature", "interpretation")
  )
  stopifnot(all(table(
    deployment_sister_contrasts$anchor_id,
    deployment_sister_contrasts$tissue,
    deployment_sister_contrasts$feature
  ) == 4L))
}
append_timing(
  "life_history_descriptive_analysis",
  proc.time()[["elapsed"]] - life_history_start
)

clique_sizes <- as.data.frame(table(candidate_cliques$n_species))
names(clique_sizes) <- c("n_species", "n_candidate_cliques")
clique_sizes$n_species <- as.integer(as.character(clique_sizes$n_species))

eligible_by_anchor <- if (nrow(topology_pairs) == 0L) {
  setNames(logical(), character())
} else {
  tapply(
    topology_pairs$enough_edges_for_correlation,
    topology_pairs$anchor_id,
    any
  )
}
modules$has_pairwise_edge_correlation <- unname(
  eligible_by_anchor[modules$anchor_id]
)
modules$has_pairwise_edge_correlation[
  is.na(modules$has_pairwise_edge_correlation)
] <- FALSE

topology_anchor_summary <- aggregate(
  cbind(
    recurrent_coverage,
    internal_edge_density,
    largest_component_fraction
  ) ~ anchor_id,
  topology,
  function(x) mean(x, na.rm = TRUE)
)
representative <- merge(
  selected[c(
    "anchor_id", "hog", "n_species", "max_q", "mean_effect_size",
    "n_multicopy_species", "max_reduced_copies", "selection_stratum"
  )],
  modules,
  by = "anchor_id",
  all.x = TRUE
)
representative <- merge(
  representative,
  topology_anchor_summary,
  by = "anchor_id",
  all.x = TRUE
)
names(representative)[names(representative) == "n_species.x"] <- "n_species"
representative$n_species.y <- NULL
representative <- representative[
  order(match(representative$anchor_id, selected$anchor_id)),
  ,
  drop = FALSE
]

write_table(clique_sizes, "candidate-clique-sizes.csv")
write_table(selected, "selected-anchors.csv")
write_table(modules, "recurrent-module-summary.csv")
write_table(module_hogs, "module-hog-support.csv")
write_table(module_membership, "module-membership.csv")
write_table(topology, "topology-profiles.csv")
write_table(topology_pairs, "topology-pairwise.csv")
write_table(
  topology_trait_exact,
  "life-history-topology-exact.csv"
)
write_table(
  topology_phylo_pairs,
  "life-history-topology-phylogenetic-pairs.csv"
)
write_table(
  topology_sister_contrasts,
  "life-history-topology-sister-contrasts.csv"
)
write_table(
  topology_sister_summary,
  "life-history-topology-sister-summary.csv"
)
write_table(context_diagnostic, "context-diagnostic.csv")
if (context$ok) {
  write_table(context$replicate_counts, "context-replicate-coverage.csv")
  write_table(deployment, "deployment-profiles.csv")
  write_table(deployment_pairs, "deployment-pairwise.csv")
  write_table(
    deployment_species,
    "life-history-deployment-species.csv"
  )
  write_table(
    deployment_sister_contrasts,
    "life-history-deployment-sister-contrasts.csv"
  )
  write_table(
    deployment_sister_summary,
    "life-history-deployment-sister-summary.csv"
  )
  write_table(
    background_summary,
    "background-expression-diagnostic.csv"
  )
  write_table(background_kernel, "background-expression-kernel.csv")
  write_table(background_pcs, "background-expression-pcs.csv")
  write_table(
    background_eigenvalues,
    "background-expression-eigenvalues.csv"
  )
}
write_table(representative, "representative-anchors.csv")
append_timing("total", proc.time()[["elapsed"]] - script_start)
timings <- timing_state$timings
write_table(timings, "timings.csv")

cat("\nClique-resolved module/deployment feasibility probe\n")
cat("===================================================\n")
cat(
  "Reduced genes:",
  paste(
    sprintf(
      "%s %d->%d",
      species,
      vapply(reductions, `[[`, numeric(1), "n_original"),
      vapply(reductions, `[[`, numeric(1), "n_reduced")
    ),
    collapse = "; "
  ),
  "\n"
)
cat(
  sprintf(
    "Analytical edges: %d tested, %d conserved at q < %.2f\n",
    nrow(edges),
    sum(edges$type == "conserved", na.rm = TRUE),
    config$alpha
  )
)
cat(
  sprintf(
    "Candidate cliques: %d spanning >= %d species\n",
    nrow(candidate_cliques),
    config$min_clique_species
  )
)
print(clique_sizes, row.names = FALSE)
cat(
  sprintf(
    "Selected anchors: %d (%d resolve >=1 multi-copy species)\n",
    nrow(selected),
    sum(selected$n_multicopy_species > 0L)
  )
)
cat(
  sprintf(
    paste0(
      "Recurrent local modules: %d/%d non-empty; ",
      "median size %.1f HOGs (range %d-%d)\n"
    ),
    sum(modules$n_recurrent_hogs > 0L),
    nrow(modules),
    stats::median(modules$n_recurrent_hogs),
    min(modules$n_recurrent_hogs),
    max(modules$n_recurrent_hogs)
  )
)
cat(
  sprintf(
    paste0(
      "Modules with >=1 species pair having enough matched module edges ",
      "for a strength correlation: %d/%d (%.1f%%)\n"
    ),
    sum(modules$has_pairwise_edge_correlation),
    nrow(modules),
    100 * mean(modules$has_pairwise_edge_correlation)
  )
)
cat("Context decision:", context$reason, "\n")
if (context$ok) {
  cat(
    sprintf(
      paste0(
        "Matched contexts: %d tissue x time-point cells; ",
        "replicates per species/context %d-%d\n"
      ),
      length(context$common_contexts),
      min(context$replicate_counts$n_replicates),
      max(context$replicate_counts$n_replicates)
    )
  )
  for (scope in c(
    sort(unique(context$metadata$tissue)),
    "within_tissue_centered"
  )) {
    correlations <- deployment_pairs$deployment_correlation[
      deployment_pairs$scope == scope
    ]
    cat(
      sprintf(
        "Median deployment correlation (%s): %.3f (IQR %.3f-%.3f)\n",
        scope,
        stats::median(correlations, na.rm = TRUE),
        stats::quantile(correlations, 0.25, na.rm = TRUE),
        stats::quantile(correlations, 0.75, na.rm = TRUE)
      )
    )
  }
  cat(
    sprintf(
      paste0(
        "Background expression diagnostic: %d shared HOGs before ",
        "leave-one-module-family-out; retained %d-%d HOGs, kernel rank ",
        "%d-%d\n"
      ),
      unique(background_summary$n_common_profiled_hogs),
      min(background_summary$n_background_hogs),
      max(background_summary$n_background_hogs),
      min(background_summary$kernel_rank),
      max(background_summary$kernel_rank)
    )
  )
  cat(
    sprintf(
      paste0(
        "Median PC1/PC2 variance fractions: %.3f / %.3f; median kernel ",
        "correlation with same life cycle %.3f and phylogenetic pair %.3f\n"
      ),
      stats::median(background_summary$variance_fraction_pc1),
      stats::median(background_summary$variance_fraction_pc2),
      stats::median(
        background_summary$similarity_trait_correlation
      ),
      stats::median(
        background_summary$similarity_phylogenetic_pair_correlation
      )
    )
  )
}

multicopy_membership <- module_membership$n_reduced_genes > 1L
cat(
  sprintf(
    paste0(
      "Multi-copy recurrent-HOG memberships: %d; anchor-edge-resolved %d; ",
      "explicitly unresolved %d\n"
    ),
    sum(multicopy_membership),
    sum(
      multicopy_membership &
        module_membership$status == "anchor_neighbor"
    ),
    sum(
      multicopy_membership &
        module_membership$status == "present_not_neighbor"
    )
  )
)

cat(
  sprintf(
    paste0(
      "Selected-anchor scope: %d/%d span all eight species; this strong ",
      "complete-conservation screen is biased against gross disruption.\n"
    ),
    sum(selected$n_species == length(species)),
    nrow(selected)
  )
)
cat(
  sprintf(
    paste0(
      "Exact free topology relabelling: 70 balanced assignments, naive ",
      "p_min %.4f, global-flip attainable floor at least %.4f; ",
      "%d/%d rows have p <= 0.05 (%d after BH).\n"
    ),
    1 / 70,
    2 / 70,
    sum(topology_trait_exact$p_value <= 0.05),
    nrow(topology_trait_exact),
    sum(topology_trait_exact$q_bh <= 0.05)
  )
)
cat(
  paste0(
    "Within-pair relabelling is not tested: four sister pairs give 16 ",
    "labellings and a global-flip floor of at least 2/16 = 0.125.\n"
  )
)

exact_overview <- safe_rbind(lapply(
  topology_metric_info$metric,
  function(metric) {
    rows <- topology_trait_exact$metric == metric
    data.frame(
      metric = metric,
      median_signed_statistic = stats::median(
        topology_trait_exact$signed_statistic[rows]
      ),
      minimum_exact_p = min(topology_trait_exact$p_value[rows]),
      median_exact_rank = stats::median(
        topology_trait_exact$rank_from_top[rows]
      )
    )
  }
))
cat("\nTrait-concordant versus trait-discordant topology\n")
print(exact_overview, row.names = FALSE)

topology_sister_overview <- aggregate(
  median_annual_minus_perennial ~ feature,
  topology_sister_summary,
  stats::median
)
cat("\nMedian annual-minus-perennial sister contrasts across anchors\n")
print(topology_sister_overview, row.names = FALSE)
if (context$ok) {
  deployment_sister_overview <- aggregate(
    median_annual_minus_perennial ~ tissue + feature,
    deployment_sister_summary,
    stats::median
  )
  cat("\nDeployment sister contrasts across anchors\n")
  print(deployment_sister_overview, row.names = FALSE)
}

cat("\nRepresentative anchors\n")
print(
  utils::head(
    representative[c(
      "anchor_id", "hog", "n_species", "n_multicopy_species",
      "n_recurrent_hogs", "recurrent_coverage",
      "internal_edge_density", "has_pairwise_edge_correlation"
    )],
    6L
  ),
  row.names = FALSE
)

cat("\nTimings (seconds)\n")
print(timings, row.names = FALSE)
cat(
  "\nInterpretation guardrails:\n",
  "- Modules were defined without life-cycle labels.\n",
  "- Tissue/time-point cells were matched exactly; replicates quantify",
  " within-cell uncertainty and are not contexts.\n",
  "- Deployment correlations are tissue-specific; the only combined value",
  " is explicitly centered within tissue before correlation.\n",
  "- Recurrent-HOG copies are selected by strongest anchor edge with a",
  " stable gene-ID tie break; present non-neighbors stay unresolved.\n",
  "- Topology values are thresholded MR association strengths, not signed",
  " correlation-reversal measurements.\n",
  "- Species, not genes/edges/contexts/replicates, are the evolutionary",
  " units; all outputs are descriptive.\n",
  "- The background-expression kernel excludes each anchor module family",
  " and is a sensitivity covariate only; it does not replace phylogeny.\n",
  "- Kernel correlations use the 28 species pairs descriptively; the four",
  " named phylogenetic pairs each contain opposite life cycles.\n",
  "- Exact topology relabelling keeps each 28-pair matrix intact. Sister",
  " contrasts summarize four paired differences, not independent rows.\n",
  sep = ""
)
