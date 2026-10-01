#!/usr/bin/env Rscript

# Exploratory heterochrony screen for the clique-resolved rewiring candidates.
# Timing, amplitude, and shape are kept separate. Results are descriptive
# because the candidate set and trajectories originate from the same data.

required_packages <- c("SummarizedExperiment")
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

config <- list(
  seed = 20260914L,
  n_boot = 499L,
  min_pair_agreement = 3L,
  min_bootstrap_support = 0.80,
  min_effect_percentile = 0.75
)

input_dir <- file.path("analysis", "cache", "rewiring-screen")
output_dir <- file.path("analysis", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- function(name) {
  path <- file.path(input_dir, name)
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  path
}

output_file <- function(name) file.path(output_dir, name)

write_result <- function(x, name) {
  utils::write.csv(x, output_file(name), row.names = FALSE, na = "")
}

write_cache <- function(x, name) {
  utils::write.csv(x, file.path(input_dir, name), row.names = FALSE, na = "")
}

safe_rbind <- function(x) {
  x <- x[lengths(x) > 0L]
  if (length(x) == 0L) {
    return(data.frame())
  }
  result <- do.call(rbind, x)
  rownames(result) <- NULL
  result
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

phase_centroid <- function(values, position) {
  scale <- stats::sd(values)
  if (!is.finite(scale) || scale <= .Machine$double.eps) {
    return(0.5)
  }
  shape <- (values - mean(values)) / scale
  weights <- exp(shape - max(shape))
  sum(position * weights) / sum(weights)
}

trajectory_metrics <- function(values, se, position) {
  variance <- stats::var(values)
  noise <- mean(se^2, na.rm = TRUE)
  reliability <- if (
    !is.finite(variance) || !is.finite(noise) || variance + noise == 0
  ) {
    NA_real_
  } else {
    variance / (variance + noise)
  }
  fit <- stats::lm(values ~ position)
  max_value <- max(values)
  min_value <- min(values)
  c(
    mean_activity = mean(values),
    amplitude = max_value - min_value,
    early_to_late = values[length(values)] - values[1L],
    linear_trend = unname(stats::coef(fit)[["position"]]),
    phase_centroid = phase_centroid(values, position),
    peak_position = mean(position[values == max_value]),
    trough_position = mean(position[values == min_value]),
    trajectory_reliability = reliability
  )
}

standardize_shape <- function(values) {
  scale <- stats::sd(values)
  if (!is.finite(scale) || scale <= .Machine$double.eps) {
    return(rep(0, length(values)))
  }
  (values - mean(values)) / scale
}

shape_phase_shift <- function(annual, perennial, position) {
  annual <- standardize_shape(annual)
  perennial <- standardize_shape(perennial)
  annual_fun <- stats::approxfun(position, annual, rule = 1)
  perennial_fun <- stats::approxfun(position, perennial, rule = 1)
  shifts <- seq(-0.25, 0.25, by = 0.01)
  errors <- vapply(shifts, function(shift) {
    lower <- max(0, -shift)
    upper <- min(1, 1 - shift)
    grid <- seq(lower, upper, length.out = 101L)
    mean((annual_fun(grid + shift) - perennial_fun(grid))^2)
  }, numeric(1))
  best <- which(errors == min(errors))
  best <- best[which.min(abs(shifts[best]))]
  shifts[best]
}

summarize_contrasts <- function(data) {
  groups <- split(
    seq_len(nrow(data)),
    paste(data$anchor_hog, data$tissue, data$layer, data$metric,
      sep = "\x1f"
    )
  )
  rows <- lapply(groups, function(index) {
    all_values <- data$annual_minus_perennial[index]
    values <- all_values[is.finite(all_values)]
    estimate <- if (length(values) == 0L) {
      NA_real_
    } else {
      stats::median(values)
    }
    direction <- if (!is.finite(estimate)) {
      "undefined"
    } else if (estimate > 0) {
      "annual_later_or_higher"
    } else if (estimate < 0) {
      "annual_earlier_or_lower"
    } else {
      "tied"
    }
    leave_one_out <- if (length(values) > 1L) {
      vapply(seq_along(values), function(drop) {
        stats::median(values[-drop])
      }, numeric(1))
    } else {
      NA_real_
    }
    sign_value <- sign(estimate)
    stable <- is.finite(sign_value) && sign_value != 0 &&
      all(sign(leave_one_out) == sign_value)
    row <- data[index[1L], c(
      "anchor_hog", "tissue", "layer", "metric"
    )]
    row$median_annual_minus_perennial <- estimate
    row$n_annual_later_or_higher <- sum(values > 0)
    row$n_annual_earlier_or_lower <- sum(values < 0)
    row$n_equal <- sum(values == 0)
    row$n_pairs_tested <- length(values)
    row$pair_agreement <- max(sum(values > 0), sum(values < 0))
    row$direction <- direction
    row$leave_one_pair_out_stable <- stable
    row$leave_one_pair_out_min <- if (all(is.na(leave_one_out))) {
      NA_real_
    } else {
      min(leave_one_out, na.rm = TRUE)
    }
    row$leave_one_pair_out_max <- if (all(is.na(leave_one_out))) {
      NA_real_
    } else {
      max(leave_one_out, na.rm = TRUE)
    }
    correlations <- data$shape_correlation[index]
    row$median_shape_correlation <- if (
      all(!is.finite(correlations))
    ) {
      NA_real_
    } else {
      stats::median(correlations[is.finite(correlations)])
    }
    row
  })
  safe_rbind(rows)
}

selected <- utils::read.csv(
  input_file("selected-anchors.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
deployment <- utils::read.csv(
  input_file("deployment-profiles.csv"),
  stringsAsFactors = FALSE
)
membership <- utils::read.csv(
  input_file("module-membership.csv"),
  stringsAsFactors = FALSE
)
hotspots <- utils::read.csv(
  output_file("rewiring-hot-hogs.csv"),
  stringsAsFactors = FALSE
)
regulators <- utils::read.csv(
  output_file("rewiring-regulator-cooccurrence.csv"),
  stringsAsFactors = FALSE
)

species <- c(
  "BDIS", "BSYL", "HVUL", "HJUB", "BMAX", "BMED", "VBRO", "FPRA"
)
stopifnot(
  all(species %in% names(selected)),
  all(table(
    deployment$anchor_hog,
    deployment$species,
    deployment$tissue
  ) == 5L)
)

regulators <- regulators[
  regulators$q_edge_swap_contrast < 0.05,
  ,
  drop = FALSE
]
regulator_hogs <- unique(regulators$hog)
extra_regulator_hogs <- setdiff(regulator_hogs, selected$hog)

copy_candidates <- membership[
  membership$module_hog %in% extra_regulator_hogs &
    membership$status == "anchor_neighbor" &
    !is.na(membership$selected_gene),
  ,
  drop = FALSE
]
copy_counts <- aggregate(
  anchor_id ~ module_hog + species + selected_gene,
  copy_candidates,
  length
)
names(copy_counts)[names(copy_counts) == "anchor_id"] <- "selection_count"
copy_strength <- aggregate(
  anchor_edge_strength ~ module_hog + species + selected_gene,
  copy_candidates,
  mean
)
names(copy_strength)[names(copy_strength) == "anchor_edge_strength"] <-
  "mean_anchor_edge_strength"
copy_choice <- merge(
  copy_counts,
  copy_strength,
  by = c("module_hog", "species", "selected_gene"),
  all = TRUE,
  sort = FALSE
)
copy_choice <- copy_choice[order(
  copy_choice$module_hog,
  copy_choice$species,
  -copy_choice$selection_count,
  -copy_choice$mean_anchor_edge_strength,
  copy_choice$selected_gene
), ]
copy_choice$selected <- !duplicated(
  copy_choice[c("module_hog", "species")]
)

extra_profiles <- data.frame()
if (length(extra_regulator_hogs) > 0L) {
  cache_path <- file.path(
    "analysis",
    "cache",
    "pooideae-probe-v1-density-030-alpha-050.rds"
  )
  if (!file.exists(cache_path)) {
    stop("Missing input: ", cache_path)
  }
  cache <- readRDS(cache_path)
  se_list <- readRDS(file.path("inst", "extdata", "pooideae_vignette.rds"))
  profile_rows <- list()
  chosen <- copy_choice[copy_choice$selected, , drop = FALSE]

  for (sp in species) {
    expression <- robust_scale_rows(cache$reductions[[sp]]$expr_matrix)
    metadata <- as.data.frame(
      SummarizedExperiment::colData(se_list[[sp]])
    )
    metadata$sample <- rownames(metadata)
    metadata <- metadata[
      match(colnames(expression), metadata$sample),
      ,
      drop = FALSE
    ]
    stopifnot(identical(metadata$sample, colnames(expression)))
    species_choices <- chosen[chosen$species == sp, , drop = FALSE]

    for (index in seq_len(nrow(species_choices))) {
      choice <- species_choices[index, , drop = FALSE]
      gene <- choice$selected_gene
      if (!gene %in% rownames(expression)) {
        next
      }
      values <- as.numeric(expression[gene, ])
      groups <- split(
        seq_along(values),
        paste(metadata$tissue, metadata$time.point, sep = "::")
      )
      for (group_name in names(groups)) {
        sample_index <- groups[[group_name]]
        profile_rows[[length(profile_rows) + 1L]] <- data.frame(
          anchor_id = paste0("regulator_", choice$module_hog),
          anchor_hog = choice$module_hog,
          species = sp,
          tissue = metadata$tissue[sample_index[1L]],
          time_point = metadata$time.point[sample_index[1L]],
          day = metadata$day[sample_index[1L]],
          context = group_name,
          n_replicates = length(sample_index),
          anchor_score_mean = mean(values[sample_index]),
          anchor_score_sd = stats::sd(values[sample_index]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  extra_profiles <- safe_rbind(profile_rows)
}

base_anchor <- deployment[c(
  "anchor_id", "anchor_hog", "species", "tissue", "time_point",
  "day", "context", "n_replicates", "anchor_score_mean",
  "anchor_score_sd"
)]
if (nrow(extra_profiles) > 0L) {
  base_anchor <- rbind(base_anchor, extra_profiles)
}
names(base_anchor)[names(base_anchor) == "anchor_score_mean"] <- "value"
base_anchor$se <- base_anchor$anchor_score_sd /
  sqrt(base_anchor$n_replicates)
base_anchor$layer <- "anchor"

module_profiles <- deployment[c(
  "anchor_id", "anchor_hog", "species", "tissue", "time_point",
  "day", "context", "n_replicates", "module_score_mean",
  "module_score_se"
)]
names(module_profiles)[names(module_profiles) == "module_score_mean"] <-
  "value"
names(module_profiles)[names(module_profiles) == "module_score_se"] <- "se"
module_profiles$layer <- "module"

profiles <- rbind(
  base_anchor[c(
    "anchor_id", "anchor_hog", "species", "tissue", "time_point",
    "day", "context", "n_replicates", "value", "se", "layer"
  )],
  module_profiles[c(
    "anchor_id", "anchor_hog", "species", "tissue", "time_point",
    "day", "context", "n_replicates", "value", "se", "layer"
  )]
)
profiles <- profiles[order(
  profiles$anchor_hog,
  profiles$layer,
  profiles$species,
  profiles$tissue,
  profiles$day
), ]
profiles$stage_position <- ave(
  profiles$day,
  profiles$anchor_hog,
  profiles$layer,
  profiles$species,
  profiles$tissue,
  FUN = function(x) {
    ranks <- rank(x, ties.method = "first")
    (ranks - 1) / (length(ranks) - 1)
  }
)

profile_groups <- split(
  seq_len(nrow(profiles)),
  paste(
    profiles$anchor_hog,
    profiles$species,
    profiles$tissue,
    profiles$layer,
    sep = "\x1f"
  )
)
metric_rows <- lapply(profile_groups, function(index) {
  rows <- profiles[index, , drop = FALSE]
  rows <- rows[order(rows$stage_position), , drop = FALSE]
  metrics <- trajectory_metrics(
    rows$value,
    rows$se,
    rows$stage_position
  )
  data.frame(
    anchor_hog = rows$anchor_hog[1L],
    species = rows$species[1L],
    tissue = rows$tissue[1L],
    layer = rows$layer[1L],
    n_time_points = nrow(rows),
    as.list(metrics),
    stringsAsFactors = FALSE
  )
})
trajectory <- safe_rbind(metric_rows)

eligible <- trajectory$anchor_hog %in% selected$hog
background_phase <- aggregate(
  phase_centroid ~ species + tissue + layer,
  trajectory[eligible, , drop = FALSE],
  stats::median
)
names(background_phase)[names(background_phase) == "phase_centroid"] <-
  "background_phase"
trajectory <- merge(
  trajectory,
  background_phase,
  by = c("species", "tissue", "layer"),
  all.x = TRUE,
  sort = FALSE
)
trajectory$background_centered_phase <- trajectory$phase_centroid -
  trajectory$background_phase

species_metadata <- unique(deployment[c("species")])
se_list <- if (exists("se_list", inherits = FALSE)) {
  se_list
} else {
  readRDS(file.path("inst", "extdata", "pooideae_vignette.rds"))
}
species_metadata <- safe_rbind(lapply(species, function(sp) {
  metadata <- as.data.frame(SummarizedExperiment::colData(se_list[[sp]]))
  data.frame(
    species = sp,
    life_cycle = as.character(unique(metadata$life_cycle)),
    phylogenetic_pair = as.character(unique(metadata$pair)),
    stringsAsFactors = FALSE
  )
}))

metric_names <- c(
  "mean_activity", "amplitude", "early_to_late", "linear_trend",
  "phase_centroid", "background_centered_phase", "peak_position",
  "trough_position", "trajectory_reliability"
)
contrast_rows <- list()
for (hog in unique(trajectory$anchor_hog)) {
  for (tissue in unique(trajectory$tissue)) {
    for (layer in unique(trajectory$layer)) {
      subset <- trajectory[
        trajectory$anchor_hog == hog &
          trajectory$tissue == tissue &
          trajectory$layer == layer,
        ,
        drop = FALSE
      ]
      if (nrow(subset) == 0L) {
        next
      }
      for (pair_name in unique(species_metadata$phylogenetic_pair)) {
        pair <- species_metadata[
          species_metadata$phylogenetic_pair == pair_name,
          ,
          drop = FALSE
        ]
        annual <- pair$species[pair$life_cycle == "annual"]
        perennial <- pair$species[pair$life_cycle == "perennial"]
        annual_row <- subset[subset$species == annual, , drop = FALSE]
        perennial_row <- subset[
          subset$species == perennial,
          ,
          drop = FALSE
        ]
        if (nrow(annual_row) != 1L || nrow(perennial_row) != 1L) {
          next
        }
        for (metric in metric_names) {
          contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
            anchor_hog = hog,
            tissue = tissue,
            layer = layer,
            metric = metric,
            phylogenetic_pair = pair_name,
            annual_species = annual,
            perennial_species = perennial,
            annual_value = annual_row[[metric]],
            perennial_value = perennial_row[[metric]],
            annual_minus_perennial =
              annual_row[[metric]] - perennial_row[[metric]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
contrasts <- safe_rbind(contrast_rows)

shape_rows <- list()
profile_keys <- unique(profiles[c("anchor_hog", "tissue", "layer")])
for (index in seq_len(nrow(profile_keys))) {
  key <- profile_keys[index, , drop = FALSE]
  subset <- profiles[
    profiles$anchor_hog == key$anchor_hog &
      profiles$tissue == key$tissue &
      profiles$layer == key$layer,
    ,
    drop = FALSE
  ]
  for (pair_name in unique(species_metadata$phylogenetic_pair)) {
    pair <- species_metadata[
      species_metadata$phylogenetic_pair == pair_name,
      ,
      drop = FALSE
    ]
    annual <- pair$species[pair$life_cycle == "annual"]
    perennial <- pair$species[pair$life_cycle == "perennial"]
    annual_rows <- subset[subset$species == annual, , drop = FALSE]
    perennial_rows <- subset[subset$species == perennial, , drop = FALSE]
    if (nrow(annual_rows) != 5L || nrow(perennial_rows) != 5L) {
      next
    }
    annual_rows <- annual_rows[order(annual_rows$stage_position), ]
    perennial_rows <- perennial_rows[
      order(perennial_rows$stage_position),
    ]
    correlation <- if (
      stats::sd(annual_rows$value) > 0 &&
        stats::sd(perennial_rows$value) > 0
    ) {
      stats::cor(annual_rows$value, perennial_rows$value)
    } else {
      NA_real_
    }
    shift <- shape_phase_shift(
      annual_rows$value,
      perennial_rows$value,
      annual_rows$stage_position
    )
    shape_rows[[length(shape_rows) + 1L]] <- data.frame(
      anchor_hog = key$anchor_hog,
      tissue = key$tissue,
      layer = key$layer,
      metric = "shape_phase_shift",
      phylogenetic_pair = pair_name,
      annual_species = annual,
      perennial_species = perennial,
      annual_value = NA_real_,
      perennial_value = NA_real_,
      annual_minus_perennial = shift,
      shape_correlation = correlation,
      stringsAsFactors = FALSE
    )
  }
}
shape_contrasts <- safe_rbind(shape_rows)
contrasts$shape_correlation <- NA_real_
contrasts <- rbind(contrasts, shape_contrasts)

lag <- reshape(
  trajectory[c(
    "anchor_hog", "species", "tissue", "layer", "phase_centroid"
  )],
  idvar = c("anchor_hog", "species", "tissue"),
  timevar = "layer",
  direction = "wide"
)
lag <- lag[
  is.finite(lag$phase_centroid.anchor) &
    is.finite(lag$phase_centroid.module),
  ,
  drop = FALSE
]
lag$anchor_minus_module_phase <- lag$phase_centroid.anchor -
  lag$phase_centroid.module
lag_rows <- list()
for (hog in unique(lag$anchor_hog)) {
  for (tissue in unique(lag$tissue)) {
    subset <- lag[
      lag$anchor_hog == hog & lag$tissue == tissue,
      ,
      drop = FALSE
    ]
    for (pair_name in unique(species_metadata$phylogenetic_pair)) {
      pair <- species_metadata[
        species_metadata$phylogenetic_pair == pair_name,
        ,
        drop = FALSE
      ]
      annual <- pair$species[pair$life_cycle == "annual"]
      perennial <- pair$species[pair$life_cycle == "perennial"]
      annual_value <- subset$anchor_minus_module_phase[
        subset$species == annual
      ]
      perennial_value <- subset$anchor_minus_module_phase[
        subset$species == perennial
      ]
      if (length(annual_value) != 1L || length(perennial_value) != 1L) {
        next
      }
      lag_rows[[length(lag_rows) + 1L]] <- data.frame(
        anchor_hog = hog,
        tissue = tissue,
        layer = "cross_layer",
        metric = "anchor_minus_module_phase",
        phylogenetic_pair = pair_name,
        annual_species = annual,
        perennial_species = perennial,
        annual_value = annual_value,
        perennial_value = perennial_value,
        annual_minus_perennial = annual_value - perennial_value,
        shape_correlation = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
contrasts <- rbind(contrasts, safe_rbind(lag_rows))
summary <- summarize_contrasts(contrasts)

percentile_groups <- split(
  seq_len(nrow(summary)),
  paste(summary$tissue, summary$layer, summary$metric, sep = "\x1f")
)
summary$absolute_effect_percentile <- NA_real_
for (index in percentile_groups) {
  values <- abs(summary$median_annual_minus_perennial[index])
  finite <- is.finite(values)
  summary$absolute_effect_percentile[index[finite]] <-
    rank(values[finite], ties.method = "average") / sum(finite)
}

set.seed(config$seed)
phase_groups <- unique(profiles[c("anchor_hog", "tissue", "layer")])
bootstrap_rows <- list()
for (index in seq_len(nrow(phase_groups))) {
  key <- phase_groups[index, , drop = FALSE]
  subset <- profiles[
    profiles$anchor_hog == key$anchor_hog &
      profiles$tissue == key$tissue &
      profiles$layer == key$layer,
    ,
    drop = FALSE
  ]
  pair_draws <- list()
  for (pair_name in unique(species_metadata$phylogenetic_pair)) {
    pair <- species_metadata[
      species_metadata$phylogenetic_pair == pair_name,
      ,
      drop = FALSE
    ]
    annual <- pair$species[pair$life_cycle == "annual"]
    perennial <- pair$species[pair$life_cycle == "perennial"]
    annual_rows <- subset[subset$species == annual, , drop = FALSE]
    perennial_rows <- subset[subset$species == perennial, , drop = FALSE]
    if (nrow(annual_rows) != 5L || nrow(perennial_rows) != 5L) {
      next
    }
    annual_rows <- annual_rows[order(annual_rows$stage_position), ]
    perennial_rows <- perennial_rows[
      order(perennial_rows$stage_position),
    ]
    annual_draws <- matrix(
      stats::rnorm(
        config$n_boot * 5L,
        rep(annual_rows$value, each = config$n_boot),
        rep(annual_rows$se, each = config$n_boot)
      ),
      nrow = config$n_boot
    )
    perennial_draws <- matrix(
      stats::rnorm(
        config$n_boot * 5L,
        rep(perennial_rows$value, each = config$n_boot),
        rep(perennial_rows$se, each = config$n_boot)
      ),
      nrow = config$n_boot
    )
    position <- annual_rows$stage_position
    annual_phase <- apply(
      annual_draws,
      1L,
      phase_centroid,
      position = position
    )
    perennial_phase <- apply(
      perennial_draws,
      1L,
      phase_centroid,
      position = position
    )
    pair_draws[[pair_name]] <- annual_phase - perennial_phase
  }
  if (length(pair_draws) != 4L) {
    next
  }
  draws <- do.call(cbind, pair_draws)
  median_draw <- apply(draws, 1L, stats::median)
  bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.frame(
    anchor_hog = key$anchor_hog,
    tissue = key$tissue,
    layer = key$layer,
    metric = "phase_centroid",
    bootstrap_probability_positive = mean(median_draw > 0),
    bootstrap_median_low = stats::quantile(
      median_draw,
      0.025,
      names = FALSE
    ),
    bootstrap_median_high = stats::quantile(
      median_draw,
      0.975,
      names = FALSE
    ),
    stringsAsFactors = FALSE
  )
}
bootstrap <- safe_rbind(bootstrap_rows)
summary <- merge(
  summary,
  bootstrap,
  by = c("anchor_hog", "tissue", "layer", "metric"),
  all.x = TRUE,
  sort = FALSE
)
summary$bootstrap_direction_support <- ifelse(
  summary$direction == "annual_later_or_higher",
  summary$bootstrap_probability_positive,
  ifelse(
    summary$direction == "annual_earlier_or_lower",
    1 - summary$bootstrap_probability_positive,
    NA_real_
  )
)

hotspot_meta <- unique(hotspots[
  hotspots$topology_supported,
  c(
    "anchor_hog", "direction", "topology_supported",
    "best_arabi_defline", "best_rice_defline"
  ),
  drop = FALSE
])
names(hotspot_meta)[names(hotspot_meta) == "direction"] <-
  "hotspot_direction"
regulator_meta <- unique(regulators[c(
  "hog", "direction", "regulator_classes", "best_arabi_defline",
  "best_rice_defline"
)])
names(regulator_meta) <- c(
  "anchor_hog", "regulator_direction", "regulator_classes",
  "regulator_arabi_defline", "regulator_rice_defline"
)
candidate_meta <- data.frame(
  anchor_hog = union(unique(hotspot_meta$anchor_hog), regulator_hogs),
  stringsAsFactors = FALSE
)
candidate_meta$is_eligible_anchor <-
  candidate_meta$anchor_hog %in% selected$hog
candidate_meta <- merge(
  candidate_meta,
  hotspot_meta,
  by = "anchor_hog",
  all.x = TRUE,
  sort = FALSE
)
candidate_meta <- merge(
  candidate_meta,
  regulator_meta,
  by = "anchor_hog",
  all.x = TRUE,
  sort = FALSE
)
candidate_meta$is_hotspot <- !is.na(candidate_meta$hotspot_direction)
candidate_meta$is_recurrent_regulator <-
  !is.na(candidate_meta$regulator_direction)

phase <- summary[
  summary$metric == "phase_centroid" &
    summary$anchor_hog %in% candidate_meta$anchor_hog,
  ,
  drop = FALSE
]
phase$strong_phase_evidence <-
  phase$pair_agreement >= config$min_pair_agreement &
  phase$leave_one_pair_out_stable &
  phase$bootstrap_direction_support >= config$min_bootstrap_support &
  phase$absolute_effect_percentile >= config$min_effect_percentile

candidate_rows <- list()
candidate_keys <- unique(phase[c("anchor_hog", "tissue")])
for (index in seq_len(nrow(candidate_keys))) {
  key <- candidate_keys[index, , drop = FALSE]
  rows <- phase[
    phase$anchor_hog == key$anchor_hog &
      phase$tissue == key$tissue,
    ,
    drop = FALSE
  ]
  anchor <- rows[rows$layer == "anchor", , drop = FALSE]
  module <- rows[rows$layer == "module", , drop = FALSE]
  anchor_strong <- nrow(anchor) == 1L && anchor$strong_phase_evidence
  module_strong <- nrow(module) == 1L && module$strong_phase_evidence
  pattern <- if (anchor_strong && module_strong) {
    if (anchor$direction == module$direction) {
      "coordinated_phase_shift"
    } else {
      "anchor_module_decoupling"
    }
  } else if (anchor_strong) {
    "anchor_only_phase_shift"
  } else if (module_strong) {
    "module_only_phase_shift"
  } else {
    "no_recurrent_phase_shift"
  }
  candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
    anchor_hog = key$anchor_hog,
    tissue = key$tissue,
    heterochrony_pattern = pattern,
    anchor_phase_difference = if (nrow(anchor) == 1L) {
      anchor$median_annual_minus_perennial
    } else {
      NA_real_
    },
    anchor_phase_direction = if (nrow(anchor) == 1L) {
      anchor$direction
    } else {
      NA_character_
    },
    anchor_pair_agreement = if (nrow(anchor) == 1L) {
      anchor$pair_agreement
    } else {
      NA_integer_
    },
    anchor_bootstrap_support = if (nrow(anchor) == 1L) {
      anchor$bootstrap_direction_support
    } else {
      NA_real_
    },
    anchor_effect_percentile = if (nrow(anchor) == 1L) {
      anchor$absolute_effect_percentile
    } else {
      NA_real_
    },
    module_phase_difference = if (nrow(module) == 1L) {
      module$median_annual_minus_perennial
    } else {
      NA_real_
    },
    module_phase_direction = if (nrow(module) == 1L) {
      module$direction
    } else {
      NA_character_
    },
    module_pair_agreement = if (nrow(module) == 1L) {
      module$pair_agreement
    } else {
      NA_integer_
    },
    module_bootstrap_support = if (nrow(module) == 1L) {
      module$bootstrap_direction_support
    } else {
      NA_real_
    },
    module_effect_percentile = if (nrow(module) == 1L) {
      module$absolute_effect_percentile
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}
candidates <- safe_rbind(candidate_rows)
candidates <- merge(
  candidates,
  candidate_meta,
  by = "anchor_hog",
  all.x = TRUE,
  sort = FALSE
)

add_metric_columns <- function(data, metric, layer, prefix) {
  rows <- summary[
    summary$metric == metric & summary$layer == layer,
    c(
      "anchor_hog", "tissue", "median_annual_minus_perennial",
      "pair_agreement", "direction", "leave_one_pair_out_stable",
      "absolute_effect_percentile", "median_shape_correlation"
    ),
    drop = FALSE
  ]
  names(rows)[-(1:2)] <- paste0(prefix, names(rows)[-(1:2)])
  merge(
    data,
    rows,
    by = c("anchor_hog", "tissue"),
    all.x = TRUE,
    sort = FALSE
  )
}

additional_metrics <- data.frame(
  metric = c(
    "background_centered_phase", "background_centered_phase",
    "shape_phase_shift", "shape_phase_shift",
    "mean_activity", "mean_activity",
    "amplitude", "amplitude",
    "anchor_minus_module_phase"
  ),
  layer = c(
    "anchor", "module", "anchor", "module", "anchor", "module",
    "anchor", "module", "cross_layer"
  ),
  prefix = c(
    "anchor_adjusted_phase_", "module_adjusted_phase_",
    "anchor_shape_", "module_shape_", "anchor_mean_", "module_mean_",
    "anchor_amplitude_", "module_amplitude_", "phase_lag_"
  ),
  stringsAsFactors = FALSE
)
for (index in seq_len(nrow(additional_metrics))) {
  candidates <- add_metric_columns(
    candidates,
    additional_metrics$metric[index],
    additional_metrics$layer[index],
    additional_metrics$prefix[index]
  )
}

candidates <- candidates[order(
  candidates$heterochrony_pattern == "no_recurrent_phase_shift",
  -pmax(
    candidates$anchor_effect_percentile,
    candidates$module_effect_percentile,
    na.rm = TRUE
  ),
  candidates$anchor_hog,
  candidates$tissue
), ]

copy_output <- copy_choice[
  copy_choice$module_hog %in% regulator_hogs,
  ,
  drop = FALSE
]
names(copy_output)[names(copy_output) == "module_hog"] <- "anchor_hog"

stopifnot(
  nrow(trajectory) > 0L,
  all(trajectory$n_time_points == 5L),
  all(trajectory$phase_centroid >= 0 & trajectory$phase_centroid <= 1),
  all(table(
    contrasts$anchor_hog[
      contrasts$metric == "phase_centroid" &
        contrasts$anchor_hog %in% selected$hog
    ],
    contrasts$tissue[
      contrasts$metric == "phase_centroid" &
        contrasts$anchor_hog %in% selected$hog
    ],
    contrasts$layer[
      contrasts$metric == "phase_centroid" &
        contrasts$anchor_hog %in% selected$hog
    ]
  ) == 4L),
  all(candidate_meta$anchor_hog %in% trajectory$anchor_hog)
)

write_cache(
  trajectory,
  "rewiring-heterochrony-species-metrics.csv"
)
write_cache(
  contrasts,
  "rewiring-heterochrony-sister-contrasts.csv"
)
write_cache(
  summary,
  "rewiring-heterochrony-summary.csv"
)
write_result(
  candidates,
  "rewiring-heterochrony-candidates.csv"
)
write_result(
  profiles[profiles$anchor_hog %in% candidate_meta$anchor_hog, ],
  "rewiring-heterochrony-candidate-profiles.csv"
)
write_result(
  copy_output,
  "rewiring-heterochrony-regulator-copies.csv"
)
write_result(
  data.frame(
    parameter = c(
      "n_boot", "seed", "min_pair_agreement",
      "min_bootstrap_support", "min_effect_percentile",
      "phase_interpretation", "selection_caveat"
    ),
    value = c(
      config$n_boot,
      config$seed,
      config$min_pair_agreement,
      config$min_bootstrap_support,
      config$min_effect_percentile,
      paste(
        "positive annual-minus-perennial phase means later annual",
        "deployment on the matched five-point axis"
      ),
      paste(
        "exploratory conditional analysis; hotspot selection used",
        "deployment from the same expression data"
      )
    )
  ),
  "rewiring-heterochrony-parameters.csv"
)

cat("Rewiring heterochrony screen\n")
cat("============================\n")
cat("Eligible anchor universe:", length(unique(selected$hog)), "\n")
cat("Detected candidate HOGs:", nrow(candidate_meta), "\n")
cat(
  "Candidates with a recurrent phase pattern:",
  sum(candidates$heterochrony_pattern != "no_recurrent_phase_shift"),
  "of",
  nrow(candidates),
  "candidate-tissue rows\n"
)
cat("\nCHR4\n")
print(
  candidates[
    candidates$anchor_hog == "HOG:0025149",
    c(
      "tissue", "heterochrony_pattern", "anchor_phase_difference",
      "anchor_pair_agreement", "module_phase_difference",
      "module_pair_agreement"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)
