# Probe E2 harness: joint-module engines on split halves of a data set.
# DS wood: expression and gene ids as prepare_data/cache/wood_041/run.R;
#   samples of each species split at random into halves A and B (seed 1).
# DS leaf: Pooideae leaf, prepare_data/data/<SP>_se.rds, 8 species; genes
#   as load_pooideae() (hog set, variance > 0, top 20,000 by variance);
#   within each time point reps split at random into A and B (seed 1).
# One network per half (density 0.03, stored at the analysis threshold).
#
# Engines (registry ENGINES): two_step (detect_modules() per species),
# match (joint_louvain(init, freeze = TRUE)), E2 (joint_louvain()) at
# lambda in {0.25, 1, 4, 16, 64}, seeds 1-3; E1 (joint_modules(), weight
# 32 x median, resolution 6, seed 1). Hooks: `hleiden` and `bethe` source
# dev/probes/<name>.R when it exists (skipped with a message otherwise).
# A hook file defines `halves_engines`: a named list of entries
#   list(lambdas = <numeric>, seeds = <integer>,
#        run = function(ctx, lambda, seed) list(lab = <labels>, ...))
# ctx is an environment: ctx$nets, ctx$ortho, ctx$dat() (jl_data(),
# built once per half), ctx$half. `lab` is either an integer vector over
# the jl_data() nodes (1-based) or a joint_lab() vector. `secs` is filled
# in when the entry does not return it. Every entry's name is an engine.
#
# Every run is cached as <OUT>/<ds>_runs_<engine>_<lambda>_<seed>_<half>
# .rds, so a later engine only adds its own runs. Wood runs of section 11
# (wood_runs_A/B.rds) are split into that cache once, not recomputed.
# Scores (modules of 20 or more genes) read every cached run with both
# halves: per-species split-half ARI, joint split-half ARI, hog purity,
# hog pairs in one module, ARI against two-step, held-out AUROC z (half-A
# modules on half-B networks, every species pair, 100 nulls). Writes
# <ds>.tsv, <ds>_auroc.rds and <ds>_summary.tsv (seed mean, min, max).
#
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 \
#     Rscript dev/probes/joint_halves.R leaf two_step,match,E2,E1
#   ... joint_halves.R wood -          (no engines: score cached runs)
#   ... joint_halves.R leaf E2 runs    (stop after the engine runs)
#   ... joint_halves.R leaf - scores   (stop before the AUROC)
suppressPackageStartupMessages({
  library(rcomplex)
  library(Matrix)
})
args <- commandArgs(TRUE)
DS <- args[1]
stopifnot(DS %in% c("wood", "leaf"))
ENG <- if (length(args) < 2L || args[2] == "-") {
  character(0)
} else {
  strsplit(args[2], ",", fixed = TRUE)[[1]]
}
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
jl <- new.env(parent = asNamespace("rcomplex"))
jl$JL_DIR <- "dev/probes"
sys.source("dev/probes/joint_louvain.R", jl)
LAMBDA <- c(0.25, 1, 4, 16, 64)
SEEDS <- 1:3
MIN_SIZE <- 20L
N_CORES <- 8L
log_t <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
}
secs_since <- function(t0) {
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# ---- Networks per half (cached) ----
half_nets <- function(expr, split_of) {
  nets <- list()
  for (half in c("A", "B")) {
    nets[[half]] <- lapply(
      stats::setNames(names(expr), names(expr)),
      function(s) {
        x <- expr[[s]][, split_of[[s]] == half]
        t0 <- Sys.time()
        n <- compute_network(x,
          density = 0.03, sparse = TRUE, store_density = 0.03,
          n_cores = N_CORES
        )
        log_t(
          "network %s %s: %d genes, %d samples, %.0f s", half, s,
          n$n_genes, ncol(x), secs_since(t0)
        )
        n
      }
    )
  }
  nets
}

load_wood <- function() {
  E <- file.path(REPO, "prepare_data/evotree/data")
  sp <- c("Asp", "Birch", "Cher", "Nor", "Scots", "Lodge")
  col_of <- c(
    Asp = "Populus_tremula", Birch = "Betula_pendula",
    Cher = "Prunus_avium", Nor = "Picea_abies",
    Scots = "Pinus_sylvestris", Lodge = "Pinus_sylvestris"
  )
  clean_id <- function(g, s) {
    if (s == "Asp") {
      sub("\\.\\d\\.p\\d$", "", g)
    } else if (s %in% c("Nor", "Cher")) {
      sub("\\.p\\d$", "", g)
    } else {
      g
    }
  }
  read_expr <- function(s) {
    f <- file.path(
      E, paste0(s, "Wood"), paste0(s, "Wood_transcriptomics.txt.gz")
    )
    d <- data.table::fread(f, showProgress = FALSE)
    x <- as.matrix(d[, -1])
    rownames(x) <- d[[1]]
    stopifnot(!anyDuplicated(rownames(x)))
    x
  }
  expr <- stats::setNames(lapply(sp, read_expr), sp)
  og <- data.table::fread(
    file.path(E, "DATA/OrthoGroups", "Orthogroups_20240823_clean_N1.tsv.gz"),
    showProgress = FALSE, data.table = FALSE
  )
  ortho <- do.call(rbind, lapply(sp, function(s) {
    v <- og[[col_of[[s]]]]
    ok <- !is.na(v) & v != ""
    g <- strsplit(v[ok], ", ", fixed = TRUE)
    d <- data.frame(
      species = s, gene = clean_id(unlist(g), s),
      hog = rep(og$OrthoGroup[ok], lengths(g))
    )
    unique(d[d$gene %in% rownames(expr[[s]]), ])
  }))
  stopifnot(!anyDuplicated(ortho[, c("species", "gene")]))
  set.seed(1)
  split_of <- lapply(expr, function(x) {
    sample(rep(c("A", "B"), length.out = ncol(x)))
  })
  list(nets = half_nets(expr, split_of), ortho = ortho, split = split_of)
}

load_leaf <- function() {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  sp <- c("BDIS", "BMAX", "BMED", "BSYL", "FPRA", "HJUB", "HVUL", "VBRO")
  set.seed(1)
  expr <- list()
  split_of <- list()
  ortho <- list()
  for (s in sp) {
    se <- readRDS(file.path(REPO, "prepare_data/data", paste0(s, "_se.rds")))
    se <- se[, colData(se)$tissue == "leaf"]
    hog <- as.character(rowData(se)$hog)
    keep <- !is.na(hog) & nzchar(hog)
    x <- as.matrix(assay(se, "vst"))[keep, , drop = FALSE]
    hog <- hog[keep]
    v <- apply(x, 1L, stats::var)
    ok <- is.finite(v) & v > 0
    x <- x[ok, , drop = FALSE]
    hog <- hog[ok]
    v <- v[ok]
    if (nrow(x) > 20000L) {
      top <- order(v, decreasing = TRUE)[seq_len(20000L)]
      x <- x[top, , drop = FALSE]
      hog <- hog[top]
    }
    stopifnot(!anyDuplicated(rownames(x)))
    # within each time point, reps alternate A, B at random
    tp <- as.character(colData(se)$time_point)
    h <- character(length(tp))
    for (t in unique(tp)) {
      i <- which(tp == t)
      h[i] <- sample(rep(c("A", "B"), length.out = length(i)))
    }
    expr[[s]] <- x
    split_of[[s]] <- h
    ortho[[s]] <- data.frame(species = s, gene = rownames(x), hog = hog)
  }
  ortho <- do.call(rbind, ortho)
  rownames(ortho) <- NULL
  list(nets = half_nets(expr, split_of), ortho = ortho, split = split_of)
}

HALVES <- file.path(OUT, if (DS == "wood") "halves.rds" else "leaf_halves.rds")
if (!file.exists(HALVES)) {
  saveRDS(if (DS == "wood") load_wood() else load_leaf(), HALVES)
  invisible(gc())
}
hv <- readRDS(HALVES)
SP <- names(hv$nets$A)
ortho <- hv$ortho

joint_lab <- function(mods) {
  unlist(lapply(names(mods), function(s) {
    m <- mods[[s]]
    stats::setNames(as.character(m), paste(s, names(m), sep = "\t"))
  }))
}
sp_of <- function(lab) sub("\t.*", "", names(lab))
# joint_lab() vector back to a per-species list of label vectors
per_species <- function(lab) {
  lapply(split(lab, factor(sp_of(lab), levels = SP)), function(l) {
    stats::setNames(l, sub("^[^\t]*\t", "", names(l)))
  })
}
# integer labels over jl_data() nodes to a joint_lab() vector
node_lab <- function(d, lab) {
  stats::setNames(
    as.character(lab),
    paste(rep(d$species, lengths(d$genes)), unlist(d$genes), sep = "\t")
  )
}

# ---- Per-run cache ----
run_file <- function(ds, engine, lambda, seed, half) {
  file.path(OUT, sprintf(
    "%s_runs_%s_%g_%d_%s.rds", ds, engine, lambda, seed, half
  ))
}

# Section 11 wood runs into the per-run cache (two_step_match -> match).
# The old files stay; the AUROC cache is copied under the new ids.
if (DS == "wood") {
  for (half in c("A", "B")) {
    f_old <- file.path(OUT, sprintf("wood_runs_%s.rds", half))
    if (!file.exists(f_old)) next
    n_new <- 0L
    for (r in readRDS(f_old)) {
      old_id <- sprintf("%s_l%g_s%d", r$engine, r$lambda, r$seed)
      r$engine <- sub("^two_step_match$", "match", r$engine)
      f <- run_file("wood", r$engine, r$lambda, r$seed, half)
      if (file.exists(f)) next
      r$params <- list(source = basename(f_old))
      saveRDS(r, f)
      n_new <- n_new + 1L
      a_old <- file.path(OUT, "wood_auroc", paste0(old_id, ".rds"))
      a_new <- file.path(OUT, "wood_auroc", sprintf(
        "%s_l%g_s%d.rds", r$engine, r$lambda, r$seed
      ))
      if (file.exists(a_old) && !file.exists(a_new)) {
        file.copy(a_old, a_new)
      }
    }
    if (n_new) log_t("migrated %d wood runs of half %s", n_new, half)
  }
}

# ---- Engine registry ----
ENGINES <- list(
  two_step = list(lambdas = 0, seeds = SEEDS, run = function(ctx, lam, s) {
    t0 <- Sys.time()
    two <- lapply(ctx$nets, function(n) {
      detect_modules(n, seed = s)$modules
    })
    list(lab = joint_lab(two), secs = secs_since(t0))
  }),
  match = list(lambdas = LAMBDA, seeds = SEEDS, run = function(ctx, lam, s) {
    two <- get_run("two_step", 0, s, ctx)
    f2 <- jl$joint_louvain(
      lambda = lam, seed = s, init = per_species(two$lab), freeze = TRUE,
      data = ctx$dat()
    )
    list(
      lab = joint_lab(lapply(f2$modules, `[[`, "modules")),
      secs = f2$params$secs + two$secs, quality = f2$quality,
      params = f2$params
    )
  }),
  E2 = list(lambdas = LAMBDA, seeds = SEEDS, run = function(ctx, lam, s) {
    fe <- jl$joint_louvain(lambda = lam, seed = s, data = ctx$dat())
    list(
      lab = joint_lab(lapply(fe$modules, `[[`, "modules")),
      secs = fe$params$secs, quality = fe$quality,
      quality_parts = fe$quality_parts, params = fe$params
    )
  }),
  E1 = list(lambdas = 32, seeds = 1L, run = function(ctx, lam, s) {
    t0 <- Sys.time()
    w0 <- stats::median(unlist(lapply(names(ctx$nets), function(sp) {
      e <- rcomplex:::.joint_edges(ctx$nets[[sp]], sp)
      e$w / max(e$w)
    })))
    f1 <- joint_modules(ctx$nets, ctx$ortho,
      weight = lam * w0, resolution = 6, seed = s
    )
    list(
      lab = joint_lab(lapply(f1$modules, `[[`, "modules")),
      secs = secs_since(t0), params = f1$params
    )
  })
)
for (hook in c("hleiden", "bethe")) {
  if (!hook %in% ENG) next
  f <- file.path("dev/probes", paste0(hook, ".R"))
  if (!file.exists(f)) {
    message(sprintf("engine %s: %s not found, skipped", hook, f))
    ENG <- setdiff(ENG, hook)
    next
  }
  he <- new.env(parent = jl)
  sys.source(f, he)
  if (!is.list(he$halves_engines)) {
    stop(sprintf("%s must define `halves_engines` (a named list).", f))
  }
  ENGINES[names(he$halves_engines)] <- he$halves_engines
  ENG <- c(setdiff(ENG, hook), names(he$halves_engines))
}
bad <- setdiff(ENG, names(ENGINES))
if (length(bad)) {
  stop(sprintf(
    "engines: got %s; need some of %s.", toString(bad),
    toString(c(names(ENGINES), "hleiden", "bethe"))
  ))
}

# Cached run of one engine, lambda, seed and half; computed when missing
get_run <- function(engine, lambda, seed, ctx) {
  f <- run_file(DS, engine, lambda, seed, ctx$half)
  if (file.exists(f)) {
    return(readRDS(f))
  }
  t0 <- Sys.time()
  r <- ENGINES[[engine]]$run(ctx, lambda, seed)
  if (is.null(r)) {
    return(log_t(
      "%s half %s: %s lambda %g seed %d returned no labels",
      DS, ctx$half, engine, lambda, seed
    ))
  }
  if (is.null(r$secs)) r$secs <- secs_since(t0)
  if (is.null(names(r$lab))) r$lab <- node_lab(ctx$dat(), r$lab)
  r$engine <- engine
  r$lambda <- lambda
  r$seed <- seed
  saveRDS(r, f)
  log_t(
    "%s half %s: %s lambda %g seed %d, %.0f s", DS, ctx$half, engine,
    lambda, seed, r$secs
  )
  r
}

for (half in c("A", "B")) {
  ctx <- new.env()
  ctx$half <- half
  ctx$nets <- hv$nets[[half]]
  ctx$ortho <- ortho
  ctx$dat <- local({
    d <- NULL
    function() {
      if (is.null(d)) {
        t0 <- Sys.time()
        d <<- jl$jl_data(ctx$nets, ctx$ortho)
        log_t(
          "%s half %s: jl_data %d nodes, %.0f s", DS, half, d$n,
          secs_since(t0)
        )
      }
      d
    }
  })
  for (e in ENG) {
    for (s in ENGINES[[e]]$seeds) {
      for (lam in ENGINES[[e]]$lambdas) get_run(e, lam, s, ctx)
    }
  }
  rm(ctx)
  invisible(gc())
}
if ("runs" %in% args[-(1:2)]) quit(save = "no")

# ---- Scores ----
read_half <- function(half) {
  fs <- list.files(OUT, sprintf("^%s_runs_.*_%s\\.rds$", DS, half),
    full.names = TRUE
  )
  r <- lapply(fs, readRDS)
  names(r) <- vapply(r, id_of, "")
  r
}
id_of <- function(r) sprintf("%s_l%g_s%d", r$engine, r$lambda, r$seed)
runs <- list(A = read_half("A"), B = read_half("B"))
ids <- intersect(names(runs$A), names(runs$B))
ari <- mclust::adjustedRandIndex
for (h in names(runs)) {
  # two-step labels are per species; make them distinct across species
  for (id in grep("^two_step_l0", names(runs[[h]]), value = TRUE)) {
    l <- runs[[h]][[id]]$lab
    runs[[h]][[id]]$lab[] <- paste(sp_of(l), l)
  }
}

# Labels of genes in modules of MIN_SIZE or more; by species when per_sp
big <- function(lab, per_sp = FALSE) {
  g <- if (per_sp) paste(sp_of(lab), lab) else lab
  lab[ave(seq_along(g), g, FUN = length) >= MIN_SIZE]
}
split_ari <- function(a, b) {
  g <- intersect(names(a), names(b))
  if (length(g) < 2L) NA else ari(a[g], b[g])
}

# Multi-species hogs of a half: copy keys and cross-species copy pairs
hog_index <- function(nets) {
  key <- paste(ortho$species, ortho$gene, sep = "\t")
  have <- unlist(lapply(names(nets), function(s) {
    paste(s, rownames(nets[[s]]$network), sep = "\t")
  }))
  ok <- key %in% have
  o <- data.frame(key = key[ok], hog = ortho$hog[ok], s = ortho$species[ok])
  nsp <- tapply(o$s, o$hog, function(x) length(unique(x)))
  o <- o[nsp[o$hog] >= 2L, ]
  dt <- data.table::as.data.table(o)
  pr <- merge(dt, dt, by = "hog", allow.cartesian = TRUE)
  pr <- pr[pr$s.x < pr$s.y, ]
  list(o = o, a = pr$key.x, b = pr$key.y)
}
hix <- lapply(hv$nets, hog_index)
hog_scores <- function(lab, hx) {
  l <- lab[hx$o$key]
  c(
    purity = mean(tapply(l, hx$o$hog, function(x) length(unique(x)) == 1)),
    pairs = mean(lab[hx$a] == lab[hx$b])
  )
}

rows <- list()
for (id in ids) {
  ra <- runs$A[[id]]
  rb <- runs$B[[id]]
  two_id <- sprintf("two_step_l0_s%d", ra$seed)
  per_sp <- vapply(SP, function(s) {
    split_ari(
      big(ra$lab[sp_of(ra$lab) == s], TRUE),
      big(rb$lab[sp_of(rb$lab) == s], TRUE)
    )
  }, 0)
  vs_two <- function(h) {
    r <- runs[[h]][[id]]
    t2 <- runs[[h]][[two_id]]
    if (is.null(t2)) {
      return(NA_real_)
    }
    stats::median(vapply(SP, function(s) {
      i <- sp_of(r$lab) == s
      ari(r$lab[i], t2$lab[names(r$lab)[i]])
    }, 0))
  }
  hs_a <- hog_scores(ra$lab, hix$A)
  hs_b <- hog_scores(rb$lab, hix$B)
  nmod <- function(lab) sum(table(lab) >= MIN_SIZE)
  rows[[id]] <- data.frame(
    arm = DS, engine = ra$engine, lambda = ra$lambda, delta = NA,
    d = NA, seed = ra$seed, half = "AB",
    n_modules_A = nmod(ra$lab), n_modules_B = nmod(rb$lab),
    ari_sp_median = stats::median(per_sp), ari_sp_min = min(per_sp),
    ari_sp_max = max(per_sp),
    ari_joint = split_ari(big(ra$lab), big(rb$lab)),
    purity_A = hs_a[["purity"]], purity_B = hs_b[["purity"]],
    pairs_A = hs_a[["pairs"]], pairs_B = hs_b[["pairs"]],
    ari_two_A = vs_two("A"), ari_two_B = vs_two("B"),
    secs_A = ra$secs, secs_B = rb$secs
  )
}
tab <- do.call(rbind, rows)
f_tab <- file.path(OUT, paste0(DS, ".tsv"))
utils::write.table(tab, f_tab, sep = "\t", quote = FALSE, row.names = FALSE)
log_t("%s split-half scores written: %d rows", DS, nrow(tab))
if ("scores" %in% args[-(1:2)]) quit(save = "no")

# ---- Held-out AUROC z: half-A modules on half-B networks ----
# The probe module_auroc() reads Species1 / Species2 columns; 0.4.1's
# .ortholog_pair_index() wants gene1 / gene2. A shim renames both ways.
ma <- new.env(parent = asNamespace("rcomplex"))
sys.source("dev/probes/module_auroc.R", ma)
Rcpp::sourceCpp("dev/probes/module_auroc.cpp", env = ma)
ma$.ortholog_pair_index <- function(net1, net2, orthologs) {
  o <- data.frame(
    gene1 = orthologs$Species1, gene2 = orthologs$Species2,
    hog = orthologs$hog
  )
  r <- rcomplex:::.ortholog_pair_index(net1, net2, o)
  r$orthologs$Species1 <- r$orthologs$gene1
  r$orthologs$Species2 <- r$orthologs$gene2
  r
}
pair_ortho <- lapply(stats::setNames(SP, SP), function(s) {
  lapply(stats::setNames(SP, SP), function(t) {
    a <- ortho[ortho$species == s, c("gene", "hog")]
    b <- ortho[ortho$species == t, c("gene", "hog")]
    m <- merge(a, b, by = "hog")
    data.frame(Species1 = m$gene.x, Species2 = m$gene.y, hog = m$hog)
  })
})
N_NULL <- 100L
A_DIR <- file.path(OUT, paste0(DS, "_auroc"))
auroc_run <- function(id) {
  f <- file.path(A_DIR, paste0(id, ".rds"))
  if (file.exists(f)) {
    out <- readRDS(f)
    if (!is.null(out)) out$id <- id
    return(out)
  }
  lab <- runs$A[[id]]$lab
  out <- list()
  for (s in SP) {
    l <- lab[sp_of(lab) == s]
    l <- big(l, TRUE)
    if (!length(l)) next
    mods <- stats::setNames(l, sub("^[^\t]*\t", "", names(l)))
    for (t in SP) {
      r <- tryCatch(suppressMessages(if (s == t) {
        ma$module_replication(mods, hv$nets$B[[t]],
          n_null = N_NULL, max_draws = N_NULL, seed = 1
        )
      } else {
        ma$module_auroc(mods, hv$nets$A[[s]], hv$nets$B[[t]],
          pair_ortho[[s]][[t]],
          n_null = N_NULL, max_draws = N_NULL, seed = 1
        )
      }), error = function(e) NULL)
      if (is.null(r)) next
      out[[length(out) + 1L]] <- data.frame(
        id = id, s = s, t = t, module = r$module, n_hogs = r$n_hogs,
        auroc = r$auroc, z = r$z
      )
    }
  }
  out <- do.call(rbind, out)
  saveRDS(out, f)
  out
}
dir.create(A_DIR, showWarnings = FALSE)
t0 <- Sys.time()
z <- parallel::mclapply(ids, auroc_run,
  mc.cores = N_CORES, mc.preschedule = FALSE
)
bad <- vapply(z, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad)) log_t("auroc failed for: %s", toString(ids[bad]))
z <- do.call(rbind, z[!bad])
log_t(
  "auroc done: %d rows, %.1f min", nrow(z), secs_since(t0) / 60
)
saveRDS(z, file.path(OUT, paste0(DS, "_auroc.rds")))
zs <- do.call(rbind, lapply(split(z, z$id), function(x) {
  cr <- x[x$s != x$t, ]
  pm <- tapply(cr$z, paste(cr$s, cr$t), stats::median)
  data.frame(
    id = x$id[1], z_cross_median = stats::median(cr$z),
    z_cross_pair_median = stats::median(pm),
    z_within_median = stats::median(x$z[x$s == x$t]),
    n_module_tests = nrow(cr)
  )
}))
tab$id <- rownames(tab)
tab <- merge(tab, zs, by = "id", all.x = TRUE, sort = FALSE)
utils::write.table(tab, f_tab, sep = "\t", quote = FALSE, row.names = FALSE)

# Seed summary per engine and lambda: mean, min, max over seeds
m <- data.frame(
  engine = tab$engine, lambda = tab$lambda,
  modules = (tab$n_modules_A + tab$n_modules_B) / 2,
  ari_sp_median = tab$ari_sp_median, ari_sp_min = tab$ari_sp_min,
  ari_sp_max = tab$ari_sp_max, ari_joint = tab$ari_joint,
  purity = (tab$purity_A + tab$purity_B) / 2,
  pairs = (tab$pairs_A + tab$pairs_B) / 2,
  ari_two = (tab$ari_two_A + tab$ari_two_B) / 2,
  z_cross = tab$z_cross_median, z_within = tab$z_within_median,
  secs = (tab$secs_A + tab$secs_B) / 2
)
key <- m[, c("engine", "lambda")]
val <- m[, -(1:2)]
agg <- function(f, suffix) {
  a <- stats::aggregate(val, key, f)
  names(a)[-(1:2)] <- paste0(names(a)[-(1:2)], suffix)
  a
}
smry <- Reduce(
  function(x, y) merge(x, y, by = c("engine", "lambda")),
  list(
    stats::aggregate(list(n_seeds = m$secs), key, length),
    agg(mean, ""), agg(min, "_lo"), agg(max, "_hi")
  )
)
utils::write.table(smry, file.path(OUT, paste0(DS, "_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
log_t("%s done", DS)
