# Probe E2, data arm 2: wood split halves.
# Expression and gene ids as prepare_data/cache/wood_041/run.R. Samples
# of each species split at random into halves A and B (seed 1); one
# network per half (density 0.03, store at the analysis threshold).
# Engines per half: two-step (detect_modules() per species), two-step +
# match (joint_louvain(init, freeze = TRUE)), E2 (joint_louvain()) at
# lambda in {0.25, 1, 4, 16, 64}, seeds 1-3; E1 (joint_modules(), weight
# 32 x median, resolution 6, seed 1).
# Scores on modules of 20 or more genes: per-species split-half ARI,
# joint split-half ARI, hog agreement, ARI against two-step, held-out
# AUROC z (half-A modules on half-B networks, every species pair).
#
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/joint_wood_halves.R
#   ... joint_wood_halves.R runs     (stop after the engine runs)
suppressPackageStartupMessages({
  library(rcomplex)
  library(Matrix)
})
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
E <- file.path(REPO, "prepare_data/evotree/data")
jl <- new.env(parent = asNamespace("rcomplex"))
jl$JL_DIR <- "dev/probes"
sys.source("dev/probes/joint_louvain.R", jl)
LAMBDA <- c(0.25, 1, 4, 16, 64)
SEEDS <- 1:3
MIN_SIZE <- 20L
log_t <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
}

# ---- Networks per half (cached) ----
HALVES <- file.path(OUT, "halves.rds")
if (!file.exists(HALVES)) {
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
  expr <- setNames(lapply(sp, read_expr), sp)
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
  rm(og)
  stopifnot(!anyDuplicated(ortho[, c("species", "gene")]))
  set.seed(1)
  split_of <- lapply(expr, function(x) {
    sample(rep(c("A", "B"), length.out = ncol(x)))
  })
  nets <- list()
  for (half in c("A", "B")) {
    nets[[half]] <- lapply(setNames(sp, sp), function(s) {
      x <- expr[[s]][, split_of[[s]] == half]
      t0 <- Sys.time()
      n <- compute_network(x,
        density = 0.03, sparse = TRUE, store_density = 0.03,
        n_cores = 8
      )
      log_t(
        "network %s %s: %d genes, %d samples, %.0f s", half, s,
        n$n_genes, ncol(x), as.numeric(difftime(Sys.time(), t0,
          units = "secs"
        ))
      )
      n
    })
  }
  saveRDS(list(nets = nets, ortho = ortho, split = split_of), HALVES)
  rm(expr)
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

# ---- Engine runs per half (cached per half) ----
for (half in c("A", "B")) {
  f_runs <- file.path(OUT, sprintf("wood_runs_%s.rds", half))
  if (file.exists(f_runs)) next
  nets <- hv$nets[[half]]
  t0 <- Sys.time()
  dat <- jl$jl_data(nets, ortho)
  log_t(
    "half %s: jl_data %d nodes, %.0f s", half, dat$n,
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  runs <- list()
  for (seed in SEEDS) {
    t0 <- Sys.time()
    two <- lapply(nets, function(n) detect_modules(n, seed = seed)$modules)
    t_two <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    runs[[length(runs) + 1L]] <- list(
      engine = "two_step", lambda = 0, seed = seed, lab = joint_lab(two),
      secs = t_two
    )
    log_t("half %s seed %d two-step %.0f s", half, seed, t_two)
    for (lam in LAMBDA) {
      f2 <- jl$joint_louvain(
        lambda = lam, seed = seed, init = two, freeze = TRUE, data = dat
      )
      runs[[length(runs) + 1L]] <- list(
        engine = "two_step_match", lambda = lam, seed = seed,
        lab = joint_lab(lapply(f2$modules, `[[`, "modules")),
        secs = f2$params$secs + t_two, quality = f2$quality
      )
      fe <- jl$joint_louvain(lambda = lam, seed = seed, data = dat)
      runs[[length(runs) + 1L]] <- list(
        engine = "E2", lambda = lam, seed = seed,
        lab = joint_lab(lapply(fe$modules, `[[`, "modules")),
        secs = fe$params$secs, quality = fe$quality,
        quality_parts = fe$quality_parts
      )
      log_t(
        "half %s seed %d lambda %g: match %.0f s Q %.3f; E2 %.0f s Q %.3f",
        half, seed, lam, f2$params$secs, f2$quality, fe$params$secs,
        fe$quality
      )
    }
  }
  rm(dat)
  invisible(gc())
  t0 <- Sys.time()
  w0 <- stats::median(unlist(lapply(SP, function(s) {
    e <- rcomplex:::.joint_edges(nets[[s]], s)
    e$w / max(e$w)
  })))
  f1 <- joint_modules(nets, ortho,
    weight = 32 * w0, resolution = 6,
    seed = 1
  )
  runs[[length(runs) + 1L]] <- list(
    engine = "E1", lambda = 32, seed = 1,
    lab = joint_lab(lapply(f1$modules, `[[`, "modules")),
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  log_t("half %s E1 %.0f s", half, runs[[length(runs)]]$secs)
  saveRDS(runs, f_runs)
}
if ("runs" %in% commandArgs(TRUE)) quit(save = "no")

# ---- Scores ----
runs <- list(
  A = readRDS(file.path(OUT, "wood_runs_A.rds")),
  B = readRDS(file.path(OUT, "wood_runs_B.rds"))
)
id_of <- function(r) sprintf("%s_l%g_s%d", r$engine, r$lambda, r$seed)
ari <- mclust::adjustedRandIndex
sp_of <- function(lab) sub("\t.*", "", names(lab))
for (h in names(runs)) {
  names(runs[[h]]) <- vapply(runs[[h]], id_of, "")
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
  ok <- ortho$gene %in% unlist(lapply(nets, function(n) {
    rownames(n$network)
  }))
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
for (id in names(runs$A)) {
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
    stats::median(vapply(SP, function(s) {
      i <- sp_of(r$lab) == s
      ari(r$lab[i], t2$lab[names(r$lab)[i]])
    }, 0))
  }
  hs_a <- hog_scores(ra$lab, hix$A)
  hs_b <- hog_scores(rb$lab, hix$B)
  nmod <- function(lab) sum(table(lab) >= MIN_SIZE)
  rows[[id]] <- data.frame(
    arm = "wood", engine = ra$engine, lambda = ra$lambda, delta = NA,
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
utils::write.table(tab, file.path(OUT, "wood.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
log_t("split-half scores written: %d rows", nrow(tab))
if ("scores" %in% commandArgs(TRUE)) quit(save = "no")

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
pair_ortho <- lapply(setNames(SP, SP), function(s) {
  lapply(setNames(SP, SP), function(t) {
    a <- ortho[ortho$species == s, c("gene", "hog")]
    b <- ortho[ortho$species == t, c("gene", "hog")]
    m <- merge(a, b, by = "hog")
    data.frame(Species1 = m$gene.x, Species2 = m$gene.y, hog = m$hog)
  })
})
N_NULL <- 100L
auroc_run <- function(id) {
  f <- file.path(OUT, "wood_auroc", paste0(id, ".rds"))
  if (file.exists(f)) return(readRDS(f))
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
dir.create(file.path(OUT, "wood_auroc"), showWarnings = FALSE)
t0 <- Sys.time()
z <- parallel::mclapply(names(runs$A), auroc_run, mc.cores = 8,
  mc.preschedule = FALSE
)
bad <- vapply(z, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad)) log_t("auroc failed for: %s", toString(names(runs$A)[bad]))
z <- do.call(rbind, z[!bad])
log_t("auroc done: %d rows, %.1f min", nrow(z),
  as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(z, file.path(OUT, "wood_auroc.rds"))
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
utils::write.table(tab, file.path(OUT, "wood.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
log_t("wood done")
