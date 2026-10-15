# Probe E2, data arm 1: species-blocked planted model.
# Six species, 4,000 genes each, 12 planted modules present in every
# species (80-200 genes each, drawn per species), the rest background.
# Binary co-expression: a module gene has d_in expected edges inside its
# module and d_out to random genes; a background gene has d_in + d_out
# random edges. Hogs follow the wood copy-count profile; a copy leaves
# its hog's module with probability delta (a defector). Run with `calib`
# for the d_in calibration grid only.
#
# Engines: two-step + match (detect_modules() per species, then
# joint_louvain(init, freeze = TRUE)), E2 (joint_louvain()), E1
# (joint_modules(resolution = 6), seed 1 only). One seed per row: the
# seed draws the planted data and seeds the engine.
#
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/joint_planted.R
#   ... joint_planted.R summary   (tables from planted.tsv only)
suppressPackageStartupMessages(library(rcomplex))
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
dir.create(file.path(OUT, "planted_labels"), showWarnings = FALSE)
jl <- new.env(parent = asNamespace("rcomplex"))
jl$JL_DIR <- "dev/probes"
sys.source("dev/probes/joint_louvain.R", jl)

S <- 6L
N <- 4000L
K <- 12L
SP <- paste0("S", seq_len(S))
LAMBDA <- c(0.25, 1, 4, 16, 64)
E1_W <- c(4, 32, 128)
# Spec (20, 40) and (4, 6) gave single-species AMI 0.54 and 0.13
# (seed 1); the calibration grid (`calib`) picked these for 0.90 and 0.49
DEG <- list(d60 = c(34, 26), d10 = c(6, 4))
DELTA <- c(0, 0.1, 0.3)
SEEDS <- 1:3

# Wood copy profile: copies per species of every multi-species hog
prof <- local({
  f <- file.path(OUT, "ortho_full.rds")
  o <- if (file.exists(f)) {
    readRDS(f)
  } else {
    readRDS(file.path(REPO, "prepare_data/cache/snr/data.rds"))$ortho
  }
  m <- unclass(table(o$hog, o$species))
  m[rowSums(m > 0) >= 2, , drop = FALSE]
})

# ---- Metrics ----
ami <- function(x, y) {
  t <- table(x, y)
  n <- sum(t)
  a <- rowSums(t)
  b <- colSums(t)
  ent <- function(p) {
    p <- p[p > 0] / n
    -sum(p * log(p))
  }
  nz <- which(t > 0, arr.ind = TRUE)
  v <- t[nz]
  mi <- sum(v / n * log(n * v / (a[nz[, 1]] * b[nz[, 2]])))
  lf <- lfactorial
  emi <- 0
  for (i in seq_along(a)) {
    for (j in seq_along(b)) {
      lo <- max(1, a[i] + b[j] - n)
      hi <- min(a[i], b[j])
      if (lo > hi) next
      k <- lo:hi
      lp <- lf(a[i]) + lf(b[j]) + lf(n - a[i]) + lf(n - b[j]) - lf(n) -
        lf(k) - lf(a[i] - k) - lf(b[j] - k) - lf(n - a[i] - b[j] + k)
      emi <- emi + sum(k / n * log(n * k / (a[i] * b[j])) * exp(lp))
    }
  }
  (mi - emi) / ((ent(a) + ent(b)) / 2 - emi)
}

# Module of most copies per hog; a tie goes to the smallest label
main_of <- function(lab, hog) {
  m <- tapply(lab, hog, function(x) {
    u <- sort(unique(x))
    u[which.max(tabulate(match(x, u)))]
  })
  m[hog]
}

scores <- function(lab, pl) {
  tr <- pl$truth
  lab <- lab[names(tr)]
  mods <- split(names(tr), tr)[as.character(seq_len(K))]
  found <- split(names(lab), lab)
  jac <- vapply(mods, function(g) {
    l <- unique(lab[g])
    max(vapply(found[as.character(l)], function(f) {
      length(intersect(g, f)) / length(union(g, f))
    }, 0))
  }, 0)
  shared <- vapply(mods, function(g) {
    sp <- sub("\t.*", "", g)
    maj <- tapply(lab[g], sp, function(x) names(which.max(table(x))))
    length(unique(maj)) == 1L
  }, TRUE)
  o <- pl$ortho
  key <- paste(o$species, o$gene, sep = "\t")
  pred <- lab[key] != main_of(lab[key], o$hog)
  tp <- sum(pred & o$defector)
  c(
    ami = ami(tr, lab), ari = mclust::adjustedRandIndex(tr, lab),
    n_modules = length(unique(lab)),
    jac_median = stats::median(jac), jac_frac50 = mean(jac > 0.5),
    shared = mean(shared),
    def_precision = if (sum(pred)) tp / sum(pred) else NA,
    def_recall = if (sum(o$defector)) tp / sum(o$defector) else NA,
    n_defectors = sum(o$defector), n_pred = sum(pred)
  )
}

# ---- Planted data ----
planted <- function(seed, d_in, d_out, delta) {
  set.seed(seed)
  # class of every gene: 1..K module, K + 1 background
  cls <- lapply(SP, function(s) {
    size <- sample(80:200, K, replace = TRUE)
    sample(c(rep(seq_len(K), size), rep(K + 1L, N - sum(size))))
  })
  names(cls) <- SP
  nets <- lapply(SP, function(s) {
    cl <- cls[[s]]
    ij <- NULL
    for (k in seq_len(K)) {
      g <- which(cl == k)
      m <- length(g)
      pr <- utils::combn(m, 2)
      keep <- stats::runif(ncol(pr)) < d_in / (m - 1)
      ij <- rbind(ij, cbind(g[pr[1, keep]], g[pr[2, keep]]))
    }
    r <- ifelse(cl <= K, d_out, d_in + d_out)
    ne <- round(sum(r) / 2)
    a <- sample.int(N, ne, replace = TRUE, prob = r)
    b <- sample.int(N, ne, replace = TRUE, prob = r)
    ij <- rbind(ij, cbind(a, b))
    ij <- ij[ij[, 1] != ij[, 2], ]
    ij <- unique(cbind(pmin(ij[, 1], ij[, 2]), pmax(ij[, 1], ij[, 2])))
    nm <- sprintf("%s_g%04d", s, seq_len(N))
    x <- Matrix::sparseMatrix(
      i = c(ij[, 1], ij[, 2]), j = c(ij[, 2], ij[, 1]), x = 1,
      dims = c(N, N), dimnames = list(nm, nm)
    )
    list(network = x, threshold = 1, store_threshold = 1, n_genes = N)
  })
  names(nets) <- SP

  # Hogs: wood profiles until half the genes are covered
  free <- lapply(cls, function(cl) split(seq_len(N), cl))
  pick <- function(s, k) {
    f <- free[[s]][[k]]
    if (!length(f)) {
      return(NA_integer_)
    }
    i <- if (length(f) == 1L) f else sample(f, 1L)
    free[[s]][[k]] <<- setdiff(f, i)
    i
  }
  tot <- table(factor(unlist(cls), levels = seq_len(K + 1L)))
  rows <- list()
  covered <- 0
  h <- 0L
  while (covered < S * N / 2) {
    h <- h + 1L
    cp <- prof[sample.int(nrow(prof), 1L), ]
    home <- sample.int(K + 1L, 1L, prob = tot)
    for (k in seq_len(S)) {
      for (cc in seq_len(cp[[k]])) {
        cl <- if (stats::runif(1) < delta) {
          sample(setdiff(seq_len(K + 1L), home), 1L,
            prob = tot[-home]
          )
        } else {
          home
        }
        i <- pick(SP[k], cl)
        # a full class drops the copy, so delta 0 has no defector
        if (is.na(i)) next
        rows[[length(rows) + 1L]] <- data.frame(
          species = SP[k], gene = rownames(nets[[k]]$network)[i],
          hog = sprintf("H%05d", h), home = home, cls = cls[[k]][i]
        )
        covered <- covered + 1
      }
    }
  }
  ortho <- do.call(rbind, rows)
  ortho$defector <- ortho$cls != ortho$home
  truth <- unlist(lapply(SP, function(s) {
    cl <- cls[[s]]
    lab <- ifelse(cl <= K, as.character(cl), paste0("bg_", s))
    stats::setNames(lab, paste(s, rownames(nets[[s]]$network), sep = "\t"))
  }))
  list(nets = nets, ortho = ortho, truth = truth, cls = cls)
}

joint_lab <- function(mods) {
  unlist(lapply(names(mods), function(s) {
    m <- mods[[s]]
    stats::setNames(as.character(m), paste(s, names(m), sep = "\t"))
  }))
}

# ---- Calibration: single-species detect_modules() ----
# AMI per species against module id + background, at delta 0.
calibrate <- function(deg, seeds) {
  out <- list()
  for (dn in names(deg)) {
    for (seed in seeds) {
      pl <- planted(seed, deg[[dn]][1], deg[[dn]][2], 0)
      a <- vapply(SP, function(s) {
        m <- detect_modules(pl$nets[[s]], seed = seed)$modules
        ami(pl$cls[[s]], m[rownames(pl$nets[[s]]$network)])
      }, 0)
      out[[length(out) + 1L]] <- data.frame(
        d = dn, d_in = deg[[dn]][1], d_out = deg[[dn]][2], seed = seed,
        ami_mean = mean(a), ami_min = min(a), ami_max = max(a),
        mean_degree = mean(Matrix::colSums(pl$nets[[1]]$network))
      )
    }
  }
  do.call(rbind, out)
}
if ("calib" %in% commandArgs(TRUE)) {
  grid <- list(
    d60_20 = c(20, 40), d60_30 = c(30, 30), d60_32 = c(32, 28),
    d60_34 = c(34, 26), d60_36 = c(36, 24), d60_40 = c(40, 20),
    d60_50 = c(50, 10), d10_4 = c(4, 6), d10_5 = c(5, 5),
    d10_6 = c(6, 4), d10_7 = c(7, 3), d10_8 = c(8, 2)
  )
  cal <- calibrate(grid, 1L)
  print(cal)
  utils::write.table(cal, file.path(OUT, "planted_calibration_grid.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  quit(save = "no")
}
summarize <- function() {
  x <- utils::read.delim(file.path(OUT, "planted.tsv"))
  m <- c(
    "ami", "jac_median", "jac_frac50", "shared", "def_precision",
    "def_recall", "n_modules", "secs"
  )
  key <- x[, c("d", "delta", "engine", "lambda")]
  sm <- stats::aggregate(x[, m], key, mean)
  sd <- stats::aggregate(x[, c("ami", "def_recall")], key, stats::sd)
  names(sd)[5:6] <- c("ami_sd", "def_recall_sd")
  sm <- merge(sm, sd)
  sm <- sm[order(sm$d, sm$delta, sm$engine, sm$lambda), ]
  utils::write.table(sm, file.path(OUT, "planted_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  # best lambda (highest mean AMI) per engine, d and delta
  best <- do.call(rbind, lapply(
    split(sm, sm[, c("d", "delta", "engine")], drop = TRUE),
    function(b) b[which.max(b$ami), ]
  ))
  best <- best[order(best$d, best$delta, best$engine), ]
  utils::write.table(best, file.path(OUT, "planted_best.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(best)
}
if ("summary" %in% commandArgs(TRUE)) {
  print(summarize(), digits = 3)
  quit(save = "no")
}

calib <- calibrate(DEG, SEEDS)
print(calib)
utils::write.table(calib, file.path(OUT, "planted_calibration.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- Sweep ----
rows <- list()
add <- function(engine, lambda, fitlab, secs, pl, key) {
  sc <- scores(fitlab, pl)
  rows[[length(rows) + 1L]] <<- data.frame(
    arm = "planted", engine = engine, lambda = lambda, delta = key$delta,
    d = key$d, seed = key$seed, half = NA, t(sc), secs = secs
  )
  saveRDS(fitlab, file.path(OUT, "planted_labels", sprintf(
    "%s_%s_l%g_delta%g_seed%d.rds", key$d, engine, lambda, key$delta,
    key$seed
  )))
}
t_all <- Sys.time()
for (dn in names(DEG)) {
  for (delta in DELTA) {
    for (seed in SEEDS) {
      key <- list(d = dn, delta = delta, seed = seed)
      pl <- planted(seed, DEG[[dn]][1], DEG[[dn]][2], delta)
      o <- pl$ortho[, c("species", "gene", "hog")]
      dat <- jl$jl_data(pl$nets, o)
      t0 <- Sys.time()
      two <- lapply(pl$nets, function(n) {
        detect_modules(n, seed = seed)$modules
      })
      t_two <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      two_lab <- joint_lab(two)
      two_lab[] <- paste(sub("\t.*", "", names(two_lab)), two_lab)
      add("two_step", 0, two_lab, t_two, pl, key)
      for (lam in LAMBDA) {
        f2 <- jl$joint_louvain(
          lambda = lam, seed = seed, init = two, freeze = TRUE,
          data = dat
        )
        add(
          "two_step_match", lam,
          joint_lab(lapply(f2$modules, `[[`, "modules")),
          f2$params$secs + t_two, pl, key
        )
        fe <- jl$joint_louvain(lambda = lam, seed = seed, data = dat)
        add(
          "E2", lam, joint_lab(lapply(fe$modules, `[[`, "modules")),
          fe$params$secs, pl, key
        )
      }
      if (seed == 1L) {
        for (w in E1_W) {
          t0 <- Sys.time()
          f1 <- joint_modules(pl$nets, o,
            weight = w, resolution = 6,
            seed = 1
          )
          add(
            "E1", w, joint_lab(lapply(f1$modules, `[[`, "modules")),
            as.numeric(difftime(Sys.time(), t0, units = "secs")), pl, key
          )
        }
      }
      tab <- do.call(rbind, rows)
      utils::write.table(tab, file.path(OUT, "planted.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      cat(sprintf(
        "%s delta %g seed %d done, %d rows, %.1f min\n", dn, delta, seed,
        nrow(tab), as.numeric(difftime(Sys.time(), t_all, units = "mins"))
      ))
    }
  }
}
print(summarize(), digits = 3)
cat("planted done\n")
