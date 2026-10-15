# Probe WP15b acceptance (spec arm B checks 1-4) and the planted sweep.
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/hleiden_check.R
#   ... hleiden_check.R planted   (checks 1-2, then the planted sweep)
#   ... hleiden_check.R planted a0   (sweep at alpha 0 only)
#   ... hleiden_check.R summary   (tables from planted_hleiden.tsv only)
# Check 2 (every module connected in the 2-section) is asserted on every
# fixture run and every planted run.
suppressPackageStartupMessages(library(rcomplex))
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
hl <- new.env(parent = asNamespace("rcomplex"))
hl$JL_DIR <- "dev/probes"
sys.source("dev/probes/hleiden.R", hl)
ari <- mclust::adjustedRandIndex
args <- commandArgs(TRUE)
SCHED <- c(0, 0.5, 1)

# Check 2 on one fit: every module connected in the 2-section
connected <- function(fit, d, what) {
  nd <- hl$hl_disconnected(fit$labels, d)
  if (nd > 0) stop(sprintf("check 2 FAILED: %s, %d disconnected", what, nd))
  invisible(TRUE)
}

# Plain-R strict term T (per unit lambda), the reference for check 1b
strict_r <- function(lab, d) {
  l <- d$link
  l$s <- d$spi[l$node]
  l$c <- stats::ave(l$node, l$hog, l$s, FUN = length)
  k <- tapply(l$s, l$hog, function(x) length(unique(x)))
  l <- l[k[l$hog] >= 2, ]
  l$k <- k[l$hog]
  l$C <- lab[l$node]
  f <- stats::aggregate(1 / l$c, l[c("hog", "C", "s")], sum)
  key <- paste(f$hog, f$C, sep = "\r")
  pr <- tapply(f$x, key, prod)
  ns <- tapply(f$x, key, length)
  kh <- k[sub("\r.*", "", names(pr))]
  obs <- sum((kh * (kh - 1) * pr)[ns == kh])
  dg <- tapply(1 / l$c, l$node, sum)
  node <- as.integer(names(dg))
  s <- d$spi[node] + 1L
  Hs <- tapply(dg, s, sum)[as.character(s)]
  ul <- unique(lab)
  P <- as.matrix(Matrix::sparseMatrix(
    i = match(lab[node], ul), j = s, x = as.vector(dg / Hs),
    dims = c(length(ul), length(d$species))
  ))
  sets <- tapply(l$s, l$hog, function(x) sort(unique(x)) + 1L)
  sets <- sets[!vapply(sets, is.null, TRUE)]
  sk <- vapply(sets, paste, "", collapse = ",")
  nul <- 0
  for (u in unique(sk)) {
    sv <- sets[[match(u, sk)]]
    w <- sum(sk == u) * length(sv) * (length(sv) - 1)
    nul <- nul + w * sum(apply(P[, sv, drop = FALSE], 1, prod))
  }
  bar <- mean(tapply(Matrix::colSums(d$co), d$spi, sum))
  (obs - nul) / bar
}

q_r <- function(lab, d, lambda, alpha, gamma = 1) {
  q0 <- hl$joint_quality(lab, d, lambda = lambda, gamma = gamma)
  sum(q0$species) + lambda * ((1 - alpha) * q0$hog_raw +
    alpha * strict_r(lab, d))
}

# ---- Fixtures of joint_louvain_check.R ----
set.seed(1)
toy_net <- function(prefix) {
  m <- matrix(0, 8, 8)
  m[upper.tri(m)] <- stats::runif(28) * stats::rbinom(28, 1, 0.6)
  m <- m + t(m)
  dimnames(m) <- list(paste0(prefix, 1:8), paste0(prefix, 1:8))
  list(network = m, threshold = 1e-9)
}
toy <- list(A = toy_net("A"), B = toy_net("B"))
toy_o <- data.frame(
  species = c("A", "A", "B", "A", "B", "B", "A", "B"),
  gene = c("A1", "A2", "B1", "A3", "B2", "B3", "A4", "B4"),
  hog = c("H1", "H1", "H1", "H2", "H2", "H2", "H3", "H3")
)
d <- hl$jl_data(toy, toy_o)
# Three species, paralogs and partial hogs: q(0) and q(1) differ here
toy3 <- list(A = toy_net("A"), B = toy_net("B"), C = toy_net("C"))
toy3_o <- data.frame(
  species = c("A", "A", "B", "C", "A", "B", "C", "C", "A", "B", "B", "C"),
  gene = c(
    "A1", "A2", "B1", "C1", "A3", "B2", "C2", "C3", "A4", "B3",
    "B4", "C4"
  ),
  hog = c(
    "H1", "H1", "H1", "H1", "H2", "H2", "H2", "H2", "H3", "H3",
    "H4", "H4"
  )
)
d3 <- hl$jl_data(toy3, toy3_o)

helper <- new.env(parent = asNamespace("rcomplex"))
sys.source("tests/testthat/helper-preservation.R", helper)
for (ex in parse("tests/testthat/test-joint-modules.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
    as.character(ex[[2]]) %in% c("jm_fixture", "long_ortho")) {
    eval(ex, helper)
  }
}
f <- helper$jm_fixture()
jm_nets <- list(A = f$netA, B = f$netB)
long2 <- helper$long_ortho(f$ortho)
long2$hog[long2$species == "B" & long2$gene == "B0041"] <- "H1"
dj <- hl$jl_data(jm_nets, long2)

if (!"summary" %in% args) {
  # ---- 1. q(0) equals joint_quality(); q(1) strict term equals plain R ----
  fx <- list(toy = d, toy3 = d3, jm = dj)
  for (nm in names(fx)) {
    dd <- fx[[nm]]
    r <- t(vapply(1:20, function(i) {
      lab <- sample.int(if (i %% 2) 4 else 40, dd$n, replace = TRUE)
      lam <- c(0, 0.3, 2.5, 16)[i %% 4 + 1]
      a <- hl$hl_quality(lab, dd, lambda = lam, gamma = 1.3, alpha = 0)
      b <- hl$joint_quality(lab, dd, lambda = lam, gamma = 1.3)
      c(
        q0 = abs(a$quality - b$quality),
        strict = abs(a$strict_raw - strict_r(lab, dd))
      )
    }, numeric(2)))
    cat(sprintf(paste0(
      "check 1: %s max |q(0) - joint_quality()| %.1e; ",
      "max |T C++ - T plain R| %.1e (20 labelings)\n"
    ), nm, max(r[, "q0"]), max(r[, "strict"])))
    stopifnot(max(r) < 1e-10)
  }

  # 1b. No single-gene move improves q(alpha) at the returned partition
  d <- d3
  for (lam in c(1, 5)) {
    for (a in list(0, 0.5, 1, SCHED)) {
      fit <- hl$hleiden(lambda = lam, alpha = a, seed = 1, data = d)
      connected(fit, d, "toy3")
      al <- a[length(a)]
      q0 <- q_r(fit$labels, d, lam, al)
      best <- -Inf
      for (i in seq_len(d$n)) {
        for (cc in c(unique(fit$labels), max(fit$labels) + 1L)) {
          lab <- fit$labels
          lab[i] <- cc
          best <- max(best, q_r(lab, d, lam, al) - q0)
        }
      }
      cat(sprintf(
        paste0(
          "check 1b: toy3 lambda %g alpha %s q %.6f (kernel %.6f, ",
          "q(0) %.6f) best single move %+.1e, iterations %d\n"
        ), lam, paste(a, collapse = "/"), q0, fit$quality,
        q_r(fit$labels, d, lam, 0), best, fit$params$iterations
      ))
      stopifnot(
        best < 1e-12, abs(q0 - fit$quality) < 1e-10,
        fit$params$converged
      )
    }
  }

  # ---- 2. Connected modules on jm_fixture (plus planted recovery) ----
  truth <- function(s) helper$true_modules(jm_nets[[s]], f$mods)$modules
  for (a in list(0, SCHED, 1)) {
    for (sd in 1:3) {
      fit <- hl$hleiden(lambda = 4, alpha = a, seed = sd, data = dj)
      connected(fit, dj, "jm_fixture")
      ar <- vapply(c("A", "B"), function(s) {
        tr <- truth(s)
        ari(tr, fit$modules[[s]]$modules[names(tr)])
      }, 0)
      h1 <- fit$hogs[fit$hogs$hog == "H1", ]
      cat(sprintf(
        paste0(
          "check 2: jm_fixture lambda 4 alpha %s seed %d connected; ",
          "ARI vs planted A %.3f B %.3f; B0041 defects %s; q %.4f\n"
        ), paste(a, collapse = "/"), sd, ar[1], ar[2],
        fit$modules$B$modules[["B0041"]] != h1$module_main, fit$quality
      ))
    }
  }
  cat("checks 1-2 passed\n")
}

# ---- 3-4. Planted sweep: hleiden at alpha 0 and the schedule ----
pe <- new.env()
for (ex in parse("dev/probes/joint_planted.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
    is.name(ex[[2]]) && as.character(ex[[2]]) %in% c(
    "S", "N", "K", "SP", "LAMBDA", "DEG", "DELTA", "SEEDS", "prof",
    "ami", "main_of", "scores", "planted", "joint_lab"
  )) {
    eval(ex, pe)
  }
}
TSV <- file.path(OUT, "planted_hleiden.tsv")
ENG <- list(hleiden_a0 = 0, hleiden_sched = SCHED)
# `a0`: alpha 0 only (the schedule arm failed check 2, see the log)
if ("a0" %in% args) ENG <- ENG["hleiden_a0"]

summarize <- function() {
  x <- utils::read.delim(TSV)
  e2 <- utils::read.delim(file.path(OUT, "planted.tsv"))
  e2 <- e2[e2$engine %in% c("E2", "two_step_match"), ]
  e2$engine[e2$engine == "two_step_match"] <- "match"
  m <- c(
    "ami", "jac_median", "jac_frac50", "def_precision", "def_recall",
    "n_modules", "secs"
  )
  x <- rbind(
    x[c("d", "delta", "engine", "lambda", "seed", m)],
    e2[c("d", "delta", "engine", "lambda", "seed", m)]
  )
  key <- x[, c("d", "delta", "engine", "lambda")]
  sm <- stats::aggregate(x[, m], key, mean)
  sd <- stats::aggregate(x[, c("ami", "def_recall")], key, stats::sd)
  names(sd)[5:6] <- c("ami_sd", "def_recall_sd")
  sm <- merge(sm, sd)
  best <- do.call(rbind, lapply(
    split(sm, sm[, c("d", "delta", "engine")], drop = TRUE),
    function(b) b[which.max(b$ami), ]
  ))
  best <- best[order(best$d, best$delta, best$engine), ]
  utils::write.table(best, file.path(OUT, "planted_hleiden_best.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  best
}
if ("summary" %in% args) {
  print(summarize(), digits = 3, row.names = FALSE)
  quit(save = "no")
}
if (!"planted" %in% args) quit(save = "no")

rows <- list()
t_all <- Sys.time()
for (dn in names(pe$DEG)) {
  for (delta in pe$DELTA) {
    for (seed in pe$SEEDS) {
      pl <- pe$planted(seed, pe$DEG[[dn]][1], pe$DEG[[dn]][2], delta)
      dat <- hl$jl_data(pl$nets, pl$ortho[, c("species", "gene", "hog")])
      for (en in names(ENG)) {
        for (lam in pe$LAMBDA) {
          fit <- hl$hleiden(
            lambda = lam, alpha = ENG[[en]], seed = seed,
            data = dat
          )
          connected(fit, dat, sprintf(
            "%s %s delta %g seed %d lambda %g",
            en, dn, delta, seed, lam
          ))
          sc <- pe$scores(
            pe$joint_lab(lapply(fit$modules, `[[`, "modules")), pl
          )
          rows[[length(rows) + 1L]] <- data.frame(
            engine = en, lambda = lam, delta = delta, d = dn, seed = seed,
            t(sc), secs = fit$params$secs,
            iterations = fit$params$iterations,
            converged = fit$params$converged, stalls = fit$params$stalls,
            disconnected = 0L
          )
        }
      }
      tab <- do.call(rbind, rows)
      utils::write.table(tab, TSV,
        sep = "\t", quote = FALSE,
        row.names = FALSE
      )
      cat(sprintf(
        "%s delta %g seed %d done, %d rows, %.1f min\n", dn, delta, seed,
        nrow(tab), as.numeric(difftime(Sys.time(), t_all, units = "mins"))
      ))
    }
  }
}
print(summarize(), digits = 3, row.names = FALSE)
d60 <- tab[tab$d == "d60", ]
cat(sprintf(
  "check 4: d60 wall time per run, median %.2f s, max %.2f s\n",
  stats::median(d60$secs), max(d60$secs)
))
cat("planted done\n")
