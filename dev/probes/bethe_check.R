# Probe E3 acceptance: species-blocked Bethe Hessian (dev/probes/bethe.R).
# Check 1: K (negative eigenvalues of H(r)) on the planted fixtures of
# dev/probes/joint_planted.R at delta 0, d 60 and 10, seeds 1-3, both
# variants, over the lambda grid. The truth has 12 joint modules plus one
# background per species: 18 classes (13 if the backgrounds join).
# Labels and scores where K <= 160 (k-means is pointless above).
#
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/bethe_check.R
suppressPackageStartupMessages(library(rcomplex))
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")

# joint_planted.R: definitions only (stop before its calibration run)
pe <- new.env()
for (e in parse("dev/probes/joint_planted.R")) {
  if (is.call(e) && identical(e[[1]], as.name("if"))) next
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("calib"))) {
    break
  }
  eval(e, pe)
}
jl <- pe$jl
sys.source("dev/probes/bethe.R", jl)
K_TRUE <- c(classes = pe$K + pe$S, joint = pe$K + 1L)

# `fine`: K on a fine lambda grid, d 60, delta 0, seed 1 (k up to 200)
if ("fine" %in% commandArgs(TRUE)) {
  pl <- pe$planted(1L, pe$DEG$d60[1], pe$DEG$d60[2], 0)
  dat <- jl$jl_data(pl$nets, pl$ortho[, c("species", "gene", "hog")])
  for (v in c("global", "blocked")) {
    for (lam in c(0, 0.25, 0.5, 1, 1.5, 1.75, 2, 2.5, 3)) {
      bh <- jl$bethe_hessian(dat$co + lam * dat$hog, dat$spi, v)
      ev <- jl$bethe_negative(bh$H, k0 = 200L, k_max = 200L)
      cat(sprintf(
        "fine %s lambda %g: K %d%s\n", v, lam, ev$K,
        if (ev$saturated) "+" else ""
      ))
    }
  }
  quit(save = "no")
}

rows <- list()
for (dn in names(pe$DEG)) {
  for (seed in pe$SEEDS) {
    pl <- pe$planted(seed, pe$DEG[[dn]][1], pe$DEG[[dn]][2], 0)
    dat <- jl$jl_data(pl$nets, pl$ortho[, c("species", "gene", "hog")])
    for (v in c("global", "blocked")) {
      for (lam in pe$LAMBDA) {
        t0 <- Sys.time()
        bh <- jl$bethe_hessian(dat$co + lam * dat$hog, dat$spi, v)
        ev <- jl$bethe_negative(bh$H, k_max = 640L)
        secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        sc <- c(ami = NA, jac_median = NA)
        if (ev$K <= 160L) {
          f <- jl$bethe(lambda = lam, variant = v, seed = seed, data = dat)
          lab <- pe$joint_lab(lapply(f$modules, `[[`, "modules"))
          sc <- pe$scores(lab, pl)[c("ami", "jac_median")]
          secs <- f$params$secs
        }
        rows[[length(rows) + 1L]] <- data.frame(
          d = dn, seed = seed, variant = v, lambda = lam, K = ev$K,
          saturated = ev$saturated, t(sc), secs = secs
        )
        print(utils::tail(rows[[length(rows)]], 1), digits = 3)
      }
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.table(tab, file.path(OUT, "bethe_check1.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
d60 <- tab[tab$d == "d60", ]
ok <- vapply(split(d60, d60$variant), function(x) {
  any(vapply(split(x, x$lambda), function(y) {
    all(abs(y$K - K_TRUE[["classes"]]) <= 1L) ||
      all(abs(y$K - K_TRUE[["joint"]]) <= 1L)
  }, TRUE))
}, TRUE)
cat(sprintf(
  "check 1: true K %d (classes) or %d (joint); pass %s\n",
  K_TRUE[["classes"]], K_TRUE[["joint"]],
  paste(names(ok), ok, collapse = ", ")
))
if (!all(ok)) stop("check 1 failed: K not recovered at d = 60.")
