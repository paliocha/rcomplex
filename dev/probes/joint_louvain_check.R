# Probe E2 kernel acceptance (spec checks 1-4).
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-e2 Rscript dev/probes/joint_louvain_check.R
#   ... joint_louvain_check.R wood   (adds check 4 on the full wood nets)
suppressPackageStartupMessages(library(rcomplex))
REPO <- "/Users/martin/Documents/R/rcomplex"
OUT <- file.path(REPO, "prepare_data/cache/probe_e2")
jl <- new.env(parent = asNamespace("rcomplex"))
jl$JL_DIR <- "dev/probes"
sys.source("dev/probes/joint_louvain.R", jl)
ari <- mclust::adjustedRandIndex

# ---- 1. Toy: C++ quality equals plain R; Louvain is a local optimum ----
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
d <- jl$jl_data(toy, toy_o)
diffs <- vapply(1:20, function(r) {
  lab <- sample.int(4, d$n, replace = TRUE)
  cq <- do.call(jl$jl_quality_cpp, c(
    list(lab = lab - 1L), jl$.jl_args(d), list(lambda = 2.5, gamma = 1.3)
  ))
  rq <- jl$joint_quality(lab, d, lambda = 2.5, gamma = 1.3)
  max(abs(c(
    cq$quality - rq$quality, cq$species - rq$species,
    cq$hog_raw - rq$hog_raw
  )))
}, 0)
cat(sprintf("check 1: max |C++ - R| over 20 labelings = %.2e\n", max(diffs)))
stopifnot(max(diffs) < 1e-12)

# A single-gene move from the Louvain result must not raise Q
for (lam in c(0, 1, 5)) {
  fit <- jl$joint_louvain(lambda = lam, seed = 1, data = d)
  q0 <- jl$joint_quality(fit, d)$quality
  best <- -Inf
  for (i in seq_len(d$n)) {
    for (c in c(unique(fit$labels), max(fit$labels) + 1L)) {
      lab <- fit$labels
      lab[i] <- c
      best <- max(best, jl$joint_quality(lab, d, lam)$quality - q0)
    }
  }
  cat(sprintf(
    "check 1b: lambda %g Q %.6f (C++ %.6f) best single move %+.2e\n",
    lam, q0, fit$quality, best
  ))
  stopifnot(best < 1e-12, abs(q0 - fit$quality) < 1e-12)
}

# ---- jm_fixture() from the test file, pres_expr() from the helper ----
helper <- new.env(parent = asNamespace("rcomplex"))
sys.source("tests/testthat/helper-preservation.R", helper)
for (ex in parse("tests/testthat/test-joint-modules.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
    as.character(ex[[2]]) %in% c("jm_fixture", "long_ortho")) {
    eval(ex, helper)
  }
}
f <- helper$jm_fixture()
long <- helper$long_ortho(f$ortho)
nets <- list(A = f$netA, B = f$netB)
d <- jl$jl_data(nets, long)
truth <- function(s) {
  helper$true_modules(nets[[s]], f$mods)$modules
}

# ---- 2. lambda = 0: per-species modularity, against igraph ----
# All-gene ARI is reported but not gated: the background of jm_fixture()
# is noise, and on B igraph's own Louvain and Leiden agree at ARI 0.4
# there. The gate is the planted genes. Q_s (median over seeds 1-5) must
# be within 1 % of detect_modules() (Leiden); igraph Louvain is reported.
for (s in c("A", "B")) {
  dm <- detect_modules(nets[[s]], seed = 1)
  g <- dm$graph
  tr <- truth(s)
  r <- t(vapply(1:5, function(sd) {
    fit0 <- jl$joint_louvain(lambda = 0, seed = sd, data = d)
    ours <- fit0$modules[[s]]$modules[igraph::V(g)$name]
    set.seed(sd)
    lv <- igraph::membership(igraph::cluster_louvain(g))
    c(
      q = fit0$quality_parts[[s]], q_ig = igraph::modularity(g, ours),
      q_lv = igraph::modularity(g, lv), q_dm = dm$modularity,
      ari_dm = ari(ours, dm$modules), ari_lv = ari(ours, lv),
      ari_pl = ari(tr, ours[names(tr)]),
      ari_pl_dm = ari(tr, dm$modules[names(tr)])
    )
  }, numeric(8)))
  m <- apply(r, 2, stats::median)
  cat(sprintf(
    paste0(
      "check 2: %s (median of 5 seeds) ARI all genes vs detect_modules ",
      "%.3f, vs cluster_louvain %.3f; planted genes vs planted %.3f ",
      "(detect_modules %.3f); Q_s %.4f, igraph Louvain %.4f (rel %+.4f), ",
      "Leiden %.4f (rel %+.4f); max |Q_s - igraph modularity()| %.1e\n"
    ), s, m[["ari_dm"]], m[["ari_lv"]], m[["ari_pl"]], m[["ari_pl_dm"]],
    m[["q"]], m[["q_lv"]], m[["q"]] / m[["q_lv"]] - 1, m[["q_dm"]],
    m[["q"]] / m[["q_dm"]] - 1, max(abs(r[, "q"] - r[, "q_ig"]))
  ))
  stopifnot(
    m[["ari_pl"]] > 0.9, abs(m[["q"]] / m[["q_dm"]] - 1) < 0.01,
    max(abs(r[, "q"] - r[, "q_ig"])) < 1e-10
  )
}

# ---- 3. lambda = 4: planted recovery, shared labels, a defector ----
long2 <- long
long2$hog[long2$species == "B" & long2$gene == "B0041"] <- "H1"
d2 <- jl$jl_data(nets, long2)
fit4 <- jl$joint_louvain(lambda = 4, seed = 1, data = d2)
for (s in c("A", "B")) {
  tr <- truth(s)
  a <- ari(tr, fit4$modules[[s]]$modules[names(tr)])
  cat(sprintf("check 3: %s ARI vs planted %.3f\n", s, a))
  stopifnot(a > 0.9)
}
hg <- fit4$hogs
planted <- hg$hog %in% paste0("H", 2:160)
h1 <- hg[hg$hog == "H1", ]
cat(sprintf(
  paste0(
    "check 3: hog purity planted %.3f, all %.3f; H1 n_modules %d, ",
    "B0041 label %d, module_main %d\n"
  ), mean(hg$n_modules[planted] == 1), mean(hg$n_modules == 1),
  h1$n_modules, fit4$modules$B$modules[["B0041"]], h1$module_main
))
stopifnot(
  mean(hg$n_modules[planted] == 1) > 0.9,
  fit4$modules$B$modules[["B0041"]] != h1$module_main
)

# ---- 4. Wood, full networks of the 0.4.1 run, one run at lambda = 1 ----
if (identical(commandArgs(TRUE)[1], "wood")) {
  w <- readRDS(file.path(REPO, "prepare_data/cache/snr/data.rds"))
  t0 <- Sys.time()
  dw <- jl$jl_data(w$nets, w$ortho[, c("species", "gene", "hog")])
  t1 <- Sys.time()
  fw <- jl$joint_louvain(lambda = 1, seed = 1, data = dw)
  t2 <- Sys.time()
  cat(sprintf(
    paste0(
      "check 4: wood n %d, co entries %d, hog entries %d; build %.0f s, ",
      "Louvain %.0f s (%d levels); modules %d; Q %.4f; hog purity %.3f\n"
    ), dw$n, length(dw$co@x), length(dw$hog@x),
    as.numeric(difftime(t1, t0, units = "secs")),
    as.numeric(difftime(t2, t1, units = "secs")), fw$params$levels,
    length(unique(fw$labels)), fw$quality,
    mean(fw$hogs$n_modules == 1)
  ))
  print(fw$quality_parts)
}
cat("all checks passed\n")
