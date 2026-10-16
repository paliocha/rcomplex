# One-off: the trait pass on the null edges of a finished run, written
# as runs/<run>/trait_pass_null.tsv (run.R now does this itself).
# Usage: Rscript trait_null.R <exp> <run>
suppressPackageStartupMessages(library(rcomplex))
args <- commandArgs(TRUE)
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
E <- file.path(ROOT, paste0("complex_", args[1]))
OUT <- file.path(E, "runs", args[2])
info <- read.delim(file.path(E, "inputs", "species.tsv"))
info <- info[info$run == args[2], ]
groups <- split(info$species, info$trait)
res <- readRDS(file.path(OUT, "res.rds"))
n_ts <- function(x) {
  length(unique(x$hog[x$hog_class == "trait_specific"]))
}
tp0 <- rcomplex:::.driver_cliques(res$edges_null, info$species,
  clades = groups, alpha = 0.1, min_power = 0.8
)$classification
data.table::fwrite(tp0, file.path(OUT, "trait_pass_null.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)
tp <- data.table::fread(file.path(OUT, "trait_pass.tsv"), data.table = FALSE)
cat("trait_specific hogs: real", n_ts(tp), "null", n_ts(tp0), "\n")
print(table(tp0$hog_class[!duplicated(tp0$hog)]))
