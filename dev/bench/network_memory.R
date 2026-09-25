# Time and peak RSS of compute_network(): dense versus blockwise.
# Usage (from the repo root):
#   Rscript dev/bench/network_memory.R [n_values=5000,10000,20000]
#     [block_sizes=256,1024] [out=dev/bench/network_memory.tsv] [lib=<path>]
# Without lib=, each child loads rcomplex with devtools::load_all(".").
# Each configuration runs in its own Rscript child under /usr/bin/time, so
# peak RSS is that of one compute_network() call plus R and package startup.

args <- c(
  n_values = "5000,10000,20000", block_sizes = "256,1024",
  out = "dev/bench/network_memory.tsv", lib = ""
)
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (!kv[1] %in% names(args)) stop("unknown argument: ", kv[1])
  args[[kv[1]]] <- paste(kv[-1], collapse = "=")
}
n_values <- as.integer(strsplit(args[["n_values"]], ",")[[1]])
block_sizes <- as.integer(strsplit(args[["block_sizes"]], ",")[[1]])
linux <- Sys.info()[["sysname"]] == "Linux"

se_path <- "prepare_data/data/BDIS_se.rds"
if (file.exists(se_path)) {
  se <- readRDS(se_path)
  pool <- SummarizedExperiment::assay(se)[
    , SummarizedExperiment::colData(se)$tissue == "leaf"
  ]
  v <- apply(pool, 1, stats::var)
  pool <- pool[v > 0, ][order(v[v > 0], decreasing = TRUE), ]
  cat(
    "data: BDIS leaf,", nrow(pool), "non-constant genes x",
    ncol(pool), "samples\n"
  )
} else {
  pool <- NULL
  cat("data: ", se_path, " absent, using seeded rnorm (100 samples)\n",
    sep = ""
  )
}

make_x <- function(n, pool) {
  set.seed(1L)
  if (is.null(pool)) {
    x <- matrix(stats::rnorm(n * 100), n)
  } else {
    x <- pool[seq_len(min(n, nrow(pool))), ]
    if (n > nrow(x)) {
      cat(
        "n =", n, "exceeds", nrow(x), "genes: padding",
        n - nrow(x), "seeded rnorm rows\n"
      )
      x <- rbind(x, matrix(
        stats::rnorm((n - nrow(x)) * ncol(x)),
        n - nrow(x)
      ))
    }
  }
  rownames(x) <- sprintf("g%06d", seq_len(n))
  x
}

loader <- if (nzchar(args[["lib"]])) {
  sprintf("library(rcomplex, lib.loc = %s)", deparse(args[["lib"]]))
} else {
  sprintf("devtools::load_all(%s, quiet = TRUE)", deparse(normalizePath(".")))
}

run_child <- function(x_file, block) {
  res_file <- tempfile(fileext = ".rds")
  extra <- if (is.na(block)) "" else sprintf(", block_size = %dL", block)
  code <- paste0(
    "suppressMessages(", loader, "); x <- readRDS('", x_file, "'); ",
    "t0 <- proc.time()[['elapsed']]; ",
    "r <- tryCatch(compute_network(x, density = 0.03, ",
    "store_density = 0.05, sparse = TRUE, n_cores = 8L", extra, "), ",
    "error = function(e) conditionMessage(e)); ",
    "s <- proc.time()[['elapsed']] - t0; ",
    "saveRDS(if (is.character(r)) list(err = r) else list(seconds = s, ",
    "threshold = r$threshold, store_threshold = r$store_threshold, ",
    "nnz = length(r$network@x)), '", res_file, "')"
  )
  out <- suppressWarnings(system2(
    "/usr/bin/time", c(
      if (linux) "-v" else "-l", "Rscript", "-e",
      shQuote(code)
    ),
    stdout = TRUE, stderr = TRUE
  ))
  rss_line <- grep("maximum resident set size", out,
    ignore.case = TRUE,
    value = TRUE
  )
  rss <- as.numeric(gsub("[^0-9]", "", rss_line))
  rss_mb <- if (length(rss) == 1) {
    rss * (if (linux) 1024 else 1) / 2^20
  } else {
    NA_real_
  }
  if (!file.exists(res_file)) {
    stop("child failed:\n", paste(out, collapse = "\n"))
  }
  r <- readRDS(res_file)
  if (!is.null(r$err)) {
    if (!grepl("unused argument", r$err)) stop("child error: ", r$err)
    cat("  blockwise not supported yet (", r$err, "): recording NA\n")
    r <- list(
      seconds = NA_real_, threshold = NA_real_,
      store_threshold = NA_real_, nnz = NA_integer_
    )
    rss_mb <- NA_real_
  }
  c(r, peak_rss_mb = rss_mb)
}

rows <- list()
for (n in n_values) {
  x_file <- tempfile(fileext = ".rds")
  saveRDS(make_x(n, pool), x_file)
  for (b in c(NA_integer_, block_sizes)) {
    r <- run_child(x_file, b)
    row <- data.frame(
      n = n, mode = if (is.na(b)) "dense" else "blockwise",
      block_size = b, seconds = r$seconds,
      peak_rss_mb = r$peak_rss_mb, threshold = r$threshold,
      store_threshold = r$store_threshold, nnz = r$nnz
    )
    msg <- sprintf(
      "n=%d %s block=%s: %.1f s, %.0f MB peak, nnz=%s",
      n, row$mode, b, row$seconds, row$peak_rss_mb, row$nnz
    )
    if (!is.na(b)) {
      d <- rows[[paste(n, NA)]]
      same <- identical(
        c(d$threshold, d$store_threshold, d$nnz),
        c(row$threshold, row$store_threshold, row$nnz)
      )
      msg <- paste0(
        msg, ", equals dense: ",
        if (is.na(row$nnz)) NA else same
      )
    }
    cat(msg, "\n")
    rows[[paste(n, b)]] <- row
  }
  unlink(x_file)
}
utils::write.table(do.call(rbind, rows), args[["out"]],
  sep = "\t",
  quote = FALSE, row.names = FALSE
)
cat("wrote", args[["out"]], "\n")
