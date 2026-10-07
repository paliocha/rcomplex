# Pure-R reference implementations for validating C++ output
# Direct port of RComPlEx.Rmd logic

#' Reference MR normalization (ascending ranks, raw mutual rank)
#' Matches RComPlEx.Rmd lines 167-172
reference_mr_raw <- function(net) {
  row_ranks <- t(apply(net, 1, rank))
  mr <- sqrt(row_ranks * t(row_ranks))
  diag(mr) <- 0
  mr
}

#' Reference CLR normalization
#' Matches RComPlEx.Rmd lines 157-165
reference_clr <- function(net) {
  z <- scale(net)
  z[z < 0] <- 0
  clr <- sqrt(t(z)^2 + z^2)
  diag(clr) <- 0
  clr
}

#' Reference density threshold
reference_density_threshold <- function(net, density) {
  vals <- sort(net[upper.tri(net, diag = FALSE)], decreasing = TRUE)
  vals[round(density * length(vals))]
}

#' Reference neighborhood comparison for a single pair
#' Direct port of RComPlEx.Rmd lines 217-276, with the self-excluded urn
#' (D5): the anchor gene is never its own neighbour, leaves the
