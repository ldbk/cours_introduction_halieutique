#' Explore r-k combinations compatible with a catch time series (Catch-MSY)
#'
#' Simplified implementation of the Catch-MSY method (Martell & Froese, 2013).
#' Pairs of intrinsic growth rate \code{r} and carrying capacity \code{k} are
#' drawn from log-uniform priors. For each pair, biomass is projected with a
#' deterministic Schaefer surplus-production model driven by the observed
#' catches. A pair is kept if the stock never collapses (B > 0), never exceeds
#' \code{k}, and ends within the assumed final depletion range. The retained
#' pairs give a distribution of MSY = r * k / 4.
#'
#' @param ct Numeric vector of annual catches (in chronological order).
#' @param r_range Lower and upper bounds for \code{r}, e.g. \code{c(0.2, 1)}
#'   (can be set from resilience).
#' @param k_range Lower and upper bounds for \code{k}, in the same unit as
#'   \code{ct}, e.g. \code{c(max(ct), 50 * max(ct))}.
#' @param n Number of random r-k draws.
#' @param startbio Range of relative biomass (B/k) at the start of the series.
#' @param finalbio Range of relative biomass (B/k) after the last catch.
#'
#' @return A data frame of the accepted combinations, with columns \code{r},
#'   \code{k} and \code{msy}.
#'
#' @references Martell, S. and Froese, R. (2013). A simple method for
#'   estimating MSY from catch and resilience. \emph{Fish and Fisheries},
#'   14(4), 504-514. \doi{10.1111/j.1467-2979.2012.00485.x}
#'
#' @examples
#' ct  <- c(10, 12, 15, 20, 25, 30, 28, 25, 22, 20, 18, 15)
#' out <- catch_msy(ct, r_range = c(0.2, 1), k_range = c(max(ct), 50 * max(ct)))
#' exp(mean(log(out$msy)))        # MSY (geometric mean)
#' plot(log(out$r), log(out$k))   # accepted r-k pairs
#'
#' @export
catch_msy <- function(ct, r_range, k_range, n = 10000,
                      startbio = c(0.5, 0.9), finalbio = c(0.01, 0.4)) {
  r  <- exp(runif(n, log(r_range[1]), log(r_range[2])))  # log-uniform draws
  k  <- exp(runif(n, log(k_range[1]), log(k_range[2])))
  B  <- runif(n, startbio[1], startbio[2]) * k           # initial biomass
  ok <- rep(TRUE, n)

  for (t in seq_along(ct)) {                             # Schaefer model
    B  <- B + r * B * (1 - B / k) - ct[t]
    ok <- ok & B > 0 & B <= k
  }
  ok <- ok & B / k >= finalbio[1] & B / k <= finalbio[2]

  res <- data.frame(r = r[ok], k = k[ok], msy = r[ok] * k[ok] / 4)
  if (nrow(res) < 10) warning("Few valid combinations (", nrow(res), ")")
  res
}
