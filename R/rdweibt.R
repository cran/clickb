#' Truncated Weibull distribution to generate sequence length
#'
#' @param n is the number of observations in one sequence
#' @param T.range is a two element vector for minimum and maximum length
#' @param shapewb is the Discrete Weibull shape parameter
#' @param scalewb is the Discrete Weibull scale parameter
#' @importFrom stats runif
#' @keywords internal
rdweibt <- function(n, T.range, shapewb, scalewb) {

  F.a <- DiscreteWeibull::pdweibull(min(T.range), q = shapewb, beta = scalewb)
  F.b <- DiscreteWeibull::pdweibull(max(T.range), q = shapewb, beta = scalewb)

  u <- runif(n, min = F.a, max = F.b)
  DiscreteWeibull::qdweibull(u, q = shapewb, beta = scalewb)
}
