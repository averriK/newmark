#' Mixture-like aggregator from multiple scenario means plus merged quantiles
#'
#' @param meanValue numeric vector of means (one per scenario). Averaged.
#' @param p numeric vector of probabilities (strictly between 0 and 1).
#' @param q numeric vector of quantiles (same length as p).
#' @return data.table with columns: mu, sd
#' @export
fitMixtureQ <- function(meanValue, p, q) {
  stopifnot(length(p) == length(q))
  if (any(p <= 0 | p >= 1)) stop("Probabilities must be in (0,1)")

  # 1. Combined mean
  mu <- mean(meanValue)

  # 2. Sort p,q
  df <- data.table(p = p, q = q)
  setorder(df, p)

  # 3. Estimate variance via piecewise trapezoidal integration
  dp <- diff(df$p)
  mid <- 0.5 * (df$q[-.N] + df$q[-1])
  var <- sum(dp * (mid - mu)^2)

  # 4. Tail corrections
  if (df$p[1] > 0)  var <- var + df$p[1] * (df$q[1] - mu)^2
  if (df$p[.N] < 1) var <- var + (1 - df$p[.N]) * (df$q[.N] - mu)^2

  data.table(mu = mu, sd = sqrt(var))
}
