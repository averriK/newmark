#' Interpolate Sa(p) at arbitrary period using log-log interpolation
#'
#' @param uhs data.table with columns Tn, p, Sa.
#' @param Tn numeric scalar — target period for interpolation.
#'
#' @return data.table with interpolated Tn, p, Sa
#' @export
interpolateSaTable <- function(uhs, Tn) {
  stopifnot(all(c("Tn", "p", "Sa") %in% names(uhs)))
  stopifnot(is.numeric(Tn), length(Tn) == 1L)

  p <- sort(unique(uhs[p != "mean"]$p))
  if (length(p) == 0L) {
    stop("interpolateSaTable(): missing quantiles (p != 'mean')")
  }

  SaTable <- data.table::rbindlist(lapply(p, function(p_i) {
    rows <- uhs[p == p_i]
    if (nrow(rows) < 2L) {
      stop("interpolateSaTable(): need >=2 Tn values for p = ", p_i)
    }

    y <- stats::approx(
      x = log(rows$Tn + 1e-8),
      y = log(rows$Sa),
      xout = log(Tn),
      rule = 2
    )$y

    data.table::data.table(Tn = Tn, p = p_i, Sa = exp(y))
  }))

  mean_row <- uhs[p == "mean"]
  if (nrow(mean_row) >= 2L) {
    y_mean <- stats::approx(
      x = log(mean_row$Tn + 1e-8),
      y = log(mean_row$Sa),
      xout = log(Tn),
      rule = 2
    )$y

    SaTable <- data.table::rbindlist(list(
      SaTable,
      data.table::data.table(Tn = Tn, p = "mean", Sa = exp(y_mean))
    ))
  } else if (nrow(mean_row) == 1L) {
    SaTable <- data.table::rbindlist(list(
      SaTable,
      data.table::data.table(Tn = Tn, p = "mean", Sa = mean_row$Sa)
    ))
  }

  data.table::setorder(SaTable, Tn, p)
  SaTable[]
}
