# nolint start
#' Standard deviation from partial quantiles
#'
#' @param meanValue numeric scalar, the known mean of X
#' @param p,q numeric vectors of equal length: probabilities 0<p<1 and their quantiles
#' @param nuStart starting df for Student-t (method 4)
#' @param deltaD acceptance band for Step 7 aggregation
#' @return numeric scalar sd estimate
#' @export sdQ
#' 
sdQ <- function(meanValue, p, q,nuStart = 10, deltaD = 0.01) {
  checkInputs(meanValue, p, q)
  
  OUT <- fitAllMethodsQ(
    meanValue=meanValue, 
    p=p, 
    q=q,
    nuStart = nuStart
  )
  
  agg <- aggregateSigmaQ(OUT, deltaD)
  sdValue <-  agg$sdValue
  if (is.na(sdValue)) {
    sdValue <- OUT$sigma[OUT$method == "6"]
  }
  return(sdValue)
}
