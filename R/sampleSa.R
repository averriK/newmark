#' @title Sample spectral acceleration with optional piecewise approach,
#'        enforcing monotonic ln(Sa).
#'
#' @param UHS data.table with Tn, Sa, p (prob labels + "mean").
#' @param Td numeric scalar – target period (s).
#' @param n integer – number of samples (default 1).
#' @param piecewise logical; if TRUE, do piecewise interpolation in (p, ln(Sa)).
#' @param forceMonotonic logical; if TRUE, clamp ln(Sa) so it's non-decreasing.
#' @param pLimits length-2 numeric, e.g. c(0,1) or c(0.1,0.9); probability range to keep.
#' @param checkAscending logical; if TRUE, we still check for descending ln(Sa) but fix it anyway if forceMonotonic=TRUE.
#' @return list(SaTable, muLnSa, sdLnSa, Sa) same as before
#' @export
sampleSa <- function(
    UHS,
    Td,
    n = 1,
    piecewise = TRUE,
    forceMonotonic = TRUE,
    pLimits = c(0,1),
    checkAscending = TRUE
) {
  ## 1) Interpolate Sa if Td not in UHS$Tn
  if (Td %in% UHS$Tn) {
    SaTable <- UHS[Tn == Td, .(Sa, p)]
  } else {
    SaTable <- UHS[
      , .(Sa = stats::approx(
        x = log(Tn),
        y = log(Sa),
        xout = log(Td),
        rule = 2
      )$y |> exp()),
      by = .(p)
    ]
  }
  
  ## 2) Extract mean row => muLnSa
  meanRow <- SaTable[p == "mean"]
  if (nrow(meanRow) != 1L) {
    stop("Expected exactly one 'mean' row for Tn=", Td,
         ". Got ", nrow(meanRow), " rows.")
  }
  muLnSa <- log(meanRow$Sa)
  
  ## Probability rows
  probRows <- SaTable[p != "mean"]
  if (nrow(probRows) < 1) {
    stop("No probability rows at Tn=", Td, " => cannot sampleSa().")
  }
  
  # Convert p to numeric
  probRows[, pnum := as.numeric(p)]
  probRows <- probRows[!is.na(pnum) & pnum > 0 & pnum < 1]
  
  data.table::setorder(probRows, pnum)
  probRows[, lnSa := log(Sa)]
  
  if (!piecewise) {
    # -------------------------------------------
    # A) Original aggregator approach
    # -------------------------------------------
    p_vec <- probRows$pnum
    q_vec <- probRows$lnSa
    
    sdLnSa <- sdQ(meanValue = muLnSa, p = p_vec, q = q_vec)
    Sa_draws <- rnormQ(meanValue = muLnSa, p = p_vec, q = q_vec, n = n) |> exp()
    
    return(list(
      SaTable = SaTable,
      muLnSa  = muLnSa,
      sdLnSa  = sdLnSa,
      Sa      = Sa_draws
    ))
    
  } else {
    # -------------------------------------------
    # B) Piecewise approach in (p, ln(Sa)) space
    # -------------------------------------------
    # 1) Clip p to pLimits
    lowerP <- pLimits[1]
    upperP <- pLimits[2]
    if (lowerP < 0 || upperP > 1 || lowerP >= upperP) {
      stop("Invalid pLimits: must be c(a,b) with 0<=a<b<=1.")
    }
    probRows <- probRows[pnum >= lowerP & pnum <= upperP]
    if (nrow(probRows) < 2) {
      stop("After pLimits=[", lowerP, ",", upperP,
           "], fewer than 2 p-rows remain => can't piecewise sample.")
    }
    
    # Optionally check for ascending lnSa
    lnVec <- probRows$lnSa
    if (checkAscending) {
      # If there's a descending step
      descIdx <- which(diff(lnVec) < 0)
      if (length(descIdx) > 0) {
        msg <- paste("lnSa is not strictly ascending at p in rows:",
                     paste(descIdx, collapse=","), 
                     "=> forcing monotonic if forceMonotonic=TRUE.")
        if (!forceMonotonic) {
          warning(msg)
        } else {
          message(msg)
        }
      }
    }
    
    # 2) Enforce monotonic if desired
    if (forceMonotonic && nrow(probRows) >= 2) {
      for (i in 2:nrow(probRows)) {
        if (probRows$lnSa[i] < probRows$lnSa[i-1]) {
          probRows$lnSa[i] <- probRows$lnSa[i-1]  # clamp up
        }
      }
    }
    
    # piecewise sampling
    p_vec  <- probRows$pnum
    ln_vec <- probRows$lnSa
    pMin <- p_vec[1]
    pMax <- p_vec[length(p_vec)]
    
    piecewiseSample <- function(nDraws) {
      U <- runif(nDraws)
      out <- numeric(nDraws)
      
      for (i in seq_len(nDraws)) {
        u <- U[i]
        if (u <= pMin) {
          out[i] <- ln_vec[1]
        } else if (u >= pMax) {
          out[i] <- ln_vec[length(ln_vec)]
        } else {
          j <- findInterval(u, p_vec)
          p1 <- p_vec[j]
          p2 <- p_vec[j+1]
          s1 <- ln_vec[j]
          s2 <- ln_vec[j+1]
          alpha <- (u - p1)/(p2 - p1)
          out[i] <- s1 + alpha*(s2 - s1)
        }
      }
      exp(out)
    }
    
    Sa_draws <- piecewiseSample(n)
    
    # Return sdLnSa=NA in piecewise mode
    return(list(
      SaTable = SaTable,
      muLnSa  = muLnSa,
      sdLnSa  = NA_real_,
      Sa      = Sa_draws
    ))
  }
}
