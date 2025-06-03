#' @title Site-Factor Amplification with Uncertainty
#'
#' @description Similar to \code{fitDn()} but for site-amplified spectral acceleration (SaF).
#'   For each oscillator period Tn in the input table, we sample hazard (if "Sa" in uncertainty),
#'   compute or sample site factor (if "F" in uncertainty), then produce draws of \eqn{SaF = Sa * F}.
#'   Finally, we aggregate them via weighted quantiles in the p dimension.
#'
#' @param uhs data.table with columns (Tn, p, Sa). Must include p="mean" plus
#'   0<p<1 rows if you do hazard sampling. Must also have Tn=0 or Tn=0.01 for PGA
#'   if your site-factor model depends on that.
#' @param vs30 numeric scalar – target site velocity.
#' @param vref numeric – reference velocity (760 by default).
#' @param NS integer – number of Monte Carlo draws (default 30).
#' @param models character – site-factor models, default "ST17".
#' @param score numeric – weighting for the model(s).
#' @param uncertainty character vector – any combination of c("Sa","F"):
#'   * "none" => no hazard fractiles, no site-factor random
#'   * "Sa"   => hazard fractiles only
#'   * "F"    => site-factor random only
#'   * "both" => hazard + site-factor random
#'
#' @return data.table with columns:
#'   \itemize{
#'     \item Tn  => oscillator period
#'     \item p   => probabilities from the original p set (excluding p="mean")
#'     \item SaF => weighted quantiles of final \eqn{Sa * F}
#'   }
#'   plus a "mean" row for p="mean".
#' @import data.table
#' @importFrom stats rnorm
#' @export
fitSaF <- function(
    uhs,
    vs30,
    vref        = 760,
    NS          = 30,
    models      = c("ST17"),
    score       = rep(1, length(models)),
    uncertainty = c("none","Sa","F","both")
) {
  # 1) parse `uncertainty`
  unc_set  <- tolower(uncertainty)
  hasSa <- any(unc_set %in% "sa")
  hasF  <- any(unc_set %in% "f")
  unc_mode <- "none"
  if (hasSa && hasF) {
    unc_mode <- "both"
  } else if (hasSa && !hasF) {
    unc_mode <- "sa"
  } else if (!hasSa && hasF) {
    unc_mode <- "f"
  }
  
  # 2) model weights
  wTable <- data.table::data.table(ID = models, weight = score/sum(score))
  
  # 3) pre-process UHS
  UHS <- data.table::copy(uhs)
  UHS[Tn == 0, Tn := 0.01, by = .(p)]
  
  # small helper: sample hazard draws if "Sa", else replicate "mean"
  getHazardDraws <- function(UHS, Td, n=1) {
    if (unc_mode %in% c("sa","both")) {
      # aggregator
      out <- sampleSa(UHS, Td=Td, n=n)$Sa
    } else {
      # deterministic => replicate mean
      meanRow <- UHS[Tn==Td & p=="mean"]
      if (nrow(meanRow)!=1L)
        stop("Need exactly one 'mean' row at Tn=",Td," for hazard=mean.")
      out <- rep(meanRow$Sa, times=n)
    }
    out
  }
  
  # We'll build a final table
  bigTable <- data.table()
  add <- function(dt, new) data.table::rbindlist(list(dt,new), use.names=TRUE)
  
  # gather unique Tn>0.01
  TnList <- sort(unique(UHS[Tn>0.009, Tn]))
  
  # sample once for PGA
  pgaDraws <- getHazardDraws(UHS, Td=0.01, n=NS)
  
  for (thisTn in TnList) {
    # hazard draws for Tn
    SaDraws <- getHazardDraws(UHS, Td=thisTn, n=NS)
    
    # We'll store all draws from all site-factor models
    # *and* we keep the original SaRock in the same row
    drawDT <- data.table()
    
    if (vs30 == vref) {
      # no amplification => ID="Rock" => SaF=SaDraws
      # store (SaRock, SaF=SaRock)
      tmpDT <- data.table(
        sample=1:NS,
        ID="Rock",
        SaRock=SaDraws,
        SaF   =SaDraws
      )
      drawDT <- add(drawDT, tmpDT)
    } else {
      # vs30 != vref => site-factor modeling
      for (modID in models) {
        # row-by-row aggregator
        rowDT <- data.table(
          SaRock=SaDraws,
          pga   = pgaDraws,
          TnVal = thisTn
        )
        
        rowDT2 <- rowDT[, {
          # call e.g. SaF_ST17
          tmp <- SaF_ST17(Sa=SaRock, pga=pga, Tn=TnVal, vs30=vs30, vref=vref)
          # returns muLnSaF, sdLnSaF, ID
          
          # If "F" in unc_mode => random => else deterministic
          if (unc_mode %in% c("f","both")) {
            LnSaF <- rnorm(1, tmp$muLnSaF, tmp$sdLnSaF)
            SaFval<- exp(LnSaF)
          } else {
            SaFval<- exp(tmp$muLnSaF)
          }
          .(SaRock=SaRock, SaF=SaFval)
        }, by=.(sample=.I)]
        
        rowDT2[, ID := modID]
        
        drawDT <- add(drawDT, rowDT2)
      }
    }
    
    # Now `drawDT` has columns: (sample, ID, SaRock, SaF).
    # We do a weighted aggregator => we want quantiles of both SaRock & SaF.
    
    wDT <- wTable[drawDT, on="ID"]  # merges "weight"
    
    # define p in (UHS)
    probs <- UHS[p!="mean", unique(as.numeric(p))]
    
    # Weighted quantiles for SaRock
    SaRock_q <- wDT[, .(
      p=probs,
      SaRock=Hmisc::wtd.quantile(
        x=SaRock,
        weights=weight,
        probs=probs,
        type="quantile",
        na.rm=TRUE
      )
    )]
    # Weighted quantiles for SaF
    SaF_q <- wDT[, .(
      p=probs,
      SaF=Hmisc::wtd.quantile(
        x=SaF,
        weights=weight,
        probs=probs,
        type="quantile",
        na.rm=TRUE
      )
    )]
    
    # mean row for both
    meanRow_rock <- data.table(p="mean", SaRock=mean(wDT$SaRock))
    meanRow_f    <- data.table(p="mean", SaF   =mean(wDT$SaF))
    
    # combine them: note that SaRock_q and SaF_q have the same length & same p order => we can cbind
    outQ <- cbind(SaRock_q[, .(p, SaRock)], SaF=SaF_q$SaF)
    
    # add the mean row => do it carefully:
    outMean <- data.table(p="mean", SaRock=meanRow_rock$SaRock, SaF=meanRow_f$SaF)
    outQ <- rbindlist(list(outQ, outMean))
    
    # store Tn
    outQ[, Tn := thisTn]
    
    # accumulate
    bigTable <- add(bigTable, outQ)
  }
  
  # reorder
  data.table::setcolorder(bigTable, c("Tn","p","SaRock","SaF"))
  bigTable[]
}

