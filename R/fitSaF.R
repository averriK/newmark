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
  ## 1) parse user input for uncertainty
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
  
  ## 2) normalise model weights
  wTable <- data.table::data.table(ID=models, weight=score/sum(score))
  
  ## 3) pre-process UHS
  UHS <- data.table::copy(uhs)
  UHS[Tn == 0, Tn := 0.01, by=.(p)]  # if Tn=0 => set Tn=0.01
  # NOTE: your site factor might rely on Tn=0.01 as "PGA" or Tn=0 => unify your approach
  
  # We'll define a small helper to get hazard draws if "Sa" is in uncertainty:
  getHazardDraws <- function(UHS, Td, n=1) {
    if (unc_mode %in% c("sa","both")) {
      out <- sampleSa(UHS, Td=Td, n=n)$Sa
    } else {
      # no fractile => replicate mean
      # find the "mean" row at Tn=Td
      meanRow <- UHS[Tn==Td & p=="mean"]
      if (nrow(meanRow)!=1L)
        stop("Need exactly one 'mean' row at Tn=",Td," for hazard = mean.")
      val <- meanRow$Sa
      out <- rep(val, times=n)
    }
    out
  }
  
  ## 4) We'll build up a big table of random draws for each Tn>0.
  # Each Tn => we get hazard draws for that Tn => call them SaRock(Tn).
  # We also get hazard draws for Tn=0.01 => pga. Then we combine with site factor if vs30!=vref
  
  bigTable <- data.table()
  add <- function(dt, new) rbindlist(list(dt,new), use.names=TRUE, fill=TRUE)
  
  # gather unique Tn>0.01 from the input
  TnList <- sort(unique(UHS[Tn>0.009, Tn]))  # e.g. skip Tn=0.01 if you treat that as PGA
  
  # sample PGA draws once => we'll reuse them for each Tn
  # or you can do it Tn=0.01 each time, but typically one set of PGA draws is enough
  pgaDraws <- getHazardDraws(UHS, Td=0.01, n=NS)  # length=NS
  
  for (thisTn in TnList) {
    # hazard draw for SaRock at Tn=thisTn
    SaDraws <- getHazardDraws(UHS, Td=thisTn, n=NS)  # length=NS
    
    # site factor logic
    # if vs30==vref => factor=1 => final SaF= SaDraws
    # else => we do a row-by-row aggregator with site-factor model(s)
    DnTable <- data.table()  # we'll store all random draws from each site-factor model
    if (vs30==vref) {
      # no amplification
      # produce one "model" => ID="NoAmpl"
      outDT <- data.table(
        ID="Rock",
        sample=1:NS,
        SaF=SaDraws  # direct
      )
      DnTable <- add(DnTable, outDT)
    } else {
      # vs30!=vref => we call each model. E.g. "ST17"
      for (mod in models) {
        # 1) get muLnF, sdLnF => depends on pgaDraws, thisTn, vs30
        # but we likely want row-by-row approach => each of the NS hazard draws => unique pga, Sa
        # => we do something akin to a "getFUncertainty()" aggregator
        # We'll define a mini table for this Tn
        rowDT <- data.table(
          SaRock=SaDraws,
          pga   = pgaDraws,
          TnVal = thisTn
        )
        # for each row, call e.g. SaF_ST17(Sa=SaRock, pga=pga, Tn=TnVal,...)
        # then if "F" in unc_mode => random draws => else deterministic
        
        rowDT2 <- rowDT[, {
          tmp <- SaF_ST17(
            Sa=SaRock, pga=pga, Tn=TnVal,
            vs30=vs30, vref=vref
          ) 
          # returns data.table with muLnSaF, sdLnSaF, ID
          
          if (unc_mode %in% c("f","both")) {
            # random draw
            LnSaF <- rnorm(1, tmp$muLnSaF, tmp$sdLnSaF)
            SaFval<- exp(LnSaF)
          } else {
            # deterministic
            SaFval<- exp(tmp$muLnSaF)
          }
          .(SaF=SaFval)
        }, by=.(sample=.I)]  # sample index => 1..NS
        rowDT2[, ID:=mod]
        
        DnTable <- add(DnTable, rowDT2)
      } # end for each model
    }
    
    # now DnTable has columns: sample, ID, SaF
    # aggregator => weighted quantiles => produce p in (p from UHS) => etc.
    # same as fitDn: define "p" from UHS if p != "mean"
    # We'll do a line like:
    wDT <- wTable[DnTable, on="ID"]
    
    # define the hazard fractiles from UHS => exclude p=="mean"
    # typically the same p used for Tn=thisTn
    # or we do all p in UHS[p!="mean", unique(as.numeric(p))] => same approach
    # for simplicity:
    probs <- UHS[p!="mean", unique(as.numeric(p))]
    
    # Weighted aggregator => Hmisc::wtd.quantile
    SaF_q <- wDT[, .(
      p  = probs,
      SaF=Hmisc::wtd.quantile(
        x=SaF,
        weights=weight,
        probs=probs,
        type="quantile",
        na.rm=TRUE
      )
    )]
    # mean row
    SaF_mean <- data.table(p="mean", SaF=mean(wDT$SaF))
    
    outQ <- rbindlist(list(SaF_q, SaF_mean))
    outQ[, Tn:=thisTn]
    
    # store
    bigTable <- add(bigTable, outQ)
  } # end for each Tn
  
  # reorder columns
  data.table::setcolorder(bigTable, c("Tn","p","SaF"))
  bigTable[]
}
