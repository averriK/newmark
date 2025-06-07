# nolint start
# R/utils.R
#' @keywords internal
#' @noRd
fitAllMethodsQ <- function(meanValue, p, q, nuStart = 10) {
    out1 <- .fitMethod1(p, q)
    out2 <- .fitMethod2(p, q)
    out3 <- .fitMethod3(p, q)
    out4 <- .fitMethod4(p, q, nuStart)
    out5 <- .fitMethod5(p, q, meanValue)
    out6 <- .fitMethod6(p, q, meanValue)
    # combine param methods 1..5
    OUT <- list(out1, out2, out3, out4, out5, out6) |> data.table::rbindlist(use.names = TRUE)


    return(OUT)
}

#' @keywords internal
#' @noRd

aggregateSigmaQ <- function(df, deltaD = 0.01) {
    # Extract only param methods 1..5
    core <- df[df$method %in% c("1", "2", "3", "4", "5"), ]

    # If no good data from methods 1..5, fall back to method 6 immediately
    if (nrow(core) == 0 || all(is.na(core$d))) {
        fallback6 <- df[df$method == "6", ][1, ] # method 6 row
        return(list(sdValue = fallback6$sigma, bestMethod = "6"))
    }

    # Normal aggregator logic among 1..5
    dStar <- min(core$d, na.rm = TRUE)
    inSet <- (core$d <= dStar + deltaD)

    chosenSSE <- core$sse[inSet]
    if (all(is.na(chosenSSE))) {
        # e.g. everything is NA => fallback to method 6
        fallback6 <- df[df$method == "6", ][1, ]
        return(list(sdValue = fallback6$sigma, bestMethod = "6"))
    }

    w <- 1 / chosenSSE
    sumW <- sum(w, na.rm = TRUE)
    wNorm <- w / sumW

    chosenSigma <- core$sigma[inSet]
    sdValue <- sum(wNorm * chosenSigma, na.rm = TRUE)

    # if aggregator yields NA/Inf, fallback to method 6
    if (!is.finite(sdValue)) {
        fallback6 <- df[df$method == "6", ][1, ]
        return(list(sdValue = fallback6$sigma, bestMethod = "6"))
    }

    # Otherwise, aggregator picks best method among 1..5
    bestIdx <- which.min(core$d)
    bestMethod <- core$method[bestIdx]

    list(sdValue = sdValue, bestMethod = bestMethod)
}



# -- Method 1: Weighted linear model -------------------------------------------
#' @keywords internal
#' @noRd
.fitMethod1 <- function(p, q) {
    z <- stats::qnorm(p)
    w <- wgt(p)
    fit <- safeLm(q, z, w)
    alpha <- fit$coefficients[1]
    beta <- fit$coefficients[2]

    # Evaluate in p-space
    CDF <- stats::pnorm((q - alpha) / beta)
    dVal <- max(abs(CDF - p))
    sseVal <- sum(w * (CDF - p)^2)

    templateDF(
        method = "1",
        sigma  = beta,
        d      = dVal,
        sse    = sseVal,
        alpha  = alpha,
        beta   = beta
    )
}

# -- Method 2: Central-range ---------------------------------------------------
#' @keywords internal
#' @noRd
.fitMethod2 <- function(p, q) {
    idx16 <- which.min(abs(p - 0.16))
    idx84 <- which.min(abs(p - 0.84))
    if (length(idx16) == 0 || length(idx84) == 0) {
        # fallback
        return(templateDF("2", NA, NA, NA))
    }
    y16 <- q[idx16]
    y84 <- q[idx84]
    sigma <- (y84 - y16) / 2
    center <- 0.5 * (y16 + y84)

    w <- wgt(p)
    CDF <- stats::pnorm((q - center) / sigma)
    dVal <- max(abs(CDF - p))
    sseVal <- sum(w * (CDF - p)^2)

    templateDF("2", sigma, dVal, sseVal)
}

# -- Method 3: MAD -------------------------------------------------------------
#' @keywords internal
#' @noRd
.fitMethod3 <- function(p, q) {
    idx25 <- which.min(abs(p - 0.25))
    idx50 <- which.min(abs(p - 0.50))
    idx75 <- which.min(abs(p - 0.75))

    if (any(c(idx25, idx50, idx75) == 0)) {
        return(templateDF("3", NA, NA, NA))
    }
    y25 <- q[idx25]
    y50 <- q[idx50]
    y75 <- q[idx75]

    madVal <- (y75 - y25) / 2
    sigma <- 1.4826 * madVal

    w <- wgt(p)
    CDF <- stats::pnorm((q - y50) / sigma)
    dVal <- max(abs(CDF - p))
    sseVal <- sum(w * (CDF - p)^2)

    templateDF("3", sigma, dVal, sseVal)
}

# -- Method 4: Student-t -------------------------------------------------------
#' @keywords internal
#' @noRd
# -- Method 4: Student-t -------------------------------------------------------
#' @keywords internal
#' @noRd
.fitMethod4 <- function(p, q, nuStart = 10) {
    idx50 <- which.min(abs(p - 0.50))
    if (length(idx50) == 0) {
        return(templateDF("4", NA, Inf, Inf))
    }
    mu <- q[idx50]

    idx02 <- which.min(abs(p - 0.02))
    idx98 <- which.min(abs(p - 0.98))
    if (idx02 == 0 || idx98 == 0) {
        return(templateDF("4", NA, Inf, Inf))
    }
    y02 <- q[idx02]
    y98 <- q[idx98]

    t98 <- stats::qt(0.98, df = nuStart)
    s0 <- (y98 - y02) / (2 * t98)

    w <- wgt(p)
    objFun <- function(par) {
        nu <- par[1]
        s <- par[2]
        if (nu <= 2 || s <= 0) {
            return(rep(1e10, length(q)))
        }
        cdf <- stats::pt((q - mu) / s, df = nu)
        cdf - p # unweighted residual
    }

    if (!requireNamespace("minpack.lm", quietly = TRUE)) {
        # fallback
        return(templateDF("4", NA, Inf, Inf))
    }

    ## ---- PATCH: increase maxiter to 200 and optionally suppress solver warnings
    ctrl <- minpack.lm::nls.lm.control(maxiter = 200)

    fitTry <- try(
        suppressWarnings(
            minpack.lm::nls.lm(
                par = c(nuStart, s0),
                fn = objFun,
                lower = c(2.001, 1e-8),
                control = ctrl
            )
        ),
        silent = TRUE
    )
    if (inherits(fitTry, "try-error")) {
        return(templateDF("4", NA, Inf, Inf))
    }

    nuHat <- fitTry$par[1]
    sHat <- fitTry$par[2]
    if (nuHat <= 2) {
        return(templateDF("4", NA, Inf, Inf))
    }

    cdf <- stats::pt((q - mu) / sHat, df = nuHat)
    resid <- cdf - p
    sseVal <- sum(w * (resid^2))
    dVal <- max(abs(resid))

    sigmaEst <- sHat * sqrt(nuHat / (nuHat - 2))

    templateDF("4", sigmaEst, dVal, sseVal,
        nuHat = nuHat, sHat = sHat, muHat = mu
    )
}


# -- Method 5: Johnson SU ------------------------------------------------------
#' @keywords internal
#' @noRd
# ---------------------------------------------------------------------------
#  Method 5  – Johnson SU fit  (same logic, higher maxiter, no new arguments)
.fitMethod5 <- function(p, q, meanValue) {
    w <- wgt(p)

    resFun <- function(par) {
        gamma <- par[1]
        delta <- par[2]
        lambda <- par[3]
        if (delta <= 0 || lambda <= 0) {
            return(rep(1e10, length(q)))
        }
        xi <- meanValue - lambda * sinh(gamma / delta) * exp(1 / (2 * delta^2))
        Z <- gamma + delta * asinh((q - xi) / lambda)
        stats::pnorm(Z) - p
    }

    if (!requireNamespace("minpack.lm", quietly = TRUE)) {
        return(templateDF("5", NA, NA, NA))
    }

    start_par <- c(
        gamma = 0,
        delta = 1,
        lambda = (max(q) - min(q)) / 4
    )

    ## ---- ONLY CHANGE: raise maxiter so the optimiser stops warning -------
    ctrl <- minpack.lm::nls.lm.control(maxiter = 200) # default was 50

    fit_try <- try(
        minpack.lm::nls.lm(
            par     = start_par,
            fn      = resFun,
            lower   = c(-10, 1e-6, 1e-6),
            upper   = c(10, 1e6, 1e6),
            control = ctrl # <-- here
        ),
        silent = TRUE
    )
    if (inherits(fit_try, "try-error")) {
        return(templateDF("5", NA, NA, NA))
    }

    par_hat <- fit_try$par
    gamma_hat <- par_hat[1]
    delta_hat <- par_hat[2]
    lambda_hat <- par_hat[3]
    xi_hat <- meanValue -
        lambda_hat * sinh(gamma_hat / delta_hat) *
            exp(1 / (2 * delta_hat^2))

    Zvec <- gamma_hat + delta_hat * asinh((q - xi_hat) / lambda_hat)
    resid <- stats::pnorm(Zvec) - p
    sseVal <- sum(w * resid^2)
    dVal <- max(abs(resid))

    term1 <- exp(1 / delta_hat^2) * cosh(2 * gamma_hat / delta_hat)
    term2 <- exp(1 / delta_hat^2) * sinh(gamma_hat / delta_hat)^2
    varEst <- lambda_hat^2 * (term1 - term2 - 1)
    if (is.na(varEst) || varEst <= 0) {
        return(templateDF("5", NA, dVal, sseVal))
    }

    sigmaEst <- sqrt(varEst)
    templateDF("5", sigmaEst, dVal, sseVal,
        gamma  = gamma_hat,
        delta  = delta_hat,
        lambda = lambda_hat,
        xi     = xi_hat
    )
}


# -- Method 6: Piecewise stepwise CDF ------------------------------------------
#' @keywords internal
#' @noRd
.fitMethod6 <- function(p, q, meanValue) {
    # 1) compute piecewise variance from midpoints
    K <- length(p)
    dp <- diff(p)
    midpoints <- 0.5 * (q[-K] + q[-1])
    varEst <- sum(dp * (midpoints - meanValue)^2)
    # tails
    if (p[1] > 0) {
        varEst <- varEst + p[1] * (q[1] - meanValue)^2
    }
    if (p[K] < 1) {
        varEst <- varEst + (1 - p[K]) * (q[K] - meanValue)^2
    }
    sigmaEst <- sqrt(varEst)

    templateDF("6", sigmaEst, d = NA, sse = NA)
}

# Validate user input ------------------------------------------------
#' @keywords internal
#' @noRd
checkInputs <- function(meanValue, p, q) {
    if (!is.numeric(meanValue) || length(meanValue) != 1L || !is.finite(meanValue)) {
        stop("`meanValue` must be one finite numeric scalar.", call. = FALSE)
    }

    if (length(p) != length(q)) {
        stop("`p` and `q` must have the same length.", call. = FALSE)
    }

    if (any(!is.finite(p)) || any(!is.finite(q))) {
        stop("`p` and `q` must not contain NA, NaN, or Inf.", call. = FALSE)
    }

    if (any(p <= 0 | p >= 1)) {
        stop("All probabilities in `p` must lie strictly between 0 and 1.", call. = FALSE)
    }

    # Re-order if not already increasing
    if (is.unsorted(p)) {
        ord <- order(p)
        p <<- p[ord]
        q <<- q[ord]
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
wgt <- function(p) 1 / (p * (1 - p))



#' @keywords internal
#' @noRd
templateDF <- function(method, sigma, d, sse,
                       alpha = NA_real_, beta = NA_real_,
                       nuHat = NA_real_, sHat = NA_real_, muHat = NA_real_,
                       gamma = NA_real_, delta = NA_real_,
                       lambda = NA_real_, xi = NA_real_) {
    data.frame(
        method = method,
        sigma  = sigma,
        d      = d,
        sse    = sse,
        alpha  = alpha,
        beta   = beta,
        nuHat  = nuHat,
        sHat   = sHat,
        muHat  = muHat,
        gamma  = gamma,
        delta  = delta,
        lambda = lambda,
        xi     = xi
    )
}


#' @keywords internal
#' @noRd
safeLm <- function(y, x, w) {
    # Use lm.wfit to handle weights
    stats::lm.wfit(x = cbind(1, x), y = y, w = w)
}

# Draw ε given the winning method -----------------------------------
#' @keywords internal
#' @noRd
simulateEps <- function(n, method, df) {
    if (method %in% c("1", "2", "3")) {
        return(stats::rnorm(n))
    }
    if (method == "4") {
        nu <- df$nuHat[df$method == "4"][1]
        if (!is.na(nu)) {
            return(stats::rt(n, df = nu))
        }
        return(stats::rnorm(n))
    }
    ## Johnson SU (5) – no closed-form ε, default to N(0,1)
    stats::rnorm(n)
}

# Vectorised piecewise linear sampler --------------------------------
#' @keywords internal
#' @noRd
samplePiecewiseY <- function(n, p, q, meanValue, sigmaFallback) {
    ## Include 0 and 1 so approx() is defined everywhere, then interpolate
    xx <- c(0, p, 1)
    yy <- c(q[1], q, q[length(q)])

    y <- stats::approx(
        x = xx,
        y = yy,
        xout = stats::runif(n),
        method = "linear",
        ties = "ordered",
        rule = 2
    )$y

    ## Mean-shift to honour the supplied meanValue (and preserve variance)
    y + (meanValue - mean(y))
}
