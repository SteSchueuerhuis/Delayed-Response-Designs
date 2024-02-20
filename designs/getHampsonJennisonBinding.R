
getHampsonJennisonBinding <- function(alpha = 0.025,
                                      beta = 0.2,
                                      a_spend = c("asP", "asOF"),
                                      b_spend = c("bsP", "bsOF"),
                                      I1 = 0.5,
                                      I_d = 0.2,
                                      tol = 1e-8) {
  design <- list()
  class(design) <- c("BindingHampsonJennison", class(design))

  a_spend <- match.arg(a_spend)
  b_spend <- match.arg(b_spend)

  I1_d <- I1 + I_d
  information <- c(I1, I1_d, 1)

  sigma <- getCovarianceFromInformation(c(I1, I1_d, 1))

  gsd <- rpact::getDesignGroupSequential(
    kMax = 2,
    alpha = alpha,
    beta = beta,
    typeOfDesign = a_spend,
    typeBetaSpending = b_spend,
    informationRates = c(I1, 1),
    bindingFutility = T
  )


  lower <- gsd$futilityBounds
  upper <- gsd$criticalValues

  l1 <- lower[1]
  u1 <- upper[1]
  d2 <- upper[2]

  d1_d <- rootSearch(f = function(d1_d) {
    bounds <- matrix(c(u1, d1_d, d2, l1, d1_d, d2), ncol = 3, byrow = T)
    P <- getDelayedResponseProbabilities(
      method = class(design)[1],
      information = information,
      bounds = bounds,
      delta = 0,
      stDev = 1,
      nMax = 2
    )

    return(P$reversalProbsUpper - P$reversalProbsLower)
  }, a = 1e-4, b = 10, tol = tol, method = "Bisection")$root


  design$alpha <- alpha
  design$beta <- beta
  design$alphaSpending <- ifelse(a_spend == "asOF", "O'Brien-Fleming", "Pocock")
  design$betaSpending <- ifelse(b_spend == "bsOF", "O'Brien-Fleming", "Pocock")
  design$recStopInfo <- I1
  design$decInfo <- I1_d
  design$u1 <- u1
  design$l1 <- l1
  design$d <- c(d1_d, d2)
  design$cp <- CP(u1, I1, I_d, d1_d)

  return(design)
}


summary.BindingHampsonJennison <- function(x, ...) {
  obj <- getHampsonJennisonBinding(
    alpha = x$alpha,
    beta = x$beta,
    a_spend = ifelse(x$alphaSpending == "O'Brien-Fleming", "asOF", "asP"),
    b_spend = ifelse(x$betaSpending == "O'Brien-Fleming", "bsOF", "bsP"),
    I1 = x$recStopInfo,
    I_d = x$decInfo - x$recStopInfo,
    tol = 1e-8
  )

  obj <- lapply(obj, roundList, digits = 3)

  alpha <- obj$alpha
  beta <- obj$beta
  alphaSpending <- obj$alphaSpending
  betaSpending <- obj$betaSpending
  recStopInfo <- obj$recStopInfo
  decInfo <- obj$decInfo
  u1 <- obj$u1
  l1 <- obj$l1
  d <- obj$d
  cp <- obj$cp

  label <- "Two-Stage Binding Hampson & Jennison (2013)-Design"

  cat("---------------------------------------------------", "\n")
  cat(bold(green("Design:")), label, "\n")
  cat("---------------------------------------------------", "\n")
  cat(green("Design parameter"), "\n")
  cat(greeks("alpha"), ":", bold(alpha), "--", greeks("beta"), ":", bold(beta), "\n")
  cat("---------------------------------------------------", "\n")
  cat(green("Interim information values"), "\n")
  cat("Information available at analysis for recruitment stopping:", bold(recStopInfo), "\n")
  cat("Information available at interim significance test:", bold(decInfo), "\n")
  cat("---------------------------------------------------", "\n")
  cat(green("Spending functions:"), "\n")
  cat("Alpha-spending:", bold(alphaSpending), "\n")
  cat("Beta-spending: ", bold(betaSpending), "\n")
  cat("---------------------------------------------------", "\n")
  cat(green("Boundaries:"), "\n")
  cat("Continuation region: (", bold(l1), ",", bold(u1), ")", "\n")
  cat("Decision boundaries: {", bold(d[1]), ",", bold(d[2]), "}", "\n")
  cat("CP at Interim: ", bold(cp))
}
