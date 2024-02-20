
getRecruitmentPauseDesign <- function(alpha = 0.025,
                                      beta = 0.2,
                                      I1 = 0.5,
                                      I_d = 0.2,
                                      l1.tilde = 0.5,
                                      tol = 1e-8,
                                      a_spend = c("asP", "asOF"),
                                      gamma_l = 0.8,
                                      gamma_u = 0.8) {
  design <- list()
  class(design) <- c("RecPauseDesign", class(design))

  a_spend <- match.arg(a_spend)
  aS <- ifelse(a_spend == "asP", P, OBF)

  I1_d <- I1 + I_d
  information <- c(I1, I1_d, 1)

  l1 <- rootSearch(f = function(l1) {
    epsilon <- CP(l1, I1, I_d, d1.tilde = l1.tilde, assumed = T, delta = 0, sd = 1, nMax = 0, upper = F) - gamma_l
    return(epsilon)
  }, a = -15, b = l1.tilde + 5, tol = 1e-8, method = "Bisection")$root

  to_be_solved <- function(x) {
    eq1 <- CP(x[1], I1, I_d, x[2], assumed = F, upper = T) - gamma_u
    bounds <- matrix(c(x[1], x[2], NA, l1, l1.tilde, NA), ncol = 3, byrow = T)

    P <- getDelayedResponseProbabilities(
      method = class(design)[1],
      information = information,
      bounds = bounds,
      delta = 0,
      stDev = 1,
      nMax = 10
    )

    eq2 <- P$efficacyStopInterim - aS(alpha, I1_d)
    return(c(eq1, eq2))
  }

  init <- c(2, 2)
  upperBoundsInterim <- nleqslv(init, to_be_solved)$x
  u1 <- upperBoundsInterim[1]
  d1_d <- upperBoundsInterim[2]

  d2 <- rootSearch(f = function(d2) {
    bounds <- matrix(c(u1, d1_d, d2, l1, l1.tilde, d2), ncol = 3, byrow = T)

    P <- getDelayedResponseProbabilities(
      method = class(design)[1],
      information = information,
      bounds = bounds,
      delta = 0,
      stDev = 1,
      nMax = 10
    )


    alpha_left <- alpha - aS(alpha, I1_d)

    return(P$efficacyFinal - alpha_left)
  }, a = 0, b = 5, method = "Bisection")$root

  design$alpha <- alpha
  design$beta <- beta
  design$alphaSpending <- ifelse(a_spend == "asOF", "O'Brien-Fleming", "Pocock")
  design$recStopInfo <- I1
  design$decInfo <- I1_d
  design$upper <- c(u1, d1_d, d2)
  design$lower <- c(l1, l1.tilde)
  design$gamma_l <- gamma_l
  design$gamma_u <- gamma_u

  return(design)
}

summary.RecPauseDesign <- function(x, ...) {
  obj <- getRecruitmentPauseDesign(
    alpha = x$alpha,
    beta = x$beta,
    a_spend = ifelse(x$alphaSpending == "O'Brien-Fleming", "asOF", "asP"),
    I1 = x$recStopInfo,
    I_d = x$decInfo - x$recStopInfo,
    tol = 1e-8,
    l1.tilde = x$lower[2],
    gamma_l = x$gamma_l,
    gamma_u = x$gamma_u
  )

  obj <- lapply(obj, roundList, digits = 3)

  alpha <- obj$alpha
  beta <- obj$beta
  alphaSpending <- obj$alphaSpending
  recStopInfo <- obj$recStopInfo
  decInfo <- obj$decInfo
  u1 <- obj$upper[1]
  l1 <- obj$lower[1]
  l1.tilde <- obj$lower[2]
  d <- obj$upper[2:3]

  label <- "Two-Stage Recuitment Pausing Method"

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
  cat("---------------------------------------------------", "\n")
  cat(green("Boundaries:"), "\n")
  cat("Continuation region: (", bold(l1), ",", bold(u1), ")", "\n")
  cat("Decision boundaries: {", bold(d[1]), ",", bold(d[2]), "}", "\n")
  cat("Futility boundary: ", bold(l1.tilde), "\n")
  cat("---------------------------------------------------", "\n")
}
