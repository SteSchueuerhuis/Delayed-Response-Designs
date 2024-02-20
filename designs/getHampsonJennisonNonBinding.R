
getHampsonJennisonNonBinding <- function(alpha = 0.025,
                                         beta = 0.2,
                                         a_spend = c("asP", "asOF"),
                                         b_spend = c("bsP", "bsOF"),
                                         I1 = 0.5,
                                         I_d = 0.2,
                                         tol = 1e-8,
                                         d1_d = qnorm(1 - 0.025)) {
  design <- list()
  class(design) <- c("NonBindingHampsonJennison", class(design))

  a_spend <- match.arg(a_spend)
  b_spend <- match.arg(b_spend)

  aS <- ifelse(a_spend == "asP", P, OBF)
  bS <- ifelse(b_spend == "bsP", P, OBF)

  I1_d <- I1 + I_d
  information <- c(I1, I1_d, 1)

  sigma <- getCovarianceFromInformation(information)

  u1 <- rootSearch(f = function(u1) {
    pmvnorm(
      lower = c(u1, d1_d),
      upper = c(Inf, Inf),
      mean = c(0, 0),
      sigma = sigma[1:2, 1:2]
    ) - aS(alpha, I1)
  }, a = 1e-4, b = 10, tol = tol, method = "Bisection")$root

  d2 <- rootSearch(f = function(d2) {
    pmvnorm(
      lower = c(-Inf, d2),
      upper = c(u1, Inf),
      mean = c(0, 0),
      sigma = sigma[c(1, 3), c(1, 3)]
    ) - (alpha - aS(alpha, I1))
  }, a = 1e-4, b = 10, tol = tol, method = "Bisection")$root

  constraints <- function(a, b) {
    e1 <- pnorm(a, mean = sqrt(I1) * b) +
      pmvnorm(
        lower = c(u1, -Inf),
        upper = c(Inf, d1_d),
        mean = c(sqrt(I1) * b, sqrt(I1_d) * b),
        sigma = sigma[1:2, 1:2]
      ) - bS(beta, I1)

    e2 <- pmvnorm(
      lower = c(a, -Inf),
      upper = c(u1, d2),
      mean = c(sqrt(I1) * b, b),
      sigma = sigma[c(1, 3), c(1, 3)]
    ) - (beta - bS(beta, I1))

    return(c(e1, e2) - c(0, 0))
  }

  to_be_solved <- function(x) crossprod(constraints(x[1], x[2]))

  sol <- optim(c(0.2, 3), to_be_solved)$par
  l1 <- sol[1]
  shift <- sol[2]


  design$alpha <- alpha
  design$beta <- beta
  design$alphaSpending <- ifelse(a_spend == "asOF", "O'Brien-Fleming", "Pocock")
  design$betaSpending <- ifelse(b_spend == "bsOF", "O'Brien-Fleming", "Pocock")
  design$recStopInfo <- I1
  design$decInfo <- I1_d
  design$u1 <- u1
  design$l1 <- l1
  design$d <- c(d1_d, d2)
  design$shift <- shift
  design$cp <- CP(u1, I1, I_d, d1_d)
  # message(cat(red(
  #   "Beware: the lower boundary", l1, "does only precisely comply with a power of", 1 - beta,
  #   "if the sample size n, the unknown effect", greeks("delta"),
  #   "as well as the true, but unkown variance", greeks("sigma^2"),
  #   "are in a ratio that equals", shift
  # )))

  return(design)
}

summary.NonBindingHampsonJennison <- function(x, ...) {
  obj <- getHampsonJennisonNonBinding(
    alpha = x$alpha,
    beta = x$beta,
    a_spend = ifelse(x$alphaSpending == "O'Brien-Fleming", "asOF", "asP"),
    b_spend = ifelse(x$betaSpending == "O'Brien-Fleming", "bsOF", "bsP"),
    I1 = x$recStopInfo,
    I_d = x$decInfo - x$recStopInfo,
    tol = 1e-8,
    d1_d = x$d[1]
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
  shift <- obj$shift
  cp <- obj$cp

  label <- "Two-Stage Non-Binding Hampson & Jennison (2013)-Design"

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
  cat("---------------------------------------------------", "\n")
  cat("Shift value:", bold(shift), "\n")
  cat("CP at Interim: ", bold(cp))
}
