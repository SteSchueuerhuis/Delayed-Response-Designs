
# FUNCTION THAT CALCULATES POWER, EXPECTED SAMPLE SIZE AND EXPECTETED TRIAL DURATION --------

analyseDesign <- function(method = c("BindingHampsonJennison", "NonBindingHampsonJennison", "RecPauseDesign"),
                          alpha = 0.025,
                          beta = 0.2,
                          a_spend = c("asP", "asOF"),
                          b_spend = c("bsP", "bsOF"),
                          I1 = 0.5,
                          I_d = 0.2,
                          delta = 0.35,
                          stDev = 1,
                          nMax = 200,
                          rr = 10,
                          gamma_l = 0.8,
                          gamma_u = 0.8) {
  a_spend <- match.arg(a_spend)
  b_spend <- match.arg(b_spend)
  method <- match.arg(method)

  I1_d <- I1 + I_d
  information <- c(I1, I1_d, 1)

  tFinal <- nMax / rr + (I_d * nMax) / rr
  t1.tilde <- nMax / rr * I1_d
  tFinal.restart <- tFinal + (I_d * nMax) / rr

  des <- switch(method,
    "BindingHampsonJennison" = getHampsonJennisonBinding,
    "NonBindingHampsonJennison" = getHampsonJennisonNonBinding,
    "RecPauseDesign" = getRecruitmentPauseDesign
  )

  if (method != "RecPauseDesign") {
    des_args <- list(
      alpha = alpha,
      beta = beta,
      a_spend = a_spend,
      b_spend = b_spend,
      I1 = I1,
      I_d = I_d,
      tol = 1e-8
    )
  } else {
    des_args <- list(
      alpha = alpha,
      beta = beta,
      a_spend = a_spend,
      I1 = I1,
      I_d = I_d,
      tol = 1e-8,
      gamma_l = gamma_l,
      gamma_u = gamma_u
    )
  }

  des_result <- do.call(des, des_args)

  l1 <- des_result$l1
  u1 <- des_result$u1
  d <- des_result$d

  if (method != "RecPauseDesign") {
    bounds <- matrix(c(u1, d[1], d[2], l1, d[1], d[2]), ncol = 3, byrow = TRUE)
  } else {
    lower <- des_result$lower
    upper <- des_result$upper

    l1 <- lower[1]
    l1.tilde <- lower[2]

    u1 <- upper[1]
    d1_d <- upper[2]
    d2 <- upper[3]

    bounds <- matrix(c(u1, d1_d, d2, l1, l1.tilde, d2), ncol = 3, byrow = TRUE)
  }

  P <- getDelayedResponseProbabilities(
    method = method,
    information = information,
    bounds = bounds,
    delta = delta,
    stDev = stDev,
    nMax = nMax
  )

  power <- P$efficacyStopInterim + P$efficacyFinal

  if (method != "RecPauseDesign") {
    ASNH1 <- (I1_d * nMax) * (1 - P$continuation) + nMax * P$continuation
    ASTH1 <- t1.tilde * (1 - P$continuation) + tFinal * P$continuation
  } else {
    ASNH1 <- (I1_d * nMax) * (P$efficacyStopInterim + P$futilityStopInterim) + nMax * (P$continuation + P$restart)
    ASTH1 <- t1.tilde * (P$efficacyStopInterim + P$futilityStopInterim) + tFinal * P$continuation + tFinal.restart * P$restart
  }

  result <- list(
    effect = delta,
    stDev = stDev,
    nMax = nMax,
    information = c(I1, I1_d, 1),
    power = ifelse(power <= 1, power, NA),
    ASNH1 = ASNH1,
    ASTH1 = ASTH1
  )

  return(result)
}
