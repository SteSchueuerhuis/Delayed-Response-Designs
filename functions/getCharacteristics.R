
# FUNCTION TO CALCULATE CHARACTERISTICS  ----------------------------------

getCharacteristics <- function(method = c("BindingHampsonJennison", "NonBindingHampsonJennison", "RecPauseDesign"),
                               alpha = 0.025,
                               power = 0.8,
                               a_spend = c("asP", "asOF"),
                               b_spend = c("bsP", "bsOF"),
                               I1 = 0.5,
                               I_d = 0.2,
                               tol = 1e-8,
                               delta = 0.3,
                               stDev = 1,
                               nMax = 200,
                               gamma_l = 0.9,
                               gamma_u = 0.9,
                               print = F) {
  parameters <- list(power, delta, nMax, I1)
  parameterNames <- c("power", "delta", "nMax", "I1")

  reqParIndex <- which(unlist(lapply(parameters, is.na)))

  if (length(reqParIndex) > 1) {
    stop("Only one of 'power', 'delta', nMax', 'I1' shall be NA.")
  }

  reqPar <- parameterNames[reqParIndex]

  pwr <- quote({
    analyseDesign(
      method = method,
      alpha = alpha,
      beta = 1 - power,
      a_spend = a_spend,
      b_spend = b_spend,
      I1 = I1,
      I_d = I_d,
      delta = delta,
      stDev = stDev,
      nMax = nMax,
      gamma_l = gamma_l,
      gamma_u = gamma_u
    )$power
  })

  if (is.na(power)) {
    power <- eval(pwr)
  }
  if (is.na(nMax)) {
    nMax <- rootSearch(function(nMax) {
      eval(pwr) - power
    }, a = 2 + 1e-10, b = 1e+07, method = "Bisection", tol = tol)$root
  }
  if (is.na(I1)) {
    I1 <- rootSearch(function(I1) {
      eval(pwr) - power
    }, a = 1e-07, b = 1 - 1e-07, method = "Bisection", tol = tol)$root
  }
  if (is.na(delta)) {
    delta <- rootSearch(function(delta) {
      eval(pwr) - power
    }, a = 1e-07, b = 6, method = "Bisection", tol = tol)$root
  }

  if (print) {
    structure(get(reqPar), names = reqPar)
    cat("The omitted variable was:", green(parameterNames[reqParIndex]), "\n")
    cat("_______________________________________________________\n\n")
    cat("Given the input of power = ", power, "the calculated required value for", parameterNames[reqParIndex], "is", green(get(reqPar)))
  }
  return(get(reqPar))
}
