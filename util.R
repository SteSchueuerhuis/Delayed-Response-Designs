# ROOT SEARCH FUNCTION ----------------------------------------

rootSearch <- function(f, a, b,
                       tol = 1e-8,
                       maxIter = 10000,
                       method = c(
                         "Bisection",
                         "Regula Falsi",
                         "Newton-Raphson",
                         "Uniroot"
                       )) {
  if (!("tibble" %in% installed.packages()) ||
    !("numDeriv") %in% installed.packages()) {
    print("Required Packages will be installed")
  }

  if (!require(numDeriv)) {
    install.packages("numDeriv")
    require(numDeriv)
  }
  if (!require(tibble)) {
    install.packages("tibble")
    require(tibble)
  }

  k <- vector(length = maxIter)

  # if (f(a) == 0) {
  #   return(a)
  #   print("lower boundary is root")
  # }
  # if (f(b) == 0) {
  #   return(b)
  #   print("upper boundary is root")
  # }

  method <- match.arg(method)
  iter <- 1

  if (method == "Bisection") {
    tryCatch(
      {
        while ((abs(a - b) > tol) && (iter <= maxIter)) {
          root.guess1 <- (a + b) / 2
          f.guess <- f(root.guess1)
          ifelse(f.guess * f(b) < 0, a <- root.guess1, b <- root.guess1)
          iter <- iter + 1
          # cat("At iteration", iter, "the root guess is:", root.guess1, "\n")
        }
      },
      error = function(e) {
        message("Root Bisection Search did not work, stats::uniroot used instead")
        return(stats::uniroot(
          f = f,
          interval = c(a, b)
        )$root)
      }
    )

    root <- root.guess1
  }

  if (method == "Regula Falsi") {
    tryCatch(
      {
        while ((abs(b - a) > tol) && (iter < maxIter)) {
          root.guess2 <- b - f(b) * (b - a) / (f(b) - f(a))
          a <- b
          b <- root.guess2
          iter <- iter + 1
          # cat("At iteration", iter, "value of x.n is:", root.guess2, "\n")
        }
      },
      error = function(e) {
        message("Regula Falsi Method did not find root, stats::uniroot used instead")
        return(stats::uniroot(
          f = f,
          interval = c(a, b)
        )$root)
      }
    )

    root <- b
  }

  if (method == "Newton-Raphson") {
    tryCatch(
      {
        while (iter <= maxIter) {
          da <- genD(func = f, x = a)$D[1]
          b <- a - (f(a) / da)
          k[iter] <- b
          iter <- iter + 1
          if (abs(b - a) < tol) {
            root <- tail(k, n = 1)
          } else {
            a <- b
          }
        }
      },
      error = function(e) {
        message("Newton-Raphson Method did not find root, stats::uniroot used instead")
        return(stats::uniroot(
          f = f,
          interval = c(a, b)
        )$root)
      }
    )
  }

  if (method == "Uniroot") {
    root <- stats::uniroot(
      f = f,
      lower = a,
      upper = b
    )$root
  }

  return(tibble(
    method = method,
    root = root
  ))
}

# UTIL FUNCTIONS ----------------------------------------


# O'Brien-Fleming spending function, X is either alpha or beta, I1 <= 1
OBF <- function(X, I1) {
  res <- ifelse(I1 == 0, 0, min(2 * (1 - pnorm(qnorm(1 - X / 2) / sqrt(I1))), X))

  return(res)
}

# Pocock spending funtion, X is either alpha or beta, I1 <= 1
P <- function(X, I1) {
  return(min(X * (log(1 + (exp(1) - 1) * I1)), X))
}

# List rounding
roundList <- function(x, digits = 3) {
  if (is.numeric(x)) {
    x <- round(x, digits = digits)
  }

  return(x)
}

# Simulation via Cholesky-Decomposition
simulateNorm <- function(n, mean, sigma) {
  sigma <- as.matrix(sigma)
  mean <- as.matrix(mean)

  if (!(ncol(sigma) == 1 && nrow(sigma) == 1)) {
    if (any(eigen(sigma)$values < 0)) {
      stop("sigma is no covariancematrix")
    }
  }

  if (nrow(sigma) != ncol(sigma)) {
    stop("sigma is no square matrix")
  }

  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma need to be of same dimension")
  }

  ncol <- ncol(sigma)
  mu <- rep(mean, each = n)
  return(mu + matrix(rnorm(n * ncol), ncol = ncol) %*% chol(sigma))
}

# Conditional power
CP <- function(z1, I1, I_d, d1.tilde, assumed = F, delta = NA, sd = NA, nMax = NA, upper = T) {
  if (assumed) {
    if (any(is.na(c(delta, sd, nMax)))) {
      stop("Effect, standard deviation and maximum sample size need to be assumed")
    }
  }

  I1.tilde <- I1 + I_d

  if (assumed) {
    # last term is = 0 if delta = 0 <=> independent of nMax
    cp <- 1 - pnorm(d1.tilde * sqrt(I1.tilde / I_d) - z1 * sqrt(I1 / I_d) - sqrt(nMax * I_d) * delta / sqrt(2 * sd^2))
  } else {
    cp <- 1 - pnorm(((d1.tilde * sqrt(I1.tilde) - z1 * sqrt(I1)) / sqrt(I_d)) - z1 * sqrt(I_d / I1))
  }

  if (!upper) {
    cp <- 1 - cp
  }

  return(cp)
}

# Probabilities
getDelayedResponseProbabilities <- function(method, information, bounds, delta, stDev, nMax) {
  effect <- delta / sqrt(2 * stDev^2) * sqrt((information * nMax))
  sigma <- getCovarianceFromInformation(information)

  efficacyStopInterim <- tryCatch(
    {
      pmvnorm(
        lower = c(bounds[1, 1], bounds[1, 2]),
        upper = c(Inf, Inf),
        mean = effect[1:2],
        sigma = sigma[1:2, 1:2]
      )[1]
    },
    error = function(e) {
      return(NA)
    }
  )

  futilityStopInterim <- tryCatch(
    {
      pmvnorm(
        lower = c(-Inf, -Inf),
        upper = c(bounds[2, 1], bounds[2, 2]),
        mean = effect[1:2],
        sigma = sigma[1:2, 1:2]
      )[1]
    },
    error = function(e) {
      return(NA)
    }
  )

  continuation <- tryCatch(
    {
      pmvnorm(
        lower = c(bounds[2, 1]),
        upper = c(bounds[1, 1]),
        mean = effect[1],
        sigma = sigma[1, 1]
      )[1]
    },
    error = function(e) {
      return(NA)
    }
  )

  efficacyFinal <- tryCatch(
    {
      pmvnorm(
        lower = c(bounds[2, 1], bounds[1, 3]),
        upper = c(bounds[1, 1], Inf),
        mean = effect[c(1, 3)],
        sigma = sigma[c(1, 3), c(1, 3)]
      )[1]
    },
    error = function(e) {
      return(NA)
    }
  )

  reversalProbsUpper <- reversalProbsLower <- NA
  restart <- NA

  if (method == "BindingHampsonJennison") {
    efficacyStopInterim <- efficacyStopInterim + tryCatch(
      {
        pmvnorm(
          lower = c(-Inf, bounds[1, 2]),
          upper = c(bounds[2, 1], Inf),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )

    futilityStopInterim <- futilityStopInterim + tryCatch(
      {
        pmvnorm(
          lower = c(bounds[1, 1], -Inf),
          upper = c(Inf, bounds[1, 2]),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )

    reversalProbsLower <- tryCatch(
      {
        pmvnorm(
          lower = c(-Inf, bounds[1, 2]),
          upper = c(bounds[2, 1], Inf),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )

    reversalProbsUpper <- tryCatch(
      {
        pmvnorm(
          lower = c(bounds[1, 1], -Inf),
          upper = c(Inf, bounds[1, 2]),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )
  }

  if (method == "RecPauseDesign") {
    efficacyStopInterim <- efficacyStopInterim + tryCatch(
      {
        pmvnorm(
          lower = c(-Inf, bounds[1, 2]),
          upper = c(bounds[2, 1], Inf),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )

    futilityStopInterim <- futilityStopInterim + tryCatch(
      {
        pmvnorm(
          lower = c(bounds[1, 1], -Inf),
          upper = c(Inf, bounds[2, 1]),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    )

    restart <- tryCatch(
      {
        pmvnorm(
          lower = c(bounds[1, 1], bounds[2, 2]),
          upper = c(Inf, bounds[1, 2]),
          mean = effect[1:2],
          sigma = sigma[1:2, 1:2]
        )[1]
      },
      error = function(e) {
        return(0)
      }
    ) +
      tryCatch(
        {
          pmvnorm(
            lower = c(-Inf, bounds[2, 2]),
            upper = c(bounds[2, 1], bounds[1, 2]),
            mean = effect[1:2],
            sigma = sigma[1:2, 1:2]
          )[1]
        },
        error = function(e) {
          return(0)
        }
      )
    efficacyFinal <- efficacyFinal + tryCatch(
      {
        pmvnorm(
          lower = c(bounds[1, 1], bounds[2, 2], bounds[1, 3]),
          upper = c(Inf, bounds[1, 2], Inf),
          mean = effect,
          sigma = sigma
        )[1]
      },
      error = function(e) {
        return(0)
      }
    ) +
      tryCatch(
        {
          pmvnorm(
            lower = c(-Inf, bounds[2, 2], bounds[1, 3]),
            upper = c(bounds[2, 1], bounds[1, 2], Inf),
            mean = effect,
            sigma = sigma
          )[1]
        },
        error = function(e) {
          return(0)
        }
      )
  }

  result <- list(
    efficacyStopInterim = efficacyStopInterim,
    futilityStopInterim = futilityStopInterim,
    continuation = continuation,
    efficacyFinal = efficacyFinal,
    reversalProbsUpper = reversalProbsUpper,
    reversalProbsLower = reversalProbsLower,
    restart = restart
  )

  return(result)
}

# Covariance
getCovarianceFromInformation <- function(information) {
  sigma <- matrix(NA, ncol = length(information), nrow = length(information))

  for (k in seq_along(information)) {
    for (h in seq_along(information)) {
      sigma[k, h] <- sqrt(information[k] / information[h])

      if (sigma[k, h] > 1) sigma[k, h] <- 1 / sigma[k, h]
    }
  }

  return(sigma)
}

# Data preparation
mergeData <- function(X, Y, Z, longformat = F) {
  by <- colnames(X$data)[1]

  if (!.isPlotObject(X)) stop("mergePlots requires objects from class 'plotObject")
  if (!.isPlotObject(Y)) stop("mergePlots requires objects from class 'plotObject")
  if (!.isPlotObject(Z)) stop("mergePlots requires objects from class 'plotObject")

  if (any(X$data[, by] != Y$data[, by] || X$data[, by] != Z$data[, by] || Y$data[, by] != Z$data[, by])) {
    stop("If methods are to be jointly plotted, considered effects need to be of same scale")
  }


  data <- data.frame(
    scale = c(X$data[, by], Y$data[, by], Z$data[, by]),
    power = c(X$data$power, Y$data$power, Z$data$power),
    ASNH1 = c(X$data$ASNH1, Y$data$ASNH1, Z$data$ASNH1),
    ASTH1 = c(X$data$ASTH1, Y$data$ASTH1, Z$data$ASTH1),
    method = c(X$data$design, Y$data$design, Z$data$design)
  )

  if (longformat) {
    data <- data.frame(
      scale = rep(data$by, 3),
      value = c(data$power, data$ASNH1, data$ASTH1),
      measure = c(rep("Power", length(data$power)), rep("ASNH1", length(data$power)), rep("ASTH1", length(data$power))),
      method = rep(data$method, 3)
    )
  }

  return(data)
}
