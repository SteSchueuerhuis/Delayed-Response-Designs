
# WRAPPER AROUND GGPLOT2 --------------------------------------------------

getPerformancePlotsByEffect <-
  function(method = c(
             "BindingHampsonJennison",
             "NonBindingHampsonJennison",
             "RecPauseDesign"
           ),
           alpha = 0.025,
           beta = 0.2,
           a_spend = c("asP", "asOF"),
           b_spend = c("bsP", "bsOF"),
           I1 = 0.5,
           I_d = 0.2,
           delta = seq(-1, 1, 0.05),
           stDev = 1,
           nMax = 200,
           rr = 10,
           # recruitment rate
           gamma_l = 0.8,
           gamma_u = 0.8,
           parallel = F,
           n.cores = detectCores() - 1) {
    I1_d <- I1 + I_d

    performance <- data.frame(
      effect = delta,
      power = rep(NA, length(delta)),
      ASNH1 = rep(NA, length(delta)),
      ASTH1 = rep(NA, length(delta)),
      design = rep(method, length(delta))
    )

    if (parallel) {
      cl <- makeCluster(n.cores)
      clusterEvalQ(cl, source("./init.R"))

      RES <- parSapply(
        cl = cl,
        X = delta,
        FUN = function(x) {
          analyseDesign(
            method = method,
            alpha = alpha,
            beta = beta,
            a_spend = a_spend,
            b_spend = b_spend,
            I1 = I1,
            I_d = I_d,
            delta = x,
            stDev = stDev,
            nMax = nMax,
            rr = rr,
            gamma_l = gamma_l,
            gamma_u = gamma_u
          )
        }
      )

      stopCluster(cl)
    } else {
      RES <- sapply(
        X = delta,
        FUN = function(x) {
          analyseDesign(
            method = method,
            alpha = alpha,
            beta = beta,
            a_spend = a_spend,
            b_spend = b_spend,
            I1 = I1,
            I_d = I_d,
            delta = x,
            stDev = stDev,
            nMax = nMax,
            rr = rr,
            gamma_l = gamma_l,
            gamma_u = gamma_u
          )
        }
      )
    }


    performance$power <- unlist(RES["power", ])
    performance$ASNH1 <- unlist(RES["ASNH1", ])
    performance$ASTH1 <- unlist(RES["ASTH1", ])

    p1 <- ggplot(performance, aes(x = delta, y = power)) +
      geom_line(size = 0.7, linetype = "dashed") +
      theme +
      # facet_wrap(~bounds, nrow = 2) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(
        limits = c(min(delta), max(delta)),
        breaks = seq(min(delta), max(delta), 0.25)
      ) +
      xlab(expression(delta)) +
      ylab("Rej. Prob.") +
      ggtitle(expression(paste("Power by ", delta))) +
      geom_hline(
        yintercept = alpha,
        linetype = "longdash",
        col = "black",
        size = 0.9
      ) +
      geom_hline(
        yintercept = 1,
        linetype = "longdash",
        size = 0.9,
        col = "black"
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "longdash",
        col = "black",
        size = 0.9
      ) +
      annotate(
        geom = "text",
        x = 0.5,
        y = alpha + 0.025,
        label = expression(paste(alpha, " = 0.025")),
        color = "red",
        cex = 4
      )

    p2 <- ggplot(performance, aes(x = effect, y = ASNH1)) +
      geom_line(size = 0.6, linetype = "dashed") +
      theme +
      scale_x_continuous(
        limits = c(min(delta), max(delta)),
        breaks = seq(min(delta), max(delta), 0.25)
      ) +
      scale_y_continuous(
        limits = c(I1 * nMax, nMax + 20),
        breaks = seq(floor(I1 * nMax) - 10, nMax + 10, 25)
      ) +
      xlab(expression(delta)) +
      ylab(expression(E[delta](N))) +
      ggtitle(expression(paste(E[delta](N), " by ", delta))) +
      geom_hline(
        yintercept = nMax,
        linetype = "longdash",
        size = 0.9,
        col = "black"
      ) +
      geom_hline(
        yintercept = (I1 + I_d) * nMax,
        linetype = "longdash",
        size = 0.9,
        col = "black"
      )

    p3 <- ggplot(performance, aes(x = effect, y = ASTH1)) +
      geom_line(size = 0.6, linetype = "dashed") +
      theme +
      scale_x_continuous(
        limits = c(min(delta), max(delta)),
        breaks = seq(min(delta), max(delta), 0.25)
      ) +
      scale_y_continuous(
        limits = c(floor(I1_d * nMax / rr) - 3, (nMax / rr) + 2 * (I_d * nMax / rr)),
        breaks = seq(floor(I1_d * nMax / rr) - 2, ceiling((nMax / rr) + 1 * (I_d * nMax / rr)), 2)
      ) +
      xlab(expression(delta)) +
      ylab(expression(E[delta](T))) +
      ggtitle(expression(paste(E[delta](T), " by ", delta))) +
      geom_hline(
        yintercept = floor(I1_d * nMax / rr),
        linetype = "longdash",
        size = 0.9,
        col = "black"
      )

    rowPlot <- annotate_figure(
      ggarrange(p1, p2, p3,
        ncol = 3, nrow = 1
      ),
      top = text_grob(
        bquote("Operating Characteristics by" ~ delta),
        face = "bold",
        size = 0.93
      )
    )

    output <- list(
      data = performance,
      powerPlot = p1,
      eNplot = p2,
      eTplot = p3,
      rowPlot = rowPlot
    )

    class(output) <- c("plotObject", class(output))

    return(output)
  }

.isPlotObject <- function(X) {
  return(class(X)[1] == "plotObject")
}


mergePlot <- function(X, Y, Z) {
  data <- mergeData(
    X = X,
    Y = Y,
    Z = Z,
    longformat = F
  )

  labels <-
    c(
      "Hampson & Jennison (Binding)",
      "Hampson & Jennison (Non-Binding)",
      "Delayed Response Pausing Design"
    )

  power <-
    ggplot(data, aes(x = scale, y = power, linetype = method, color = method)) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotted"),
      labels = labels
    ) +
    geom_line(size = 0.9) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(
      limits = c(min(data$scale), max(data$scale)),
      breaks = seq(min(data$scale), max(data$scale), 0.25)
    ) +
    scale_color_manual(values = c("seagreen4", "orange2", "blue"), labels = labels) +
    xlab(expression(delta)) +
    ylab("Rej. Prob.") +
    ggtitle("Power") +
    geom_hline(
      yintercept = 0.025,
      linetype = "longdash",
      col = "black",
      size = 0.7
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "longdash",
      size = 0.7,
      col = "black"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "longdash",
      col = "black",
      size = 0.7
    ) +
    annotate(
      geom = "text",
      x = 0.5,
      y = 0.025 + 0.03,
      label = expression(paste(alpha, " = 0.025")),
      color = "black",
      cex = 4
    ) +
    theme

  ASNH1 <-
    ggplot(data, aes(x = scale, y = ASNH1, linetype = method, color = method)) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotted"),
      labels = labels
    ) +
    geom_line(size = 0.9) +
    scale_x_continuous(
      limits = c(min(data$scale), max(data$scale)),
      breaks = seq(min(data$scale), max(data$scale), 0.25)
    ) +
    scale_y_continuous(
      limits = c(round(min(data$ASNH1)), 200 + 10),
      breaks = seq(round(min(data$ASNH1)) - 30, 200 + 10, 10)
    ) +
    scale_color_manual(values = c("seagreen4", "orange2", "blue"), labels = labels) +
    xlab(expression(delta)) +
    ylab(expression(E[delta](N))) +
    ggtitle(expression(paste(E[delta](N)))) +
    geom_hline(
      yintercept = 200,
      linetype = "longdash",
      size = 0.7,
      col = "black"
    ) +
    geom_hline(
      yintercept = min(round(data$ASNH1)),
      linetype = "longdash",
      size = 0.7,
      col = "black"
    ) +
    theme

  ASTH1 <-
    ggplot(data, aes(x = scale, y = ASTH1, linetype = method, color = method)) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotted"),
      labels = labels
    ) +
    xlab(expression(delta)) +
    geom_line(size = 0.9) +
    scale_x_continuous(
      limits = c(min(data$scale), max(data$scale)),
      breaks = seq(min(data$scale), max(data$scale), 0.25)
    ) +
    scale_y_continuous(
      limits = c(round(min(data$ASTH1)) - 1, ceiling(max(data$ASTH1)) + 2),
      breaks = seq(round(min(data$ASTH1)) - 1, ceiling(max(data$ASTH1)) + 2, 1)
    ) +
    scale_color_manual(values = c("seagreen4", "orange2", "blue"), labels = labels) +
    xlab(expression(delta)) +
    ylab(expression(E[delta](T))) +
    ggtitle(expression(E[delta](T))) +
    geom_hline(
      yintercept = min(round(data$ASTH1)),
      linetype = "longdash",
      size = 0.7,
      col = "black"
    ) +
    theme

  output <- list(
    power = power,
    ASNH1 = ASNH1,
    ASTH1 = ASTH1,
    mergedPlot = ggarrange(
      power,
      ASNH1,
      ASTH1,
      ncol = 3,
      common.legend = T,
      legend = "bottom"
    )
  )

  return(output)
}
