rm(list = ls())

# RUN THIS BEFORE WORKING -------------------------------------------------

source("util.R")

# required packages
packages <- c(
  # "tidyverse",
  "MASS",
  "rpact",
  "mvtnorm",
  "greekLetters",
  "numDeriv",
  "crayon",
  "nleqslv",
  "parallel",
  "ggplot2",
  "tictoc",
  "ggpubr",
  "kableExtra"
)

# wrapper around require-install.packages call
.reqORinst <- function(package) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    require(package)
  }
}

lapply(packages, .reqORinst)

theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13, face = "italic"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.2, "cm"),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold")
  )


files1 <- paste0("functions/", list.files("./functions/", pattern = "*.R$"))
files2 <- paste0("designs/", list.files("./designs/", pattern = "*.R$"))

sapply(files1, source)
sapply(files2, source)
