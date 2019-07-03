## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
papv <- packageVersion("PersonAlytics")



## ----eval=FALSE----------------------------------------------------------
#  install.deps <- function(pkg)
#  {
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    old.pkg <- pkg[ (pkg %in% installed.packages()[, "Package"])]
#    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE,
#                                          repos = "https://cran.rstudio.com/")
#    if (length(old.pkg)) update.packages(old.pkg, dependencies = TRUE)
#  }
#  packages <- c("gamlss", "ggplot2", "multilevel", "nlme", "R6")
#  install.deps(packages)

## ---- message=FALSE------------------------------------------------------
library(PersonAlytics)

## ------------------------------------------------------------------------
?PersonAlytic

