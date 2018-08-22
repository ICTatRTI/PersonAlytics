## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(nlme)

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
#library(PersonAlytics)

## ------------------------------------------------------------------------
#' t1 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  time="Time")
#' # individual level analyses
#' w <- OvaryICT$Mare==1
#' t2 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  subgroup=w)
#' summary(t2)
#' # individual level analyses with standardized data
#' t3 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  subgroup=w,
#'                  standardize=TRUE)
#' summary(t3)
#' # interaction terms with nlme estimation
#' set.seed(123)
#' OvaryICT$iv1 <- factor(cut(runif(nrow(OvaryICT)), breaks=2), labels=c(0,1))
#' t4 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  interactions=list(c('Time', 'Phase'),
#'                                    c('Time', 'iv1'),
#'                                    c('Phase',   'iv1')),
#'                  time="Time",
#'                  package="nlme")
#' summary(t4)

