## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(PersonAlytics)

doEval <- FALSE

depends <- function()
{
  d <- scan('../DESCRIPTION', what = 'character', sep='\n')
  wd <- which( grepl("*Depends*", d) )
  wr <- which( grepl("RoxygenNote*", d) )
  pkgs <- c(d[(wd+1):(wr-1)])
  pkgs <- gsub( " ", "", pkgs)
  pkgs <- gsub( "\t", "", pkgs)
  pkgs <- gsub( ",", "", pkgs)
  pkgs
}
pkgs <- depends()
pkgss <-  matrix( unlist( strsplit(pkgs, "\\(") ), ncol=2, byrow=TRUE)[,1]

## ----eval=doEval---------------------------------------------------------
#  install.packages('devtools')

## ----eval=doEval---------------------------------------------------------
#  devtools::install_github("https://github.com/ICTatRTI/PersonAlytics")

## ---- message=FALSE------------------------------------------------------
library(PersonAlytics)

## ---- eval=doEval--------------------------------------------------------
#  ?PersonAlytic

## ---- echo=FALSE---------------------------------------------------------
kable(head(PersonAlytics::OvaryICT))

## ---- eval=doEval, message=FALSE-----------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  time="Time")

## ---- eval=doEval--------------------------------------------------------
#  eg_nlme <- PersonAlytic(output="nlme_example",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  autoDetect=list())

## ---- eval=doEval--------------------------------------------------------
#  eg_gamlss <- PersonAlytic(output="gamlss_example",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  package="gamlss",
#                  autoDetect=list(),
#                  sigma.formula = ~ Targe1)

## ---- eval=doEval--------------------------------------------------------
#  class(eg_nlme)
#  class(eg_gamlss)

## ---- eval=doEval--------------------------------------------------------
#  summary(eg_nlme)

## ---- eval=doEval--------------------------------------------------------
#  summary(eg_gamlss)

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time")

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  ivs=list("Target1", "Target2", "Target3"))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  target_ivs=list("Target4", "Target5", "Target6"))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  ivs=list("Target1", "Target2", "Target3"),
#                  target_ivs=list("Target4", "Target5", "Target6"))

## ---- eval=FALSE---------------------------------------------------------
#  OvaryICT$Target1 <- factor(OvaryICT$Target1)

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  interactions=list(c("Target1", "Target2"), c("Target1", "Target3"),
#                           c("Target2", "Target3")))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  time_power=2)

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(output="MyResults",
#                  data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  correlation="ARMA(p=3,q=4)")

## ---- eval=doEval--------------------------------------------------------
#  t1.normal <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  autoDetect=NULL,
#                  family=NO(),
#                  package="gamlss")
#  plot(t1.normal)

## ---- eval=doEval--------------------------------------------------------
#  t1.skew <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  autoDetect=NULL,
#                  family=WEI2(),
#                  package="gamlss")
#  plot(t1.skew)
#  

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  subgroup=OvaryICT$Mare<6)

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  target_ivs=list("Target1", "Target2", "Target3"),
#                  standardize=list(dv=FALSE, ivs=TRUE, byids=FALSE))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  method="ML")

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  method="ML",
#                  package="gamlss",
#                  autoDetect = list(DIST=list()))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  method="ML",
#                  package="gamlss",
#                  autoDetect = list(TO=list(polyMax=4)))

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  ivs=list("Target1", "Target2", "Target3"),
#                  method="ML",
#                  package="gamlss",
#                  autoDetect = list(),
#                  sigma.formula = formula(.~Target1)
#                  )

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  target_ivs=list("Target1", "Target2", "Target3"),
#                  p.method='BY',
#                  alpha=.1
#                  )

## ---- eval=doEval--------------------------------------------------------
#  t1 <- PersonAlytic(data=OvaryICT,
#                  ids="Mare",
#                  dvs="follicles",
#                  phase="Phase",
#                  time="Time",
#                  ivs=list("Target1", "Target2", "Target3"),
#                  fpc=6000
#                  )

## ------------------------------------------------------------------------
parallel::detectCores()

