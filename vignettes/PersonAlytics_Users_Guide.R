## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(PersonAlytics)

OvaryICT <- PersonAlytics::OvaryICT

doEval <<- TRUE

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

## ----eval=FALSE----------------------------------------------------------
#  install.packages('devtools')

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("https://github.com/ICTatRTI/PersonAlytics")

## ---- message=FALSE------------------------------------------------------
library(PersonAlytics)

## ---- eval=FALSE---------------------------------------------------------
#  ?PersonAlytic

## ---- echo=FALSE---------------------------------------------------------
kable(head(PersonAlytics::OvaryICT), digits = 2)

## ---- eval=doEval, message=FALSE-----------------------------------------
eg_required <- PersonAlytic(data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                time="Time",
                autoSelect=NULL)

## ---- eval=doEval--------------------------------------------------------
eg_nlme <- PersonAlytic(output="nlme_example",
                data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                phase="Phase",
                time="Time",
                autoSelect=NULL)

## ---- eval=doEval--------------------------------------------------------
eg_gamlss <- PersonAlytic(output="gamlss_example",
                data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                phase="Phase",
                time="Time",
                package="gamlss",
                autoSelect=NULL)

## ---- eval=doEval--------------------------------------------------------
class(eg_nlme)
class(eg_gamlss)

## ---- eval=doEval, comment=NA--------------------------------------------
summary(eg_nlme)

## ---- eval=doEval, comment=NA--------------------------------------------
summary(eg_gamlss)

## ---- eval=doEval--------------------------------------------------------
OvaryICT$follicles2 <- OvaryICT$follicles^2
OvaryICT$folliclesr <- sqrt(OvaryICT$follicles2)
eg_htp <- PersonAlytic(output="htp_example",
                data=OvaryICT,
                ids="Mare",
                dvs=list("follicles", "follicles2", "folliclesr"),
                phase="Phase",
                time="Time",
                autoSelect=NULL)

## ---- eval=doEval, echo=FALSE, comment=NA--------------------------------
names(eg_htp)

