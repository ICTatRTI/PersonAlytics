## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(PersonAlytics)

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
#  install.deps <- function(pkg)
#  {
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    old.pkg <- pkg[ (pkg %in% installed.packages()[, "Package"])]
#    if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE,
#                                          repos = "https://cran.rstudio.com/")
#    if (length(old.pkg)) update.packages(old.pkg, dependencies = TRUE)
#  }

## ---- echo=FALSE, comment=NA---------------------------------------------
cat(paste("install.deps(c(", paste("'", pkgss[-1], "'", collapse=', ', sep=""), "))", 
          sep=""), "\n")

## ---- message=FALSE------------------------------------------------------
library(PersonAlytics)

## ---- eval=FALSE---------------------------------------------------------
#  ?PersonAlytic

## ---- echo=FALSE---------------------------------------------------------
kable(head(PersonAlytics::OvaryICT))

## ------------------------------------------------------------------------
t1 <- PersonAlytic(data=PersonAlytics::OvaryICT,
                ids="Mare",
                dvs="follicles",
                phase="Phase",
                time="Time")

## ------------------------------------------------------------------------
class(t1)

