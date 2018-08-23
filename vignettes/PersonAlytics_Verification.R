## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(nlme)

## ---- message=FALSE------------------------------------------------------
library(PersonAlytics)

## ---- warning=FALSE------------------------------------------------------
t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                  time='Time', phase='Phase')

## ------------------------------------------------------------------------
all.equal(t1$fixed, formula("follicles ~ Time * Phase"))
all.equal(t1$random, formula("~Time | Mare"))
all.equal(t1$formula, formula(paste("follicles ~ Time * Phase +",
                                       "re(random = ~Time | Mare,",
                                        "method = 'REML', correlation = NULL)")))
all.equal(t1$correlation, NULL)

## ---- warning=FALSE, message=FALSE---------------------------------------
  t1.gamlss <- t1$gamlss() # Palytic gamlss
  t1.lme    <- t1$lme() # Palytic lme

## ---- message=FALSE, warning=FALSE---------------------------------------
  library(gamlss)
  t0.gamlss <- gamlss(follicles ~ Time * Phase +
                        re(random = ~Time | Mare, method = "REML"),
                      data = OvaryICT) # 'manual' gamlss
  t0.lme <- nlme::lme(t1$fixed, OvaryICT, t1$random) # 'manual' lme
  capture.output( t0.gamlss.table <- summary( t0.gamlss ), file='NUL')
  capture.output( t1.gamlss.table <- summary( t1.gamlss ), file='NUL')

## ------------------------------------------------------------------------
  all.equal(t0.gamlss.table, t1.gamlss.table, tolerance = .002)

## ------------------------------------------------------------------------
  all.equal(summary(t0.lme)$tTable, summary(t1.lme)$tTable)

## ------------------------------------------------------------------------
  all.equal(t1.gamlss.table[1:4,1], summary(t1.lme)$tTable[,1], tolerance = .005)

## ------------------------------------------------------------------------
  all.equal(as.vector( t0.gamlss.table[,1] ), 
            c(10.6898310,  -0.8637036,  10.9724586,  -8.6439015,   1.1633329))
  all.equal(as.vector( summary(t1.lme)$tTable[,1] ), 
            c(10.6616306, -0.8689801, 10.9720945, -8.6438502 ))

## ---- warning=FALSE------------------------------------------------------
  t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase')

## ------------------------------------------------------------------------
t1$correlation
t1$formula

## ------------------------------------------------------------------------
  t1$correlation <- "corARMA(p=1, q=0)"
  all.equal(t1$correlation, "nlme::corARMA(p=1, q=0)")

## ------------------------------------------------------------------------
  t1$formula

## ---- warning=FALSE------------------------------------------------------
t1simple <- PersonAlytic(data=OvaryICT,
                 ids="Mare",
                 dv="follicles",
                 phase="Phase",
                 time="Time")
capture.output(t1simple.summary <- summary(t1simple), file='NUL')

## ---- warning=FALSE------------------------------------------------------
t1palytic <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase')
t1palytic.summary <- summary(t1palytic$lme())
all.equal(t1simple.summary$tTable, t1palytic.summary$tTable)

