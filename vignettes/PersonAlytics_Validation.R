## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(nlme)

## ---- warning=FALSE------------------------------------------------------
t1palytic <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase')
t1palytic.mare1 <- t1palytic$lme(t1palytic$data$Mare==1)
summary(t1palytic.mare1)$tTable

## ---- warning=FALSE------------------------------------------------------
t1palytic.mare1 <- t1palytic$lme(t1palytic$data$Mare==1)
summary(t1palytic.mare1)$tTable

