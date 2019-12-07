## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(nlme)
library(gamlss)
library(PersonAlytics)

## ------------------------------------------------------------------------
egaov <- aov(follicles ~ factor(Mare), data = Ovary)

eglme <- lme(follicles ~ 1, data = Ovary, random = ~ 1 | Mare,
                   method = 'ML')

eggamlss1 <- gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
                            data = Ovary,
                            sigma.formula = ~ 1)

eggamlss2 <- gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
                            data = Ovary,
                            sigma.formula = ~ 0)

# compare variance components
VarCorr(eglme)
VarCorr(getSmo(eggamlss1))
VarCorr(getSmo(eggamlss2))

# compare iccs
ICC(eglme)
ICC(eggamlss1)
ICC(eggamlss2)
multilevel::ICC1(egaov)

## ---- error=TRUE---------------------------------------------------------
# show the error
summary(eggamlss2)
# compare estimates
summary(eglme)$tTable
capture.output(eggamlss1.tt <- summary(eggamlss1), file='NUL'); eggamlss1.tt
capture.output(eggamlss2.tt <- summary(eggamlss2), file='NUL'); eggamlss2.tt

