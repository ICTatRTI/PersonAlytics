## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(nlme)

## ----mods----------------------------------------------------------------
library(nlme)

# select the 1st mare
wm   <- Ovary$Mare==1

# mean and standard error of the mean, ignoring autocorrelation
f.mn <- mean(Ovary$follicles[wm])
f.se.mn <- sd(Ovary$follicles[wm])/sqrt(length(Ovary$follicles[wm]))

# naive linear model
f.lm <- lm(follicles ~ 1, data=Ovary[wm,])

# Autoregressive model
f.ts <- ar(ts(Ovary$follicles[wm]))
f.ts.mn <- f.ts$x.mean
f.ts.se.mn <- f.ts$var.pred/sqrt(length(Ovary$follicles[wm]))

# the best AR order is 1
f.ts$order
# but do we need ARIMA? No, AR(1) is still adequate
(f.arima <- forecast::auto.arima(Ovary$follicles[wm]))

# random intercept model
f.re.i <-lme(follicles ~ 1, data=Ovary[wm,], random = ~ 1 | Mare)

# random intercept model, AR(1)
f.re.i.ar1 <-lme(follicles ~ 1, data=Ovary[wm,], random = ~ 1 | Mare,
                 correlation = corAR1())

# random intercept and slope model
f.re.is <-lme(follicles ~ 1, data=Ovary[wm,], random = ~ Time | Mare)

# random intercept and slope model, AR(1)
f.re.is.ar1 <-lme(follicles ~ 1, data=Ovary[wm,], random = ~ Time | Mare,
                  correlation = corAR1())


# random intercept model, fixed effect for time
f.ret.i <-lme(follicles ~ Time, data=Ovary[wm,], random = ~ 1 | Mare)

# random intercept model, AR(1), fixed effect for time
f.ret.i.ar1 <-lme(follicles ~ Time, data=Ovary[wm,], random = ~ 1 | Mare,
                 correlation = corAR1())

# random intercept and slope model, fixed effect for time
f.ret.is <-lme(follicles ~ Time, data=Ovary[wm,], random = ~ Time | Mare)

# random intercept and slope model, AR(1), fixed effect for time
f.ret.is.ar1 <-lme(follicles ~ Time, data=Ovary[wm,], random = ~ Time | Mare,
                  correlation = corAR1())

# table the results
rownms <- c('Descr.', 'LM', 'AR1', 'RE.i', 'RE.i.AR1', 'RE.is', 'RE.is.AR',
            'RE.t.i', 'RE.t.i.AR1', 'RE.t.is', 'RE.t.is.AR1')
outmat <- 
data.frame(Model=rownms,
           matrix(c(f.mn, f.se.mn,
                    summary(f.lm)$coef[1:2],
                    f.ts.mn, f.ts.se.mn,
                    summary(f.re.i)$tTable[1:2],
                    summary(f.re.i.ar1)$tTable[1:2],
                    summary(f.re.is)$tTable[1:2],
                    summary(f.re.is.ar1)$tTable[1:2],
                    summary(f.ret.i)$tTable[1,1:2],
                    summary(f.ret.i.ar1)$tTable[1,1:2],
                    summary(f.ret.is)$tTable[1,1:2],
                    summary(f.ret.is.ar1)$tTable[1,1:2]), ncol=2, byrow=TRUE))
names(outmat) <- c('Model', 'Mean', 'se')


## ----echo=FALSE, results='asis'------------------------------------------
kable(outmat, caption='Table 1. Comparisons of estimates of the mean of a time series')

## ----echo=FALSE, results='asis'------------------------------------------
f.resids <- data.frame(f.ts = f.ts$resid,
                       f.lm = resid(f.lm),
                       f.re.i       = resid(f.re.i      ),
                       f.re.i.ar1   = resid(f.re.i.ar1  ),
                       f.re.is      = resid(f.re.is     ),
                       f.re.is.ar1  = resid(f.re.is.ar1 ),
                       f.ret.i      = resid(f.ret.i     ),
                       f.ret.i.ar1  = resid(f.ret.i.ar1 ),
                       f.ret.is     = resid(f.ret.is    ),
                       f.ret.is.ar1 = resid(f.ret.is.ar1))
pairs(f.resids)

## ----echo=FALSE, results='asis'------------------------------------------
kable(round(cor(f.resids, use = 'pairwise.complete.obs'),2), caption='Table 2. Correlations of residuals from the models in Table 1')

## ----echo=FALSE, results='asis'------------------------------------------
plot(Ovary$Time[wm], Ovary$follicles[wm], type='l')

## ----cubic---------------------------------------------------------------
f.ret3.is3.ar1 <-lme(follicles ~ Time + I(Time^2) + I(Time^3), 
                    data=Ovary[wm,], 
                    random = ~ Time + I(Time^2) + I(Time^3) | Mare,
                  correlation = corAR1())
summary(f.ret3.is3.ar1)$tTable

## ----comps---------------------------------------------------------------
# random slopes do not improve model fit beyond random intercepts
anova(f.re.i, f.re.is)
# 
f.ret3.is.ar1 <-lme(follicles ~ Time + I(Time^2) + I(Time^3), 
                    data=Ovary[wm,], 
                    random = ~ 1 | Mare,
                  correlation = corAR1())
anova(f.ret3.is.ar1, f.ret3.is3.ar1 )

