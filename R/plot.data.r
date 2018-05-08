# this needs to be given and a class and a plot method (use s3), which also
# work with traj.plot() or turn them into one plot method

plot.data <- function(id, dv, time, n, fit.nlme, fit.lme4, phase, which.fit="lme4")
{
  ### why am i using lme4? does nlme work?
  if(!exists("which.fit")) which.fit="lme4"

  ### data
  # extract raw and fitted data
  dat      <- dat.nlme <- fit.nlme$data
  dat.lme4 <- model.frame(fit.lme4)
  # get fitted values
  fitted.nlme <- fit.nlme$"Individual Model Based"[,"fixed"] # options are "fixed" and id
  fitted.lme4 <- fitted(fit.lme4)     # all(fitted(fit.lme4) == predict(fit.lme4))
  # choose which fitted value to use
  if(which.fit=="nlme") fitted.val <- fitted.nlme
  if(which.fit=="lme4") fitted.val <- fitted.lme4

  ### create new data across the range of the time variable
  x.range <- range(dat[,time], na.rm=TRUE)
  t.dat   <- 0:(x.range[2]-1) # you need data on all x in X
  x.new   <- data.frame(t.dat, rep(NA, length(t.dat)), rep(NA, length(t.dat)))

  ### create the design matrix
  Designmat <- model.matrix(eval(fit.nlme$call$fixed)[-2], data=dat) #!# do I need the new dat line above??
  # remove reference columns to get conformable arguments for diag
  if(ncol(Designmat)!=dim(fit.nlme$varFix)[1])
  {
    Designmat <- Designmat[, colnames(Designmat)%in%colnames(fit.nlme$varFix) ]
  }

  ### get fitted data and standard errors
  #   http://glmm.wikidot.com/faq
  predvar <- diag(Designmat %*% fit.nlme$varFix %*% t(Designmat))
  SE  <- sqrt(predvar)                  # confidence interval
  SE2 <- sqrt(predvar+fit.nlme$sigma^2) # prediction interval
  dat <- data.frame(dat, SE=SE, SE2=SE2, predvar=predvar)

  ### set up list of matrices to store plotting data
  plotdat  <- vector("list", 0)

  ### make ID character
  dat[,id] <- as.character(dat[,id])

  ### populate plotdat with default colors, lwds, etc.
  #   Observed
  nr.dat <- nrow(dat)
  plotdat$"Individual Observed" <- data.frame(x=dat[,time],
                                              y=dat[,dv],
                                              group=dat[,id],
                                              leg=rep("Individual Observed", nr.dat),
                                              phase=dat[,phase])
  #   Fitted
  plotdat$"Individual Model Based"   <- data.frame(x=dat[,time],
                                                   y=fitted.val,
                                                   group=dat[,id],
                                                   leg=rep("Individual Model Based", nr.dat),
                                                   phase=dat[,phase])
  save(dat, plotdat, id, dv, time, n, phase,
       file="plotData.Rdata")
}
