### function to refit lme model in lme4
lmeTOlmer <- function(lme.mod, data, method='REML')
{
  lmer.frm <- as.formula( paste( c(
    paste(deparse(lme.mod$call$fixed), collapse=''),
    paste("(",substr(deparse(lme.mod$call$random), 2, 10000), ")")),
    collapse="+"))
  lmer.mod <- lmer(lmer.frm, data, REML=(method=='REML'))
  #!# need a warings catch here
  #warnings()
  return(lmer.mod)
}
