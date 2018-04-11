### function to convert ML model from pact.auto to REML
MLtoREML <- function(mod, data)
{
  lmee(fixed=mod$call$fixed,
       data=data,
       random=mod$call$random,
       correlation=mod$call$correlation,
       method="REML")
}
