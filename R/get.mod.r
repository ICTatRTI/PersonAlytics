### function to extract model from lme object
get.mod <- function(mod)
{
  correlation=attr(summary(mod$call$correlation), "structName")
  if(is.null(correlation)) correlation = "NULL"
  c(
    fixed=paste(deparse(mod$call$fixed), collapse="") ,
    random=deparse(mod$call$random) ,
    correlation=correlation
  )
  ### put warnings not captured by lmee in log
  #pact.err(warnings(), save_model, "get.mod()")
}
