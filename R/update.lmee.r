# this needs to be changed to a method for the s3 generic update

### update lmee one variable x at a time
update.lmee <- function(mod, data, x, z=NULL)
{
  if( is.null(z)) formula.phase <- as.formula( paste(paste(deparse(mod$call$fixed),collapse=""), "+", x))
  isx <- any(attr(terms(mod$call$fixed), which="term.labels")==x)
  if(!is.null(z))
  {
    isz <- any(attr(terms(mod$call$fixed), which="term.labels")==z)
    if(!isx & !isz) formula.phase <- as.formula( paste(c(paste(deparse(mod$call$fixed), collapse=""), x, z, paste(x, "*", z, sep="")), collapse=" + "))
    if( isx & !isz) formula.phase <- as.formula( paste(c(paste(deparse(mod$call$fixed), collapse=""),    z, paste(x, "*", z, sep="")), collapse=" + "))
    if(!isx &  isz) formula.phase <- as.formula( paste(c(paste(deparse(mod$call$fixed), collapse=""), x,    paste(x, "*", z, sep="")), collapse=" + "))
    if( isx &  isz) formula.phase <- as.formula( paste(c(paste(deparse(mod$call$fixed), collapse=""),       paste(x, "*", z, sep="")), collapse=" + "))
  }
  mod.update <- lmee(fixed=formula.phase,
                     data=data,
                     random=mod$call$random,
                     correlation=mod$call$correlation,
                     method="ML")

  # put warnings not captured by lmee in log
  #pact.err(warnings(), save_model, "update.lmee()")

  # return updated model
  return(mod.update)
}
