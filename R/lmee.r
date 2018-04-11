#' Function to fit lme() models with error handling and automatically refitting using optim if nlminb fails to converge.
#'
#'
#' @param fixed Formula. See \code{\link{lme}}.
#' @param data data.frame. See \code{\link{lme}}.
#' @param random Formula. See \code{\link{lme}}.
#' @param correlation Correlation structure class. \code{\link{lme}} and \code{\link{corClasses}}.
#' @param method Logical See ?\code{\link{lme}}.
#'
#' @return Either a \code{lme} object or a character message that the model did not converge.
#'
#' @author Stephen Tueller \email{stueller@rti.org}
#'
#'

# TODO
# - consider, when no converge, try dropping AR() or other model aspects

lmee <- function(fixed, data, random, correlation = NULL, method = "REML", ...)
{
  m1 <- try(nlme::lme(fixed=fixed,
                data=data,
                random=random,
                correlation=correlation,
                method=method,
                keep.data=FALSE,
                na.action=na.omit,
                ...), silent = TRUE)

  if( 'try-error'%in%class(m1) | !eds(m1) )
  {
    ctrl <- nlme::lmeControl(opt='optim')
    m1 <- try(nlme::lme(fixed=fixed,
                  data=data,
                  random=random,
                  correlation=correlation,
                  method=method,
                  control=ctrl,
                  keep.data=FALSE,
                  na.action=na.omit,
                  ...), silent = TRUE)
  }
  if( 'lme'%in%class(m1) )
  {
    m1$call$fixed <- fixed
    m1$call$random <- random
    return(m1)
  }
  else
  {
    return('Model did not converge')
  }
}

