#' Personalytic, a simplified user interface for a linear mixed effects model
#' via a Palytic object
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#' @importFrom gamlss re
#'
#' @description A simplified user interface for initializing a Palytic object
#' and fitting a linear mixed effects model via \code{gamlss}.
#'
#' @param data A \code{\link{data.frame}} with the variables assinged to \code{ids},
#' \code{dv}, \code{time}, and optionally, \code{phase} and \code{ivs}. If left as
#' \code{NULL}, an example will run using the Ovary data from the \code{\link{nlme}}
#' package, see \code{\link{OvaryITC}}.
#' @param ids Character. Name of the ID variable. The ID variable must be numeric to ensure accurate
#' sorting.
#' @param dv Character. Name of dependent variable.
#' @param time Character. Name of the time variable.
#' @param phase Charcter. Name of the phase variable.
#' @param ivs Character list of covariates, e.g., \code{list('iv2', 'iv2')}.
#' @param interactions List of vector pairs of variable names for which interaction
#' terms should be specified, e.g., \code{list(c('time', 'phase'), c('time', 'iv1'),
#' c('iv1', 'iv2'))} where \code{'iv1'} is the name of a variable in the liste \code{ivs}.
#' @param time_power Numeric. Power of the time variable (e.g., \code{time^time_power}).
#' A quadratic or cubic growth model would be specified using
#' \code{time_power=2} or \code{time_power=3}, respectively. See also \code{maxOrder}.
#' @param correlation See \code{\link{corStruct}} in \code{\link{nlme}}. Must be passed as a
#' character, e.g. \code{"corARMA(p=1)"}.
#' @param family See \code{\link{gamlss.family}}. The default is normal.
#' @param subgroup Logical vector where \code{length(subgroup)==nrow(data)} indicating
#' which subset of the data should be used for analysis.
#' @param standardize Should the dependent and independent variables be standardized?
#' Does not apply to factor variables.
#' @param package Which package should be used? Options are \code{\link{gamlss}}
#' (the default) and \code{\link{nlme}} passed as character strings, e.g., \code{"gamlss"}.
#' @param maxAR Integer, the largest AR(p) value that should be considered in the
#' automated residual correlation structure search. See \code{\link{nlme::corARMA}}.
#' If \code{maxAR=>0} or \code{maxMA>0}, the automated residual correlation
#' structure search will be used. Otherwise, the value for \code{correlation} will
#' be used.
#' @parame maxMA Integer, the largest MA(q) value that should be considered in the
#' automated residual correlation structure search. See \code{\link{nlme::corARMA}}.
#' If \code{maxAR=>0} or \code{maxMA>0}, the automated residual correlation
#' structure search will be used. Otherwise, the value for \code{correlation} will
#' be used.
#' @param IC Wich information criterion should be used for selecting among ARMA(p,q)
#' models in the automated residual correlation structure search.
#' @param lrt logical. Should likelihood ration tests (lrt) be used? If \code{FALSE},
#' the smallest information criterion will be used (see \code{IC}.) The should be
#' false when both AR(p) and MA(q) models are considered since they are not nested
#' and the lrt requires nested models.
#' @param alpha value greater than 0 and less than 1. The Type I error rate for the lrt.
#' @maxOrder integer. The highest order of \code{time_power} for automatically
#' detecting the functional form of the relationship between time and the outcome.
#' If \code{maxOrder>0}, the automated search is conducted using likelihood ratio
#' tests, otherwise the specified value for \code{time_power} is used.
#'
#' @return \code{PersonAlytic} returns an object of class \code{\link{lme}} or
#' \code{\link{gamlss}}.
#'
#'
#' @examples
#' #
#' t1 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  time="Time")
#'
#' # for more examples, see the user's guide:
#' \dontrun{
#' vignette('PersonAlytics_Users_Guide', package='PersonAlytics')
#' }
#'

PersonAlytic <- function(data=NULL,
                         ids,
                         dv,
                         time,
                         phase=NULL,
                         ivs=NULL,
                         interactions=NULL,
                         time_power=1,
                         correlation=NULL,
                         family=gamlss.dist::NO(),
                         subgroup=NULL,
                         standardize=FALSE,
                         package=c('nlme', 'gamlss'),
                         maxAR=0,
                         maxMA=0,
                         IC=c("BIC", "AIC"),
                         lrt=FALSE,
                         alpha=.05,
                         maxOrder=0)
{
  # if no data are given, use a test data set
  if(is.null(data))
  {
    data   <- OvaryICT
    dv     <- "follicles"
    ids    <- "Mare"
    time   <- "Time"
    phase  <- "Phase"
    ivs    <- NULL
    interactions<- NULL
    time_power  <- 1
    correlation <- "nlme::corARMA(p=1)"
    subgroup    <- NULL
    standardize <- FALSE
    package     <- 'gamlss'
    maxAR=3
    maxMA=3
    IC=c("BIC", "AIC")
    lrt=FALSE
    alpha=.05
    maxOrder=7
  }

  if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))

  t1 <- Palytic$new(data=data[subgroup,]      ,
                    dv=dv                     ,
                    ids=ids                   ,
                    time=time                 ,
                    phase=phase               ,
                    ivs=ivs                   ,
                    interactions=interactions ,
                    time_power=time_power     ,
                    correlation=correlation   ,
                    standardize=standardize   )
  if(maxAR > 0 | maxMA > 0) t1$GroupAR_order(dv, maxAR, maxMA, IC, lrt, alpha)
  if(maxOrder>0) t1$GroupTime_Power(maxOrder)
  if(package=="gamlss") Grp.out <- t1$gamlss()
  if(package=="nlme")   Grp.out <- t1$lme()

  return(Grp.out)

}
