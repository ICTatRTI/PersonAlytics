#' PersonAlytics: Analytics for single-case, small N, and Idiographic Clinical Trials
#'
#' The PersonAlytics package provides the simplified user interface
#' for implementing linear mixed effects models for idiographic clinical trials
#' (ICT) data, single case studies, and small N studies with intensive longitudinal
#' designs. Contact us via https://personalytics.rti.org/ for
#' licensing options.
#'
#' @details
#' The basic mixed effects model is \eqn{dv=time+phase+phase*time}
#' with random intercepts and random slopes for time. The phase variable is optional.
#' Additional independent variables (or covariates) can be included.
#' The \code{PersonAlytics} package provides the simplified user interface
#' for implementing this model using \code{\link{gamlss}} or \code{\link{lme}}. The
#' primary function of \code{PersonAlytics} is \code{\link{PersonAlytic}}.
#'
#' Key features of the \code{PersonAlytics} package include:
#'
#' \strong{Automated detection of the residual covariance structure.}
#' \code{PersonAlytics} automates model comparisons for determining autocorrelation
#' structure for all patients or for each patient.
#'
#' \strong{Automated detection of the function form for the time variable.}
#' \code{PersonAlytics} automates model comparisons for determining the functional
#' form of the relationship between time and the outcome
#' (i.e., linear vs. quadratic vs. cubic growth models) for all patients or for each patient.
#'
#' \strong{Estimation.} The automated covariance structure and function form for time is
#' done using maximum likelihood (ML) estimators. Final results are estimated using
#' restricted maximum likelihood (REML).
#'
#' \strong{High Throughput.} When users have a list of outcomes (dependent variables),
#' a list of target covariates, and/or
#' or desire the analyses to be repeated for each individual in the data set, high
#' throughput options automate the model fitting process.
#'
#' \strong{False Discovery Rate Adjustment.} When high throughput options are requested,
#' Type I error correction and false discovery rate adjustments are implemented
#' post-implementation across target covariates (and individuals if requested) within
#' each outcome.
#'
#' \strong{Linear and Generalized Linear Mixed Effects Models.} Linear mixed effects models
#' can be fit in either the \code{\link{nlme}} framework or the \code{\link{gamlss}}
#' approach. The two approaches give nearly identical fixed effects estimates but differ
#' in their computation of standard errors and random effects. Generalized linear
#' mixed effects models can be fit using the \code{\link{gamlss}} option
#' (see \code{\link{gamlss.family}}). The \code{\link{gamlss}}  appoarch also allows models
#' for dealing with heteroscedasticity implemented by including mixed effects models for
#' the variance.
#'
#' @section PersonAlytics functions
#' pending...
#'
#' @section The Palytic class
#' See \code{\link{Palytic}}
#'
#' @docType package
#' @name PersonAlytics
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
NULL

#' Curelator test data
#'
#' @name CurEx
#' @docType data
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' Example data from Curelator
#'
#' @examples
#' # no predictor model
#' CurExMod <- PersonAlytic(data=CurEx,
#'                  ids="id",
#'                  dvs="severity.migraine",
#'                  time="time",
#'                  time_power=3,
#'                  standardize=list(dv=FALSE, iv=TRUE, byids=TRUE),
#'                  package="nlme")
#' summary(CurExMod)
#'
#' # angriness predictor model
#' CurExModAng <- PersonAlytic(data=CurEx,
#'                  ids="id",
#'                  dvs="severity.migraine",
#'                  time="time",
#'                  ivs=list('angriness'),
#'                  time_power=3,
#'                  standardize=list(dv=FALSE, iv=TRUE, byids=TRUE),
#'                  package="nlme")
#' summary(CurExModAng)
#'
#' # angriness predictor model, add correlation
#' CurExModAngAR1 <- PersonAlytic(data=CurEx,
#'                  ids="id",
#'                  dvs="severity.migraine",
#'                  time="time",
#'                  ivs=list('angriness'),
#'                  time_power=3,
#'                  correlation="corARMA(p=1,q=0)",
#'                  standardize=list(dv=FALSE, iv=TRUE, byids=TRUE),
#'                  package="nlme")
#' summary(CurExModAngAR1)
#'
NULL

#' Ovary data from nlme modified for ideographic clinical trial analysis
#'
#' @name OvaryICT
#' @docType data
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
#'
#' The Ovary data set from the nlme package, modified as shown in the example.
#'
#' \itemize{
#'   \item Mare. The mare ID.
#'   \item Time. The time variable.
#'   \item follicles. See \code{\link{Ovary}}.
#'   \item Phase. A faux phase variable created at \code{Time > .5}.
#'   \item Time. Equal to \code{sin(2*pi*Ovary$Time)}. See Pinheiro & Bates (2000).
#' }
#'
#' @references Pinero, J., & Bates, D. (2000). Mixed-effects models in S and S-PLUS
#' (statistics and computing).
#'
#' @examples
#' # A simple mixed effects models using PersonAlytic and lme
#' t1 <- PersonAlytic(data=OvaryICT,
#'                    ids="Mare",
#'                    dvs="follicles",
#'                    phase="Phase",
#'                    time="Time",
#'                    package="nlme",
#'                    autoDetect=list(),
#'                    standardize = list(dvs=FALSE,ivs=FALSE,byids=FALSE),
#'                    alignPhase = 'none')
#' summary(t1)
#'
#' # Verify the PersonAlytic results to a direct call to lme
#' t2 <- lme(follicles ~ Time * Phase,
#'           data = OvaryICT,
#'           random = ~Time | Mare,
#'           method = "REML",
#'           control = t1$call$control)
#' summary(t2)
#'
#' # verification tests - are the results the same?
#' message( '\n\nAre fixed effects equal?\n', all.equal(summary(t1)$tTable, summary(t2)$tTable) )
#' message( '\n\nAre variance components equal?\n',all.equal(VarCorr(t1), VarCorr(t2)) )
#'
#' # verify estimates against known values
#' message('\n\nAre stored values replicated?\n',
#' all.equal( summary(t1)$tTable[,1],
#' c(`(Intercept)`  = 10.6616306,
#'    Time          = -0.8689801,
#'    Phase         = 10.9720945,
#'    `Time:Phase`  = -8.6438502) ) )
#'
#' # Illustrate a multiphase study, Phase alignment (To do: add visualization)
#' t3 <- PersonAlytic(data=OvaryICT,
#'                    ids="Mare",
#'                    dvs="follicles",
#'                    phase="Phase2",
#'                    time="Time",
#'                    package="nlme",
#'                    autoDetect=list(),
#'                    standardize = list(dvs=FALSE,ivs=FALSE,byids=FALSE),
#'                    alignPhase = 'align')
#' summary(t3)
#'
#' # definition of OvaryICT
#' \dontrun{
#' OvaryICT <- as.data.frame(nlme::Ovary)
#' OvaryICT$Mare <- as.numeric(OvaryICT$Mare)
#' OvaryICT$Phase <- as.numeric(OvaryICT$Time > .5)
#' OvaryICT$Phase2 <- factor( cut(OvaryICT$Time,
#'                                breaks = c(-Inf, .2, .6, 1, Inf)),
#'                            labels = 1:4 )
#' #OvaryICT$Time <- sin(2*pi*OvaryICT$Time)
#' set.seed(1)
#' t1 <- matrix( sample(1:4, 3*nrow(OvaryICT), replace = TRUE), ncol = 3)
#' t1 <- data.frame(t1)
#' names(t1) <- paste('Target', 1:3, sep='')
#' set.seed(2)
#' t2 <- data.frame( matrix( rnorm(3*nrow(OvaryICT)), ncol = 3 ) )
#' names(t2) <- paste('Target', 4:6, sep='')
#' OvaryICT <- data.frame(OvaryICT, t1, t2)
#' OvaryICT <- OvaryICT[order(OvaryICT$Mare),]
#' }
#'


if(1==2)
{
  OvaryICT <- as.data.frame(nlme::Ovary)
  OvaryICT$Mare <- as.numeric(OvaryICT$Mare)
  OvaryICT$Phase <- as.numeric(OvaryICT$Time > .5)
  OvaryICT$Phase2 <- factor( cut(OvaryICT$Time,
                                 breaks = c(-Inf, .2, .6, 1, Inf)),
                             labels = 1:4 )
  #OvaryICT$Time <- sin(2*pi*OvaryICT$Time)
  set.seed(1)
  t1 <- matrix( sample(1:4, 3*nrow(OvaryICT), replace = TRUE), ncol = 3)
  t1 <- data.frame(t1)
  names(t1) <- paste('Target', 1:3, sep='')
  t1 <- do.call(data.frame, lapply(t1, factor))
  set.seed(2)
  t2 <- data.frame( matrix( rnorm(3*nrow(OvaryICT)), ncol = 3 ) )
  names(t2) <- paste('Target', 4:6, sep='')
  OvaryICT <- data.frame(OvaryICT, t1, t2)
  OvaryICT <- OvaryICT[order(OvaryICT$Mare),]
  rm(t1, t2)
  devtools::use_data(OvaryICT, overwrite = TRUE)
}
