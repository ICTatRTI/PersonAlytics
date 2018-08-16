#' PersonAlytics: Analytics for single-case, small N, and Idiographic Clinical Trials
#'
#' The PersonAlytics package provides the simplified user interface
#' for implementing linear mixed effects models for idiographic clinical trials
#' (ICT) data, single case studies, and small N studies with intensive longitudinal
#' designs. Contact us via \link{https://personalytics.rti.org/} for
#' licensing options.
#'
#' @details
#' The basic mixed effects model is \eqn{dv=time+phase+phase*time}
#' with random intercepts and random slopes for time. The phase variable is optional.
#' Additional independent variables (or covariates) can be included.
#' The \code{PersonAlytics} package provides the simplified user interface
#' for implementing this model using \code{\link{gamlss}} or \code{\link{lme}}.
#'
#' Key features of the \code{PersonAlytics} package include:
#'
#' \strong{Automated detection of the residual covariance structure.}
#' \code{PersonAlytics} automates model comparisons for determining autocorrelation
#' structure for all patients or for each patient.
#'
#' \strong{Automated detection of the function form for the time variable.}
#' \code{PersonAlytics} automates model comparisons for determining the the functional
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
#' \strong{False Discovery Rate Adjustment.} When hight throughput options are requested,
#' Type I error correction and false discovery rate adjustments are implemented
#' post-implementation across target covariates (and individuals if requested) within
#' each outcome.
#'
#' \strong{Linear and Generalized Linear Mixed Effects Models.} Linear mixed effects models
#' can be fit in either the \code{\link{nlme}} framework or the \code{\link{gamlss}}
#' approach. The two approaches give nearly identicle fixed effects estimates but differ
#' in their computation of standard errors and random effects. Generalized linear
#' mixed effects models can be fit using the \code{\link{gamlss}} option
#' (see \code{\link{gamlss.family}}). The \code{\link{gamlss}} appoarch also allows models
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
#'                  dv="severity.migraine",
#'                  time="time",
#'                  time_power=3,
#'                  standardize=TRUE,
#'                  package="nlme")
#' summary(CurExMod)
#' # angriness predictor model
#' CurExModAng <- PersonAlytic(data=CurEx,
#'                  ids="id",
#'                  dv="severity.migraine",
#'                  time="time",
#'                  ivs=list('angriness'),
#'                  time_power=3,
#'                  standardize=TRUE,
#'                  package="nlme")
#' summary(CurExModAng)
#' # angriness predictor model, add correlation
#' CurExModAngAR1 <- PersonAlytic(data=CurEx,
#'                  ids="id",
#'                  dv="severity.migraine",
#'                  time="time",
#'                  ivs=list('angriness'),
#'                  time_power=3,
#'                  correlation="corARMA(p=1)",
#'                  standardize=TRUE,
#'                  package="nlme")
#' summary(CurExModAngAR1)
#'
#' # demonstrate how random effects modeling components affect bias and precision
#' # in eastimating the mean of a repeated measures variable with n=1
#' CurEx.mn <- mean(CurEx$severity.migraine)
#' CurEx.sd <- sd(CurEx$severity.migraine)/sqrt(nrow(CurEx))
#' CurEx.lm <- lm(severity.migraine~1, data=CurEx)
#' CurEx.arm <- ar(CurEx$severity.migraine)$x.mean
#' CurEx.ars <- ar(CurEx$severity.migraine)$var.pred
#' CurEx.ri <- nlme::lme(severity.migraine~1, data=CurEx, random=~1|id)
#' CurEx.rs <- nlme::lme(severity.migraine~1, data=CurEx, random=~time|id)
#' CurEx.rsar1 <- nlme::lme(severity.migraine~1, data=CurEx,
#'                          random=~time|id, correlation=nlme::corAR1())
#' CurEx.rst <- nlme::lme(severity.migraine~time, data=CurEx, random=~time|id)
#' CurEx.rstar1 <- nlme::lme(severity.migraine~time, data=CurEx, random=~time|id,
#'                           correlation=nlme::corAR1())
#' rownms <- c('Descriptive', 'Linear Model', 'Autoregressive Model',
#'             'Random Intercept', 'Random Intercept & Slope',
#'             'Random Intercept & Slope  & AR(1)',
#'             'Random Intercept & Slope & Fixed Effects for Time',
#'             'Random Intercept & Slope & AR(1) & Fixed Effects for Time')
#' outdf <- data.frame(rownms,
#'                     matrix(c(CurEx.mn, CurEx.sd,
#'                              summary(CurEx.lm)$coef[1:2],
#'                              CurEx.arm, CurEx.ars,
#'                              summary(CurEx.ri)$tTable[1:2],
#'                              summary(CurEx.rs)$tTable[1:2],
#'                              summary(CurEx.rsar1)$tTable[1:2],
#'                              summary(CurEx.rst)$tTable[1,1:2],
#'                              summary(CurEx.rstar1)$tTable[1,1:2]), ncol=2, byrow = TRUE)
#' )
#' names(outdf) <- c('Model', 'Mean', 'se')
#' outdf

#' Ovary data from nlme modified for ideographic clinical trial analysis
#'
#' @name OvaryICT
#' @docType data
#' @author Stephen Tueller \email{stueller@@rti.org}
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
#' # this example
#' t1 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dv="follicles",
#'                  phase="Phase",
#'                  time="Time")
#' summary(t1)
#' # is equivalent to
#' t2 <- PersonAlytic()
#' summary(t2)
#' # check
#' identical(t1,t2)
#' # NOTE. Time is a trigonometric function of Time, warnings are produced
#' # because Time is not ordered.
#' # definition of OvaryICT
#' \dontrun{
#' OvaryICT <- as.data.frame(nlme::Ovary)
#' OvaryICT$Mare <- as.numeric(OvaryICT$Mare)
#' OvaryICT$Phase <- as.numeric(OvaryICT$Time > .5)
#' OvaryICT$Time <- sin(2*pi*OvaryICT$Time)
#' OvaryICT <- OvaryICT[order(OvaryICT$Mare),]
#' }
#'
#' # output verification of t1
#' capture.output(t1.coef <- summary(t1), file = 'NUL')
#' all.equal(as.vector(t1.coef[,1]),
#' c(11.2421402, -1.9730995, 0.9126069, -1.9746526, 1.0967642))

if(1==2)
{
  OvaryICT <- as.data.frame(nlme::Ovary)
  OvaryICT$Mare <- as.numeric(OvaryICT$Mare)
  OvaryICT$Phase <- as.numeric(OvaryICT$Time > .5)
  OvaryICT$Time <- sin(2*pi*OvaryICT$Time)
  OvaryICT <- OvaryICT[order(OvaryICT$Mare),]
  devtools::use_data(OvaryICT, overwrite = TRUE)
}
