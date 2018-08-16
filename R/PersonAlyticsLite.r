#' PersonAlyticsLite: A package for Ideographic Clinical Trial analyses.
#'
#' The PersonAlyticsLite package provides the simplified user interface
#' \code{\link{PersonAlytic}}
#' for implementing linear mixed effects models for ideagraphic clinical trial
#' (ICT) data using \code{\link{gamlss}} fit using restricted maximum likelihood
#' (REML; ML is also an option).
#'
#' @details
#' The basic model for ICTs is \eqn{dv=time+phase+phase*time} with random intercepts
#' and random slopes for time. The phase variable is optional. Additional covariates
#' can be added.
#'
#' This is done via the \code{\link{Palytic}} object, which
#' offers additional options for advanced users including generalized linear mixed effects
#' models (see \code{\link{gamlss.family}}), user specific correlation structures, and
#' user specified random effects structures.
#'
#' The companion package \code{PersonAlyticsPro} implements automated section of
#' autoregressive structures and time polynomials via maximum likelihood (ML) model
#' comparisons. High throughput estimation of individual level models, iterations through
#' multiple dependent variables, and iterations through multiple independent variables
#' is also implemented. Contact us via \link{https://personalytics.rti.org/} for
#' licensing options.
#'
#' @section PersonAlytics functions
#' pending...
#'
#' @section The Palytic class
#' See \code{\link{Palytic}}
#'
#' @docType package
#' @name PersonAlyticsLite
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
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
NULL

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
NULL

makeData <- FALSE
if(makeData)
{
  OvaryICT <- as.data.frame(nlme::Ovary)
  OvaryICT$Mare <- as.numeric(OvaryICT$Mare)
  OvaryICT$Phase <- as.numeric(OvaryICT$Time > .5)
  OvaryICT$Time <- sin(2*pi*OvaryICT$Time)
  OvaryICT <- OvaryICT[order(OvaryICT$Mare),]
  devtools::use_data(OvaryICT, overwrite = TRUE)
}



