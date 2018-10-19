#' \code{Palytic} class generator.
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @import R6
#' @import nlme
#' @importFrom nlme corAR1
#' @importFrom nlme corARMA
#' @importFrom nlme corCompSymm
#' @importFrom nlme corSymm
#' @import gamlss
#' @importFrom gamlss re
#'
#' @export
#' @format A \code{Palytic} generator
#' @keywords data
#'
#' @usage Palytic(data, ids, dv, time)
#'
#' @details
#' The fields \code{data}, \code{ids}, \code{dv}, and \code{time} are required.
#' Using these, the default model \eqn{dv=time} with random intercepts for \code{ids}
#' and random intercepts for \code{time} is constructed. See the example. If
#' \code{phase} is provided, the default model is \eqn{dv=time+phase+phase*time}, and if
#' \code{ivs} are provided they are included in the model.
#'
#' @field data A \code{\link{data.frame}} that contains as variables \code{ids},
#' \code{dv}, \code{phase}, and \code{time}. Optionally, additional independent
#' variables can be included in \code{ivs}. \code{fixed} and \code{random} formulae
#' for \code{\link{lme}} models and \code{formula} for \code{\link{gamlss}} models are
#' automatically generated when a \code{Palytic} object is created if these fields
#' are left \code{NULL}.
#' @field ids A character string giving the name of the id variable in \code{data}.
#' @field dv A character string giving the name of the dependent variable in \code{data}.
#' @field time A character string giving the name of the time variable in \code{data}.
#' Random slopes for time are inclued by default. This can be overridden by specifying
#' \code{fixed} and \code{random} formula for \code{\link{lme}} models or by specifying
#' the \code{formula} for \code{\link{gamlss}} models.
#' @field phase A character string giving the name of the phase variable in \code{data}.
#' The \code{phase*time} interaction is included by default. This can be overridden by
#' specifying \code{fixed} and \code{random} formula for \code{\link{lme}} models or by
#' specifying the \code{formula} for \code{\link{gamlss}} models.
#' @field ivs A \code{\link{list}} of one or more character strings giving the names
#' of additional variables in \code{data}, e.g., \code{list('iv2', 'iv2')}.
#' @field interactions List of vector pairs of variable names for which interaction
#' terms should be specified, e.g., \code{list(c('time', 'phase'), c('time', 'iv1'),
#' c('iv1', 'iv2'))} where \code{'iv1'} is the name of a variable in the liste \code{ivs}.
#' @field time_power The polynomial for \code{time}, e.g., \code{time^time_power}. Fixed
#' effects for \code{time^1...time^time_power} will be included in models. Future
#' releases will allow for other functions of time such as \code{\link{sin}}, but these
#' can be applied directly by transforming the \code{time} variable.
#' @field correlation See \code{\link{corStruct}}. Defaults to \code{NULL}, see
#' \code{\link{lme}}. Used by both \code{\link{lme}} and \code{\link{gamlss}} models.
#' @field family The \code{\link{gamlss.family}} distribution.
#' @field fixed The \code{fixed} effects model for \code{\link{lme}} models.
#' @field random The \code{random} effects model for \code{\link{lme}} models.
#' @field formula The \code{formula} effects model for \code{\link{gamlss}} models.
#' \code{sigma.formula}, \code{nu.formula}, and \code{tau.formula} will be implemented in
#' a future release.
#' @field method See \code{method} in \code{\link{lme}}. Is usef for both \code{\link{lme}}
#' and \code{\link{gamlss}} models.
#' @field standardize Named logical list. Which variables should be standardized? The default
#' is \code{list(dv=FALSE, ivs=FALSE, byids=FALSE)}. See \code{dv} and \code{ivs}. The option
#' \code{byids} controls whether standardization is done by individuals or by group. Any time
#' variables are changed (e.g., \code{ivs}), the data are subset, or the options in
#' \code{standardize} are changed, the raw data will be restandardized (see \code{datac}).
#' @field corStructs Vector. A \code{correlation} structure for each case in \code{ids}. Not
#' user accesible. Populated by \code{\link{PersonAlytic}}.
#' @field time_powers Vector. A \code{time_order} for each case in \code{ids}. Not
#' user accesible. Populated by \code{\link{PersonAlytic}}.
#' @field datac data.frame. Cleaned data. Cleaning involves the following steps:
#' 1. Check that the variables in \code{ids}, \code{dv}, \code{time}, \code{phase},
#' \code{ivs}, and \code{interactions} are in \code{data}.
#' 2. The \code{ids} variable is forced to be numeric.
#' 3. The validity of the \code{correlation} structure is checked, see \code{\link{corStruct}}.
#' 4. check that variables have non-zero variance.
#' 5. If standardization is requested, standardize the data (see \code{standardize}).
#' 6. Sort the data on \code{ids} and \code{time}.
#' 7. If patients have < 2 observations, they are dropped from the data set.
#' 8. If phase alignment is requested, align phases (see \code{alignPhases}).
#' @field warnings A list of warnings that will be populated as methods are called on a
#' \code{Palytic} object.
#' @field errors A list of errors that will be populated as methods are called on a
#' \code{Palytic} object.
#' @field try_silent Logical flag for testing error handling in \code{Palytic} methods.
#'
#' @section Methods:
#' \describe{
#'
#'   \item{\code{summary}}{This method provides a summary of the inputs, the cleaned data,
#'   and the raw data.}
#'
#'   \item{\code{describe}}{This method gives the correlation between \code{dv} and
#'   each continuous variable in \code{ivs} (as well as the \code{time} variable), or,
#'   if variablse are factors (including \code{phase}), the mean of \code{dv} is given
#'   for each factor level of each variable.}
#'
#'   \item{\code{lme(subgroup = NULL, dropVars = NULL)}}{This method fits the
#'   linear mixed effects
#'    \code{lme} model implied by the \code{Palytic} fields \code{ids},
#'    \code{dv}, \code{phase}, \code{time}, and optionally \code{ivs},
#'    \code{time_power}
#'    \code{correlation}. The default formula can be overridden by specifying
#'    \code{fixed} and \code{random} and optionally \code{method} and
#'    \code{correlation}. The see \code{\link{lme}} for the parameter
#'    \code{subgroup}. The \code{dropVars} parameter indicates which fixed
#'    effects should be dropped for a likelihood ration test (LRT). This is
#'    used by \code{\link{PersonAlytic}} to test \code{target_ivs}.}
#'
#'   \item{\code{gamlss(subgroup = NULL, sigma.formula = ~1, family = NULL,
#'   dropVars = NULL)}}{This method fits the \code{lme} model implied
#'   by the \code{Palytic} fields \code{ids}, \code{dv}, \code{phase}, \code{time}
#'   and optionally \code{ivs}, \code{time_power}, \code{correlation}. The default
#'   formula can be overridden by specifying \code{formula} and
#'   optionally \code{method}, \code{correlation}, and \code{family} (which can be
#'   used to specify generalized linear mixed effects models,
#'   see \code{\link{gamlss.family}}).
#'   The parameter \code{subgroup} operates as in \code{\link{lme}}. The parameter
#'   \code{sigma.formula} and \code{family} are desribed in \code{\link{gammlss}}.
#'    The \code{dropVars} parameter indicates which fixed
#'    effects should be dropped for a likelihood ration test (LRT). This is
#'    used by \code{\link{PersonAlytic}} to test \code{target_ivs}.}
#'
#'   \item{\code{arma(subgroup=NULL, max.p=3, max.q=3, dropVars=NULL,
#'   max.P=0, max.Q=0, max.d=0, max.D=0)}}{For individual level models,
#'   random intercepts and
#'   random slopes are not defined. In this situatation, an \code{ARMA(p,q)}
#'   should be used. This is implemented using the \code{xreg}
#'   option of \code{\link{arima}}. The residual correlation search is achieved using
#'   \code{\link{auto.arima}}.  The \code{dropVars} parameter indicates which fixed
#'    effects (in \code{xreg}) should be dropped for a likelihood ration test
#'    (LRT). This is used by \code{\link{PersonAlytic}} to test \code{target_ivs}.
#'   The other parameters are as in \code{\link{auto.arima}}.}
#'
#'   \item{\code{getAR_order(P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05)}}{
#'   --- slated for deprication, superceded by \code{arma}. ---
#'   This method automates the task of determining the correlation structure for
#'   each case in \code{ids} (see \code{\link{PersonAlytic} or
#'   \code{\link{PersonAlytic}}}). \code{P} and \code{Q} set the highest
#'   autoregressive and moving average parameters to be tested. If the time
#'   variable is approximatetly equally spaced, \code{whichIC} is the criterion
#'   used for determining the correlation structure for each \code{ids} using
#'   the \code{\link{auto.arima}} function. If the time variable is unequally
#'   spaced, \code{whichIC} as also the criterion for model selection via mixed
#'   effects models using \code{\link{lme}} if \code{lrt=FALSE}. If
#'   \code{lrt=TRUE} likelihood ratios are used via the \code{\link{anova}}
#'   method for \code{\link{lme}} objects. This is NOT reccomended unless Q=0
#'   and only AR models are considered since AR and MA models are not nested.
#'   Calling \code{getAR_order} populates the \code{corStructs} field of a
#'   \code{Palytic} object. For usage, see the examples.}
#'
#'   \item{\code{GroupAR_order(dv, P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05)}}{The
#'   same as \code{getAR_order} when the ARMA order is desired for the full sample.}
#'   \item{\code{getTime_Power(subset, maxOrder)}}{This method automates the task of
#'   determining  \code{time_power} for each case in \code{ids}
#'   (see \code{\link{PersonAlytic} or \code{\link{PersonAlytic}}}). For example,
#'   if \code{getTime_Power} returns \code{time_power=3},
#'   then \code{time + time^2 + time^3}
#'   will be added to the fixed effects of the model.
#'   Calling \code{getTime_Power} populates the
#'   \code{GroupTime_power} field of a \code{Palytic} object. For usage, see the examples.}
#'
#'   \item{\code{groupTime_Power(subset, maxOrder)}}{The same as \code{getTime_power} when
#'   the polynomial of time is desired for the full sample.}
#' }
#'
#' @examples
#' # construct a new Payltic object and examine the default formulae
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
#'                   time='Time', phase='Phase')
#' t1$fixed
#' t1$random
#' t1$formula
#' t1$correlation
#'
#' # summary function
#' t1$summary()
#'
#' # descriptive statistics
#' t1$describe()
#'
#' # compare gamlss and lme output, in which the models of the default formulae are fit
#' (t1.gamlss <- summary(t1$gamlss()))
#' (t1.lme    <- summary(t1$lme()))
#'
#' # parameter estimates are equal within 0.01
#' message('\n\nAre `lme` and `gamlss` results equivalent?\n',
#' all.equal(t1.gamlss[1:4,1], t1.lme$tTable[,1], tolerance = 0.01) )
#'
#' # now change the correlation structure and compare gamlss and lme output,
#' # noting that the intercepts are very different now
#' t1$correlation <- "corARMA(p=1, q=0)"
#' summary(t1$gamlss())
#' summary(t1$lme())
#'
#' # fit the model only to the first mare with ML instead of REML
#' t1$method <- 'ML'
#' summary(t1$gamlss(OvaryICT$Mare==1))
#'
#' # change the formula
#' t2 <- t1$clone()
#' t2$formula <- formula(follicles ~ Time * Phase +
#'                       re(random = ~Time + I(Time^2) | Mare, method = "ML",
#'                       correlation = corARMA(p=1,q=0)))
#' t2$formula
#'
#' # random intercept only model
#' t2 <- t1$clone()
#' t2$random <- formula(~1|Mare)
#' t2$random
#' t2$formula
#'
#' # random slope only model
#' t2$random <- formula(~0+Time|Mare)
#' t2$random
#' t2$formula
#'
#' # automatically select the polynomial order of time with getTime_Power
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare',
#'                   dv='follicles', time='Time', phase='Phase')
#' t1$getTime_Power()
#' t1$time_powers
#'
#' # automatically select the ARMA model for residual correlation getAR_order
#' # this runs slow because a model is selected for each individual
#' \dontrun{
#' t1$getAR_order()
#' t1$corStructs
#' }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# start of Palytic function ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Palytic <- R6::R6Class("Palytic",
                       private = list(
                         .data        = NULL, # consider pass by reference environment
                         .ids         = NULL,
                         .dv          = NULL,
                         .time        = NULL,
                         .phase       = NULL,
                         .ivs         = NULL,
                         .interactions= NULL,
                         .time_power  = NULL,
                         .correlation = NULL,
                         .family      = NULL,
                         .fixed       = NULL,
                         .random      = NULL,
                         .formula     = NULL,
                         .method      = NULL,
                         .standardize = list(dv=FALSE, iv=FALSE, byids=FALSE),
                         .corStructs  = NULL,
                         .time_powers = NULL,
                         .alignPhase  = FALSE,
                         .ismonotone  = NULL,
                         .datac       = NULL,
                         .warnings    = list(),
                         .errors      = list(),
                         .try_silent  = TRUE
                       ),


                       #active = .active(), # depricate until we can find work around, low priority
                       # this may be achievable by adding a utility package that lodas first
                       # and using @importFrom in the preamble
                       active =   list(
                         data = function(value)
                         {
                           if( missing(value) ){ private$.data }
                           else
                           {
                             stop("`$data` is read only", call. = FALSE)
                           }
                         },

                         ids = function(value)
                         {
                           if( missing(value) ) private$.ids
                           else
                           {
                             if(! is.character(value) )
                             {
                               stop("`ids` must be a character variable name in data")
                             }
                             if( is.null(private$.data[[value]]) )
                             {
                               stop( paste(value, "is not in the data") )
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = value,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         # use functions in the way of iscorStruct to update all formula any time
                         # any formula related object is changed
                         dv = function(value)
                         {
                           if( missing(value) ) private$.dv
                           else
                           {
                             if(! is.character(value) )
                             {
                               stop("`dv` must be a character variable name in data")
                             }
                             if( is.null(private$.data[[value]]) )
                             {
                               stop( paste(value, "is not in the data") )
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = value,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         time = function(value)
                         {
                           if( missing(value) ) private$.time
                           else
                           {
                             if(! is.character(value) )
                             {
                               stop("`time` must be a character variable name in the data")
                             }
                             if( is.null(private$.data[[value]]) )
                             {
                               stop( paste(value, "is not in the data") )
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = value,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             private$.ismonotone   <- monotone(private$.ids,
                                                               private$.time,
                                                               private$.data)
                             self
                           }
                         },

                         phase = function(value)
                         {
                           if( missing(value) ) private$.phase
                           else
                           {
                             if(! is.null(value))
                             {
                               if(! is.character(value) )
                               {
                                 stop("`phase` must be a character variable name in the data")
                               }
                               if( is.null(private$.data[[value]]) )
                               {
                                 stop( paste(value, "is not in the data") )
                               }
                               if( ! var(self$datac[[value]], na.rm = TRUE) > 0 )
                               {
                                 stop( paste(value, "has zero variance") )
                               }
                             }

                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = value,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         ivs = function(value)
                         {
                           if( missing(value) ) private$.ivs
                           else
                           {
                             if( ! any( c('list', 'character', "NULL") %in%
                                      is(value) ) )
                             {
                               stop("`ivs` must be a character list of ",
                                    "variables in the data or `NULL`.")
                             }
                             if( ! all(value %in% names(private$.data) ) )
                             {
                               nov <- value[ which(! value %in% names(private$.data )) ]
                               if(length(nov)==1)
                               {
                                 stop( paste(nov, "is not in the data") )
                               }
                               if(length(nov)>=2)
                               {
                                 stop( paste(paste(nov, collapse=", "), "are not in the data") )
                               }
                             }
                             frms <- forms(private$.data        ,
                                           PalyticObj   = self  ,
                                           ids          = NULL  ,
                                           dv           = NULL  ,
                                           time         = NULL  ,
                                           phase        = NULL  ,
                                           ivs          = value ,
                                           interactions = NULL  ,
                                           time_power   = NULL  ,
                                           correlation  = NULL  ,
                                           family       = NULL  ,
                                           fixed        = NULL  ,
                                           random       = NULL  ,
                                           formula      = NULL  ,
                                           method       = NULL  )
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         interactions = function(value)
                         {
                           if( missing(value) ) private$.interactions
                           else
                           {
                             if(! is.list(value) )
                             {
                               stop("`interactions` must be a list of vector pairs of variables in the data")
                             }
                             if( ! all(unlist(value) %in% names(private$.data)) & !is.null(value) )
                             {
                               nov <- value[ which(! value %in% names(private$.data )) ]
                               if(length(nov)==1) stop( paste(nov, "is not in the data") )
                               if(length(nov)>=2)
                               {
                                 stop( paste(paste(nov, collapse=", "), "are not in the data") )
                               }
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = value,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         time_power = function(value)
                         {
                           if( missing(value) ) private$.time_power
                           else
                           {
                             if(! is.numeric(value) ) stop("`time_power` must be numeric")
                             if( round(value, 0) != value | value < 1 )
                             {
                               stop("`time_power` must be a positive whole number")
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = value,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         correlation = function(value)
                         {
                           if( missing(value) ) private$.correlation
                           else
                           {
                             fixcor <- function(x)
                             {
                               if(!is.null(x))
                               {
                                 if( x != "NULL" )
                                 {
                                   x <- unlist( strsplit(x, "::") )
                                   return( paste("nlme", x[length(x)], sep="::" ) )
                                 }
                                 if( x == "NULL" )
                                 {
                                   return(NULL)
                                 }
                               }
                               if( is.null(x)  )
                               {
                                 return(x)
                               }
                             }
                             #iscor <- try( iscorStruct(fixcor(value)), TRUE )
                             value <- fixcor(value)
                             if(! iscorStruct(value) )
                             {
                               stop( paste("`correlation` must be of class `corStruct`.",
                                           "See `?nlme::corS.truct`",
                                           "For example, for AR(1) use 'corAR1()' or 'corARMA(p=1)'.",
                                           "For are 2 or higher use 'corARMA(p=2)'.",
                                           "For ARMA(3,1) use 'corARMA(p=3,q=1)'.",
                                           sep="\n") )
                             }
                             # NULL is a valid `value` for correlation, so we need one more
                             # parameter `corFromPalyticObj` to get it right
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = value,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL,
                                           corFromPalyticObj = FALSE)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         family = function(value)
                         {
                           if( missing(value) ) private$.family
                           else
                           {
                             if(! "gamlss.family" %in% class(value) )
                             {
                               stop("`family` is not in gamlss.family, see `?gamlss.dist::gamlss.family`")
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = value,
                                           fixed        = NULL,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         fixed = function(value)
                         {
                           if( missing(value) ) private$.fixed
                           else
                           {
                             if(! "formula" %in% class(value) )
                             {
                               stop("`fixed` must be a formula, see `?formula` and `??nlme::lme`")
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = value,
                                           random       = NULL,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         random = function(value)
                         {
                           if( missing(value) ) private$.random
                           else
                           {
                             if(! "formula" %in% class(value) )
                             {
                               stop("`random` must be a formula, see `?formula` and `??nlme::lme`")
                             }
                             frms <- forms(private$.data,
                                           PalyticObj   = self,
                                           ids          = NULL,
                                           dv           = NULL,
                                           time         = NULL,
                                           phase        = NULL,
                                           ivs          = NULL,
                                           interactions = NULL,
                                           time_power   = NULL,
                                           correlation  = NULL,
                                           family       = NULL,
                                           fixed        = NULL,
                                           random       = value,
                                           formula      = NULL,
                                           method       = NULL)
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         formula = function(value)
                         {
                           if( missing(value) ) private$.formula
                           else
                           {
                             if(! "formula" %in% class(value) )
                             {
                               stop("`formula` must be a formula, see `?formula` and `??gamlss::gamlss`")
                             }
                             suppressWarnings(
                               frms <- forms(private$.data,
                                             PalyticObj   = NULL,
                                             ids          = NULL,
                                             dv           = NULL,
                                             time         = NULL,
                                             phase        = NULL,
                                             ivs          = NULL,
                                             interactions = NULL,
                                             time_power   = NULL,
                                             correlation  = NULL,
                                             family       = NULL,
                                             fixed        = NULL,
                                             random       = NULL,
                                             formula      = value,
                                             method       = NULL))
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         method = function(value)
                         {
                           if( missing(value) ) private$.method
                           else
                           {
                             if(!any( c('ML', 'REML') %in% value ))
                             {
                               stop("`method` should be `ML` or `REML`")
                             }
                             suppressWarnings(
                               frms <- forms(private$.data,
                                             PalyticObj   = self,
                                             ids          = NULL,
                                             dv           = NULL,
                                             time         = NULL,
                                             phase        = NULL,
                                             ivs          = NULL,
                                             interactions = NULL,
                                             time_power   = NULL,
                                             correlation  = NULL,
                                             family       = NULL,
                                             fixed        = NULL,
                                             random       = NULL,
                                             formula      = NULL,
                                             method       = value))
                             private$.ids          <- frms$ids
                             private$.dv           <- frms$dv
                             private$.time         <- frms$time
                             private$.phase        <- frms$phase
                             private$.ivs          <- frms$ivs
                             private$.interactions <- frms$interactions
                             private$.time_power   <- frms$time_power
                             private$.correlation  <- frms$correlation
                             private$.family       <- frms$family
                             private$.fixed        <- frms$fixed
                             private$.random       <- frms$random
                             private$.formula      <- frms$formula
                             private$.method       <- frms$method
                             self
                           }
                         },

                         standardize = function(value)
                         {
                           if( missing(value) ) private$.standardize
                           else
                           {
                             if(! is.logical( value ) )
                             {
                               stop("`standardize` should be `TRUE` or `FALSE`")
                             }
                             private$.standardize <- value
                             self
                           }
                         },

                         corStructs = function(value)
                         {
                           if( missing(value) ) private$.corStructs
                           else
                           {
                             private$.corStructs <- value
                             self
                           }
                         },

                         time_powers = function(value)
                         {
                           if( missing(value) ) private$.time_powers
                           else
                           {
                             private$.time_powers <- value
                             self
                           }
                         },

                         monotone = function(value)
                         {
                           if( missing(value) ) private$.monotone
                           else
                           {
                             stop("`monotone` is read only", call. = FALSE)
                           }
                         },

                         alignPhase = function(value)
                         {
                           if( missing(value) ) private$.alignPhase
                           else
                           {
                             private$.alignPhase <- value
                             self
                           }
                         },

                         datac = function(value)
                         {
                           if( missing(value) ) private$.datac
                           else
                           {
                             stop("`$datac (cleaned data) is read only", call. = FALSE)
                           }
                         },

                         try_silent = function(value)
                         {
                           if( missing(value) ) private$.try_silent
                           else
                           {
                             if(! "logical" %in% class(value)) stop("`try_silent` should be logical")
                             private$.try_silent <- value
                             self
                           }
                         }

                         ### need to add active bindings for warnings, errors


                       ),

                       #########################################################
                       # initialize ####
                       #########################################################
                       public = list(
                         initialize = function
                         (
                           data        ,
                           ids         = NULL,
                           dv          = NULL,
                           time        = NULL,
                           phase       = NULL,
                           ivs         = NULL,
                           interactions= NULL,
                           time_power  = NULL,
                           correlation = NULL,
                           family      = gamlss.dist::NO(),
                           fixed       = NULL,
                           random      = NULL,
                           formula     = NULL,
                           method      = "REML",
                           standardize = list(dv=FALSE, iv=FALSE, byids=FALSE),
                           corStructs  = NULL,
                           time_powers = NULL,
                           ismonotone  = NULL,
                           alignPhase  = FALSE,
                           datac       = NULL,
                           warnings    = list(), # hide? or just make them read only?
                           errors      = list(),
                           try_silent  = TRUE
                         )
                         {
                           ### consider adding option to read a file, could autodetect file type
                           #if( is.character(data) )
                           #{
                           #  read
                           #}

                           # checks that get used multiple times
                           is.min <- !(is.null(ids) | is.null(c) | is.null(time))
                           is.lme <- !(is.null(fixed) | is.null(random))
                           is.frm <- !is.null(formula)

                           # check whether there is sufficient information
                           if( !any(is.min, is.lme, is.frm) )
                           {
                             stop( paste('You must provide one of the following\n',
                                         '1. `ids`, `dv`, and `time`\n',
                                         '2. `fixed` and `random`\n',
                                         '3. `formula`')
                             )
                           }

                           # clean the data
                           datac <- clean(data, ids, dv, time, phase, ivs,
                                         fixed, random, formula, correlation, family,
                                         dvs=NULL, target_ivs=NULL, standardize,
                                         sortData=TRUE, alignPhase)

                           # create the formulae
                           frms <- forms(datac                     ,
                                         PalyticObj=NULL           ,
                                         ids=ids                   ,
                                         dv=dv                     ,
                                         time=time                 ,
                                         phase=phase               ,
                                         ivs=ivs                   ,
                                         interactions=interactions ,
                                         time_power=time_power     ,
                                         correlation=correlation   ,
                                         family=family             ,
                                         fixed=NULL                ,
                                         random=NULL               ,
                                         formula=NULL              ,
                                         method=method             )

                           # check whether time is monotorically increasing
                           ismonotone <- monotone(ids, time, datac)

                           # populate private
                           private$.data        <- data
                           private$.ids         <- ids
                           private$.dv          <- dv
                           private$.time        <- time
                           private$.phase       <- phase
                           private$.ivs         <- ivs
                           private$.interactions<- interactions
                           private$.time_power  <- time_power
                           private$.correlation <- correlation
                           private$.family      <- family
                           private$.fixed       <- frms$fixed
                           private$.random      <- frms$random
                           private$.formula     <- frms$formula
                           private$.method      <- method
                           private$.standardize <- standardize
                           private$.corStructs  <- corStructs
                           private$.time_powers <- time_powers
                           private$.ismonotone  <- ismonotone
                           private$.alignPhase  <- alignPhase
                           private$.datac       <- datac
                           private$.warnings    <- warnings
                           private$.errors      <- errors
                           private$.try_silent  <- try_silent

                         }

                       )

)

# add methods
Palytic$set("public", "summary",
            function()
            {
              variables <- all.vars( self$formula )
              dropvars  <- grep(pattern = '\\(', variables)
              if(length(dropvars) > 0 ) variables <- variables[-dropvars]

              tempCorrelation <- ifelse(is.null(self$correlation), "NULL",
                                        self$correlation)

              # return the summary
              list(      ids          = self$ids            ,
                         dv           = self$dv             ,
                         time         = self$time           ,
                         phase        = self$phase          ,
                         ivs          = self$ivs            ,
                         interactions = self$interactions   ,
                         time_power   = self$time_power     ,
                         correlation  = tempCorrelation     ,
                         family       = self$family         ,
                         package      = self$package        ,
                         fixed        = self$fixed          ,
                         random       = self$random         ,
                         formula      = self$formula        ,
                         method       = self$method         ,
                         standardize  = self$standardize    ,
                         corStructs   = self$corStructs     ,
                         time_powers  = self$time_powers    ,
                         ismonotone   = self$ismonotone     ,
                         warnings     = self$warnings       ,
                         errors       = self$errors         ,
                         try_silent   = self$try_silent     ,
                         datac        = summary(self$datac[,variables]) ,
                         data         = summary(self$data[,variables])  )
            })

Palytic$set("public", "describe",
            function(subgroup=NULL)
            {
              # concatenate all the ivs with double checks for dropped terms
              ivall <- all.vars(self$fixed)[-1]

              if(! self$time[1] %in% ivall ) ivall <- c(ivall, self$time[1])
              if(length(self$phase) > 0)
              {
                if(! self$phase %in% ivall ) ivall <- c(ivall, self$phase)
              }

              # subset the data
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- self$data[subgroup, c(self$dv, ivall)]

              # get descriptive statistics for each variable from the raw data
              ivstats <- list()
              for(i in ivall)
              {
                if(is.factor(tempData[[i]]) | length(table(tempData[[i]]))==2)
                {
                  mns <- aggregate(tempData[[self$dv]], list(tempData[[i]]),
                                   mean, na.rm=TRUE)
                  nms <- paste('mean of', self$dv, 'for',
                               paste(i, mns[,1],sep='='),
                               sep=' ')
                  for(j in 1:nrow(mns)) ivstats[[nms[j]]] <- unname( mns$x[j] )
                }
                else
                {
                  des <- paste('correlation of', self$dv, 'and', i, sep=' ')

                  hasVarDV <- sd(tempData[[self$dv]], na.rm=TRUE)!=0
                  hasVarIV <- sd(tempData[[i]]      , na.rm=TRUE)!=0

                  if(is.na(hasVarDV)) hasVarDV <- FALSE
                  if(is.na(hasVarIV)) hasVarIV <- FALSE

                  if(hasVarDV & hasVarIV)
                  {
                    ivstats[[des]] <- cor(tempData[[self$dv]], tempData[[i]],
                                          use = 'pairwise.complete.obs')
                  }
                  else
                  {
                    ivstats[[des]] <- paste("No variance in",
                                            ifelse(!hasVarDV, self$dv, ""),
                                            ifelse(!hasVarIV, i      , ""),
                                            sep=": ")
                  }
                }
              }

              # resturcture the results
              statName = as.list(  names(ivstats) )
              names(statName) <- paste('statName', 1:length(statName), sep='')

              statValue = ivstats
              names(statValue) <- paste('statValue', 1:length(statValue), sep='')

              idx <- order(c(seq_along(statName), seq_along(statValue)))
              ivstats <- c(statName, statValue)[idx]

              # return
              return( ivstats )
            })

Palytic$set("public", "arma",
            function(subgroup=NULL, max.p=3, max.q=3, dropVars=NULL,
                     max.P=0, max.Q=0, max.d=0, max.D=0, ...)
            {
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- data.frame( na.omit( subset(self$datac, subgroup,
                                 all.vars(self$formula)) ) )

              # check that only one participant is in the data
              if( length(unique(tempData[[self$ids]])) != 1 & nrow(tempData) > 0 )
              {
                stop('`arma` requires n=1 data')
              }

              # if the variances are zero, set escape
              # note that this isn't needed by lme/gamlss b/c they fail
              # to converge, whereas arima will 'converge' but coeftest
              # fails
              zerovar <- FALSE
              if( any( unlist( lapply(tempData[,names(tempData) != self$ids],
                              function(x) all(duplicated(x)[-1L])) ) ) )
              {
                zerovar <- TRUE
              }

              # `xreg` requires pre-constructed interaction terms, here we
              # 1. create the model.matrix RHS from self$fixed
              # 2. drop the intercept column
              xdat <- model.matrix(self$fixed, tempData)[,-1]
              if(is.null(dim(xdat))) xdat <- matrix(xdat)
              if(ncol(xdat)==1) colnames(xdat) <- all.vars(self$fixed)[-1]

              # auto detect residual correlation structure here, time power
              # must be detected elsewhere or added here
              m1 <- try( forecast::auto.arima(y     = tempData[[self$dv]],
                                         xreg  = xdat,
                                         max.p = max.p,
                                         max.q = max.q,
                                         max.P = max.P,
                                         max.Q = max.Q,
                                         max.d = max.d,
                                         max.D = max.D,
                                         ...), silent = TRUE)
              # TODO() refitting if there are convergence issues

              # if dropVars is provided, fit the sub-model and estimate LRT
              wasLRTrun <- FALSE
              lrtp <- as.numeric(NA)
              if(!is.null(dropVars))
              {
                xdat0 <- xdat[,-(! colnames(xdat) %in% dropVars)]
                m0 <- try( forecast::auto.arima(y     = tempData[[self$dv]],
                                              xreg  = xdat0,
                                              max.p = max.p,
                                              max.q = max.q,
                                              max.P = max.P,
                                              max.Q = max.Q,
                                              max.d = max.d,
                                              max.D = max.D,
                                              ...), silent = TRUE)
                if( any("ARIMA" %in% class(m0) ) &
                    any("ARIMA" %in% class(m1) ) )
                {
                  l0 <- logLik(m0)
                  l1 <- logLik(m1)
                  df0 <- strsplit( unlist( strsplit(capture.output(l0), "=") )[2] , ")")
                  df1 <- strsplit( unlist( strsplit(capture.output(l1), "=") )[2] , ")")
                  df0 <- as.numeric( unlist(df0) )
                  df1 <- as.numeric( unlist(df1) )
                  lrtest <- as.numeric(2*(l1-l0))
                  lrtp <- pchisq(lrtest, df1-df0, lower.tail = FALSE)
                  wasLRTrun <- TRUE
                }
              }

              # output
              if( zerovar ) m1 <- ""
              if(  any("ARIMA" %in% class(m1)) ) tTable = lmtest::coeftest(m1)
              if(! any("ARIMA" %in% class(m1)) )
              {
                m1 <- "Model did not converge"
                tTable = NA
              }

              m1 <- list(arima = m1, tTable = tTable,
                         PalyticSummary = self$summary(),
                         xregs = colnames(xdat),
                         lrt = list(wasLRTrun=wasLRTrun, lrtp=lrtp) )
              return(m1)
            },
            overwrite = TRUE
)

#TODO(Stephen) you can't just add ... to lme() and get things passed, you'll
#need to add some sort of evaluation, e.g.
# > args <- list(...)
# > if("contrasts" %in% names(args)) # then pass `contrasts` to lme
# or some other method of matching arguments, see ?match.arg
Palytic$set("public", "lme",
            function(subgroup=NULL, dropVars=NULL, PQ=c(3,3), ...)
            {
			        if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- na.omit(subset(self$datac, subgroup,
                                 all.vars(self$formula)))
              # github issue #1
              cor <- eval(parse(text = ifelse(!is.null(self$correlation),
                                              self$correlation,
                                              'NULL')))
              wm <- 1
              m1 <- try(nlme::lme(fixed=self$fixed,
                                  data=tempData,
                                  random=self$random,
                                  correlation=cor,
                                  method=self$method),
                        silent = TRUE)


              ctrl <- nlme::lmeControl()
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 2
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }
              # try without correlation structure
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 3
                ctrl <- nlme::lmeControl(opt="optim")
                self$correlation <- "NULL" # not updating
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }
              # try without random slopes or correlation structure
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 4
                newformula <- forms(data       = self$datac ,
                                    PalyticObj = self      ,
                                    dropTime   = TRUE      )
                self$random <- newformula$random
                self$correlation <- NULL
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }

              # if n=1, detect the correlation structure here using getARnEQ1()
              # as it may affect the lrt (we should probably do the same for
              # group ar...);
              if( length(table(tempData[[self$ids]]))==1 )
              {
                cor  <- getARnEQ1(m1, PQ, self$dv)
                ctrl <- nlme::lmeControl(opt="optim")
                m1   <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }

			  	    cat("Palytic$lme:",
				      "\ntempData: ", names(tempData),
				      "\ndim(tempData): ", toString(dim(tempData)),
				      "\ncor: ", toString(self$correlation),
				      "\nfixed: ", toString(self$fixed),
				      "\nrandom: ", toString(self$random),
				      "\nm1: ", class(m1),
				      "\n\n", file="last_lme.txt", append=FALSE)

              # clean up the call - may not need this
              m1 <- cleanCall(modelResult=m1, PalyticObj=self)

              # lrt
              wasLRTrun <- FALSE
              lrtp <- as.numeric(NA)
              if( "lme" %in% class(m1) & !is.null(dropVars) )
              {
                frm0 <- all.vars(self$fixed)[-1]
                frm0 <- frm0[! frm0 %in% dropVars]
                frm0 <- formula( paste(self$dv, '~', paste(frm0, collapse = '+')) )
                #m0   <- update(m1, frm0) # this fails b/c the call is messed up
                m0 <- try(nlme::lme(fixed=frm0,
                                    data=tempData,
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method),
                          silent = TRUE)
                if("lme" %in% class(m1) & "lme" %in% class(m1))
                {
                  lrt <- anova(m1, cleanCall(m0, self))
                  if(nrow(lrt)==2 & "p-value" %in% names(lrt))
                  {
                    lrtp <- lrt$"p-value"[2]
                  }
                }
                wasLRTrun <- TRUE
              }

              # return
              if( "lme" %in% class(m1) )
              {
                m1$PalyticSummary  <- self$summary()
                m1$PalyticSummary$correlation <- cor
                m1$whichPalyticMod <- paste('Palytic lme model #', wm)
                m1$lrt <- list(wasLRTrun=wasLRTrun, lrtp=lrtp)
                return(m1)
              }
              else
              {
                return("Model did not converge")
                # consider replacing with
                #attr(m1, 'condition')$message
                # for the actual convergence error message
              }
            },
            overwrite = TRUE
)

#' cleanCall - clean up the call in in palytic objects
#' @keywords internal
cleanCall <- function(modelResult, PalyticObj)
{
  if("lme" %in% class(modelResult) )
  {
    modelResult$call$fixed       <- PalyticObj$fixed
    modelResult$call$random      <- PalyticObj$random
    modelResult$call$correlation <- PalyticObj$correlation
    modelResult$call$method      <- PalyticObj$method
  }
  if( "gamlss" %in% class(modelResult) )
  {
    modelResult$call$formula       <- PalyticObj$formula
    modelResult$call$sigma.formula <- PalyticObj$sigma.formula
    modelResult$call$family        <- PalyticObj$family
  }
  return(modelResult)
}

#' getARnEQ1
#' @keywords internal
getARnEQ1 <- function(m, PQ, dv)
{
  correlation <- "NULL"
  faa <- try( forecast::auto.arima(y     = m$residuals[,2],
                                   xreg  = m$data[[dv]],
                                   max.p = PQ[1],
                                   max.q = PQ[2],
                                   max.P = 0,
                                   max.Q = 0,
                                   max.d = 0,
                                   max.D = 0),
              silent = TRUE)
  if( "Arima" %in% class(faa) )
  {
    faaOrder <- forecast::arimaorder(faa)
    correlation <- paste('corARMA(p=', faaOrder[1],
                         ',q=', faaOrder[2], ')', sep='')
    correlation <- eval(parse(text = ifelse(!is.null(correlation),
                             correlation,
                             'NULL')))
  }
  return( correlation )
}

### TODO(Stephen) add residual correlation search for n=1
Palytic$set("public", "gamlss",
            function(subgroup=NULL, sigma.formula = ~1, family=NULL,
                     dropVars=NULL, ...)
            {
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- na.omit( subset(self$datac, subgroup,
                                          all.vars(self$formula)) )

              # allow for family to be changed on the fly
              currentFamily <- self$family
              if(!is.null(family)) currentFamily <- family

              wm <- 1 # default model
              ctrl <- gamlss::gamlss.control()
              m1 <- try(gamlss::gamlss(formula = self$formula,
                                       sigma.formula = sigma.formula,
                                       data = tempData,
                                       family = currentFamily),
                        silent = self$try_silent)

              # default model with increased n.cyc
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 2
                ctrl <- gamlss::gamlss.control(n.cyc=100)
                m1 <- try(refit(gamlss::gamlss(formula = self$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }
              if( "try-error" %in% class(m1) )
              {
                wm <- 3 # drop the random slope(s)
                newformula <- forms(data = self$datac ,
                                    PalyticObj = self ,
                                    dropTime = TRUE ,
                                    family = currentFamily)
                self$formula <- newformula$formula
                #self$family <- newformula$family
                ctrl <- gamlss::gamlss.control(n.cyc=100)
                m1 <- try(refit(gamlss::gamlss(formula = self$formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }

              # lrt
              # TODO: this is untested
              # TODO: this is the same in lme and gamlss, move to separate function
              wasLRTrun <- FALSE
              lrtp <- as.numeric(NA)
              if( "lme" %in% class(m1) & !is.null(dropVars) )
              {
                frm0 <- all.vars(self$fixed)[-1]
                frm0 <- frm0[! frm0 %in% dropVars]
                frm0 <- formula( paste(self$dv, '~', paste(frm0, collapse = '+')) )
                m0 <- update(m1, frm0)
                lrt <- anova(m1, m0)
                if(nrow(lrt)==2 & "p-value" %in% names(lrt))
                {
                  lrtp <- lrt$"p-value"[2]
                }
                wasLRTrun <- TRUE
              }

              # output
              if("gamlss" %in% class(m1))
              {
                m1$PalyticSummary <- self$summary()
                m1$whichPalyticMod <- paste('Palytic gamlss model #', wm)
                m1$lrt <- list(wasLRTrun=wasLRTrun, lrtp=lrtp)
                return(m1)
              }
              if("try-error" %in% class(m1))
              {
                return("Model did not converge")
              }

            },
            overwrite = TRUE
)


# whichIC can take AIC or BIC
Palytic$set("public", "GroupAR_order",
            function(P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05,
                       subgroup=NULL)
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "residual correlation structure starting...")
              start <- Sys.time()

              nullMod <- self$lme(subgroup)

              if( "lme" %in% class(nullMod) )
              {
                require(foreach)

                DIMS <- expand.grid(p=0:P, q=0:Q)
                DIMS <- DIMS[-1,]

                funcs <- c("mean")
                cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
                snow::clusterExport(cl, funcs)
                doSNOW::registerDoSNOW(cl)
                pkgs  <- c("gamlss", "nlme")
                pb <- txtProgressBar(max = P, style = 3)
                progress <- function(n) setTxtProgressBar(pb, n)
                opts <- list(progress = progress)

                corMods <- foreach(p=DIMS$p, q=DIMS$q, .packages = pkgs,
                                   .options.snow = opts,
                                   .export = c("self"))  %dopar%
                {
                  clone             <- self$clone(deep=TRUE)
                  clone$method      <- "ML"
                  clone$correlation <- "NULL"

                  cortemp <- paste("nlme::corARMA(p=", p, ",
                                 q=", q, ")", sep="")
                  cortemp <- gsub('\n', '', cortemp)
                  cortemp <- gsub(' ', '', cortemp)
                  clone$correlation <- cortemp
                  corMod <- clone$lme(subgroup)
                  if( class(corMod) != "lme" )
                  {
                    corMod <- NULL
                  }

                  cat('\n\n')

                  return(corMod)
                }
                parallel::stopCluster(cl)

                corMods <- corMods[!unlist(lapply(corMods, is.null))]

                names(corMods) <- unlist( lapply(corMods,
                                                 function(x) x$PalyticSummary$correlation) )
                corMods <- corMods[!is.na(names(corMods))]

                if(!lrt)
                {
                  if(whichIC=="AIC") ICs <- data.frame( unlist( lapply(corMods, AIC) ) )
                  if(whichIC=="BIC") ICs <- data.frame( unlist( lapply(corMods, BIC) ) )
                  else( stop( paste(whichIC, "is not a valid value for `whichIC`,",
                                    "use AIC or BIC.")))
                  ICs <- rbind("NULL" = ifelse(whichIC=="AIC", AIC(nullMod),
                                               BIC(nullMod)), ICs)
                  bestCor <- c("NULL", names(corMods))[which.min(ICs[,1])]
                }
              }
              if(! "lme" %in% class(nullMod) )
              {
                bestCor <- "NULL"
              }

              message("\nAutomatic detection of the residual\n",
                      "correlation structure took: ",
                      capture.output(Sys.time() - start), ".\n\n")
            },
            overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
Palytic$set("public", "getTime_Power",
            function(maxOrder=3, whichIC="BIC")
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              uid <- sort( as.numeric( unique(self$datac[[self$ids]]) ) )
              time_powers <- list()

              for(id in uid)
              {
                aics <- list()
                for(i in 1:maxOrder)
                {
                  clone <- self$clone(deep=TRUE)

                  clone$method     <- "ML"
                  clone$time_power <- i

                  mod0 <- clone$lme(self$datac[[self$ids]]==id)
                  if("lme" %in% class(mod0))
                  {
                    aics[[i]] <- AIC( mod0 )
                  }
                  else aics[[i]] <- NA

                }
                aics <- unlist(aics)
                time_powers[[id]] <- ifelse( all( is.na(aics) ),
                                             1,
                                             which.min( aics ) )
              }
              self$time_powers <- data.frame(ids=uid, time_power=unlist(time_powers))

              message("\nAutomatic detection of the time/outcome\n",
                      "relationship took: ",
                      capture.output(Sys.time() - start), ".\n\n")
            },
            overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
Palytic$set("public", "GroupTime_Power",
            function(subgroup=NULL, maxOrder=3, whichIC="BIC")
            {

              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))

              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              mods <- list()
              for(i in 1:maxOrder)
              {
                clone <- self$clone(deep=TRUE)

                clone$method     <- "ML"
                clone$time_power <- i

                mods[[i]] <- clone$lme(subgroup)
                if(! "lme" %in% class(mods[[i]]))
                {
                  mods[[i]] <- NULL
                }
                rm(clone)
              }
              if(whichIC=="BIC") bestMods <- unlist( lapply(mods, BIC) )
              if(whichIC=="AIC") bestMods <- unlist( lapply(mods, AIC) )
              self$time_power <- which.min( bestMods )

              message("\nAutomatic detection of the time/outcome\n",
                      "relationship took: ",
                      capture.output(Sys.time() - start), ".\n\n")

            },
            overwrite = TRUE)


