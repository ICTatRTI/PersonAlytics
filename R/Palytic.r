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
#' @field standardize Logical. Should \code{y} and \code{ivs} be standardized?
#' @field corStructs Vector. A \code{correlation} structure for each case in \code{ids}. Not
#' user accesible. Populated by \code{\link{PersonAlytic}}.
#' @field time_powers Vector. A \code{time_order} for each case in \code{ids}. Not
#' user accesible. Populated by \code{\link{PersonAlytic}}.
#' @field is_clean Logical flag for data cleaning. Not yet implemented.
#' @field warnings A list of warnings that will be populated as methods are called on a
#' \code{Palytic} object.
#' @field errors A list of errors that will be populated as methods are called on a
#' \code{Palytic} object.
#' @field try_silent Logical flag for testing error handling in \code{Palytic} methods.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{lme(subgroup = NULL)}}{This method fits the linear mixed effects
#'    \code{lme} model implied by the \code{Palytic} fields \code{ids},
#'    \code{dv}, \code{phase}, \code{time}, and optionally \code{ivs}, \code{time_power}
#'    \code{correlation}. The default formula can be overridden by specifying
#'    \code{fixed} and \code{random} and optionally \code{method} and \code{correlation}.
#'    The see \code{\link{lme}}
#'   for the parameter \code{subgroup}. }
#'   \item{\code{gamlss(subgroup = NULL)}}{This method fits the \code{lme} model implied
#'   by the \code{Palytic} fields \code{ids}, \code{dv}, \code{phase}, \code{time}
#'   and optionally \code{ivs}, \code{time_power}, \code{correlation}. The default
#'   formula can be overridden by specifying \code{formula} and
#'   optionally \code{method}, \code{correlation}, and \code{family} (which can be
#'   used to specify generalized linear mixed effects models,
#'   see \code{\link{gamlss.family}}).
#'   The parameter \code{subgroup} operates as in \code{\link{lme}}.}
#'   \item{\code{getAR_order(P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05)}}{
#'   This method automates the task of determining the correlation structure for each case in
#'   \code{ids} (see \code{\link{PersonAlytic} or \code{\link{PersonAlytic}}}).
#'   \code{P} and \code{Q} set the highest autoregressive and moving
#'   average parameters to be tested. If the time variable is approximatetly equally spaced,
#'   \code{whichIC} is the criterion used for determining the correlation structure for each
#'   \code{ids} using the \code{\link{auto.arima}} function. If the time variable is unequally
#'   spaced, \code{whichIC} as also the criterion for
#'   model selection via mixed effects models using \code{\link{lme}} if \code{lrt=FALSE}.
#'   If \code{lrt=TRUE} likelihood ratios are used via the \code{\link{anova}}
#'   method for \code{\link{lme}} objects. This is NOT reccomended unless Q=0 and only AR
#'   models are considered since AR and MA models are not nested. Calling \code{getAR_order}
#'   populates the
#'   \code{corStructs} field of a \code{Palytic} object. For usage, see the examples.}
#'   \item{\code{GroupAR_order(dv, P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05)}}{The
#'   same as \code{getAR_order} when the ARMA order is desired for the full sample.}
#'   \item{\code{getTime_Power(subset, maxOrder)}}{This method automates the task of
#'   determining  \code{time_power} for each case in \code{ids}
#'   (see \code{\link{PersonAlytic} or \code{\link{PersonAlytic}}}). For example,
#'   if \code{getTime_Power} returns \code{time_power=3}, then \code{time + time^2 + time^3}
#'   will be added to the fixed effects of the model.
#'   Calling \code{getTime_Power} populates the
#'   \code{GroupTime_power} field of a \code{Palytic} object. For usage, see the examples.}
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
#' # compare gamlss and lme output, in which the models of the default formulae are fit
#' (t1.gamlss <- summary(t1$gamlss()))
#' (t1.lme    <- summary(t1$lme()))
#' # parameter estimates are equal within 0.01
#' all.equal(t1.gamlss[1:4,1], t1.lme$tTable[,1], tolerance = 0.01)
#' # now change the correlation structure and compare gamlss and lme output,
#' # noting that the intercepts are very different now
#' t1$correlation <- "corARMA(p=1, q=0)"
#' summary(t1$gamlss())
#' summary(t1$lme())
#' # fit the model only to the first mare with ML instead of REML
#' t1$method <- 'ML'
#' summary(t1$gamlss(OvaryICT$Mare==1))
#'
#' # change the formula (note limitations with the interaciton)
#' t2 <- t1$clone()
#' t2$formula <- formula(follicles ~ Time * Phase +
#'                       re(random = ~Time + I(Time^2) | Mare, method = "ML",
#'                       correlation = corARMA(p=1,q=0)))
#' t2$formula
#'
#' # getTime_Power
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare',
#'                   dv='follicles', time='Time', phase='Phase')
#' t1$getTime_Power()
#' t1$time_powers
#'
#' # getAR_order works on one case at a time
#' t1$getAR_order()
#' t1$corStructs


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
                         .standardize = FALSE,
                         .corStructs  = NULL,
                         .time_powers = NULL,
                         .alignPhase  = FALSE,
                         .ismonotone  = NULL,
                         .is_clean    = FALSE,
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
                             if(! is.character(value) )
                             {
                               stop("`phase` must be a character variable name in the data")
                             }
                             if( is.null(private$.data[[value]]) )
                             {
                               stop( paste(value, "is not in the data") )
                             }
                             if( ! var(self$data[[value]], na.rm = TRUE) > 0 )
                             {
                               stop( paste(value, "has zero variance") )
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
                             if(! all(is.character(value)) & ! is.null(value) )
                             {
                               stop("`ivs` must be a character list of variables in the data (or `NULL`)")
                             }
                             if( ! all(value %in% names(private$.data) ) )
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
                                           ivs          = value,
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

                         is_clean = function(value)
                         {
                           if( missing(value) ) private$.is_clean
                           else
                           {
                             if(! "logical" %in% class(value)) stop("`is_clean` should be logical")
                             private$.is_clean <- value
                             self
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
                           standardize = FALSE,
                           corStructs  = NULL,
                           time_powers = NULL,
                           ismonotone  = NULL,
                           alignPhase  = FALSE,
                           is_clean    = FALSE,
                           warnings    = list(), # can we hide these? or just make them read only?
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

                           data <- clean(data, ids, dv, time, phase, ivs,
                                         fixed, random, formula, correlation, family,
                                         dvs=NULL, target_ivs=NULL, standardize,
                                         sortData=TRUE, alignPhase)

                           frms <- forms(data,
                                         PalyticObj=NULL,
                                         ids=ids,
                                         dv=dv,
                                         time=time,
                                         phase=phase,
                                         ivs=ivs,
                                         interactions=interactions,
                                         time_power=time_power,
                                         correlation=correlation,
                                         family=family,
                                         fixed=NULL,
                                         random=NULL,
                                         formula=NULL,
                                         method=method)

                           ismonotone <- monotone(ids, time, data)


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
                           private$.is_clean    <- is_clean
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

            })

selfsumm <- function(x)
{
  list(      data         = x$data         ,
             ids          = x$ids          ,
             dv           = x$dv           ,
             time         = x$time         ,
             phase        = x$phase        ,
             ivs          = x$ivs          ,
             interactions = x$interactions ,
             time_power   = x$time_power   ,
             correlation  = x$correlation  ,
             family       = x$family       ,
             fixed        = x$fixed        ,
             random       = x$random       ,
             formula      = x$formula      ,
             method       = x$method       ,
             standardize  = x$standardize  ,
             corStructs   = x$corStructs   ,
             time_powers  = x$time_powers  ,
             ismonotone   = x$ismonotone   ,
             is_clean     = x$is_clean     ,
             warnings     = x$warnings     ,
             errors       = x$errors       ,
             try_silent   = x$try_silent   )
}

Palytic$set("public", "arma",
            function(subgroup=NULL, max.p=3, max.q=3,
                     max.P=0, max.Q=0, max.d=0, max.D=0, ...)
            {
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$data))
              tempData <- subset(self$data, subgroup,
                                 all.vars(self$formula))

              # check that only one participant is in the data
              if( length(unique(tempData[[self$ids]])) != 1 )
              {
                stop('`arma` requires n=1 data')
              }

              # `xreg` requires pre-constructed interaction terms, here we
              # 1. create the model.matrix RHS from self$fixed
              # 2. drop the intercept column
              xdat <- model.matrix(self$fixed, tempData)[,-1]

              # auto detect residual correlation structure here, time power
              # must be detected elsewhere or added here
              m1 <- forecast::auto.arima(y     = tempData[[self$dv]],
                                         xreg  = xdat,
                                         max.p = max.p,
                                         max.q = max.q,
                                         max.P = max.P,
                                         max.Q = max.Q,
                                         max.d = max.d,
                                         max.D = max.D,
                                         ...)



              m1 <- list(arima = m1, tTable = lmtest::coeftest(m1),
                         PalyticSummary = selfsumm(self))

            }
              )

Palytic$set("public", "lme",
            function(subgroup=NULL, ...)
            {
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$data))
              tempData <- subset(self$data, subgroup,
                                 all.vars(self$formula))
              # github issue #1
              cor <- eval(parse(text = ifelse(!is.null(self$correlation),
                                              self$correlation,
                                              'NULL')))
              wm <- 1
              m1 <- try(nlme::lme(fixed=self$fixed,
                                  data=na.omit(tempData),
                                  random=self$random,
                                  correlation=cor,
                                  method=self$method,
                                  keep.data=FALSE,
                                  ...),
                        silent = TRUE)

              ctrl <- nlme::lmeControl()
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 2
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=na.omit(tempData),
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method,
                                    control=ctrl,
                                    ...),
                          silent = TRUE)
              }
              # try without correlation structure
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 3
                ctrl <- nlme::lmeControl(opt="optim")
                self$correlation <- "NULL" # not updating
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=na.omit(tempData),
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=ctrl,
                                    ...),
                          silent = TRUE)
              }
              # try without random slopes or correlation structure
              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 4
                newformula <- forms(data       = self$data ,
                                    PalyticObj = self      ,
                                    dropTime   = TRUE      )
                self$random <- newformula$random
                self$correlation <- NULL
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=na.omit(tempData),
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=ctrl,
                                    ...),
                          silent = TRUE)
              }
              if( "lme" %in% class(m1) )
              {
                m1$PalyticSummary  <- selfsumm(self)
                m1$whichPalyticMod <- paste('Palytic lme model #', wm)
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



### this needs to be expanded to include our list, e.g.,
#   R:\PaCCT\Process\MMTA Process and Record Keeping.docx
# -- this will be done in Palytic augmentations in .Palytic
Palytic$set("public", "gamlss",
            function(subgroup=NULL, sigma.formula = ~1, family=NULL, ...)
            {
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$data))
              tempData <- na.omit( subset(self$data, subgroup,
                                          all.vars(self$formula)) )

              # allow for family to be changed on the fly
              currentFamily <- self$family
              if(!is.null(family)) currentFamily <- family

              wm <- 1 # default model
              ctrl <- gamlss::gamlss.control()
              m1 <- try(gamlss::gamlss(formula = self$formula,
                                       sigma.formula = sigma.formula,
                                       data    = tempData,
                                       family  = currentFamily,
                                       ...),
                        silent = self$try_silent)

              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- 2 # default model with increased n.cyc
                ctrl <- gamlss::gamlss.control(n.cyc=100)
                m1 <- try(refit(gamlss::gamlss(formula = self$formula,
                                               sigma.formula = sigma.formula,
                                               data    = tempData,
                                               family  = currentFamily,
                                               control = ctrl,
                                               ...)),
                          silent = self$try_silent)
              }
              if( "try-error" %in% class(m1) )
              {
                wm <- 3 # drop the random slope(s)
                newformula <- forms(data       = self$data   ,
                                    PalyticObj = self        ,
                                    dropTime   = TRUE        ,
                                    family     = currentFamily,
                                    ... )
                self$formula <- newformula$formula
                #self$family  <- newformula$family
                ctrl <- gamlss::gamlss.control(n.cyc=100)
                m1 <- try(refit(gamlss::gamlss(formula = self$formula,
                                               data    = tempData,
                                               family  = currentFamily,
                                               control = ctrl,
                                               ...)),
                          silent = self$try_silent)
              }
              if("gamlss" %in% class(m1))
              {
                m1$PalyticSummary  <- selfsumm(self)
                m1$whichPalyticMod <- paste('Palytic gamlss model #', wm)
                return(m1)
              }
              if("try-error" %in% class(m1))
              {
                return("Model did not converge")
              }

            },
            overwrite = TRUE
)

# this will only be applied to one participant at a time
# whichIC can take on AIC or BIC
Palytic$set("public", "getAR_order",
            function(P=3    ,
                     Q=3    ,
                     whichIC="BIC" ,
                     lrt=FALSE  , # only valid for nested models, AR, MA are not nested
                     alpha=.05  )
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "residual correlation structure starting.")
              start <- Sys.time()

              saveMethod  <- self$method
              saveFormula <- self$formula

              self$method <- "ML"

              # check whether time is (approximately) equally spaced, use the
              # first participant for now; the method is the standard deviation
              # of the 1st order difference score
              uid <- sort( as.numeric( unique(self$data[[self$ids]]) ) )
              tt  <- self$data[[self$time]][self$data[[self$ids]]==uid[1]]
              ttd <- diff( tt )
              eqSpace <- sd(ttd) < .1

              # for bluedoor review 20180827 set eqSpace <- 0, 'ic' is not
              # passing correctly to auto.arima()
              eqSpace <- FALSE

              # if time is equally spaced use the auto.arima function
              # NOTE: this has not been tested agaignst the lme approach, but
              # it is much faster
              if(eqSpace)
              {
                AR_orders <- by(self$data[[self$dv]],
                                self$data[[self$ids]],
                                FUN = forecast::auto.arima,
                                ic=tolower(whichIC[1]))

                AR_orders <- lapply(AR_orders, function(x)as.data.frame(t(x)))
                AR_orders <- plyr::rbind.fill(AR_orders)
                AR_orders <- data.frame(ids=as.numeric(row.names(AR_orders)),
                                        p=AR_orders[,1],
                                        q=AR_orders[,2])
                AR_orders <- AR_orders[order(AR_orders$ids),]
                row.names(AR_orders) <- NULL

                AR_orders <- data.frame(ids=AR_orders$ids,
                                        arma=paste("corARMA(p=", AR_orders$p, ",q=",
                                                   AR_orders$q, ")", sep=""))
                AR_orders$arma <- as.character( AR_orders$arma )
                AR_orders[AR_orders=="corARMA(p=0,q=0)"] <- "NULL"

                AR_orders[,2]   <- as.character( AR_orders[,2] )
                self$corStructs <- AR_orders
              }

              # this approach is residual, auto.arima is observed
              # an option to increase speed would be to use auto.arima
              # on the residuals
              if(!eqSpace)
              {
                # this is supposedly taboo but I've been unable to work around it b/c we
                # cannot :: the %dopar% operator

                # parralelization
                require(foreach)
                funcs <- c("mean")
                cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
                snow::clusterExport(cl, funcs)
                doSNOW::registerDoSNOW(cl)
                pkgs  <- c("gamlss", "nlme")
                pb <- txtProgressBar(max = length(uid), style = 3)
                progress <- function(n) setTxtProgressBar(pb, n)
                opts <- list(progress = progress)

                #bestCors <- list()
                #for(id in uid)
                bestCors <- foreach(id=uid, .packages = pkgs,
                                    .options.snow = opts)  %dopar%
                {
                  corModsid <- list(); cc = 1
                  self$correlation <- "NULL"
                  nullMod <- self$lme(self$data[[self$ids]]==id)
                  print(id)
                  if( "lme" %in% class(nullMod) )
                  {
                    for(p in 0:P)
                    {
                      for(q in 0:Q)
                      {
                        if(p>0 | q>0)
                        {
                          # will this automatically update the fixed effects? it doesn't
                          # need to for nlme which takes `correlation` directly, but would
                          # need to be updated for gamlss; hence, lme for now (faster too)
                          cortemp <- paste("nlme::corARMA(p=", p, ",
                                         q=", q, ")", sep="")
                          cortemp <- gsub('\n', '', cortemp)
                          cortemp <- gsub(' ', '', cortemp)
                          self$correlation <- cortemp
                          corModsid[[cc]]  <- self$lme(self$data[[self$ids]]==id)
                          #print( corModsid[[cc]]$dims$N )
                          if( any(corModsid[[cc]]=="Model did not converge") )
                          {
                            corModsid[[cc]] <- NULL
                          }
                          else cc = cc + 1
                          self$formula <- saveFormula # restore formulae, incl. correlation
                        }
                      }
                    }
                    names(corModsid) <- unlist( lapply(corModsid,
                                                function(x) x$PalyticSummary$correlation) )
                    corModsid <- corModsid[!is.na(names(corModsid))]

                    if(!lrt)
                    {
                      if(whichIC=="AIC") ICs <- data.frame( Cor = names(corModsid),
                                                       unlist( lapply(corModsid, AIC) ) )
                      if(whichIC=="BIC") ICs <- data.frame( Cor = names(corModsid),
                                                       unlist( lapply(corModsid, BIC) ) )
                      else( stop( paste(whichIC, "is not a valid value for `whichIC`,",
                                        "use AIC or BIC.")))
                      names(ICs) <- c('Cor', 'IC')
                      nullModIC  <- c(Cor = "NULL",
                                       IC = ifelse(whichIC=="AIC", AIC(nullMod),
                                                  BIC(nullMod)))
                      nullModIC    <- as.data.frame(t(nullModIC))
                      nullModIC$IC <- as.numeric(as.character(nullModIC$IC))
                      ICs <- rbind(nullModIC, ICs)
                      return( as.character( ICs$Cor[which.min(ICs[,2])] ) )
                      #bestCors[[id]] <- as.character( ICs$Cor[which.min(ICs[,2])] )
                    } # end of if(!lrt)
                  } # end of if( "lme" %in% class(nullMod) )
                  if(! "nlme" %in% class(nullMod) )
                  {
                    #bestCors[[id]] <- "NULL"
                    return( "NULL" )
                  }
                } # end of dopar
                parallel::stopCluster(cl)

                self$corStructs <- data.frame(ids=uid,
                                              arma=as.character( unlist(bestCors)) )
                self$method  <- saveMethod
                self$formula <- saveFormula # restore formulae, incl. correlation
              } #oef !eqspace
              #self$corStructs$arma[self$corStructs$arma=="NULL"] <- NULL
              message("\nAutomatic detection of the residual\n",
                      "correlation structure took ",
                      round((Sys.time() - start)/60,1), " minutes.\n\n")
            },
            overwrite = TRUE)

# whichIC can take AIC or BIC
Palytic$set("public", "GroupAR_order",
            function(P=3, Q=3, whichIC="BIC", lrt=FALSE, alpha=.05,
                       subgroup=NULL)
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "residual correlation structure starting...")
              start <- Sys.time()

              saveMethod      <- self$method
              saveFormula     <- self$formula

              self$method <- "ML"
              self$correlation <- "NULL"

              nullMod <- self$lme(subgroup)

              if( "lme" %in% class(nullMod) )
              {
                require(foreach)
                funcs <- c("mean")
                cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
                snow::clusterExport(cl, funcs)
                doSNOW::registerDoSNOW(cl)
                pkgs  <- c("gamlss", "nlme")
                pb <- txtProgressBar(max = P, style = 3)
                progress <- function(n) setTxtProgressBar(pb, n)
                opts <- list(progress = progress)

                corMods <- foreach(p=0:P, .packages = pkgs,
                                   .options.snow = opts)  %dopar%
                {
                  corMod <- list()
                  for(q in 0:Q)
                  {
                    if(p>0 | q>0)
                    {
                      cortemp <- paste("nlme::corARMA(p=", p, ",
                                     q=", q, ")", sep="")
                      cortemp <- gsub('\n', '', cortemp)
                      cortemp <- gsub(' ', '', cortemp)
                      self$correlation <- cortemp
                      corMod[[q+1]] <- self$lme(subgroup)
                      if( any(corMod[[q+1]]=="Model did not converge") )
                      {
                        corMod[[q+1]] <- NULL
                      }
                      self$formula <- saveFormula # restore formulae, incl. correlation
                    }
                  }
                  return(corMod)
                }
                parallel::stopCluster(cl)

                corMods <- do.call(list, unlist(corMods, recursive=FALSE))
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

              self$correlation <- bestCor
              self$method      <- saveMethod

              message("\nAutomatic detection of the residual\n",
                      "correlation structure took ",
                      round((Sys.time() - start)/60,1), " minutes.\n\n")
            },
            overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
Palytic$set("public", "getTime_Power",
            function(maxOrder=3, whichIC="BIC")
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              self$method <- "ML"
              uid <- sort( as.numeric( unique(self$data[[self$ids]]) ) )
              time_powers <- list()
              temp <- Palytic$new(self$data, self$ids, self$dv,
                                  self$time)
              for(id in uid)
              {
                aics <- list()
                for(i in 1:maxOrder)
                {
                  temp$time_power <- i
                  mod0 <- temp$lme(self$data[[self$ids]]==id)
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
                      "relationship took ",
                      round((Sys.time() - start)/60,1), " minutes.\n\n")
            },
            overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
Palytic$set("public", "GroupTime_Power",
            function(maxOrder=3, whichIC="BIC")
            {
              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              self$method <- "ML"
              #temp <- Palytic$new(self$data, self$ids, self$dv,
              #                    self$time)
              mods <- list()
              for(i in 1:maxOrder)
              {
                self$time_power <- i
                mods[[i]] <- self$lme()
                if(! "lme" %in% class(mods[[i]]))
                {
                  mods[[i]] <- NULL
                }
              }
              if(whichIC=="BIC") bestMods <- unlist( lapply(mods, BIC) )
              if(whichIC=="AIC") bestMods <- unlist( lapply(mods, AIC) )
              self$time_power <- which.min( bestMods )

              message("\nAutomatic detection of the time/outcome\n",
                      "relationship took ",
                      round((Sys.time() - start)/60,1), " minutes.\n\n")
            },
            overwrite = TRUE)


