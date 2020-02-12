#' .active - active bindings for Palytic objects
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
.active <- function()
{
  list(
    data = function(value)
    {
      if( missing(value) ){ private$.data }
      else
      {
        stop("\n`$data` is read only", call. = FALSE)
      }
    },

    ids = function(value)
    {
      if( missing(value) ) private$.ids
      else
      {
        if(! is.character(value) )
        {
          stop("\n`ids` must be a character variable name in data")
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, "is not in the data") )
        }
        frms <- forms(self$datac,
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
          stop("\n`dv` must be a character variable name in data")
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, "is not in the data") )
        }
        frms <- forms(self$datac,
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
          stop("\n`time` must be a character variable name in the data")
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, "is not in the data") )
        }
        frms <- forms(self$datac,
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
                                          private$.time$raw,
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
            stop("\n`phase` must be a character variable name in the data")
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

        frms <- forms(self$datac,
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
          stop("\n`ivs` must be a character list of ",
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
          stop("\n`interactions` must be a list of vector pairs of variables in the data")
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
        frms <- forms(self$datac,
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
        if(! is.numeric(value) ) stop("\n`time_power` must be numeric")
        if( round(value, 0) != value | value < 1 )
        {
          stop("\n`time_power` must be a positive whole number")
        }
        frms <- forms(self$datac,
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
        #TODO: move to utils
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
        frms <- forms(self$datac,
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

    correlation0 = function(value)
    {
      if( missing(value) ) private$.correlation0
      else
      {
        stop("\n`correlation0`, the input value of `correlation`, cannot be edited.")
      }
    },

    family = function(value)
    {
      if( missing(value) ) private$.family
      else
      {
        if(! "gamlss.family" %in% class(as.gamlss.family(value)) )
        {
          stop("\n`family`=", value, " is not in ",
               "gamlss.family, see `?gamlss.dist::gamlss.family`")
        }
        frms <- forms(self$datac,
                      PalyticObj   = self,
                      ids          = NULL,
                      dv           = NULL,
                      time         = NULL,
                      phase        = NULL,
                      ivs          = NULL,
                      interactions = NULL,
                      time_power   = NULL,
                      correlation  = NULL,
                      family       = as.gamlss.family(value),
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
          stop("\n`fixed` must be a formula, see `?formula` and `??nlme::lme`")
        }
        frms <- forms(data         = self$datac  ,
                      PalyticObj   = self        ,
                      ids          = NULL        ,
                      dv           = NULL        ,
                      time         = NULL        ,
                      phase        = NULL        ,
                      ivs          = NULL        ,
                      interactions = NULL        ,
                      time_power   = NULL        ,
                      correlation  = NULL        ,
                      family       = NULL        ,
                      fixed        = value       ,
                      random       = NULL        ,
                      formula      = NULL        ,
                      method       = NULL        )
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
          stop("\n`random` must be a formula, see `?formula` and `??nlme::lme`")
        }
        frms <- forms(self$datac,
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
          stop("\n`formula` must be a formula, see `?formula` and `??gamlss::gamlss`")
        }
        suppressWarnings(
          frms <- forms(self$datac,
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
          stop("\n`method` should be `ML` or `REML`")
        }
        suppressWarnings(
          frms <- forms(self$datac,
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
          stop("\n`standardize` should be `TRUE` or `FALSE`")
        }
        private$.standardize <- value
        self
      }
    },

    autoSelect = function(value)
    {
      if( missing(value) ) private$.autoSelect
      else
      {
        opts <- c("AR", "TO", "DIST", "Done")
        if(! all(names(value) %in% opts ) )
        {
          stop("\n`autoselect` must include each of\n",
               paste(opts, collapse="\n"))
        }
        private$.autoSelect <- value
        self
      }
    },

    whichIC = function(value)
    {
      if( missing(value) ) private$.whichIC
      else
      {
        opts <- c("BIC", "AIC")
        if(! value %in% opts )
        {
          stop("\n`whichIC` must be one of ", paste(opts, collaps=", "))
        }
        private$.whichIC <- value
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

    alignPhase = function(value)
    {
      if( missing(value) ) private$.alignPhase
      else
      {
        opts <- c('piecewise', 'align', 'none')
        if(value %in% opts )
        {
          private$.alignPhase <- value
        }
        else stop('`alignPhase` must take on one of ',
                  paste('`', opts, '`', sep=''))
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

    ismonotone = function(value)
    {
      if( missing(value) ) private$.ismonotone
      else
      {
        stop("\n`ismonotone` is read only", call. = FALSE)
      }
    },

    datac = function(value)
    {
      if( missing(value) ) private$.datac
      else
      {
        stop("\n`$datac (cleaned data) is read only", call. = FALSE)
      }
    },

    try_silent = function(value)
    {
      if( missing(value) ) private$.try_silent
      else
      {
        if(! "logical" %in% class(value)) stop("\n`try_silent` should be logical")
        private$.try_silent <- value
        self
      }
    },

    debugforeach = function(value)
    {
      if( missing(value) ) private$.debugforeach
      else
      {
        stop("\n`$debugforeach is read only", call. = FALSE)
      }
    }

  )
}


#' Palytic class generator.
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
#' @import moments
#'
#' @export
#'
#' @keywords data
#'
#' @usage Palytic$new(data, ids, dv, time)
#'
#' @details
#' The fields \code{data}, \code{ids}, \code{dv}, and \code{time} are required.
#' Using these, the default model \eqn{dv=time} with random intercepts for \code{ids}
#' and random intercepts for \code{time} is constructed. See the example. If
#' \code{phase} is provided, the default model is \eqn{dv=time+phase+phase*time}, and if
#' \code{ivs} are provided they are included in the model.
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' # construct a new Payltic object and examine the default formulae#'
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
#'                   time='Time', phase='Phase', autoSelect=list())
#'
#' # summary, descriptive, and plot methods
#' t1$summary()
#' t1$describe()
#' t1$plot()
#'
#' # check the formulae creation
#' t1$fixed
#' t1$random
#' t1$formula
#' t1$correlation
#'
#' # Compare gamlss and lme output, in which the models of the default formulae
#' # are fit. Note that the estimates are the same (within rounding) but that
#' # the gamlss SE are much smaller. This is due to gamlss modeling the variance
#' # which reduces the redisudual variance
#' t1.gamlss <- t1$gamlss()
#' t1.lme    <- t1$lme()
#' t1.gamlss$tTable
#' t1.lme$tTable
#'
#' # now change the correlation structure and compare gamlss and lme output,
#' # noting that the intercepts are very different now
#' t1$correlation <- "corARMA(p=1, q=0)"
#' t1$gamlss()$tTable
#' t1$lme()$tTable
#'
#' # fit the model only to the first mare with ML instead of REML
#' t1$method <- 'ML'
#' t1$gamlss(OvaryICT$Mare==1)$tTable
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
#' # note that prior examples set
#' # `autoSelect=list()`, here we use the default, which is to autoselect
#' # the correlation structure (AR), the polynomial order of time (TO) and
#' # the distribution
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare',
#'                   dv='follicles', time='Time', phase='Phase',
#'                   autoSelect=list(AR=list(P=3, Q=3)     ,
#'                                TO=list(polyMax=3)    ,
#'                                DIST=list())  )
#'
#' # automatically select the polynomial order of time with getTO
#' t1$getTO()
#' t1$time_powers
#'
#' # automatically select the ARMA model for residual correlation getAR
#' t1$GroupAR()
#' t1$correlation
#' t1$corStructs
#'
#' # automatically select the distribution, noting that calling $dist() updates $family
#' t1$dist()
#' t1$family
#'
#' # automatically select all three which follows the order
#' # 1. DIST (which will switch package to gamlss for TO and AR)
#' # 2. TO (which can subsequently impact AR)
#' # 3. AR
#' t1$detect()
#'
#' # construct a new Payltic object with no phase variable
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
#'                   time='Time', phase=NULL)
#' t1$plot()
#'
#' # piecewise example
#' OvaryICT$TimeP <- round(30*OvaryICT$Time)
#' t1 <- Palytic$new(data = OvaryICT, ids = 'Mare',
#'                   dv = 'follicles', time = 'TimeP', phase = 'Phase',
#'                   alignPhase = 'piecewise', autoSelect=list())
#' t1$time
#' t1$lme()
#'
#' # piecewise with finite population correction for a population of N=200
#' t1$lme()$tTable
#' t1$lme(fpc = TRUE, popsize2 = 200)$FPCtTable
#'
#' }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# start of Palytic class ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Palytic <- R6::R6Class("Palytic",
                       private = list(
                         .data        = NULL, # consider pass by reference environment
                         .ids         = NULL,
                         .dv          = NULL,
                         .time        = NULL,
                         .phase        = NULL,
                         .ivs          = NULL,
                         .interactions = NULL,
                         .time_power   = NULL,
                         .correlation  = NULL,
                         .correlation0 = NULL,
                         .family       = NULL,
                         .fixed        = NULL,
                         .random       = NULL,
                         .formula      = NULL,
                         .method       = NULL,
                         .standardize  = list(dv=FALSE, iv=FALSE, byids=FALSE),
                         .autoSelect   = NULL,
                         .whichIC      = NULL,
                         .corStructs   = NULL,
                         .time_powers  = NULL,
                         .alignPhase   = NULL,
                         .ismonotone   = NULL,
                         .datac        = NULL,
                         .debugforeach = NULL,
                         .try_silent   = TRUE
                       ),
                       active = .active(),

                       #--------------------------------------------------------
                       # initialize ####
                       #--------------------------------------------------------
                       public = list(

                         #' @description
                         #' Create a new Palytic Object
                         #' @param data A \code{\link{data.frame}} that contains as variables \code{ids},
                         #' \code{dv}, \code{phase}, and \code{time}. Optionally, additional independent
                         #' variables can be included in \code{ivs}. \code{fixed} and \code{random} formulae
                         #' for \code{\link{lme}} models and \code{formula} for \code{\link{gamlss}} models are
                         #' automatically generated when a \code{Palytic} object is created if these fields
                         #' are left \code{NULL}.
                         #'
                         #' @param ids A character string giving the name of the id variable in \code{data}.
                         #'
                         #' @param dv A character string giving the name of the dependent variable in \code{data}.
                         #'
                         #' @param time A character string giving the name of the time variable in \code{data}.
                         #' Random slopes for time are included by default. This can be overridden by specifying
                         #' \code{fixed} and \code{random} formula for \code{\link{lme}} models or by specifying
                         #' the \code{formula} for \code{\link{gamlss}} models.
                         #'
                         #' @param phase A character string giving the name of the phase variable in \code{data}.
                         #' The \code{phase*time} interaction is included by default. This can be overridden by
                         #' specifying \code{fixed} and \code{random} formula for \code{\link{lme}} models or by
                         #' specifying the \code{formula} for \code{\link{gamlss}} models.
                         #'
                         #' @param ivs A \code{\link{list}} of one or more character strings giving the names
                         #' of additional variables in \code{data}, e.g., \code{list('iv2', 'iv2')}.
                         #'
                         #' @param interactions List of vector pairs of variable names for which interaction
                         #' terms should be specified, e.g., \code{list(c('time', 'phase'), c('time', 'iv1'),
                         #' c('iv1', 'iv2'))} where \code{'iv1'} is the name of a variable in the list \code{ivs}.
                         #'
                         #' @param time_power The polynomial for \code{time}, e.g., \code{time^time_power}. Fixed
                         #' effects for \code{time^1...time^time_power} will be included in models. Future
                         #' releases will allow for other functions of time such as \code{\link{sin}}, but these
                         #' can be applied directly by transforming the \code{time} variable.
                         #'
                         #' @param correlation See \code{\link{corStruct}}. Defaults to \code{NULL}, see
                         #' \code{\link{lme}}. Used by both \code{\link{lme}} and \code{\link{gamlss}} models.
                         #'
                         #' @param correlation0 Other options such as \code{autoSelect} can change \code{correlation},
                         #' \code{correlation0} retains the original value provided by the user.
                         #'
                         #' @param family The \code{\link{gamlss.family}} distribution.
                         #'
                         #' @param fixed The \code{fixed} effects model for \code{\link{lme}} models.
                         #'
                         #' @param random The \code{random} effects model for \code{\link{lme}} models.
                         #'
                         #' @param formula The \code{formula} effects model for \code{\link{gamlss}} models.
                         #' \code{sigma.formula}, \code{nu.formula}, and \code{tau.formula} will be implemented in
                         #' a future release.
                         #'
                         #' @param method See \code{method} in \code{\link{lme}}. Is usef for both \code{\link{lme}}
                         #' and \code{\link{gamlss}} models.
                         #'
                         #' @param standardize Named logical list. Which variables should be standardized? The default
                         #' is \code{list(dv=FALSE, ivs=FALSE, byids=FALSE)}. See \code{dv} and \code{ivs}. The option
                         #' \code{byids} controls whether standardization is done by individuals or by group. Any time
                         #' variables are changed (e.g., \code{ivs}), the data are subset, or the options in
                         #' \code{standardize} are changed, the raw data will be restandardized (see \code{datac}).
                         #'
                         #' @param autoSelect List. The default is
                         #' \code{
                         #' list(AR=list(P=3, Q=3)     ,
                         #'   TO=list(polyMax=3)       ,
                         #'   DIST=list()) }.
                         #'
                         #' If no automated model selection for the residual covariance structure (\code{AR}),
                         #' the polynomial order for the relationship between time and the dependent variable
                         #' (\code{TO}), or the dependent variable distribution is desired, an empty list
                         #' should be passed (e.g., \code{autoSelect=list()}).
                         #'
                         #' If \code{AR} is in the list,
                         #' the residual correlation structure will be automatically selected from
                         #' among \code{ARMA(p,q)} models. See \code{correlation}. Since these models are
                         #' not generally nested, model selection is done using information information
                         #' criterion (see \code{whichIC}). Model selection for the residual covariance
                         #' structure is searches among
                         #' \code{p=1,...,P} and \code{p=1,...,Q}, where \code{P} and \code{Q} are taken
                         #' from \code{PQ}, i.e., \code{PQ=c(P,Q)}. The values of \code{p} and \code{p}
                         #' are passed to \code{\link{corARMA}} ( e.g., \code{corARMA(p=p,q=q)}).
                         #' If \code{individual_mods=FALSE}, this done
                         #' comparing \code{lme} modes for N>1 data. If \code{individual_mods=TRUE},
                         #' this is done using the \code{\link{auto.arima}} function on the residuals for
                         #' each individual. For more detail, see the \code{$GroupAR()}
                         #' method.
                         #'
                         #' If \code{TO} is in the list, models with polynomial powers of time from 1 to
                         #' \code{polyMax} will be tested.
                         #' For example, if \code{polyMax=3} (implying a cubic growth model), the models
                         #' compared include \code{time}, \code{time + I(time^2)}, and
                         #' \code{time + I(time^2)+I(time^3)}. Since these models are nested, the best
                         #' fitting model is selected using likelihood ratio tests with mixed effects
                         #' models fit using maximum likelihood estimators in \code{\link{lme}}.
                         #' This is done separately for each individual in \code{ids} if
                         #' \code{individual_mods=TRUE}. For more detail, see the \code{$getTO()}
                         #' method.
                         #'
                         #' If \code{DIST} is in the list and \code{package='gamlss'}, each dependent
                         #' variable in \code{dvs} will utilize the \code{\link{fitDist}} function of
                         #' the gamlss package, and the best fitting distribution will be used for each
                         #' dependent variable. For more detail, see the \code{$dist()} method in
                         #' \code{\link{Palytic}}.
                         #'
                         #' @param whichIC Character. The default is \code{whichIC="BIC"}.
                         #'
                         #' Either the Akaike Information Criterion (\code{whichIC="AIC"}) or
                         #' the Bayesian Information Criterion (\code{whichIC="BIC"}).
                         #'
                         #' If the \code{time} variable is equally spaced, this is
                         #' done using the function \code{\link{forecast}}. If the \code{time} variable
                         #' is not equally spaced, this is done using comparisons of
                         #' mixed effects models using \code{\link{lme}} fit using maximum likelihood
                         #' estimators.
                         #'
                         #' Residual autocorrelation structure is detected separately for each individual
                         #' in \code{ids} if \code{individual_mods=TRUE}.
                         #'
                         #' @param corStructs Vector. A \code{correlation} structure for each case in \code{ids}. Not
                         #' user accessible. Populated by \code{\link{PersonAlytic}}.
                         #'
                         #' @param time_powers Vector. A \code{time_order} for each case in \code{ids}. Not
                         #' user accessible. Populated by \code{\link{PersonAlytic}}.
                         #'
                         #' @param alignPhase Character. Options include
                         #'    a. 'none', no changes are made to the time or phase variable.
                         #'    b. 'align', align the time variable to be zero at the transition between
                         #'       the first and second phase (see \code{\link{alignPhases}}).
                         #'    c. 'piecewise', add 'pwtime#' variables, which will replace time and
                         #'       time_power to create a piecewise linear growth curve model, and where `#`
                         #'       is the number of phases (i.e., one linear growth curve model per phase).
                         #'
                         #' @param ismonotone Logical. Is the \code{time} variable for each case monotonically
                         #' increasing (i.e., no returns to prior values). This is determining in data cleaning as
                         #' described for \code{datac}.
                         #'
                         #' @param datac data.frame. Cleaned data. Cleaning involves the following steps:
                         #' 1. Check that the variables in \code{ids}, \code{dv}, \code{time}, \code{phase},
                         #' \code{ivs}, and \code{interactions} are in \code{data}.
                         #' 2. The \code{ids} variable is forced to be numeric.
                         #' 3. The validity of the \code{correlation} structure is checked, see \code{\link{corStruct}}.
                         #' 4. check that variables have non-zero variance.
                         #' 5. If standardization is requested, standardize the data (see \code{standardize}).
                         #' 6. Sort the data on \code{ids} and \code{time}.
                         #' 7. If patients have < 2 observations, they are dropped from the data set.
                         #' 8. Phase alignment (if any, see \code{alignPhase}).
                         #'
                         #' @param debugforeach Logical flag for testing error handling in parallelized runs.
                         #'
                         #' @param try_silent Logical flag for testing error handling in \code{Palytic} methods.
                         #'
                         #' @return A new `Palytic` object
                         initialize = function
                         (
                           data                                                ,
                           ids          = NULL                                  ,
                           dv           = NULL                                  ,
                           time         = NULL                                  ,
                           phase        = NULL                                  ,
                           ivs          = NULL                                  ,
                           interactions = NULL                                  ,
                           time_power   = NULL                                  ,
                           correlation  = NULL                                  ,
                           family       = gamlss.dist::NO()                     ,
                           fixed        = NULL                                  ,
                           random       = NULL                                  ,
                           formula      = NULL                                  ,
                           method       = "REML"                                ,
                           standardize  = list(dv=FALSE, iv=FALSE, byids=FALSE) ,
                           autoSelect   = list(AR=list(P=3, Q=3)     ,
                                               TO=list(polyMax=3)    ,
                                               DIST=list())                     ,
                           whichIC      = c("BIC", "AIC")                       ,
                           corStructs   = NULL                                  ,
                           time_powers  = NULL                                  ,
                           ismonotone   = NULL                                  ,
                           alignPhase   = 'none'                                ,
                           datac        = NULL                                  ,
                           debugforeach = FALSE                                 ,
                           try_silent   = TRUE
                         )
                         {
                           ### consider adding option to read a file, could
                           #   autoselect file type

                           #if( is.character(data) )
                           #{
                           #  read
                           #}

                           if(debugforeach) message("\nDebugging is ON.\n\n")

                           if( is.null(phase))
                           {
                             message("\nThere is no phase variable, changing",
                                     "\nalignPhase to 'none'.")
                           }
                           if(!is.null(phase))
                           {
                             if(alignPhase == 'piecewise' &
                                length(table(data[[phase]])) <= 1)
                             {
                               message("\nThere are 0 or 1 phases, changing",
                                       "\nalignPhase to 'none'.")
                               alignPhase <- 'none'
                             }
                           }

                           # checks that get used multiple times
                           is.min <- !(is.null(ids) | is.null(dv) | is.null(time))
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

                           # update the time variable
                           time <- list(raw      = time        ,
                                        power    = time_power  ,
                                        analysis = time        )

                           # clean the data
                           if(is.null(datac))
                           {
                             datac <- clean(
                               data         = data            ,
                               ids          = ids             ,
                               dv           = dv              ,
                               time	        = time            ,
                               phase	      = phase           ,
                               ivs	        = ivs             ,
                               fixed        = fixed           ,
                               random       = random          ,
                               formula	    = formula         ,
                               correlation  = correlation     ,
                               family       = family          ,
                               dvs	        = NULL            ,
                               target_ivs   = NULL            ,
                               standardize  = standardize     ,
                               sortData	    = TRUE            ,
                               alignPhase   = alignPhase      ,
                               debugforeach = debugforeach           )
                           }


                           # if alignPhase == 'piecewise' update time
                           if(alignPhase == 'piecewise')
                           {
                             time$analysis <- names(datac)[grepl('pwtime', names(datac))]
                           }

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
                           ismonotone <- monotone(ids, time$raw, datac)

                           # add 'Done' to autoSelect
                           if(! any(names(autoSelect) %in% "Done"))
                           {
                             autoSelect$Done <- FALSE
                           }

                           # populate private
                           private$.data         <- data
                           private$.ids          <- ids
                           private$.dv           <- dv
                           private$.time         <- time
                           private$.phase        <- phase
                           private$.ivs          <- ivs
                           private$.interactions <- interactions
                           private$.time_power   <- time_power
                           private$.correlation  <- correlation
                           private$.correlation0 <- correlation
                           private$.family       <- family
                           private$.fixed        <- frms$fixed
                           private$.random       <- frms$random
                           private$.formula      <- frms$formula
                           private$.method       <- method
                           private$.standardize  <- standardize
                           private$.autoSelect   <- autoSelect
                           private$.whichIC      <- whichIC
                           private$.corStructs   <- corStructs
                           private$.time_powers  <- time_powers
                           private$.ismonotone   <- ismonotone
                           private$.alignPhase   <- alignPhase
                           private$.datac        <- datac
                           private$.debugforeach <- debugforeach
                           private$.try_silent   <- try_silent

                           # cleanup
                           rm(frms)

                         }

                       )

)

# add methods ####

#' @description
#' Print a summary of the inputs, the cleaned data,
#' and the raw data in a Palytic object
#' @param  wm Which model.
#' @return A list inlcuding the user's inputs and descriptive statistics for
#' the raw and cleaned data.
Palytic$set("public", "summary",
            function(wm=NULL)
            {
              variables <- all.vars( self$formula )
              dropvars  <- grep(pattern = '\\(', variables)
              if(length(dropvars) > 0 ) variables <- variables[-dropvars]

              tempCorrelation <- ifelse(is.null(self$correlation), "NULL",
                                        self$correlation)

              varsInData <- names(self$data)[names(self$data) %in% variables]

              augmentFormula <- self$formula
              if(!is.null(wm)) augmentFormula <- list(formula = self$formula,
                                                      whichPalyticMod = wm)

              # print descriptive statistics to the console
              if(!is.null(self$phase)) phase <- self$datac[[self$phase]]
              else phase <- self$phase

              # return the summary
              list(      ids          = self$ids                         ,
                         dv           = self$dv                          ,
                         time         = self$time                        ,
                         phase        = self$phase                       ,
                         ivs          = self$ivs                         ,
                         interactions = self$interactions                ,
                         time_power   = self$time_power                  ,
                         correlation  = tempCorrelation                  ,
                         family       = self$family                      ,
                         package      = self$package                     ,
                         fixed        = self$fixed                       ,
                         random       = self$random                      ,
                         formula      = augmentFormula                   ,
                         method       = self$method                      ,
                         standardize  = self$standardize                 ,
                         corStructs   = self$corStructs                  ,
                         time_powers  = self$time_powers                 ,
                         ismonotone   = self$ismonotone                  ,
                         debugforeach = self$debugforeach                ,
                         try_silent   = self$try_silent                  ,
                         datac        = summary(self$datac[,variables])  ,
                         data         = summary(self$data[,varsInData])  ,
                         skew_kurt    = dstats(self$datac[,self$dv],
                                               phase, print=FALSE)       )
            },
            overwrite = TRUE
)

#' @description
#' Conduction autoselection based on the inputs to the \code{autoSelect} parameter
#' when initializing an new Palytic object.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param model See \code{\link{fitDist}}.
#' @param parallel See \code{\link{fitDist}}.
#' @param plot Logical. Should plots of your dependent variable be displayed.
#' @param userFormula List. Formulae to override the default formulae created
#' when initializing a Palytic object. Options must be named and include either
#' \code{formula}, or both \code{fixed} and \code{random}. These are described
#' in the section documenting the \code{new()} method.
#' @param dims List. Used internally for high throughput. The user should not
#' adjust this parameter.
#' @param package Character. Options are \code{"nlme"} and \code{"gamlss"} as
#' described in the section documenting the \code{new()} method.
Palytic$set("public", "detect",
            function(subgroup    =  NULL                  ,
                     model       =  NULL                  ,
                     parallel    =  "snow"                ,
                     plot        =  TRUE                  ,
                     userFormula =  list(
                                      fixed=NULL,
                                      random=NULL,
                                      formula=NULL),
                     dims        =  list(ID="All Cases")  ,
                     package     =  "nlme"
                     )
            {

              # prevent recursive calls to autoselect
              temp <- self$autoSelect
              temp$Done <- TRUE
              self$autoSelect <- temp; rm(temp)

              # extract logicals
              detectAR=FALSE  ; if("AR" %in% names(self$autoSelect))   detectAR=TRUE
              detectTO=FALSE  ; if("TO" %in% names(self$autoSelect))   detectTO=TRUE
              detectDist=FALSE; if("DIST" %in% names(self$autoSelect)) detectDist=TRUE

              # allow for formula override so that we can test intercept only and
              # slope only models
              if( any(unlist(lapply(userFormula, function(x) !is.null(x)))) )
              {
                if( isNullOrForm(userFormula$fixed) ) self$fixed <- userFormula$fixed
                if( isNullOrForm(userFormula$random) ) self$random <- userFormula$random
                if( isNullOrForm(userFormula$formula) ) self$formula <- userFormula$formula
              }

              # individual vs. group model has no bearing on detectDist, but
              # detectDist should happen first for when detectAR and detectTO
              # are generalized for gamlss()
              if(detectDist)
              {
                self$dist(self$autoSelect$DIST$to01,
                          model=model, parallel=parallel,
                          plot=plot)
              }
              if(!detectDist)
              {
                # no change, leave default or user supplied self$family
              }

              #...........................................................................
              # get the TO and AR
              #...........................................................................
              if(dims$ID[1]!="All Cases")
              {
                if( detectTO) self$getTO()
                if(!detectTO) self$time_powers <- data.frame(ids=dims$ID,
                                                           rep(1, length(dims$ID)))

                if( detectAR)
                {
                  # this occurs inside the $lme() and $gamlss() functions
                }
                if(!detectAR)
                {
                  self$corStructs <- data.frame(ids=dims$ID,
                                              arma=rep( ifelse(is.null(self$correlation),
                                                               "NULL", self$correlation),
                                                        length(dims$ID)) )
                }

              }

              if(dims$ID[1]=="All Cases")
              {
                if(detectTO) self$GroupTO(subgroup, package)
                if(detectAR) self$GroupAR(subgroup,
                                                doForeach=parallel=="snow",
                                                package)
              }
              # no return (for now)
            },
            overwrite = TRUE
)

#' @description
#' This method plots the density of your dependent variable if \code{plot=TRUE} and
#' lets the user implement \code{gamlss} automated
#' distribution comparisons. If \code{model=NULL}, the \code{\link{fitDist}}
#' function is used to compare unconditional models for all applicable distributions.
#' If \code{model} is a \code{gamlss} model, the condition models are fit for
#' all applicable distributions using \code{\link{chooseDist}}.
#' If you want to rescale your dependent variable to the (0,1)
#' range, set \code{to01=TRUE}. If your dependent variable is multinomial,
#' automated distribution comparisons are not implemented.
#'
#' @param to01 Logical. Should the dependent variable be rescaled to the (0,1) range?
#' See \code{\link{to01}}. The option \code{squeeze} is not yet implemented.
#' @param model See \code{\link{fitDist}}.
#' @param parallel See \code{\link{fitDist}}.
#' @param plot Logical. Should plots of your dependent variable be displayed.
Palytic$set("public", "dist",
            function(to01=FALSE, model=NULL,
                     parallel="snow", plot=TRUE, type=NULL, extra=NULL)
            {
              # extract the dv for convenience
              dv <- na.omit( self$datac[[self$dv]] )

              #
              if(is.null(to01)) to01 <- FALSE

              # plot
              if(plot)
              {
                y  <- data.frame(y=dv)
                yq <- ggplot(y, aes(sample=y)) + stat_qq() + stat_qq_line()
                yd <- ggplot(y, aes(x=y)) + geom_density() + xlab(self$dv)
                yh <- ggplot(y, aes(x=y)) + geom_histogram() + xlab(self$dv)
                capture.output(
                    suppressMessages(print(grid.arrange(yq, yd, yh, nrow=3))),
                  file='NUL')
                #grid.arrange(yq, yd, yh, nrow=3)
              }

              # beta
              if(to01) dv <- to01(dv)

              # get the bounds on the dv
              ncats   <- length(table(dv))
              isInt   <- identical(dv, round(dv,0))
              isBin   <- ncats==2
              is01    <- min(dv, na.rm=TRUE) >= 0 & max(dv, na.rm=TRUE) <= 1
              min0    <- min(dv, na.rm=TRUE) >= 0
              isMult  <- isInt & !isBin & ncats <= 5
              isCount <- isInt & min0
              isCont  <- !isCount & !isInt & !isBin & !isMult

              # set the type parameter
              type <- as.character(NA)
              if(!isCount & !min0) type <- "realline"
              if(!isCount &  min0) type <- "realplus"
              if( isCount        ) type <- "counts"; extra <- c("realplus", extra)
              if(isBin)            type <- "binom"
              if(isMult)           type <- paste("MULTIN(type = '", ncat,"')", sep='')
              if(is01)             type <- "real0to1"

              if(is.na(type))
              {
                stop("\nDistribution comparison cannot be implemented for ", self$dv)
              }
              if(!is.na(type))
              {

                cat("\nThe variable", self$dv, "has the following characteristics:",
                    "\nInteger     : ", isInt  ,
                    "\nBinary      : ", isBin  ,
                    "\nProportion  : ", is01   ,
                    "\nPositive    : ", min0   ,
                    "\nMultinomial : ", isMult , " (integer with >2 & <= 5 categories)",
                    "\nCount       : ", isCount, " (positive integer)"                 ,
                    "\nContinuous  : ", isCont , " (non-integer)"                      ,
                    "\n\n")
              }

              # set the family dependent on type
              if(!is.null(model) & ! "gamlss" %in% class(model))
              {
                stop("\nThe model you provided is not a `gamlss` object.\n",
                     self$dv, " will be tested unconditionally instead of using a model.")
                model <- NULL
              }

              if(is.null(model) & !isMult)
              {
                # waring suppresion works, but error suppresion not working, I've tried
                # suppressWarnings, suppressMessages, sink, capture.output,
                # R.utils::captureOutput, try, tryCatch, invisible,

                options(warn = -1)#, error = function(){cat('')})
                capture.output(
                family <- fitDist(dv, type = type, try.gamlss = TRUE,
                                  extra = distTypes(extra)),
                file = 'NUL')
                options(warn =  0, error = NULL)

                print(family)

                message("\nTo explore this distribution install `gamlss.demo` and type\n",
                        "\ndev.new()",
                        "\ngamlss.demo::demoDist()\n",
                        "\ninto the console, find your distribution, and use the",
                        "\nslider bars to select the parameters printed above",
                        "\n(mu, sigma, nu, and tau).\n\n")

                family <- family$family[1]
              }

              if(!is.null(model) & "gamlss" %in% class(model))
              {
                family <- chooseDist(model, type = "extra",
                                     extra = c(distTypes(type), distTypes(extra)),
                                     parallel = parallel,
                                     ncpus = parallel::detectCores() - 1)
                family <- names(getOrder(family))[which.min(family)]
              }

              self$family <- as.gamlss.family(family)

              # print descriptive statistics to the console
              if(!is.null(self$phase)) phase <- self$datac[[self$phase]]
              else phase <- self$phase
              dstats(dv, phase, print=TRUE)

            },
            overwrite = TRUE
)

#' @description
#' Descriptive statistics for data in a Palytic object. This method gives the
#' correlation between \code{dv} and each continuous variable in \code{ivs}
#' (as well as the \code{time} variable), or, if variablse are factors (including
#' \code{phase}), the mean of \code{dv} is given for each factor level of each
#' variable.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
Palytic$set("public", "describe",
            function(subgroup=NULL)
            {
              subCheck(subgroup, self$datac)

              # concatenate all the ivs with double checks for dropped terms or non variable terms
              ivall <- all.vars(self$fixed)[-1]
              ivall <- ivall[! ivall %in% c("0", "1")] # 0,1 make come from random effects

              if(! self$time$raw %in% ivall ) ivall <- c(ivall, self$time$raw)
              if(length(self$phase) > 0)
              {
                if(! self$phase %in% ivall ) ivall <- c(ivall, self$phase)
              }

              ivall <- ivall[! ivall %in% c("0", "1")] # 0,1 may come from random effects

              # subset the data
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- self$datac[subgroup, c(self$dv, ivall)]

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
            },
            overwrite = TRUE
)

#' @description
#' For individual level models, random intercepts and
#' random slopes are not defined. In this situatation, an \code{ARMA(p,q)}
#' should be used. This is implemented using the \code{xreg}
#' option of \code{\link{arima}}. The residual correlation search is achieved using
#' \code{\link{auto.arima}}. The \code{dropVars} parameter indicates which fixed
#' effects (in \code{xreg}) should be dropped for a likelihood ration test
#' (LRT). This is used by \code{\link{PersonAlytic}} to test \code{target_ivs}.
#' The other parameters are as in \code{\link{auto.arima}}.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param max.p See \code{\link{arima}} and the \code{P} option in the
#' \code{autoSelect} parameter when initializing a Palytic object.
#' @param max.q See \code{\link{arima}} and the \code{Q} option in the
#' \code{autoSelect} parameter when initializing a Palytic object.
#' @param dropVars Character or character vector. Names of variables to drop
#' from the the analysis.
#' @param max.P See \code{\link{arima}}
#' @param max.Q See \code{\link{arima}}
#' @param max.d See \code{\link{arima}}
#' @param max.D See \code{\link{arima}}
Palytic$set("public", "arma",
            function(subgroup=NULL, max.p=3, max.q=3, dropVars=NULL,
                     max.P=0, max.Q=0, max.d=0, max.D=0, ...)
            {
              subCheck(subgroup, self$datac)

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
              if(  any("ARIMA" %in% class(m1)) ) tTable <- lmtest::coeftest(m1)
              if(! any("ARIMA" %in% class(m1)) )
              {
                m1 <- "Model did not converge"
                tTable <- NA
              }

              m1 <- list(arima = m1, tTable = tTable,
                         PalyticSummary = self$summary(),
                         xregs = colnames(xdat),
                         lrt = list(wasLRTrun=wasLRTrun, lrtp=lrtp) )
              return(m1)
            },
            overwrite = TRUE
)

# $lme() ####
#TODO(Stephen) you can't just add ... to lme() and get things passed, you'll
#need to add some sort of evaluation, e.g.
# > args <- list(...)
# > if("contrasts" %in% names(args)) # then pass `contrasts` to lme
# or some other method of matching arguments, see ?match.arg

#' @description
#' This method fits the linear mixed effects
#' \code{lme} model implied by the \code{Palytic} fields \code{ids},
#' \code{dv}, \code{phase}, \code{time}, and optionally \code{ivs},
#' \code{time_power}, and \code{correlation}. The default formula can be
#' overridden by specifying \code{fixed} and \code{random} and optionally
#' \code{method} and \code{correlation}. The see \code{\link{lme}} for the parameter
#' \code{subgroup}. The \code{dropVars} parameter indicates which fixed
#' effects should be dropped for a likelihood ration test (LRT). This is
#' used by \code{\link{PersonAlytic}} to test \code{target_ivs}.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param dropVars Character or character vector. Names of variables to drop
#' from the the analysis and conduct a likelihood ratio test (LRT) for the full
#' model and the reduced model without the variables in \code{dropVars}.
#' @param fpc Logical. Should the finite population correrction (FCP) be made
#' to the analysis.
#' @param popsize2 Numeric. If \code{fpc=TRUE}, what is the population size to
#' be used in the FPC.
Palytic$set("public", "lme",
            function(subgroup=NULL, dropVars=NULL,
                     fpc=FALSE, popsize2, ...)
            {
              subCheck(subgroup, self$datac)
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- na.omit(subset(self$datac, subgroup,
                                 all.vars(self$formula)))

              # detect
              if(any(names(self$autoSelect) %in% c("AR", "TO", "DIST")) &
                 !self$autoSelect$Done)
              {
                self$detect()
              }

              # check fpc inputs
              n <- length(table(tempData[[self$ids]]))
              if(fpc) fpcCheck(popsize2, n)

              # eval parse ok b/c correlation is validated
              cor <- eval(parse(text = ifelse(!is.null(self$correlation),
                                              self$correlation,
                                              'NULL')))
              wm <- "Full Model"
              m1 <- try(nlme::lme(fixed=self$fixed,
                                  data=tempData,
                                  random=self$random,
                                  correlation=cor,
                                  method=self$method),
                        silent = TRUE)


              ctrl <- nlme::lmeControl() # used later if not overwritten

              if( "try-error" %in% class(m1) )# | !eds(m1) )
              {
                wm <- "Full Model, opt='optim'"
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
                wm <- "No correlation structure, opt='optim'"
                ctrl <- nlme::lmeControl(opt="optim")
                self$correlation <- "NULL"
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
                wm <- "No correlation structure, no slopes, opt='optim'"
                newformula <- forms(data       = self$datac ,
                                    PalyticObj = self       ,
                                    dropTime   = "time"     )
                self$correlation <- NULL
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=newformula$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }
              # if the model is piecewise, try taking out the time X phase
              if( "try-error" %in% class(m1) )
              {
                wm <- "No correlation structure, no slopes, no time X phase, opt='optim'"
                newformula <- forms(data = self$datac,
                                    PalyticObj = self,
                                    dropTime = "int" )
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=newformula$fixed,
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
              if( length(table(tempData[[self$ids]]))==1 &
                  !is.null(self$autoSelect$AR) )
              {
                PQ <- c(self$autoSelect$AR$P, self$autoSelect$AR$Q)
                self$correlation  <- getARnEQ1(m1, PQ, self$dv)
                cor  <- eval(parse(text = ifelse(!is.null(self$correlation),
                                                 self$correlation,
                                                 'NULL')))
                ctrl <- nlme::lmeControl(opt="optim")
                m1   <- try(nlme::lme(fixed=self$fixed,
                                    data=tempData,
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
              }

              # clean up the call for accurate printing and model.matrix() in FPC
              if( exists("newformula") )
              {
                m1 <- cleanCall(modelResult=m1, PalyticObj=self,
                                newformula)
              }
              if(!exists("newformula") )
              {
                m1 <- cleanCall(modelResult=m1, PalyticObj=self)
              }
              #!#m1 <<- m1
              #!#print( formula(m1) )
              #!#m1 <- cleanCall(modelResult=m1, PalyticObj=self,
              #!#                newformula)

              # lrt
              # H0: m1 == m0
              # if p<.05, reject H0, retain better fitting model
              # only if lik(m1) always > lik(m0) do we retain m1, this is the assumption Ty
              # currently uses
              # if that does not hold, we need to use some other criterion
              # TODO: add a check of whether lik(m1) > lik(m0)
              wasLRTrun <- FALSE
              lrtp <- as.numeric(NA)
              if( "lme" %in% class(m1) & !is.null(dropVars) )
              {
                frm0 <- remove_terms(self$fixed, dropVars)
                m0 <- try(nlme::lme(fixed=frm0,
                                    data=tempData,
                                    random=self$random,
                                    correlation=cor,
                                    method=self$method,
                                    control=ctrl),
                          silent = TRUE)
                if("lme" %in% class(m1) & "lme" %in% class(m0))
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
                if( fpc ) m1$FPCtTable <- as.matrix(FPC(object=m1, popsize2=popsize2))
                m1$tTable <- summary(m1)$tTable # easier for simulation studies
                m1$PalyticSummary <- self$summary(wm)
                m1$whichPalyticMod <- wm
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

              cat("Palytic$lme:",
                  "\ntempData: ", names(tempData),
                  "\ndim(tempData): ", toString(dim(tempData)),
                  "\ncorrelation: ", toString(self$correlation),
                  "\nfixed: ", toString(self$fixed),
                  "\nrandom: ", toString(self$random),
                  "\nm1: ", class(m1),
                  "\n\n", file="./PAlogs/$lme.log", append=TRUE)
            },
            overwrite = TRUE
)

# cleanCall ####
#' cleanCall - clean up the call in in palytic objects
#' @keywords internal
cleanCall <- function(modelResult, PalyticObj, newformula=NULL)
{
  if("lme" %in% class(modelResult) )
  {
    if(!is.null(newformula))
    {
      modelResult$call$fixed       <- newformula$fixed
      modelResult$call$random      <- newformula$random
    }
    if( is.null(newformula))
    {
      modelResult$call$fixed       <- PalyticObj$fixed
      modelResult$call$random      <- PalyticObj$random
    }
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

# getARnEQ1 ####
#' getARnEQ1
#' @keywords internal
getARnEQ1 <- function(m, PQ, dv, debug=FALSE)
{
  if(is.null(PQ))
  {
    stop("\nIn `getARnEQ1` `PQ` is NULL.\n\n")
  }
  correlation <- "NULL"
  if( "gamlss" %in% class(m) ) resid <- m$residuals
  if( "lme"    %in% class(m) ) resid <- m$residuals[,2]
  faa <- try( forecast::auto.arima(y     = resid,
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
  }
  if(!file.exists('./PAlogs/getARnEQ1run.log') & debug)
  {
    dir.create('./PAlogs')
    cat('', file='./PAlogs/getARnEQ1run.log')
  }
  if(debug) cat( correlation, '\n\n', file = './PAlogs/getARnEQ1run.log', append=TRUE)
  return( correlation )
}

# $gamlss() ####
### TODO(Stephen) add residual correlation search for n=1

#' @description
#' This method fits the \code{lme} model implied
#' by the \code{Palytic} fields \code{ids}, \code{dv}, \code{phase}, \code{time}
#' and optionally \code{ivs}, \code{time_power}, \code{correlation}. The default
#' formula can be overridden by specifying \code{formula} and
#' optionally \code{method}, \code{correlation}, and \code{family} (which can be
#' used to specify generalized linear mixed effects models,
#' see \code{\link{gamlss.family}}).
#' The parameter \code{subgroup} operates as in \code{\link{lme}}. The parameter
#' \code{sigma.formula} and \code{family} are desribed in \code{\link{gamlss}}.
#' The \code{dropVars} parameter indicates which fixed
#' effects should be dropped for a likelihood ration test (LRT). This is
#' used by \code{\link{PersonAlytic}} to test \code{target_ivs}.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param sigma.formula Formula. A formula for the sigma part of a
#' \code{gamlss.family} distribution. Later versions will add nu.formula and
#' tau.formula.
#' @param family. A \code{\link{gamlss.family}} distribution.
#' @param dropVars Character or character vector. Names of variables to drop
#' from the the analysisand conduct a likelihood ratio test (LRT) for the full
#' model and the reduced model without the variables in \code{dropVars}.
Palytic$set("public", "gamlss",
            function(subgroup=NULL, sigma.formula = ~1, family=NULL,
                     dropVars=NULL, ...)
            {
              subCheck(subgroup, self$datac)

              options(warn = -1)

              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- na.omit( subset(self$datac, subgroup,
                                          c(all.vars(self$formula),
                                            all.vars(sigma.formula))) )

              # detect
              if(any(names(self$autoSelect) %in% c("AR", "TO", "DIST"))  &
                 !self$autoSelect$Done)
              {
                self$detect()
              }

              # allow for family to be changed on the fly in $gamlss()
              currentFamily <- self$family
              if(!is.null(family)) currentFamily <- family

              wm <- "Full Model"

              # turn the trace of to prevent printing it to the console,
              # increase c.crit from .001 to .005 to speed things up a little
              ctrl <- gamlss::gamlss.control(trace=FALSE, c.crit=.005)

              m1 <- try(gamlss::gamlss(formula = self$formula,
                                       sigma.formula = sigma.formula,
                                       data = tempData,
                                       family = currentFamily,
                                       control = ctrl),
                        silent = self$try_silent)

              # default model with increased n.cyc
              if( "try-error" %in% class(m1) )
              {
                wm <- "Full Model, n.cyc=100"
                #ctrl <- gamlss::gamlss.control(n.cyc=100, trace=FALSE)
                m1 <- try(refit(gamlss::gamlss(formula = self$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }
              # try without correlation structure
              if( "try-error" %in% class(m1) )
              {
                wm <- "No correlation structure, n.cyc=100"
                self$correlation <- "NULL"
                newformula <- forms(data = self$datac ,
                                    PalyticObj = self ,
                                    dropTime = "time" ,
                                    family = currentFamily)
                self$fixed   <- newformula$fixed
                # this works but fails saving
                m1 <- try(refit(gamlss::gamlss(formula = newformula$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }
              # try without random slopes or correlation structure
              if( "try-error" %in% class(m1) )
              {
                wm <- "No correlation structure, no slopes, n.cyc=100"
                newformula <- forms(data = self$datac ,
                                    PalyticObj = self ,
                                    dropTime = "time" ,
                                    family = currentFamily)
                self$fixed   <- newformula$fixed
                # this works but fails saving
                m1 <- try(refit(gamlss::gamlss(formula = newformula$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }
              # if the model is piecewise, try taking out the time X phase
              if( "try-error" %in% class(m1) )
              {
                wm <- "No correlation structure, no slopes, no time X phase, n.cyc=100"
                newformula <- forms(data = self$datac ,
                                    PalyticObj = self ,
                                    dropTime = "int" ,
                                    family = currentFamily)
                self$fixed   <- newformula$fixed
                # this works but fails saving
                m1 <- try(refit(gamlss::gamlss(formula = newformula$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }

              # if n=1, detect the correlation structure here using getARnEQ1()
              # as it may affect the lrt (we should probably do the same for
              # group ar...);
              if( length(table(tempData[[self$ids]]))==1 &
                  !is.null(self$autoSelect$AR) )
              {
                PQ <- c(self$autoSelect$AR$P, self$autoSelect$AR$Q)
                self$correlation <- getARnEQ1(m1, PQ, self$dv)
                cor <- eval(parse(text = ifelse(!is.null(self$correlation),
                                                self$correlation,
                                                'NULL')))
                newformula <- forms(data = self$datac ,
                                    PalyticObj = self ,
                                    correlation = cor ,
                                    family = currentFamily)
                m1 <- try(refit(gamlss::gamlss(formula = newformula$formula,
                                               sigma.formula = sigma.formula,
                                               data = tempData,
                                               family = currentFamily,
                                               control = ctrl)),
                          silent = self$try_silent)
              }

              # lrt
              # see note in $lme()
              # here the model with the smaller deviance is preferred (as opposed to
              # that with the higher likelihood) and we are assuming that
              # dev(m1) always < dev(m0), but this should be checked TODO
              wasLRTrun <- FALSE
              lrtp <- as.numeric(NA)
              if( "gamlss" %in% class(m1) & !is.null(dropVars) )
              {
                frm0 <- remove_terms(self$fixed, dropVars)
                m0   <- update(m1, frm0)

                # non-sig dropVars can lead to m1 and m0 with the same deviance which
                # caueses the LR.test error
                #
                # "The null model has smaller deviance than the alternative
                #  The models should be nested"
                #
                # do a check
                if(! m0$G.deviance <= m1$G.deviance) lrt  <- LR.test(m0, m1, print = FALSE)
                if(  m0$G.deviance <= m1$G.deviance) lrt  <- 1
                if(length(lrt)==3 & "p.val" %in% names(lrt))
                {
                  lrtp <- lrt$p.val
                }
                wasLRTrun <- TRUE
              }

              # output
              if("gamlss" %in% class(m1))
              {
                capture.output( tTable <- summary(m1), file='NUL')
                m1$tTable <- tTable
                m1$PalyticSummary <- self$summary()
                m1$whichPalyticMod <- paste('Palytic gamlss model #', wm)
                m1$lrt <- list(wasLRTrun=wasLRTrun, lrtp=lrtp)
                return(m1)
              }
              if("try-error" %in% class(m1))
              {
                return("Model did not converge")
              }

              options(warn = 0)

            },
            overwrite = TRUE
)

# this needs much editing
#' @description
#' The same as \code{getAR} when the ARMA order is desired for the full sample.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param doForeach Logical. Should the search for the best residual ARMA(p,q) order
#' be parallelized? This option can speed up searches when large values of P and Q
#' are specified in the autoSelect option of your Palytic object.
#' @param package Character. Options are \code{"nlme"} and \code{"gamlss"} as
#' described in the section documenting the \code{new()} method.
#' @param manual Logical. Should the \code{manual} option be used when making
#' clusters when \code{doForeach=TRUE}?
Palytic$set("public", "GroupAR",
            function(subgroup=NULL, doForeach=TRUE, package="nlme", manual=FALSE)
            {
              subCheck(subgroup, self$datac)

              # check for autoselect inputs
              if( is.null(self$autoSelect$AR$P) | is.null(self$autoSelect$AR$Q) )
              {
                stop("\nYour Palytic object does not contain both P and Q in",
                     "\nthe `autoSelect$AR` field.")
              }

              # if family is !NO override package
              if( self$family$family[2] != as.gamlss.family(NO)$family[2] )
              {
                package <- "gamlss"
                message("\nFamily is set as `", self$family$family[2], "` which",
                        "\nis not supported by package='nlme'.",
                        "\nSwitching to package='gamlss'.")
              }

              # prevent recursive calls to autoselect
              temp <- self$autoSelect
              temp$Done <- TRUE
              self$autoSelect <- temp; rm(temp)

              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "residual correlation structure starting...")
              start <- Sys.time()

              if(package %in% c("nlme", "arma")) nullMod <- self$lme(subgroup)
              if(package=="gamlss") nullMod <- self$gamlss(subgroup)

              if( class(nullMod) %in% c("lme", "gamlss")  )
              {
                DIMS <- expand.grid(p=0:self$autoSelect$AR$P,
                                    q=0:self$autoSelect$AR$Q)
                DIMS <- DIMS[-1,]

                capture.output( pb <- txtProgressBar(max = self$autoSelect$AR$P,
                                                     style = 3),
                                file='NUL')

                if(doForeach)
                {
                  cl       <- snow::makeCluster(parallel::detectCores(),
                                                type="SOCK", manual=manual)
                  snow::clusterExport(cl, list())
                  doSNOW::registerDoSNOW(cl)
                  pkgs     <- c("gamlss", "nlme")
                  progress <- function(n) setTxtProgressBar(pb, n)
                  opts     <- list(progress = progress)
                  t0       <- self$clone(deep=TRUE)

                  corMods <- foreach(p=DIMS$p, q=DIMS$q, .packages = pkgs,
                                     .options.snow = opts)  %dopar%
                  {
                    clone <- t0$clone(deep=TRUE)
                    ARpq(clone, p, q, subgroup, package)
                  }
                  parallel::stopCluster(cl)
                }
                if(!doForeach)
                {
                  corMods <- list()
                  for(i in seq_along(DIMS[,1]))
                  {
                    clone <- self$clone(deep=TRUE)
                    corMods[[i]] <- ARpq(clone, DIMS$p[i], DIMS$q[i],
                                         subgroup, package)
                  }
                }

                corMods <- corMods[!unlist(lapply(corMods, is.null))]

                names(corMods) <- unlist( lapply(corMods,
                                    function(x) x$PalyticSummary$correlation) )
                corMods <- corMods[!is.na(names(corMods))]

                if(self$whichIC[1]=="AIC") ICs <- data.frame( unlist( lapply(corMods, AIC) ) )
                if(self$whichIC[1]=="BIC") ICs <- data.frame( unlist( lapply(corMods, BIC) ) )
                else( stop( paste(self$whichIC, "is not a valid value for `whichIC`,",
                                  "use AIC or BIC.")))
                ICs <- rbind("NULL" = ifelse(self$whichIC[1]=="AIC", AIC(nullMod),
                                             BIC(nullMod)), ICs)
                bestCor <- c("NULL", names(corMods))[which.min(ICs[,1])]

              }
              if(! class(nullMod) %in% c("lme", "gamlss") )
              {
                bestCor <- "NULL"
              }

              message("\nAutomatic detection of the residual\n",
                      "correlation structure took: ",
                      capture.output(Sys.time() - start), ".\n\n")

              message("\nThe best correlation structure among those tested is ",
                      bestCor, "\n\n")

              self$correlation <- bestCor
            },
            overwrite = TRUE
)


#' ARpq - helper function for $GroupAR
#' @keywords internal
#' @author Stephen Tueller \email{stueller@@rti.org}
ARpq <- function(clone, p, q, subgroup, package="nlme")
{
  clone$method      <- "ML"
  clone$correlation <- "NULL"

  cortemp <- paste("nlme::corARMA(p=", p, ", q=", q, ")", sep="")
  cortemp <- gsub('\n', '', cortemp)
  cortemp <- gsub(' ', '', cortemp)
  clone$correlation <- cortemp
  if(package %in% c("nlme", "arma")) corMod <- clone$lme(subgroup)
  if(package=="gamlss") corMod <- clone$gamlss(subgroup)
  if( class(corMod) != "lme" )
  {
    corMod <- NULL
  }
  return(corMod)
}

#' @description
#' Use model comparisons to test for the polynomial order of time
#' @param package Character. Options are \code{"nlme"} and \code{"gamlss"} as
#' described in the section documenting the \code{new()} method.
Palytic$set("public", "getTO",
            function(package="nlme")
            {
              # extract polyMax and whichIC, which should give an error if
              # absent in $autoSelect
              polyMax <- self$autoSelect$TO$polyMax
              whichIC <- self$whichIC[1]

              # if family is !NO override package
              if( self$family$family[2] != as.gamlss.family(NO)$family[2] )
              {
                package <- "gamlss"
                message("\nFamily is set as `", self$family$family[2], "` which",
                        "\nis not supported by package='nlme'.",
                        "\nSwitching to package='gamlss'.")
              }

              # prevent recursive calls to autoselect
              temp <- self$autoSelect
              temp$Done <- TRUE
              self$autoSelect <- temp; rm(temp)

              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              uid <- sort( as.numeric( unique(self$datac[[self$ids]]) ) )
              time_powers <- list()

              for(id in uid)
              {
                aics <- mods <- list()
                for(i in 1:polyMax)
                {
                  clone <- self$clone(deep=TRUE)

                  clone$method     <- "ML"
                  clone$time_power <- i

                  if(package=="nlme")   mods[[i]] <- clone$lme(self$datac[[self$ids]]==id)
                  if(package=="gamlss") mods[[i]] <- clone$gamlss(self$datac[[self$ids]]==id)

                  if( class(mods[[i]]) %in% c("lme", "gamlss") )
                  {
                    aics[[i]] <- AIC( mods[[i]] )
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
                      capture.output(Sys.time() - start), ".\n",
                      "\nThe best polynomial order for time from 1 to ", polyMax,
                      " for each participant was ",
                      paste(self$time_powers[,2], collpase = ", ", sep=""))
            },
            overwrite = TRUE)

# GroupTO ####
# TODO:hard coded lme at this point, option for gamlss later
# needs much editing

#' @description
#' The same as \code{getTime_power} when the polynomial of time is desired for
#' the full sample.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param package Character. Options are \code{"nlme"} and \code{"gamlss"} as
#' described in the section documenting the \code{new()} method.
Palytic$set("public", "GroupTO",
            function(subgroup=NULL, package="nlme")
            {
              subCheck(subgroup, self$datac)

              # extract polyMax and whichIC, which should give an error if
              # absent in $autoSelect
              polyMax <- self$autoSelect$TO$polyMax
              whichIC <- self$whichIC[1]

              # if family is !NO override package
              if( self$family$family[2] != as.gamlss.family(NO)$family[2] )
              {
                package <- "gamlss"
                message("\nFamily is set as `", self$family$family[2], "` which",
                        "\nis not supported by package='nlme'.",
                        "\nSwitching to package='gamlss'.")
              }

              # prevent recursive calls to autoselect
              temp <- self$autoSelect
              temp$Done <- TRUE
              self$autoSelect <- temp; rm(temp)

              message("\n\nPersonAlytics: Automatic detection of the\n",
                      "time/outcome relationship starting...")
              start <- Sys.time()

              mods <- list()
              for(i in 1:polyMax)
              {
                clone <- self$clone(deep=TRUE)

                clone$method     <- "ML"
                clone$time_power <- i

                if(package %in% c("nlme", "arma")) mods[[i]] <- clone$lme(subgroup)
                if(package=="gamlss") mods[[i]] <- clone$gamlss(subgroup)
                if(! any( class(mods[[i]]) %in% c("lme", "gamlss") ) )
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
                      capture.output(Sys.time() - start), ".\n",
                      "\nThe best polynomial order for time from 1 to ", polyMax,
                      " was ", self$time_power)
            },
            overwrite = TRUE
            )

# This functions borrows from ICTviz() in PersonAlyticsPower, but the nature
# of the data for ICTviz is theoretical, this function is for real data

#' @description
#' Plot the data in your Palytic object with design elements including time and
#' phase.
#'
#' @param subgroup Logical vector. If \code{NULL} results are provided for the full
#' sample. If a logical vector of the same length as the nmuber of rows in the data,
#' the results are provided for the subgroup of cases for who \code{subgroup==TRUE}.
#' @param groupVar Character. The name of a grouping variable in the data that
#' will be used to stratify the plot.
#' @param type Character. See \code{\link{plotICT}}.
#' @param ylim See \code{\link{plotICT}}.
#' @param title Character. A title for you plot.
Palytic$set("public", "plot",
            function(subgroup=NULL, groupvar=NULL, type='histogram', ylim=NULL,
                     title=NULL)
            {
              subCheck(subgroup, self$datac)

              # qc input
              if(!is.null(groupvar))
              {
                if(length(groupvar)!=1 | !is.character(groupvar))
                {
                  stop('`groupvar` should be one variable name.')
                }
                if(! groupvar %in% names(self$datac))
                {
                  stop('`groupvar` is not in the data.')
                }
              }

              # subset the data if requested
              if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(self$datac))
              tempData <- subset(self$datac, subgroup,
                                 c(all.vars(self$formula), groupvar))

              if(!is.null(groupvar)) ug <- unique(tempData[[groupvar]])
              if( is.null(groupvar)) ug <- 1
              dens <- traj <- list()
              for(i in seq_along(ug))
              {
                legendName <- NULL
                if(length(ug)!=1) legendName <- paste(groupvar, ug[i], sep='=')
                if(!is.null(groupvar)) wg <- tempData[[groupvar]]==ug[i]
                if( is.null(groupvar)) wg <- rep(TRUE, nrow(tempData))
                plotICTs  <- plotICT(self, tempData[wg,],
                                     legendName = legendName,
                                     type = type,
                                     ylim = ylim)
                dens[[i]] <- plotICTs$d
                traj[[i]] <- plotICTs$s
                rm(plotICTs)
              }
              plotICTs <- c(dens, traj)

              top <- paste(title[[1]],
                           'Density and Average Trajectory with SD bars',
                           sep=": ")

              suppressMessages(
              suppressWarnings(
              gridExtra::marrangeGrob(plotICTs, nrow=length(dens),
                                      ncol=2, widths = c(3,6),
                                      top=top)
              ))
            },
            overwrite = TRUE
            )



