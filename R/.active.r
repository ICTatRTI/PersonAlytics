#'
#'
#'

# When does it make sense to use active bindings?
# - if you need to check the data each time the field is accessed
# - e.g., yes for formula b/c it may get changed
# the missing() first step just indicates to print the current value back,
# the rest implements checks if the user tries updating the values

.active <- function()
{

  list(
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
        private$.ids <- value
        self
      }
    },

    y = function(value)
    {
      if( missing(value) ) private$.y
      else
      {
        if(! is.character(value) )
        {
          stop("`y` must be a character variable name in data")
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, "is not in the data") )
        }
        private$.y <- value
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
        private$.phase <- value
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
        private$.time <- value
        self
      }
    },

    ivs = function(value)
    {
      if( missing(value) ) private$.ivs
      else
      {
        if(! is.character(value) )
        {
          stop("`ivs` must be a character vector of variables in the data")
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
        private$.ivs <- value
        self
      }
    },

    time_power = function(value)
    {
      if( missing(value) ) private$.time_power
      else
      {
        if(! is.numeric(value) ) stop("`time_power` must be numeric")
        if( round(value, 0) != value | vaue < 1 )
        {
          stop("`time_power` must be a positive whole number")
        }
        private$.time_power <- value
        self
      }
    },

    correlation = function(value)
    {
      if( missing(value) ) private$.correlation
      else
      {
        if(! "corStruct" %in% class(value) )
        {
          stop( paste("`correlation` must be of class corStruct, see `?nlme::corStruct`",
                      "For AR(1) use `corAR1()` or `corARMA(p=1)`. For are 2 or higher",
                      "use `corARMA(p=2)`. For ARMA(3,1) use `corARMA(p=3,q=1)`.",
                      sep="\n") )
        }
        private$.ar_order <- value
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
        private$.family <- value
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
        private$.fixed <- value
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
        private$.random <- value
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
          private$.formula <- value
          self
        }
      }
    },

    method = function(value)
    {
      if( missing(value) ) private$.method
      else
      {
        if(! c('ML', 'REML') %in% value )
        {
          stop("`method` should be `ML` or `REML`")
          private$.method <- value
          self
        }
      }
    },

    is_clean = function(value)
    {
      if( missing(value) ){ private$.is_clean }
      else
      {
        if(! "logical" %in% class(value)) stop("`is_clean` should be logical")
        private$.is_clean <- value
        self
      }
    }
  )

}
