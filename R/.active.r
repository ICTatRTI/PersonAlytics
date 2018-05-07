#'
#'
#'

# When does it make sense to use active bindings?
# - if you need to check the data each time the field is accessed
# - e.g., yes for formula b/c it may get changed

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

    y = function(value)
    {
      if( missing(value) ) private$.y
      else
      {
        if(! is.character(value) )
        {
          stop('y must be a character variable name in data')
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, 'is not in the data') )
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
          stop('phase must be a character variable name in the data')
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, 'is not in the data') )
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
          stop('time must be a character variable name in the data')
        }
        if( is.null(private$.data[[value]]) )
        {
          stop( paste(value, 'is not in the data') )
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
          stop('ivs must be a character vector of variables in the data')
        }
        if( ! all(value %in% names(private$.data) ) )
        {
          nov <- value[ which(! value %in% names(private$.data )) ]
          if(length(nov)==1) stop( paste(nov, 'is not in the data') )
          if(length(nov)>=2)
          {
            stop( paste(paste(nov, collapse=', '), 'are not in the data') )
          }
        }
        private$.ivs <- value
        self
      }
    },

    time_power = function(value)
    {
      if( missing(value) ){ private$.time_power }
      else
      {
        if(! is.numeric(value) ) stop('time_power must be numeric')
        if( round(value, 0) != value | vaue < 1 )
        {
          stop('time_power must be a positive whole number')
        }
        private$.time_power <- value
        self
      }
    },

    ar_order = function(value)
    {
      if( missing(value) ){ private$.ar_order }
      else
      {
        if(! is.numeric(value) ) stop('ar_order must be numeric')
        if( round(value, 0) != value | vaue < 0 )
        {
          stop('ar_order must be a positive whole number or 0')
        }
        private$.ar_order <- value
        self
      }
    },

    ma_order = function(value)
    {
      if( missing(value) ){ private$.mr_order }
      else
      {
        if(! is.numeric(value) ) stop('ma_order must be numeric')
        if( round(value, 0) != value | vaue < 0 )
        {
          stop('ma_order must be a positive whole number or 0')
        }
        private$.ma_order <- value
        self
      }
    },

    family = function(value)
    {
      if( missing(value) ){ private$.family }
      else
      {
        stopifnot("gamlss.family" %in% class(value))
        private$.family <- value
        self
      }
    },

    is_clean = function(value)
    {
      if( missing(value) ){ private$.is_clean }
      else
      {
        stopifnot("logical" %in% class(value))
        private$.is_clean <- value
        self
      }
    }
  )

}
