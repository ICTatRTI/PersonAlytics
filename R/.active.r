
.active <- function()
{

  list(
    data = function(value)
    {
      if( missing(value) ){ private$.data }
      else
      {
        stopifnot(is.data.frame(value) | is.matrix(value))
        #stop("`$fixed` is read only", call. = FALSE)
        print('jigglywiggly') # the data function in .active isn't invoked on data given at init
        private$.data <- value
        self
      }
    },

    fixed = function(value)
    {
      if( missing(value) ){ private$.fixed }
      else
      {
        stopifnot("formula" %in% class(value))
        #stop("`$fixed` is read only", call. = FALSE)
        private$.fixed <- value
        self
      }
    },

    random = function(value)
    {
      if( missing(value) ){ private$.random }
      else
      {
        stopifnot("formula" %in% class(value))
        #stop("`$random` is read only", call. = FALSE)
        private$.random <- value
        self
      }
    },

    time_power = function(value)
    {
      if( missing(value) ){ private$.time_power }
      else
      {
        stopifnot(is.numeric(value))
        private$.time_power <- value
        self
      }
    },

    ar_order = function(value)
    {
      if( missing(value) ){ private$.ar_order }
      else
      {
        stopifnot(is.numeric(value))
        private$.ar_order <- value
        self
      }
    },

    ma_order = function(value)
    {
      if( missing(value) ){ private$.mr_order }
      else
      {
        stopifnot(is.numeric(value))
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
