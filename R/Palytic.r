#'
#'
#'
#'
#'
#'

# so we need to separate public and private parts. regardless of whether we
# enscrypt parts of the code, we can use inheritance to separate free and
# private parts of the Palytic class, e.g., define Palytic in the free
# package, and in the paid package, update palytic through inheritance
# http://www.orbifold.net/default/2015/04/24/r6-classes/

library(gamlss) # the formula fails including gamlss::re
source('./r/eds.r')


Palytic <- R6::R6Class("Palytic",
  private = list(
    .data       = NULL, # consider pass by reference environment
    .fixed      = NULL,
    .random     = NULL,
    .time_power = 1,
    .ar_order   = 1,
    .ma_order   = 0,
    .family     = NO()
  ),



  # use active fields to do simple data validation
  # https://adv-r.hadley.nz/r6
  # (use and $validate() field for complex validation, do this for data)
  active = list(
    fixed = function(value)
    {
      if( missing(value) ){ private$.fixed }
      else
      {
        stopifnot(is.formula(value))
        private$.fixed <- value
        self
      }
    },

    random = function(value)
    {
      if( missing(value) ){ private$.random }
      else
      {
        stopifnot(is.formula(value))
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
    }
  ),

  public = list(
    initialize = function
    (
      data       = NULL,
      fixed      = NULL,
      random     = NULL,
      time_power = 1,
      ar_order   = 1,
      ma_order   = 0,
      family     = gamlss.dist::NO()
    )
    {
      private$.data       <- data
      private$.fixed      <- fixed
      private$.random     <- random
      private$.time_power <- time_power
      private$.ar_order   <- ar_order
      private$.ma_order   <- ma_order
      private$.family     <- family
    }
  )

)

t0 <- Palytic$new()
t0$.fixed
t0test <- Palytic$new(.fixed='juan')
t0test$.fixed

##### below this line should be migrated or depricated

Palytic <- R6::R6Class("Palytic",

             public = list
             (
                       errors     = NULL,
                       warnings   = NULL,

                       nlme0      = NULL,
                       method     = NULL,
                       family     = NULL,
                       TrySilent  = NULL,

                       initialize = function(errors     = list(),
                                             warnings   = list(),
                                             data       = NULL,
                                             fixed      = NULL,
                                             random     = NULL,
                                             time_power = 1,
                                             ar_order   = 1,
                                             nlme0      = NULL,
                                             method     = "ML",
                                             family     = NO(),
                                             TrySilent  = TRUE)
                       {
                         self$errors     = errors
                         self$warnings   = warnings
                         self$data       = data
                         self$fixed      = fixed
                         self$random     = random
                         self$time_power = time_power
                         self$ar_order   = ar_order
                         self$nlme0      = nlme0
                         self$method     = method
                         self$family     = family
                         self$TrySilent  = TrySilent
                       } # eof initialize
              ) # eof public
)


Palytic$set("public", "lme",
            function(w=NULL)
            {
              if(is.null(w)) w <- 1:nrow(self$data)
              m1 <- try(nlme::lme(fixed=self$fixed,
                                  data=self$data[w,],
                                  random=self$random,
                                  correlation=self$correlation,
                                  method=self$method,
                                  keep.data=FALSE,
                                  na.action=na.omit),
                        silent = self$TrySilent)

              if( 'try-error' %in% class(m1) | !eds(m1) )
              {
                ctrl <- nlme::lmeControl(opt='optim')
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=self$data[w,],
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=self$ctrl,
                                    keep.data=FALSE,
                                    na.action=na.omit),
                          silent = self$TrySilent)
              }
              if( 'lme' %in% class(m1) )
              {
                m1$call$fixed  <- self$fixed
                m1$call$random <- self$random
                return(m1)
              }
              else
              {
                return('Model did not converge')
              }
            },
            overwrite = TRUE
)



Palytic$set("public", "gamlss_formula",
            function()
            {
              # Note that there is no gamlss::re in frm, but there is
              # gamlss::gamlss in the model fitting. I'm still trying to figure
              # out why having the former causes a crash, while failing to have
              # the latter also causes a crash, yet the discrepancy works
              frm <- paste(deparse(self$fixed),
                           '+ re(random = ',
                           deparse(self$random),
                           ', method=',
                           deparse(self$method),
                           ', correlation =',
                           deparse(self$correlation),
                           ')')
              frm <- formula(frm)
              return(frm)
            },
            overwrite = TRUE
)


### this needs to be expanded to include our list, e.g.,
#   R:\PaCCT\Process\MMTA Process and Record Keeping.docx
Palytic$set("public", "gamlss",
            function(w=NULL)
            {
              if(is.null(w)) w <- 1:nrow(self$data)


              # I don't like having na.omit here, you need to test
              # 1. whether it does the whole data set (likely) or just the
              #    variables in the analysis
              # 2. if the former, missing data treatment needs to be augmented,
              #    e.g., you could add a function to pre-treat the data based
              #    on variables in the analysis

              m1 <- try(gamlss::gamlss(formula = frm,
                               data    = self$data[w,], # no missing data handling yet
                               family  = self$family),
                      silent = self$TrySilent)

              if('gamlss' %in% class(m1))
              {
                m1$call$formula <- frm
                m1$family <- self$family
                return(m1)
              }
              if('try-error' %in% class(m1))
              {
                return('Model did not converge')
              }


            },
            overwrite = TRUE
)

Palytic$set("public", "gamlssw",
             function(w)
             {
               self$gamlss(w)
             },
             overwrite = TRUE
)

Ovary <- as.data.frame(nlme::Ovary)
Ovary$Mare <- factor(Ovary$Mare, ordered = FALSE)
Ovary$Phase <- as.numeric(Ovary$Time > .5)
Ovary$TimeSin <- sin(2*pi*Ovary$Time)

# check with debug
#Palytic$debug("lme")
#Palytic$debug("gamlss")
# check without debug, require reinitializing object
#Palytic$undebug("lme")
#Palytic$undebug("gamlss")

t1 <- Palytic$new(data=Ovary, fixed = formula("follicles ~ TimeSin*Phase"),
                  random = formula(" ~ TimeSin | Mare"), TrySilent=FALSE)
t1$gamlss()
t1$lme()

# single subject test
w <- which(Ovary$Mare==1)
t1$gamlss(w)
t1$gamlssw(w)
t1$lme(w)







