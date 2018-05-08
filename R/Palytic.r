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

#library(gamlss) # the formula fails including gamlss::re
source("./r/eds.r")
source("./r/.clean.r")
source("./r/.active.r")


Palytic <- R6::R6Class("Palytic",
  private = list(
    .data       = NULL, # consider pass by reference environment
    .ids        = NULL,
    .y          = NULL,
    .phase      = NULL,
    .time       = NULL,
    .ivs        = NULL,
    .time_power = 1,
    .correlation= NULL,
    .family     = NULL,
    .fixed      = NULL,
    .random     = NULL,
    .formula    = NULL,
    .method     = NULL,
    .is_clean   = FALSE,
    .warnings   = list(),
    .errors     = list(),
    .try_silent = TRUE
  ),

  active = .active(),

  public = list(
    initialize = function
    (
      data        = NULL,
      ids         = NULL,
      y           = NULL,
      phase       = NULL,
      time        = NULL,
      ivs         = NULL,
      time_power  = 1,
      correlation = NULL,
      family      = gamlss.dist::NO(),
      fixed       = NULL,
      random      = NULL,
      formula     = NULL,
      method      = "ML",
      is_clean    = FALSE,
      warnings    = list(), # can we hide these? or just make them read only?
      errors      = list(),
      try_silent = TRUE
    )
    {
      if(is.null(fixed)) fixed <- formula( paste(y, "~", phase, "*", time) )
      if(is.null(random)) random <- formula( paste("~", time, "|", ids))
      if(is.null(formula))
      {
        formula <- formula( paste(deparse(fixed),
                     "+ re(random = ",
                     #"+ gamlss::re(random = ", # produces error, require(gamlss)
                     deparse(random),
                     ", method=",
                     deparse(method),
                     ", correlation =",
                     deparse(correlation),
                     ")")
                    )
      }

      # may need to add a check whether fixed, random, and formula conform,
      # e.g., if a user provides a formula not implied by y, phase, time, ids
      # but does not provide fixed and random


      # if we leave data as read only, it must me done here, otherwise
      # active applies making it read only
      #  -- na.omit should be applied here in some way
      #if(! is.data.frame(data)) stop("data must be a data.frame")
      #data$timmy <- data[,1]^2

      private$.data        <- data # add data cleaning, if any,here (e.g. phase.align)

      private$.ids         <- ids
      private$.y           <- y
      private$.phase       <- phase
      private$.time        <- time
      private$.ivs         <- ivs
      private$.time_power  <- time_power
      private$.correlation <- correlation
      private$.family      <- family
      private$.fixed       <- fixed
      private$.random      <- random
      private$.formula     <- formula
      private$.method      <- method
      private$.is_clean    <- is_clean
      private$.warnings    <- warnings
      private$.errors      <- errors
      private$.try_silent  <- TRUE

    }

  )

)


Ovary <- as.data.frame(nlme::Ovary)
Ovary$Mare <- factor(Ovary$Mare, ordered = FALSE)
Ovary$Phase <- as.numeric(Ovary$Time > .5)
Ovary$TimeSin <- sin(2*pi*Ovary$Time)

t0 <- Palytic$new(data=Ovary, ids="Mare", y="follicles", phase="Phase",
                  time="TimeSin")
t0$fixed
t0$random
t0$formula
t0$method

# add methods
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
                        silent = self$try_silent)

              if( "try-error" %in% class(m1) | !eds(m1) )
              {
                ctrl <- nlme::lmeControl(opt="optim")
                m1 <- try(nlme::lme(fixed=self$fixed,
                                    data=self$data[w,],
                                    random=self$random,
                                    correlation=self$correlation,
                                    method=self$method,
                                    control=self$ctrl,
                                    keep.data=FALSE,
                                    na.action=na.omit),
                          silent = self$try_silent)
              }
              if( "lme" %in% class(m1) )
              {
                m1$call$fixed  <- self$fixed
                m1$call$random <- self$random
                return(m1)
              }
              else
              {
                return("Model did not converge")
              }
            },
            overwrite = TRUE
)



### this needs to be expanded to include our list, e.g.,
#   R:\PaCCT\Process\MMTA Process and Record Keeping.docx
# -- this will be done in Palytic augmentations in .Palytic
Palytic$set("public", "gamlss",
            function(w=NULL)
            {
              if(is.null(w)) w <- 1:nrow(self$data)
              require(gamlss)
              m1 <- try(gamlss::gamlss(formula = self$formula,
                               data    = self$data[w,],
                               family  = self$family),
                      silent = self$try_silent)

              if("gamlss" %in% class(m1))
              {
                m1$call$formula <- self$formula
                m1$family <- self$family
                return(m1)
              }
              if("try-error" %in% class(m1))
              {
                return("Model did not converge")
              }

            },
            overwrite = TRUE
)


# check with debug
#Palytic$debug("lme")
#Palytic$debug("gamlss")
# check without debug, require reinitializing object
#Palytic$undebug("lme")
#Palytic$undebug("gamlss")


t0 <- Palytic$new(data=Ovary, id="Mare", y="follicles", phase="Phase",
                  time="TimeSin", try_silent=FALSE)

t1 <- Palytic$new(data=Ovary, fixed = formula("follicles ~ TimeSin*Phase"),
                  random = formula(" ~ TimeSin | Mare"), try_silent=FALSE)
t1$gamlss()
t1$lme()

# single subject test
w <- which(Ovary$Mare==1)
t1$gamlss(w)
t1$gamlssw(w)
t1$lme(w)







