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
source('./r/.clean.r')
source('./r/.active.r')


Palytic <- R6::R6Class("Palytic",
  private = list(
    .data       = NULL, # consider pass by reference environment
    .fixed      = NULL,
    .random     = NULL,
    .time_power = 1,
    .ar_order   = 1,
    .ma_order   = 0,
    .family     = NO(),
    .is_clean   = FALSE,
    .warnings   = list(),
    .errors     = list()
  ),

  active = .active(),

  public = list(
    initialize = function
    (
      data       = NULL,
      fixed      = NULL,
      random     = NULL,
      time_power = 1,
      ar_order   = 1,
      ma_order   = 0,
      family     = gamlss.dist::NO(),
      is_clean   = FALSE,
      warnings   = list(),
      errors     = list()
    )
    {
      # if we leave data as read only, it must me done here, otherwise
      # active applies making it read only
      if(! is.data.frame(data)) stop('data must be a data.frame')
      data$timmy <- data[,1]^2

      private$.data       <- data
      private$.fixed      <- fixed
      private$.random     <- random
      private$.time_power <- time_power
      private$.ar_order   <- ar_order
      private$.ma_order   <- ma_order
      private$.family     <- family
      private$.is_clean   <- FALSE
      private$.warnings   <- list() # can we not put this in public and still use them?
      private$.errors     <- list()

      # move simple checks and/or 1 time checks
      # I'm thinking we should make data immutable, you give it once
      # then you can do a lot of things to it (active can check each new
      # formula, especially in updates from .Palytic), but you don't change
      # the data b/c there are expensive checks that should only be done once.
      # Don't think of checking data for variables in the formula, check the
      # formula for presence in the data.

    },

    # note that data cleaning won't happen if the data or formula are
    # overwritten, can we dispatch a method whenever something gets changed?
    # it looks like active binding invokes it's checks, so we don't need them
    # here, commit and test
    clean = function()
    {
        self$data     <- .clean(self$data, self$fixed, self$random)
        self$is_clean <- TRUE
    },

    test = function()
    {
      if( !self$is_clean ) self$clean()
    },

    hist = function()
    {
      hist(self$data[,1])
    }

  )

)

# active binding tests
# https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html#active-bindings
# "Active bindings look like fields, but each time they are accessed,
# they call a function. They are always publicly visible."
t0 <- Palytic$new()
t0$ar_order
t0$ar_order <- 3
t0$ar_order
t0$ar_order <- 'test' # error as it should be
t0$fixed
t0$fixed <- 'test' # error as it should be
t0$fixed <- formula('follicles ~ TimeSin * Phase')
t0$fixed
rm(t0)

# test active binding and clean function
#Palytic$debug('clean')
t0 <- Palytic$new()
t0$hist() # this attempts to use self$data, but does access the data function
t0$clean()# this fails b/c data is not the right value
t0$is_clean

t0 <- Palytic$new(data='bob') # active binging not invoked, apparent the field isn't being 'accesssed'
t0$hist() # nor here, it fails b/c it is the wrong type of data, not b/c check was invoked

t0$data <- 'bob'
t0$hist()

t0$data <- matrix(rnorm(300), ncol=3) # this invokes the active binding test
t0 <- Palytic$new(data = data.frame(matrix(rnorm(300), ncol=3))) # this invokes the active binding test
t0$hist()
t0$clean()
t0$is_clean
t0 <- Palytic$new(data=Ovary, fixed = formula("follicles ~ TimeSin*Phase"),
                  random = formula(" ~ TimeSin | Mare")) # this does not invoke active binding
t0$clean() # this invokes the active binding test
t0$test() # but this does not


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







