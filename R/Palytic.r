#'
#'
#'
#'
#'
#'

library(gamlss) # the formula fails including gamlss::re
source('./r/eds.r')

Palytic <- R6::R6Class("Palytic",
             public = list
             (
                       errors     = NULL,
                       warnings   = NULL,
                       data       = NULL, # consider pass by reference environment
                       fixed      = NULL,
                       random     = NULL,
                       time_power = NULL,
                       ar_order   = NULL,
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

Palytic$set("public", "gamlss",
            function(w=NULL)
            {
              if(is.null(w)) w <- 1:nrow(self$data)
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
t1$lme(w)







