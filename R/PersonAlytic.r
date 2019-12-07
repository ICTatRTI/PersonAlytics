#' Personalytic, a simplified user interface for linear and generalized linear
#' mixed effects models with high throughput options.
#'
#' @export
#' @import snow
#' @import doSNOW
#' @import foreach
#' @import plyr
#'
#' @description A simplified user interface for longitudinal linear mixed
#' effects models (also known as growth models or hierarchical linear models) using
#' \code{\link{nlme}} and generalized linear mixed effects models using
#' \code{\link{gamlss}} (not currently implemented).
#' The basic mixed effects model is \eqn{dv=time+phase+phase*time}
#' with random intercepts and random slopes for time. The phase variable is optional.
#' Additional independent variables (or covariates) can be included.
#'
#' High throughput capabilities are invoked when the model
#' is to be fit to multiple individuals, multiple dependent variables, iterated
#' over multiple target independent variables, or any combination of these three.
#'
#' The residual correlation structure \code{ARMA(p,q)} is automatically selected
#' among \code{p=1,...,P} and \code{q=1,...,Q}
#' where \code{P} and \code{Q} are selected by the user. Model selection is done
#' by comparing fit indices using either the \code{BIC} or \code{AIC} (ARMA models
#' are not generally nested, precluding the use of likelihood ratio tests).
#'
#' The functional form of the relationship between time and the dependent variable
#' is automatically selected up to \code{time^time_power}  where \code{time_power} is
#' selected by the user. Model selection is done using likelihood
#' ratio tests (LRT) with maximum likelihood estimators. For example,
#' if \code{time_power=3} (implying a maximum of a cubic growth model), linear,
#' quadratic and cubic models will be fit. If the best fitting model is a quadratic
#' growth model, then there will be random and fixed effects for \code{time} and
#' \code{time^2}.
#'
#' When there are multiple dependent variables and multiple target independent
#' variables, Type I error corrections or false discovery rate corrections are
#' made across target independent variables within each dependent variable (and
#' within each person if individual models are being fit to the data). Correction
#' options are given in \code{\link{p.adjust}}.
#'
#' @param output Character. The default is \code{NULL}.
#'
#' A character string that will be used to name a file for saving
#' output. If left \code{NULL}, the default is
#' `PersonAlytic_Output`. Do not give a file extension, these will be added
#' automatically. For example, if  \code{output='MyResults'}, the output file will
#' be called \code{'MyResults.csv'} if high throughput options are used, or
#' \code{'MyResults.txt'} if only one model is fit. A full path with `/` instead
#' of `\` can also be used. For example, \code{output='C:/MyResults'} will produce
#' the files \code{'C:/MyResults.csv'} or \code{'C:/MyResults.txt'}. If a full
#' path is not given, the results will be saved in the working directory (see
#' \code{\link{getwd}}).
#'
#' @param data A \code{\link{data.frame}}. \code{data} must be provided by the user.
#'
#' \code{data}  that contains the variables \code{ids},
#' \code{dv}, \code{time}, and optionally, \code{phase}, \code{ivs},
#' \code{target_ivs}.
#'
#' @param ids Character. \code{ids} must be provided by the user.
#'
#' The name of the ID variable. The ID variable must be a variable in \code{data}
#' and must be numeric.
#'
#' @param dvs Character list. \code{dvs} must be provided by the user.
#'
#' A list of one or more character dependent variable names in \code{data}.
#' The linear mixed effects model
#' \code{dvs[d] ~ phase + time + phase*time + target_ivs[c] + ivs}
#' with random effects \code{~ time | ids[i]} will be fit to the data.
#' The iterators \code{[d]}, \code{[c]}, and \code{[i]}
#' indicate that the model will be fit for each combination of
#' dependent variable in \code{dvs},
#' independent variable in \code{target_ivs}, and
#' each unique ID in \code{ids} (which can be overridden using
#' \code{individual_mods=FALSE}) with each model
#' controlling for the independent variables in \code{ivs}.
#'
#' @param time Character. \code{time} must be provided by the user.
#'
#' The name of the time variable. The time variable must be a variable in \code{data}
#' and must be numeric.
#'
#' @param phase Character. The default value is \code{NULL}.
#'
#' Name of the phase or treatment variable. For model fitting,
#' \code{phase} is treated the same as any variable names in \code{ivs} or
#' \code{target_ivs} but is used for visualizing treatment effects.
#'
#' @param ivs Character list. The default value is \code{NULL}.
#'
#' A list of names of covariates, e.g., \code{list('iv2', 'iv2')}.
#' Note that the variables in \code{ivs} cannot also be in  \code{target_ivs}.
#'
#' @param target_ivs Character list. The default value is \code{NULL}.
#'
#' Independent variables that are iterated over  (see \code{dv} for details).
#' Effects for these variables are labeled as 'target predictor' in the output.
#'
#' @param interactions Character list. The default value is \code{NULL}.
#'
#' List of pairs of variable names for interactions
#' to include in the model, e.g., \code{list(c('iv1','phase'), c('iv1','iv2'))}.
#'
#' @param time_power Numeric. The default is \code{time_power=3}.
#'
#' Power of the time variable (e.g., \code{time^time_power}).
#' A quadratic or cubic growth model would be specified using
#' \code{time_power=2} or \code{time_power=3}, respectively.
#' If \code{autoDetect$TO$polyMax} is specified, \code{polyMax} is the largest
#' value tested for. The default value is \code{3}, testing up to a cubic growth
#' model. If a linear growth model is desired, set
#' \code{autoDetect$TO=NULL} and \code{time_power=1}.
#'
#' @param correlation Character. The default value is \code{NULL}.
#'
#' See \code{\link{corStruct}} in \code{\link{nlme}}.
#' Must be passed as a character, e.g. \code{"corARMA(p=1)"}.
#' The default value is \code{NULL}, assuming no residual
#' autocorrelation. If \code{auoDetect$AR} is specified, \code{correlation} will
#' be ignored.
#'
#' @param family See \code{\link{gamlss.family}}. The default is normal,
#' \code{family=NO()}. This option is not yet implemented.
#'
#' A list of the same length as \code{length(dv)} can be supplied. For example if
#'\code{dv=list('normal_y', 'binary_y', 'beta_y')}, then family could be
#' \code{family=list(NO(), BI(), BEINF())}. If  \code{length(dv)>1} and
#' \code{length(family)==1}, the one distribution in \code{family} will be
#' applied to all outcomes.
#' The \code{family} parameter is ignored if \code{package="nlme"}.
#'
#' @param subgroup Logical vector. The default is \code{subgroup=NULL} wich
#' results in a logical vector of the same length as the number of rows in
#' \code{data} with all values equal to \code{TRUE}.
#'
#' A vector where \code{length(subgroup)==nrow(data)} indicating
#' which subset of the data should be used for analysis. For example, if a model
#' should only be fit to females, \code{subgroup=gender=='female'} might be used.
#'
#' @param standardize Named logical vector. The default is
#' \code{list(dv=FALSE, ivs=FALSE, byids=FALSE)}.
#'
#' Which variables should be standardized?
#' (i.e., rescaled to have 0 mean and unit variance; see
#' \code{\link{scale}})? See \code{dv} and \code{ivs}. The option
#' \code{byids} controls whether standardization is done by individuals or by group.
#' Does not apply to factor variables. The default is \code{TRUE}.
#' Standardization makes parameter estimate magnitudes comparable across
#' individuals, outcomes in
#' \code{dvs}, and covariates in \code{target_ivs}. For dependent variables in
#' \code{dvs}, standardization is only applied for normal outcomes, see \code{family}.
#'
#' @param package Character. The default is \code{"nlme"}.
#'
#' Which package should be used to fit the models?
#' Options are \code{\link{nlme}}, \code{\link{gamlss}}, and \code{\link{arma}}.
#' It is passed as character strings, e.g., \code{"gamlss"} or \code{"nlme"}. If
#' there is only one participant in \code{ids}, "arma" will be used.
#'
#' @param method character. The default is \code{"REML"}.
#'
#' Which likelihood methods should be used to fit the models? Options are
#' \code{"REML"}, which should be used for final parameter estimates and
#' effect sizes, and \code{"ML"}, which should be used for model comparisons
#' (e.g., this is done for automatic residual correlation structure detection
#' and automatic time order detection).
#'
#' @param individual_mods Logical. The default is \code{individual_mods=FALSE}.
#'
#' Should individual models be fit for each ID in \code{ids}?
#'
#' @param PalyticObj See \code{\link{Palytic}}. Not currently implemented.
#'
#' If \code{PalyticObj} is submitted
#' then only \code{dvs}, \code{target_ivs}, and \code{individual_mods} will be
#' used. This allows users access to additional options including generalized linear
#' mixed effects models via the \code{family} option, user specified \code{correlation}
#' structures (in non-\code{NULL} this will override the automated correlation structure
#' search), and user specified models via \code{formula}.
#'
#' @param autoDetect List. The default is
#' \code{
#' list(AR=list(P=3, Q=3)     ,
#'   TO=list(polyMax=3)       ,
#'   DIST=list()) }.
#'
#' If no automated model selection for the residual covariance structure (\code{AR}),
#' the polynomial order for the relationship between time and the dependent variable
#' (\code{TO}), or the dependent variable distribution is desired, an empty list
#' should be passed (e.g., \code{autoDetect=list()}).
#'
#' If \code{AR} is in the list,
#' the residual correlation structure will be automatically selected from
#' among \code{ARMA(p,q)} models? See \code{correlation}. Since these models are
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
#' method in \code{\link{Palytic}}.
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
#' method in \code{\link{Palytic}}.
#'
#' If \code{DIST} is in the list and \code{package='gamlss'}, each dependent
#' variable in \code{dvs} will utilize the \code{\link{fitDist}} function of
#' the gamlss package, and the best fitting distribution will be used for each
#' depedent variable. For more detail, see the \code{$dist()} method in
#' \code{\link{Palytic}}. To narrow the distributions that will be tested,
#' the user must specify whether to
#' rescale the dependent variable to the (0,1) range with \code{to01}.
#' Currently settings
#' for \code{DIST} apply to all depedent variables in \code{dvs}. This will be
#' generalized to dependent variable specifics settings in a future release. If
#' your dependent variables require different \code{DIST} settings, use separate
#' calls to \code{PersonAlytic}.
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
#'
#' @param charSub Character list. The default in \code{charSub=NULL}.
#'
#' A list of paired character strings for character substitution
#' in the output. If the names of the target predictors
#' in \code{target_ivs} had to be edited to make valid variable names, this
#' parameter allows users put the illegal characters back in. For example,
#' if the original variable name was "17.00_832.2375m/z", a letter would need to
#' prefix the variable name and the
#' "/" would need to be replaced with another character, e.g., "X17.00_832.2375m.z".
#' To get the row names of the output back to original variable name, use
#' \code{charSub=list(c("X", ""), c("m.z", "m/z"))}. Note that inputs to charSub
#' must be in double quotes and are case sensitive. All duplicates will be substituted.
#' For example, if the variable name was "X1X23.x" and \code{charSub=list(c("X", ""))},
#' the resulting row label for this variable would be "123.x".
#'
#' @param sigma.formula Not currently implemented.
#'
#' A formula for the variance under \code{\link{gamlss}}.
#' Currently static: it will not change dynamically over iterations nor will it be
#' updated by \code{time_power} or \code{autoDetect}. If model fitting using this
#' option fails, another attempt will be made after reseting it to its defaul,
#' i.e., \code{~1}.
#'
#' @param p.method See \code{\link{p.adjust.methods}}. When \code{individual_mods=TRUE},
#' \code{length(dvs)>1}, or \code{length(target_ivs)>1}, p-value adjustments
#' are made and reported in separate output saved to the working directory.
#'
#' @param alpha Numeric value in the (0,1) interval. The Type I error rate for
#' adjusting p-values.
#'
#' @param nbest Numeric integer value. The number of best \code{target_ivs} to
#' report on before and after adjusting p-values. This is only used if
#' \code{length(target_ivs)>1} and the \code{target_ivs} are all continuous.
#'
#' @param alignPhase Character. The default is \code{'none'} which leaves time
#' as is.
#'
#' Other options are \code{'piecewise'} which leads to a piecewise growth model
#' and sets time to be 0 at the begining of each phase, see \code{\link{pwtime}}.
#' The other option is \code{'align'} which aligns time at 0 between the first
#' and second phase.
#'
#' Should the time variable be realigned at the phase?
#' If \code{TRUE} (the default), the time for first observation in the second phase
#' becomes 0 and times prior to the second phase are negative. For example, if the
#' time variable is \code{c(0,1,2,3,4,5)} and the phase variable is
#' \code{c(0,0,0,1,1,1)}, phase alignment yields a time variable of
#' \code{c{-3,-2,-1,0,1,2}}. This is useful when the timing of the transition
#' between the first and second phases varies by individual, especially for
#' graphing. This approach does not generalize to three or more phases, and
#' alignment only happens at the first phase transition. If there are three or
#' more phases, the later phase will not be aligned.
#'
#' @param fpc Numeric. The default is 0. If the value of fpc is greater than
#' the number of individuals in the data set, fpc is taken as the size of the
#' finite population from which the data were sampled, and a finite population
#' correction is made for the fixed effects standard errors, t-tests, and
#' p-values.
#'
#' @param debugforeach Logical. The default is \code{debugforeach=FALSE}. Not
#' implemented.
#'
#' @param cores Integer. The defaults is \code{parallel::detectCores()-1}, or
#' one fewer cores than what is detected on the machine.
#'
#' @param ... Not currently used.
#'
#' @examples
#'
#' \dontrun{
#'
#' # full sample model
#' t0 <- PersonAlytic(output  = 'Test0'     ,
#'                    data    = OvaryICT    ,
#'                    ids     = "Mare"      ,
#'                    dvs     = "follicles" ,
#'                    phase   = "Phase"     ,
#'                    time    = "Time"      ,
#'                    package = "nlme"      )
#'
#' # individual models
#' t1 <- PersonAlytic(output          = 'Test1'     ,
#'                    data            = OvaryICT    ,
#'                    ids             = "Mare"      ,
#'                    dvs             = "follicles" ,
#'                    phase           = "Phase"     ,
#'                    time            = "Time"      ,
#'                    package         = "arma"      ,
#'                    individual_mods = TRUE        )
#'
#' summary(t0)
#' summary(t1)
#'
#'
#' # full sample model with a finite population correction (FPC)
#' t0fpc <- PersonAlytic(output  = 'Test0FPC'   ,
#'                       data    = OvaryICT     ,
#'                       ids     = "Mare"       ,
#'                       dvs     = "follicles"  ,
#'                       phase   = "Phase"      ,
#'                       time    = "Time"       ,
#'                       package = "nlme"       ,
#'                       fpc     = 100          )
#'
#' # gamlss with two distributions - features not implemented
#' #OvaryICT$follicles01 <- to01(OvaryICT$follicles)
#' #t2 <- PersonAlytic(data=OvaryICT,
#' #                 ids="Mare",
#' #                 dvs=list("follicles", "follicles01"),
#' #                 phase="Phase",
#' #                 time="Time",
#' #                 family=c(NO(), BEINF()),
#' #                 package='gamlss')
#'
#'
#' # individual models with target variables
#' t3 <- PersonAlytic(output          = 'TargetIVStest'       ,
#'                    data            = OvaryICT              ,
#'                    ids             = "Mare"                ,
#'                    dvs             = "follicles"           ,
#'                    phase           = "Phase"               ,
#'                    time            = "Time"                ,
#'                    package         = "arma"                ,
#'                    individual_mods = TRUE                  ,
#'                    target_ivs      = names(OvaryICT)[6:11] )
#'
#' # multiple DVs with no target_ivs and group models
#' t4 <- PersonAlytic(output          = 'MultiDVnoIDnoIV'          ,
#'                    data            = OvaryICT                   ,
#'                    ids             = "Mare"                     ,
#'                    dvs             = names(OvaryICT)[c(3,9:11)] ,
#'                    phase           = "Phase"                    ,
#'                    time            = "Time"                     ,
#'                    package         = "lme"                      ,
#'                    individual_mods = FALSE                      ,
#'                    target_ivs      = NULL                       ,
#'                    autoDetect      = list()                     )
#'
#'
#' # repeat t4 with finite population correction
#' t5 <- PersonAlytic(output          = 'MultiDVnoIDnoIVFPC'       ,
#'                    data            = OvaryICT                   ,
#'                    ids             = "Mare"                     ,
#'                    dvs             = names(OvaryICT)[c(3,9:11)] ,
#'                    phase           = "Phase"                    ,
#'                    time            = "Time"                     ,
#'                    package         = "lme"                      ,
#'                    individual_mods = FALSE                      ,
#'                    target_ivs      = NULL                       ,
#'                    autoDetect      = list()                     ,
#'                    fpc             = 200                        )
#'
#' # repeat t5 with piecewise model, first making time an integer (required for piecewise)
#' OvaryICT$TimeP <- round(30*OvaryICT$Time)
#' t6 <- PersonAlytic(output          = 'PiecewiseExample'         ,
#'                    data            = OvaryICT                   ,
#'                    ids             = "Mare"                     ,
#'                    dvs             = names(OvaryICT)[c(3,9:11)] ,
#'                    phase           = "Phase"                    ,
#'                    time            = "TimeP"                    ,
#'                    package         = "lme"                      ,
#'                    individual_mods = FALSE                      ,
#'                    target_ivs      = NULL                       ,
#'                    autoDetect      = list()                     ,
#'                    fpc             = 200                        ,
#'                    alignPhase      = "piecewise"                )
#'
#' # clean up
#' file.remove( 'PiecewiseExample_PersonAlytic.csv' )
#' file.remove( 'REMLlme.txt' )
#' file.remove( 'MultiDVnoIDnoIVFPC_PersonAlytic.csv' )
#' file.remove( 'MultiDVnoIDnoIV_PersonAlytic.csv' )
#' file.remove( 'TargetIVStest_PersonAlytic.csv' )
#' file.remove( 'NotREMLarma.txt' )
#' file.remove( 'Test0FPC.txt' )
#' file.remove( 'Test1_PersonAlytic.csv' )
#' file.remove( 'Test0.txt' )
#'
#' }


# \dontrun{
# # if you wish to delete the automatically created csv file, run
# #NOT IMPLEMENTED YET
# }

# non-user options allowed in ...
# packageTest - override the package, for research purposes
# userFormula - override the formulae, named list that can include any of
#   fixed, random, formula


PersonAlytic <- function(output          = NULL                                  ,
                         data                                                    ,
                         ids                                                     ,
                         dvs                                                     ,
                         time                                                    ,
                         phase           = NULL                                  ,
                         ivs             = NULL                                  ,
                         target_ivs      = NULL                                  ,
                         interactions    = NULL                                  ,
                         time_power      = 1                                     ,
                         correlation     = NULL                                  ,
						             family          = gamlss.dist::NO()                     ,
                         subgroup        = NULL                                  ,
                         standardize     = list(dv=FALSE, iv=FALSE, byids=FALSE) ,
                         method          = 'REML'                                ,
                         package         = 'nlme'                                ,
                         individual_mods = FALSE                                 ,
                         PalyticObj      = NULL                                  ,
						             autoDetect      = list(AR=list(P=3, Q=3)     ,
						                                 TO=list(polyMax=3)       ,
						                                 DIST=list())                        ,
                         whichIC         = c("BIC", "AIC")                       ,
                         charSub         = NULL                                  ,
                         sigma.formula   = ~1                                    ,
                         p.method        = "BY"                                  ,
                         alpha           = .05                                   ,
                         nbest           = NULL                                  ,
                         alignPhase      = 'none'                                ,
                         fpc             = 0                                     ,
						             debugforeach    = FALSE                                 ,
                         cores           = parallel::detectCores()-1             ,
                         ...)
{
  #message('alignPhase=',alignPhase)

  if(length(whichIC)>1) whichIC <- whichIC[1]
  if(is.null(correlation)) correlation <- "NULL"

  # output labeling ####
  #TODO(Stephen) consider adding date/time/version info to output
  if(is.null(output))  fileLabel <- "PersonAlytics_Output"
  if(!is.null(output)) fileLabel <- paste(output, 'PersonAlytic', sep='_')

  # check the subgroup input ####
  if(!is.logical(subgroup) & !is.null(subgroup))
  {
    stop('`subgroup` must be a logical vector with TRUE/FALSE values.')
  }
  if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))

  # override user defaults if there is a data/package/individual_mods mismatch
  n <- length(unique(data[[ids]][subgroup]))
  if( n == 1 & package != "arma" )
  {
    message('\nThere is only one participant,', '
            automatically switching to `package="arma"`.\n')
    package <- 'arma'
  }
  if( n > 1 & package == 'arma'  &
      individual_mods != TRUE)
  {
    message('There is more than one participant, using `package="nlme"`.',
            '\nIf individual models are needed, switch to `individual_mods=TRUE`.')
    package <- 'nlme'
  }
  if( individual_mods == TRUE & package != 'arma' )
  {
    message('\nWhen `individual_mods=TRUE` you must use `package="arma"`.',
            '\nAutomatically switching to `package="arma".\n')
    package <- 'arma'
  }

  # for internal testing use only, allows us to compare n=1 lme models to
  # arma models
  #args <- list(packageTest='nlme')
  args <- list(...)
  if(package=='lme') package <- 'nlme'
  if("packageTest" %in% names(args))
  {
    if(args$packageTest=='lme') args$packageTest <- 'nlme'
    if(args$packageTest %in% c('nlme', 'gamlls', 'arma')) package <- args$packageTest
    message('\npackage was overridden by packageTest to be `', package, '`\n')
  }

  # for internal testing use only, allows us to override formulae, e.g., to
  # test random intercepts only or random slopes only models
  userFormula = list(
    fixed=NULL,
    random=NULL,
    formula=NULL)
  if('userFormula'  %in% names(args))
  {
    userFormula <- args$userFormula
  }

  # override individual_mods if only 1 id, 1 dv, <= 1 target_iv
  luid <- length(unique(data[[ids]]))
  if(individual_mods==TRUE & luid==1 &
     length(dvs)==1 & length(target_ivs)<=1)
  {
    individual_mods <- FALSE
  }

  # finite population correction (fpc) settings
  if( fpc > n ) fpcCheck(fpc, n)
  if( is.numeric(fpc) )
  {
    popsize2 <- fpc
    .fpc <- TRUE
  }
  if( !is.numeric(fpc) | fpc <= n )
  {
    popsize2 <- NULL
    .fpc <- FALSE
  }
  fpc <- .fpc; rm(.fpc)

  # augment time
  time <- list(raw      = time        ,
               power    = time_power  ,
               analysis = time        )


  #cat(package, file=paste(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"), 'txt', sep='.'))

  # call 'methods'
  if(individual_mods==FALSE & length(dvs)==1 & length(target_ivs)<=1)
  {
    pout <- pa1()
  }
  if(individual_mods==TRUE | length(dvs)>1 | length(target_ivs)>1)
  {
    pout <- paHTP()
  }

  while( sink.number() > 0 ) sink()

  return(pout)
}

#TODO(Stephen): individual_mods = FALSE is yielding different results than
#if running via htp (which can be forced), see the example in
#R:\PaCCT\02 Clients\curelator\Analyses\20181004\20181004 Curelator Analyses
#PersonAlytics_v0.1.2.2.r
#' pa1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
pa1 <- function(e=parent.frame())
{
  # set output filename ####
  if(is.null(e$output))
  {
    e$output <- paste(e$fileLabel, 'txt', sep='.')
  }
  if(!is.null(e$output))
  {
    e$output <- paste(e$output, '.txt', sep='')
  }

  if(is.null(e$subgroup)) e$subgroup <- rep(TRUE, nrow(e$data))

  # concatenate ivs
  ivs <- c(e$ivs, e$target_ivs)
  ivs <- ivs[!duplicated(ivs)]

  t1 <- Palytic$new(data=e$data                 ,
                    ids=e$ids                   ,
                    dv=e$dvs[[1]]               ,
                    time=e$time$raw             ,
                    phase=e$phase               ,
                    ivs=ivs                     ,
                    interactions=e$interactions ,
                    time_power=e$time_power     ,
                    alignPhase=e$alignPhase     ,
                    correlation=e$correlation   ,
                    family=e$family             ,
                    method=e$method             ,
                    standardize=e$standardize   ,
                    autoDetect=e$autoDetect     )

  # allow for formula override so that we can test intercept only and
  # slope only models
  if( any(unlist(lapply(e$userFormula, function(x) !is.null(x)))) )
  {
    isnnform <- function(x)
    {
      if( !is.null(x) )
      {
        return( is.formula(x) )
      }
      else return(FALSE)
    }
    if( isnnform(e$userFormula$fixed) ) t1$fixed <- e$userFormula$fixed
    if( isnnform(e$userFormula$random) ) t1$random <- e$userFormula$random
    if( isnnform(e$userFormula$formula) ) t1$formula <- e$userFormula$formula
  }

  # autodetection
  if(e$package=="nlme" | e$package=="arma")
  {
    temp <- t1$autoDetect
    temp$DIST <- NULL
    t1$autoDetect <- temp; rm(temp)
  }
  dims <- list(ID="All Cases")
  if(var(t1$datac[[t1$ids]][e$subgroup], na.rm=TRUE)==0)
  {
    dims <-list(ID=sort(unique(e$data[[e$ids]])))
  }
  t1$detect(model=NULL, parallel="snow", plot=FALSE, userFormula=e$userFormula,
            dims=dims)

  # t1$correlation
  # t1$formula

  # fit the models
  if(e$package=="gamlss") Grp.out <- t1$gamlss(e$subgroup, sigma.formula=e$sigma.formula)
  if(e$package=="nlme")   Grp.out <- t1$lme(e$subgroup, fpc=e$fpc,
                                            popsize2=e$popsize2)
  if(e$package=="arma")   Grp.out <- t1$arma(e$subgroup )

  sink(file=e$output)
  print( Grp.out$tTable )
  sink()

  return(Grp.out)

}

#' paHTP
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
paHTP <- function(e=parent.frame())
{
  # set output filename ####
  if(is.null(e$output))
  {
    e$output <- paste(e$fileLabel, 'csv', sep='.')
  }
  if(!is.null(e$output))
  {
    e$output <- paste(e$fileLabel, '.csv', sep='')
  }

  # check that dvs, target_ivs are lists, if not, force
  if( ! "list" %in% class(e$dvs) ) e$dvs <- as.list(e$dvs)
  if( ! "list" %in% class(e$ivs) ) e$ivs <- as.list(e$ivs)
  #if( ! "list" %in% class(e$target_ivs) ) e$target_ivs <- as.list(e$target_ivs)

  # check that inputs conform. This is also done when creating a Palytic
  # object, but we do it early on here to avoid problems after loops start.
  # Set standardize to FALSE here, standardization will be done later if
  # requested by the user (this avoids standardizing standardized variables,
  # which may differ under different subsets, e.g., individual_models = TRUE).
  e$data <- clean(
    data         = e$data                                ,
    ids          = e$ids                                 ,
    dv           = NULL                                  ,
    time         = e$time                                ,
    phase        = e$phase                               ,
    ivs          = e$ivs                                 ,
    fixed        = NULL                                  ,
    random       = NULL                                  ,
    formula      = NULL                                  ,
    correlation  = e$correlation                         ,
    family       = e$family                              ,
    dvs          = e$dvs                                 ,
    target_ivs   = e$target_ivs                          ,
    standardize  = list(dvs=FALSE,ivs=FALSE,byids=FALSE) ,
    alignPhase   = "none"                                ,
    debugforeach = e$debugforeach                        )

  # standardize target_ivs
  if(e$standardize$iv)
  {
    message('\nPersonAlytics is standardizing the variables in `targe_ivs`.\n')
    for(i in seq_along(e$target_ivs))
    {
      tiv <- e$target_ivs[[i]]
      if(is.numeric(e$data[tiv]))
      {
        e$data[[tiv]] <- scale(e$data[[tiv]])
      }
    }
  }

  # subgroup the data and delete the parameter, after this point, it is only
  # used to subgroup to unique ids
  if( is.null(e$subgroup)) e$subgroup <- rep(TRUE, nrow(e$data))
  if(!is.null(e$data))     e$data <- e$data[e$subgroup,]

  # check whether any variables in ivs are in target_ivs -
  # in the future, split them out automatically
  if(any(e$ivs %in% e$target_ivs) | any(e$target_ivs %in% e$ivs))
  {
    stop('target_ivs and ivs cannot share any variables.')
  }

  # autodetection
  if(e$package=="nlme" | e$package=="arma")
  {
    e$autoDetect$DIST <- NULL
  }

  # unique ids
  uids <- sort(unique(e$data[[e$ids]]))

  # dimensions for loops
  ID <- indID <- uids
  IV <- seq_along(e$target_ivs)
  if(is.null(e$target_ivs) | length(e$target_ivs)==0) IV <- 1
  DV <- seq_along(e$dvs)
  dims <- list(ID=ID, IV=IV, DV=DV, indID=indID)

  if( e$debugforeach )
  {
    cat('\n\n\n')
    print(head(e$data)   ); cat('\n')
    print(dims           ); cat('\n')
    print(e$dvs          ); cat('\n')
    print(e$phase        ); cat('\n')
    print(e$ids          ); cat('\n')
    print(uids           ); cat('\n')
    print(e$time         ); cat('\n')
    print(e$ivs          ); cat('\n')
    print(e$target_ivs   ); cat('\n')
    print(e$interactions ); cat('\n')
    print(e$time_power   ); cat('\n')
    print(e$alignPhase   ); cat('\n')
    print(e$correlation  ); cat('\n')
    print(e$family       ); cat('\n')
    print(e$standardize  ); cat('\n')
    print(e$fpc          ); cat('\n')
    print(e$popsize2     ); cat('\n')
    print(e$package      ); cat('\n')
    print(e$autoDetect   ); cat('\n')
    print(e$PQ           ); cat('\n')
    print(e$whichIC      ); cat('\n')
    print(e$polyMax     ); cat('\n')
    print(e$sigma.formula); cat('\n')
    print(e$debugforeach ); cat('\n')
    print(e$cores        ); cat('\n')
    cat('\n\n\n')
  }

  #
  if( e$individual_mods )
  {
    DVout <- htp(data          = e$data          ,
                 dims          = dims            ,
                 dvs           = e$dvs           ,
                 phase         = e$phase         ,
                 ids           = e$ids           ,
                 uids          = uids            ,
                 time          = e$time          ,
                 ivs           = e$ivs           ,
                 target_ivs    = e$target_ivs    ,
                 interactions  = e$interactions  ,
                 time_power    = e$time_power    ,
                 alignPhase    = e$alignPhase    ,
                 correlation   = e$correlation   ,
                 family        = e$family        ,
                 standardize   = e$standardize   ,
                 fpc           = e$fpc           ,
                 popsize2      = e$popsize2      ,
                 package       = e$package       ,
                 autoDetect    = e$autoDetect    ,
                 whichIC       = e$whichIC       ,
                 sigma.formula = e$sigma.formula ,
                 debugforeach  = e$debugforeach  ,
                 cores         = e$cores         )
  }
  if( !e$individual_mods )
  {
    grp.dims <- dims
    grp.dims$ID <- "All Cases"

    DVout <- htp(data          = e$data          ,
                 dims          = grp.dims        ,
                 dvs           = e$dvs           ,
                 phase         = e$phase         ,
                 ids           = e$ids           ,
                 uids          = uids            ,
                 time          = e$time          ,
                 ivs           = e$ivs           ,
                 target_ivs    = e$target_ivs    ,
                 interactions  = e$interactions  ,
                 time_power    = e$time_power    ,
                 alignPhase    = e$alignPhase    ,
                 correlation   = e$correlation   ,
                 family        = e$family        ,
                 standardize   = e$standardize   ,
                 fpc           = e$fpc           ,
                 popsize2      = e$popsize2      ,
                 package       = e$package       ,
                 autoDetect    = e$autoDetect    ,
                 whichIC       = e$whichIC       ,
                 sigma.formula = e$sigma.formula ,
                 debugforeach  = e$debugforeach  ,
                 cores         = e$cores         )

  }

  # clean up variable names
  if(!is.null(e$charSub))
  {
    DVout$target_iv <- as.character(DVout$target_iv)
    for(i in seq_along(e$charSub))
    {
      DVout$target_iv <- gsub(e$charSub[[i]][1], e$charSub[[i]][2], DVout$target_iv)
    }
  }

  # fix columns that are lists and reorder the columns
  DVout <- do.call(data.frame, lapply(DVout, nnull))

  # attempt adjusted p-values
  if(!is.null(e$p.method) & length(e$target_ivs) > 1)
  {
    DVpsuite <- try( psuite(DVout, e$ids,
                            output = e$output,
                     rawdata=e$data,
                     method=e$p.method,
                     nbest=e$nbest,
                     alpha=e$alpha), silent = TRUE )

    if( ncol(DVpsuite) == (ncol(DVout)+length(e$p.method)+1) )
    {
      DVout <- DVpsuite
    }
  }

  # save the output
  write.csv(DVout, file=e$output, row.names=FALSE)

  return(DVout)


}


#' nnull - fix columns that are actually lists
#' @keywords internal
nnull <- function(x)
{
  if(is.list(x))
  {
    if( any(as.character(x)=="NULL") )
    {
      return( unlist(as.character(x)) )
    }
    if(!any(as.character(x)=="NULL") )
    {
      return( unlist(x) )
    }
  }
  if(!is.list(x))
  {
    return( unlist(x) )
  }
}
