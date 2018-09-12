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
#' output. If left \code{NULL}, the date and time will prefix
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
#' If \code{detectTO=TRUE}, this is the largest value tested for \code{detectTO}.
#' The default value is \code{3}, testing up to a cubic growth model if
#' \code{detectTo=TRUE}. If a linear growth model is desired, set
#' \code{detectTO=FALSE} and \code{time_power=1}.
#'
#' @param correlation Character. The default value is \code{NULL}.
#'
#' See \code{\link{corStruct}} in \code{\link{nlme}}.
#' Must be passed as a character, e.g. \code{"corARMA(p=1)"}.
#' If \code{detectAR=TRUE}, \code{correlation} will be ignored, see \code{PQ}
#' and \code{detectAR}. The default value is \code{NULL}, assuming no residual
#' autocorrelation. If \code{detectAR=TRUE}, \code{correlation}
#' will be ignored.
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
#' @param subgroup Logical vector. The default is \code{subgroup==NULL}.
#'
#' A vector where \code{length(subgroup)==nrow(data)} indicating
#' which subset of the data should be used for analysis. For example, if a model
#' should only be fit to females, \code{subgroup=gender=='female'} might be used.
#'
#' @param standardize Logical. The default is \code{FALSE}.
#'
#' Should the dependent and independent variables
#' be standardized (i.e., rescaled to have 0 mean and unit variance; see
#' \code{\link{scale}})? Does not apply to factor variables. The default is \code{TRUE}
#' which makes parameter estimate magnitudes comparable across individuals, outcomes in
#' \code{dvs}, and covariates in \code{target_ivs}. For dependent variables in
#' \code{dvs}, standardization is only applied for normal outcomes, see \code{family}.
#'
#' @param package Character. The default is \code{"nlme"}.
#'
#' Which package should be used to fit the models?
#' Options are \code{\link{nlme}} and \code{\link{gamlss}}.
#' It is passed as character strings, e.g., \code{"gamlss"} or \code{"nlme"}.
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
#' @param detectAR Logical. The default is \code{detectAR=TRUE}.
#'
#' Should the residual autocorrelation structure be automatically selected from
#' among \code{ARMA(p,q)} models? See \code{correlation}. Since these models are
#' not generally nested, model selection is done using information information
#' criterion (see \code{whichIC}).
#'
#' @param PQ Numeric vector of length 2. The default is \code{PQ=c(3,3)}.
#'
#' If \code{detectAR=TRUE}, automatic selection of the residual covariance
#' structure is invoked initializing a search among
#' \code{p=1,...,P} and \code{p=1,...,Q}, where \code{P} and \code{Q} are taken
#' from \code{PQ}, i.e., \code{PQ=c(P,Q)}. The values of \code{p} and \code{p}
#' are passed to \code{\link{corARMA}} ( e.g., \code{corARMA(p=p,q=q)}) for
#' testing (see \code{detectAR}).
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
#' @param detectTO Logical. The default is \code{detectTO=TRUE}.
#'
#' Should the \code{time_power} value
#' be automatically selected? Values from 1 to \code{maxOrder} will be tested.
#' For example, if \code{maxOrder=3} (implying a cubic growth model), the models
#' compared include \code{time}, \code{time + I(time^2)}, and
#' \code{time + I(time^2)+I(time^3)}. Since these models are nested, the best
#' fitting model is selected using likelihood ratio tests with mixed effects
#' models fit using maximum likelihood estimators in \code{\link{lme}}.
#'
#' This is done separately for each individual in \code{ids} if
#' \code{individual_mods=TRUE}.
#'
#' @param maxOrder Numeric. The default is \code{maxOrder=3}.
#'
#' See \code{detectTO}.
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
#' updated by \code{time_power} or \code{detectTO}. If model fitting using this
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
#' @param alignPhase Logical. The default is \code{alignPhase=TRUE}.
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
#' @param debugforeach Logical. The default is \code{debugforeach=FALSE}. Not
#' implemented.
#'
#' @param ... Not currently used.
#'
#' @examples
#' # full sample model
#' t0 <- PersonAlytic(output='Test0',
#'                  data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  package="nlme")
#'
#' # individual models
#' t1 <- PersonAlytic(output='Test1',
#'                  data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  package="nlme",
#'                  individual_mods=TRUE)
#'
#' summary(t0)
#' summary(t1)
#'
#' # delete the output if this was run in the development directory
#' if(getwd()=="R:/PaCCT/Repository/PersonAlytics")
#' {
#'   file.remove( dir(getwd(), glob2rx("*.txt")) )
#'   file.remove( dir(getwd(), glob2rx("*.csv")) )
#' }
#'
#' # gamlss with two distributions - features not implemented
#' #OvaryICT$follicles01 <- to01(OvaryICT$follicles)
#' #t3 <- PersonAlytic(data=OvaryICT,
#' #                 ids="Mare",
#' #                 dvs=list("follicles", "follicles01"),
#' #                 phase="Phase",
#' #                 time="Time",
#' #                 family=c(NO(), BEINF()),
#' #                 package='gamlss')
#'

# \dontrun{
# # if you wish to delete the automatically created csv file, run
# #NOT IMPLEMENTED YET
# }

PersonAlytic <- function(output=NULL              ,
                         data                     ,
                         ids                      ,
                         dvs                      ,
                         time                     ,
                         phase=NULL               ,
                         ivs=NULL                 ,
                         target_ivs=NULL          ,
                         interactions=NULL        ,
                         time_power=1             ,
                         correlation=NULL         ,
                         family=gamlss.dist::NO() ,
                         subgroup=NULL            ,
                         standardize=TRUE         ,
                         method='REML'            ,
                         package='nlme'           ,
                         individual_mods=FALSE    ,
                         PalyticObj=NULL          ,
                         detectAR=TRUE            ,
                         PQ=c(3,3)                ,
                         whichIC=c("BIC", "AIC")  ,
                         detectTO=TRUE            ,
                         maxOrder=3               ,
                         charSub=NULL             ,
                         sigma.formula=~1         ,
                         p.method = "BY"          ,
                         alpha = .05              ,
                         alignPhase = FALSE       ,
                         debugforeach = FALSE     )
{
  if(length(whichIC)>1) whichIC <- whichIC[1]
  if(is.null(correlation)) correlation <- "NULL"
  pav <- paste("-PAv", packageVersion("PersonAlytics"), "-", sep='')

  if(individual_mods==FALSE & length(dvs)==1 & length(target_ivs)<=1)
  {
    pout <- pa1()
  }
  if(individual_mods==TRUE | length(dvs)>1 | length(target_ivs)>1)
  {
    pout <- paHTP()
  }

  return(pout)

  while( sink.number() > 0 ) sink()
}

#' pa1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
pa1 <- function(e=parent.frame())
{
  # if no data are given, use a test data set
  if(is.null(e$data))
  {
    output          <- NULL
    data            <- OvaryICT
    ids             <- "Mare"
    dvs             <- "follicles"
    time            <- "Time"
    phase           <- "Phase"
    ivs             <- NULL
    target_ivs      <- NULL
    interactions    <- NULL
    time_power      <- 3
    correlation     <- NULL
    family          <- gamlss.dist::NO()
    subgroup        <- NULL
    standardize     <- FALSE
    package         <- 'nlme'
    method          <- 'REML'
    individual_mods <- FALSE
    PalyticObj      <- NULL
    detectAR        <- FALSE
    PQ              <- c(3,3)
    whichIC         <- c("BIC", "AIC")
    detectTO        <- FALSE
    charSub         <- NULL
    sigma.formula   <- ~1
    p.method        <- "BY"
    alpha           <- .05
    alignPhase      <- FALSE
    debugforeach    <- FALSE
    maxOrder        <- 3
    e               <- parent.frame()
  }

  if(is.null(e$output))
  {
    e$output <- gsub(":", ".", paste(Sys.time(), e$pav,
                               'PersonAlyticHTP_Output.txt'))
  }
  if(!is.null(e$output))
  {
    e$output <- paste(e$output, e$pav, 'txt', sep='.')
  }

  if(is.null(e$subgroup)) e$subgroup <- rep(TRUE, nrow(e$data))

  t1 <- Palytic$new(data=e$data[e$subgroup,]    ,
                    ids=e$ids                   ,
                    dv=e$dvs                    ,
                    time=e$time                 ,
                    phase=e$phase               ,
                    ivs=e$ivs                   ,
                    interactions=e$interactions ,
                    time_power=e$time_power     ,
                    correlation=e$correlation   ,
                    family=e$family             ,
                    method=e$method             ,
                    standardize=e$standardize   ,
                    alignPhase=e$alignPhase     )

  # check time order first, that way the time order carries over to
  # the AR, which should be residuals on the fullest model

  if(e$detectTO) t1$GroupTime_Power(e$maxOrder, e$whichIC[1])
  # t1$time_power
  # t1$formula

  if(e$detectAR) t1$GroupAR_order(P = e$PQ[1]  ,
                                  Q = e$PQ[2]  ,
                                  whichIC= e$whichIC[1]  )
  # t1$correlation
  # t1$formula

  # fit the models
  if(e$package=="gamlss") Grp.out <- t1$gamlss()
  if(e$package=="nlme")   Grp.out <- t1$lme()

  sink(file=e$output)
  print( Grp.out )
  sink()

  return(Grp.out)

}

#' paHTP
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
paHTP <- function(e=parent.frame())
{

  if(is.null(e$output))
  {
    e$output <- gsub(":", ".", paste(Sys.time(), e$pav, 'PersonAlyticHTP_Output.csv'))
  }
  if(!is.null(e$output))
  {
    e$output <- paste(e$output, e$pav, 'csv', sep='.')
  }

  # check that dvs, target_ivs are lists, if not, force
  if( ! "list" %in% class(e$dvs) ) e$dvs <- as.list(e$dvs)
  if( ! "list" %in% class(e$ivs) ) e$ivs <- as.list(e$ivs)

  # check that inputs conform. this is also done when creating a Palytic
  # object, but we do it early on here to avoid problems after loops start.
  e$data <- PersonAlytics:::clean(e$data, e$ids, dv=NULL, e$time, e$phase, e$ivs,
                                  fixed=NULL, random=NULL, formula=NULL,
                                  e$correlation, e$family,
                                  e$dvs, e$target_ivs, e$standardize,
                                  e$alignPhase)

  # subgroup the data and delete the parameter, after this point, it is only
  # used to subgroup to unique ids
  if( is.null(e$subgroup)) e$subgroup <- rep(TRUE, nrow(e$data))
  if(!is.null(e$data)) e$data <- e$data[e$subgroup,]

  # check whether any variables in ivs are in target_ivs -
  # in the future, split them out automatically
  if(any(e$ivs %in% e$target_ivs) | any(e$target_ivs %in% e$ivs))
  {
    stop('target_ivs and ivs cannot share any variables.')
  }

  ## if no data are given, use a test data set
  if(is.null(e$data))
  {
    data   <- OvaryICT
    dvs    <- "follicles"
    phase  <- "Phase"
    ids    <- "Mare"
    time   <- "TimeSin"
  }

  # unique ids
  uids <- sort(as.numeric(unique(e$data[[e$ids]])))

  # dimensions for loops
  ID <- uids
  IV <- 1:length(e$target_ivs); if(is.null(e$target_ivs)) IV <- 1
  DV <- 1:length(e$dvs)
  dims <- list(ID=ID, IV=IV, DV=DV)

  #
  if( e$individual_mods )
  {
    DVout <- htp.foreach(data  = e$data                  ,
                         dims  = dims                    ,
                         dvs   = e$dvs                   ,
                         phase = e$phase                 ,
                         ids   = e$ids                   ,
                         uids  = uids                    ,
                         time  = e$time                  ,
                         ivs   = e$ivs                   ,
                         target_ivs    = e$target_ivs    ,
                         interactions  = e$interactions  ,
                         time_power    = e$time_power    ,
                         correlation   = e$correlation   ,
                         family        = e$family        ,
                         standardize   = e$standardize   ,
                         package       = e$package       ,
                         detectAR      = e$detectAR      ,
                         PQ            = e$PQ            ,
                         whichIC       = e$whichIC       ,
                         detectTO      = e$detectTO      ,
                         maxOrder      = e$maxOrder      ,
                         sigma.formula = e$sigma.formula ,
                         debugforeach  = e$debugforeach  )
  }
  if( !e$individual_mods )
  {
    grp.dims <- dims
    grp.dims$ID <- "All Cases"

    DVout <- htp.foreach(data  = e$data                  ,
                         dims  = grp.dims                ,
                         dvs   = e$dvs                   ,
                         phase = e$phase                 ,
                         ids   = e$ids                   ,
                         uids  = uids                    ,
                         time  = e$time                  ,
                         ivs   = e$ivs                   ,
                         target_ivs    = e$target_ivs    ,
                         interactions  = e$interactions  ,
                         time_power    = e$time_power    ,
                         correlation   = e$correlation   ,
                         family        = e$family        ,
                         standardize   = e$standardize   ,
                         package       = e$package       ,
                         detectAR      = e$detectAR      ,
                         PQ            = e$PQ            ,
                         whichIC       = e$whichIC       ,
                         detectTO      = e$detectTO      ,
                         maxOrder      = e$maxOrder      ,
                         sigma.formula = e$sigma.formula ,
                         debugforeach  = e$debugforeach  )

  }

  # clean up variable names
  if(!is.null(e$charSub))
  {
    DVout$target_iv <- as.character(DVout$target_iv)
    for(i in 1:length(e$charSub))
    {
      DVout$target_iv <- gsub(e$charSub[[i]][1], e$charSub[[i]][2], DVout$target_iv)
    }
  }

  # some columns are actually lists, fix that
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
  DVout <- do.call(data.frame, lapply(DVout, nnull))
  write.csv(DVout, file=paste(e$output, '.csv', sep=''), row.names=FALSE)

  return(DVout)

  if(!is.null(e$p.method) & length(e$target_ivs) > 1)
  {
    DVout <- try( psuite(DVout,
                    rawdata=e$data,
                    method=e$p.method,
                    alpha=e$alpha), silent = TRUE )
  }
}
