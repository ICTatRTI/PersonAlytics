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
#' effects models using
#' \code{\link{nlme}} and generalized linear mixed effects models using
#' \code{\link{gamlss}}.
#' The basic mixed effects model is \eqn{dv=time+phase+phase*time}
#' with random intercepts and random slopes for time. The phase variable is optional.
#' Additional independent variables (or covariates) can be included.
#'
#' High throughput capabalities are invoked When the model
#' is to be fit to multiple indidivuals, multiple dependent variables, iterated
#' over multiple target independent variables, or any combination of these three.
#'
#' The residual correlation structure is automatically selected up to \code{ARMA(p,q)}
#' where \code{p} and \code{q} are selected by the user. Selection is done by comparing
#' fit indices, and either the \code{BIC} or \code{AIC} can be used.
#'
#' The functional form of the relationship between time and the dependent variable
#' is automatically selected up to \code{time^time_power} where \code{time_power} is
#' selected by the user. Selection is done using maximum likelihood based likelihood
#' ratio tests (LRT).
#'
#' Type I error corrections are implemented using the options available in
#' \code{\link{p.adjust}}.
#'
#' @param file The file name (or full path with `/` instead of `\`) where output should
#' be saved. If left \code{NULL}, the date and time will prefix
#' `PersonAlyticHTP_Output.csv`.
#'
#' @param data A \code{\link{data.frame}} with the variables assinged to \code{ids},
#' \code{dv}, \code{time}, and optionally, \code{phase} and \code{ivs}. If left as
#' \code{NULL}, an example will run using the Ovary data from the \code{\link{nlme}}
#' package, see \code{\link{OvaryITC}}.
#'
#' @param ids Character. Name of the ID variable. The ID variable must be numeric.
#'
#' @param dvs A list of one or more character dependent variable names in \code{data}.
#' The linear mixed effects model
#' \code{dvs[d] ~ phase + time + phase*time + target_ivs[c] + ivs}
#' with random effects \code{~ time | ids[i]} will be fit to the data.
#' The iterators \code{[d]}, \code{[c]}, and \code{[i]}
#' indicate that the model will be fit for each combination of
#' dependent variable in \code{dvs},
#' independent variable in \code{target_ivs}, and
#' each unique ID in \code{ids} (which can be overridden using \code{ind.mods=FALSE})
#' controlling for one or more independent variables in \code{ivs}.
#'
#' @param time Character. Name of the time variable.
#'
#' @param phase Charcter. Name of the phase or treatment variable. For model fitting,
#' \code{phase} is treated the same as any variable names in \code{ivs} or
#' \code{target_ivs} but is used for visualizing treatment effects.
#'
#' @param ivs Character list of covariates, e.g., \code{list('iv2', 'iv2')}.
#' Note that the variables in \code{ivs} cannot also be in \code{target_ivs}.
#'
#' @param target_ivs Independent variables that are iterated over one at a time
#' (see \code{dv} for details).
#' Effects for these variables are labeled as 'target predictor' in the output.
#'
#' @param interactions List of vector pairs of variable names for interactions
#' to include in the model, e.g., \code{list(c('iv1','phase'), c('iv1','iv2'))}.
#'
#' @param time_power Numeric. Power of the time variable (e.g., \code{time^time_power}).
#' A quadratic or cubic growth model would be specified using
#' \code{time_power=2} or \code{time_power=3}, respectively.
#' If \code{detectTO=TRUE}, this is the largest value tested for \code{detectTO}.
#'
#' @param correlation See \code{\link{corStruct}} in \code{\link{nlme}}.
#' Must be passed as a character, e.g. \code{"corARMA(p=1)"}.
#' If \code{detectAR=TRUE}, \code{correlation} will be ignored, see \code{PQ}
#' and \code{detectAR}.
#'
#' @param family See \code{\link{gamlss.family}}. The default is normal. A list
#' of the same length as \code{length(dv)} can be supplied. For example if
#' \code{dv=list('normaly', 'binaryy', 'betay')}, then family could be
#' \code{family=list(NO(), BI(), BEINF())}. If \code{length(dv)>1} and
#' \code{length(family)==1}, the distribution will be applied to all outcomes.
#'
#' @param subgroup Logical vector where \code{length(subgroup)==nrow(data)} indicating
#' which subset of the data should be used for analysis. For example, if a model
#' should only be fit to females, \code{subgroup=gender=='female'}.
#'
#' @param standardize Logical. Should the dependent and independent variables
#' be standardized (i.e., rescaled to have 0 mean and unit variance; see
#' \code{\link{scale}})? Does not apply to factor variables. The default is \code{TRUE}
#' which makes parameter estimate magnitudes comparable across individuals, outcomes in
#' \code{dvs}, and covariates in \code{target_ivs}. For dependent variables in
#' \code{dvs}, standardization is only applied for normal outcomes, see \code{family}.
#'
#' @param package Which package should be used? Options are
#' \code{\link{nlme}} (the default) and \code{\link{gamlss}}
#' It is passed as character strings, e.g., \code{"gamlss"}. The \code{family}
#' parameter is ignored if \code{package="nlme"}.
#'
#' @param ind.mods Logical, defaults to \code{FALSE}. Should individual models be
#' fit for each ID?
#'
#' @param PalyticObj See \code{\link{Palytic}}. If \code{PalyticObj} is submitted
#' then only \code{dvs}, \code{target_ivs}, and \code{ind.mods} will be
#' used. This allows users access to additional options including generalized linear
#' mixed effects models via the \code{family} option, user specified \code{correlation}
#' structures (in non-\code{NULL} this will override the automated correlation structure
#' search), and user specified models via \code{formula}.
#'
#' @param detectAR Logical, defaults to \code{TRUE}. Should the residual
#' autocorrelation structure be automatically selected from among
#' \code{ARMA(p,q)} models? See \code{correlation}. Since these models are not
#' nested, model selection is done using information information criterion
#' (see \code{IC}).
#'
#' @param PQ Numeric vector of length 2, e.g., \code{PQ=c(3,3)}.
#' If \code{detectAR=TRUE}, automatic selection of the residual covariance
#' structure is invoked initializing a search among
#' \code{p=1,...,P} and \code{p=1,...,Q}, where \code{P} and \code{Q} are taken
#' from \code{PQ}, i.e., \code{PQ=c(P,Q)}. The values of \code{p} and \code{p}
#' are passed to \code{\link{nlme::corARMA}} ( e.g., \code{corARMA(p=p,q=q)}) for
#' testing (see \code{detectAR}). If \code{detectAR=TRUE}, \code{correlation}
#' will be ignored.
#'
#' @param IC Either the Akaike Information Criterion (\code{IC="AIC"}) or
#' the Bayesian Information Criterion (\code{IC="BIC"}, the default).
#'
#' If the \code{time} variable is equally spaced, this is
#' done using the function \code{\link{forecast}}. If the \code{time} variable
#' is notequally spaced, this is done using comparisons of
#' mixed effects models using \code{\link{lme}} fit using maximum likelihood
#' estimators.
#'
#' Residual autocorrelation structure is done separately for each case in
#' \code{ids} if \code{ind.mods=TRUE}.
#'
#' @param detectTO Logical, defaults to \code{TRUE}. Should the \code{time_power} value
#' be automatically selected? Values from 1 to \code{time_power} will be tested. For
#' example, if \code{time_power=3} (implying a cubic growth model), the models
#' compared include \code{time}, \code{time + I(time^2)}, and
#' \code{time + I(time^2)+I(time^3)}. Since these models are nested, the best
#' fitting model is selected using likelihood ratio tests with mixed effects
#' models fit using maximum likelihood estimators in \code{\link{lme}}.
#'
#' This is done separately for each case in \code{ids} if \code{ind.mods=TRUE}.
#'
#' @param charSub list of paired character strings for character substitution
#' in the output. If the names of the target predictors
#' in \code{target_ivs} had to be edited to make valid variable names, this
#' parameter allows users put the illegal characters back in. For example,
#' if the original variable name was "17.00_832.2375m/z", a letter would need to
#' prefix the variable name and the
#' "/" would need to be replaced with another character, e.g., "X17.00_832.2375m.z".
#' To get the row names of the output back to original varibale name, use
#' \code{charSub=list(c("X", ""), c("m.z", "m/z"))}. Note that inputs to charSub
#' must be in double quotes and are case sensitive. All duplicates will be substituted.
#' For example, if the variable name was "X1X23.x" and \code{charSub=list(c("X", ""))},
#' the resulting row label for this variable would be "123.x".
#'
#' @param sigma.formula A formula for the variance under \code{\link{gamlss}}.
#' Curently static: it will not change dynamically over iterations nor will it be
#' updated by \code{time_power} or \code{detectTO}. If model fitting using this
#' option fails, another attempt will be made after reseting it to its defaul,
#' i.e., \code{~1}.
#'
#' @param p.method See \code{\link{p.adjust.methods}}. When \code{ind.mods=TRUE},
#' \code{length(dvs)>1}, or \code{length(target_ivs)>1}, p-value adjustments
#' are made and reported in separate output saved to the working directory.
#'
#' @param alpha Numeric value in the (0,1) interval. The Type I error rate for
#' adjusting p-values.
#'
#' @param alignPhase Logical. Should the time variable be realigned at the phase?
#' If \code{TRUE} (the default), the time for first observation in the second phase
#' becomes 0 and times prior to the secord phase are negative. For example, if the
#' time variable is \code{c(0,1,2,3,4,5)} and the phase variable is
#' \code{c(0,0,0,1,1,1)}, phase alignment yields a time variable of
#' \code{c{-3,-2,-1,0,1,2}}. This is useful when the timing of the transition
#' between the first and second phases varies by individual, especially for
#' graphing. This approach does not generalize to three or more phases, and
#' alignment only happens at the first phase transition. If there are three or
#' more phases, the later phase will not be aligned.
#'
#' @examples
#' # full sample model
#' t0 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  package='nlme',
#'                  ind.mods=FALSE)
#'
#' # individual models (using defaults)
#' t1 <- PersonAlytic(data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="Time",
#'                  package='nlme')
#'
#' # gamlss with two distributions - features not implemented
#' #OvaryICT$follicles01 <- to01(OvaryICT$follicles)
#' #t1 <- PersonAlytic(data=OvaryICT,
#' #                 ids="Mare",
#' #                 dvs=list("follicles", "follicles01"),
#' #                 phase="Phase",
#' #                 time="Time",
#' #                 family=c(NO(), BEINF()),
#' #                 package='gamlss')
#'
#' summary(t0)
#' summary(t1)

# \dontrun{
# # if you wish to delete the automatically created csv file, run
# #NOT IMPLEMENTED YET
# }

PersonAlytic <- function(file=NULL                ,
                         data=NULL                ,
                         ids                      ,
                         dvs                      ,
                         time                     ,
                         phase=NULL               ,
                         ivs=NULL                 ,
                         target_ivs=NULL          ,
                         interactions=NULL        ,
                         time_power=3             ,
                         correlation=NULL         ,
                         family=gamlss.dist::NO() ,
                         subgroup=NULL            ,
                         standardize=TRUE         ,
                         package='nlme'           ,
                         ind.mods=TRUE            ,
                         PalyticObj=NULL          ,
                         detectAR=TRUE            ,
                         PQ=c(3,3)                ,
                         IC=c("BIC", "AIC")       ,
                         detectTO=TRUE            ,
                         charSub=NULL             ,
                         sigma.formula=~1         ,
                         p.method = "BY"          ,
                         alpha = .05              ,
                         alignPhase = TRUE        ,
                         ...)
{
  if(length(IC)>1) IC <- IC[1]
  maxOrder <- time_power
  time_power <- 1

  if(ind.mods==FALSE & length(dvs)==1 & length(target_ivs)<=1)
  {
    pa1(file           ,
        data           ,
        ids            ,
        dvs            ,
        time           ,
        phase          ,
        ivs            ,
        target_ivs     ,
        interactions   ,
        time_power     ,
        maxOrder       ,
        correlation    ,
        family         ,
        subgroup       ,
        standardize    ,
        package        ,
        ind.mods       ,
        PalyticObj     ,
        detectAR       ,
        PQ             ,
        IC             ,
        detectTO       ,
        charSub        ,
        sigma.formula  ,
        p.method       ,
        alpha          )
  }
  if(ind.mods==TRUE | length(dvs)>1 & length(target_ivs)>1)
  {
    paHTP(file           ,
          data           ,
          ids            ,
          dvs            ,
          time           ,
          phase          ,
          ivs            ,
          target_ivs     ,
          interactions   ,
          time_power     ,
          maxOrder       ,
          correlation    ,
          family         ,
          subgroup       ,
          standardize    ,
          package        ,
          ind.mods       ,
          PalyticObj     ,
          detectAR       ,
          PQ             ,
          IC             ,
          detectTO       ,
          charSub        ,
          sigma.formula  ,
          p.method       ,
          alpha          )
  }
}

#' paHTP
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
paHTP <- function(file           ,
                  data           ,
                  ids            ,
                  dvs            ,
                  time           ,
                  phase          ,
                  ivs            ,
                  target_ivs     ,
                  interactions   ,
                  time_power     ,
                  maxOrder       ,
                  correlation    ,
                  family         ,
                  subgroup       ,
                  standardize    ,
                  package        ,
                  ind.mods       ,
                  PalyticObj     ,
                  detectAR       ,
                  PQ             ,
                  IC             ,
                  detectTO       ,
                  charSub        ,
                  sigma.formula  ,
                  p.method       ,
                  alpha          ,
                  debugforeach = FALSE)
{
  #
  if(is.null(file))
  {
    file <- gsub(":", ".", paste(Sys.time(), 'PersonAlyticHTP_Output.csv'))
  }

  # check that dvs, target_ivs are lists, if not, force
  if( ! "list" %in% class(dvs) ) dvs <- as.list(dvs)
  if( ! "list" %in% class(ivs) ) ivs <- as.list(ivs)

  # check that inputs conform. this is also done when creating a Palytic
  # object, but we do it early on here to avoid problems after loops start.
  # This is why `clean()` has inputs that apply to PersonAlyticHTP but not to
  # PersonAlytic
  data <- clean(data, ids, dv=NULL, time, phase, ivs,
                fixed=NULL, random=NULL, formula=NULL,
                correlation, family,
                dvs, target_ivs, standardize)

  # subgroup the data and delete the parameter, after this point, it is only
  # used to subgroup to unique ids
  if( is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))
  if(!is.null(data)) data <- data[subgroup,]; rm(subgroup)

  ###########################################################################
  # 20180728 - commented out by Stephen Tueller when debugging metabolomics
  # some DSST output was missing time and phase (but not ids). This code is
  # not yet essential and may be the culprit. PalyticObj is use by utils in
  # PersonAlytics and we may have issues with scope allowing this to be
  # non-null (unlikely though, b/c this is before the loops and only some
  # loops are affected).
  ###########################################################################
  # if a PalyticObj is given, overwrite other objects, avoid if possible,
  # but we want to leave ids, dvs, phase, time required unless PalyticObj is
  # provided
  #if(!is.null(PalyticObj))
  #{
  #  if(! class(PalyticObj) %in% 'Palytic')
  #  {
  #    stop('PalyticObj is not a Palytic object. See ?Palytic')
  #  }
  #  ids=NULL
  #  phase=NULL
  #  time=NULL
  #}

  # check whether any variables in ivs are in target_ivs -
  # in the future, split them out automatically
  if(any(ivs %in% target_ivs) | any(target_ivs %in% ivs))
  {
    stop('target_ivs and ivs cannot share any variables.')
  }

  ## if no data are given, use a test data set
  if(is.null(data))
  {
    data   <- OvaryICT
    dvs    <- "follicles"
    phase  <- "Phase"
    ids    <- "Mare"
    time   <- "TimeSin"
  }

  # unique ids
  uids <- sort(as.numeric(unique(data[[ids]])))

  # dimensions for loops
  ID <- uids
  IV <- 1:length(target_ivs); if(is.null(target_ivs)) IV <- 1
  DV <- 1:length(dvs)
  dims <- list(ID=ID, IV=IV, DV=DV)

  #
  if( ind.mods )
  {
    DVout <- htp.foreach(data, dims, dvs, phase, ids, uids, time, ivs,
                         target_ivs, interactions, time_power, correlation,
                         family, standardize, package,
                         detectAR, PQ, IC,
                         detectTO, maxOrder = time_power,
                         sigma.formula, debugforeach)
  }
  if( !ind.mods )
  {
    grp.dims <- dims
    grp.dims$ID <- "All Cases"

    DVout <- htp.foreach(data, grp.dims, dvs, phase, ids, uids, time, ivs,
                         target_ivs, interactions, time_power, correlation,
                         family, standardize, package,
                         detectAR, PQ, IC,
                         detectTO, maxOrder = time_power,
                         sigma.formula, debugforeach)

  }

  # clean up variable names
  if(!is.null(charSub))
  {
    DVout$target_iv <- as.character(DVout$target_iv)
    for(i in 1:length(charSub))
    {
      DVout$target_iv <- gsub(charSub[[i]][1], charSub[[i]][2], DVout$target_iv)
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
  write.csv(DVout, file=file, row.names=FALSE)

  if(!is.null(p.method) & length(target_ivs) > 1) DVout <- psuite(DVout,
                                                                  rawdata=data,
                                                                  method=p.method,
                                                                  alpha=alpha)


  return(DVout)
}


#' pa1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
pa1 <- function(file           ,
                data           ,
                ids            ,
                dvs            ,
                time           ,
                phase          ,
                ivs            ,
                target_ivs     ,
                interactions   ,
                time_power     ,
                maxOrder       ,
                correlation    ,
                family         ,
                subgroup       ,
                standardize    ,
                package        ,
                ind.mods       ,
                PalyticObj     ,
                detectAR       ,
                PQ             ,
                IC             ,
                detectTO       ,
                charSub        ,
                sigma.formula  ,
                p.method       ,
                alpha          )
{
  # if no data are given, use a test data set
  if(is.null(data))
  {
    data   <- OvaryICT
    dvs    <- "follicles"
    ids    <- "Mare"
    time   <- "Time"
    phase  <- "Phase"
    ivs    <- NULL
    interactions<- NULL
    time_power  <- 1
    correlation <- NULL
    PQ      <- c(3,3)
    subgroup    <- NULL
    standardize <- FALSE
    package     <- 'nlme'
    IC=c("BIC", "AIC")
    lrt=FALSE
    alpha=.05
  }

  if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))

  t1 <- Palytic$new(data=data[subgroup,]      ,
                    dv=dvs                    ,
                    ids=ids                   ,
                    time=time                 ,
                    phase=phase               ,
                    ivs=ivs                   ,
                    interactions=interactions ,
                    time_power=time_power     ,
                    correlation=correlation   ,
                    standardize=standardize   )

  if(detectAR) t1$GroupAR_order(dV    = dvs    ,
                                maxAR = PQ[1]  ,
                                maxMA = PQ[2]  ,
                                IC    = IC[1]  )
  # t1$correlation
  # t1$formula

  if(detectTO) t1$GroupTime_Power(maxOrder)
  # t1$time_power
  # t1$formula

  if(package=="gamlss") Grp.out <- t1$gamlss()
  if(package=="nlme")   Grp.out <- t1$lme()

  return(Grp.out)

}
