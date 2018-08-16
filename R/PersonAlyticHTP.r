#' High through-put options for PalyticHTP objects.
#'
#' @export
#' @import PersonAlyticsLite
#' @import snow
#' @import doSNOW
#' @import foreach
#' @import plyr
#'
#' @description A user interface for creating a PalyticHTP object
#' (which inherits from \code{\link{Palytic}} and invoking
#' high through-put options including automated autocorrelation detection, automated
#' detection of a polynomial for time, and high throughput analyses of all individuals
#' in the data set, multiple dependent variables, and target independent variables.
#' Type I error corrections are implemented via [pending].
#'
#' @param file The file name (or full path with `/` instead of `\`) where output should
#' be saved. If left \code{NULL}, the date and time will prefix `PersonAlyticHTP_Output.csv`.
#' @param data See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#' @param dvs A list of one or more character dependent variable names in \code{data}.
#' The linear mixed effects model \code{dvs[d] ~ phase + time + phase*time + ivsl[c] + ivs}
#' with random effects \code{~ time | ids[i]} will be fit to the data using
#' \code{\link{gamlss}}. The iterators \code{[d]}, \code{[d]}, and \code{[d]}
#' indicates the model will be fit for each combination of dependent variable in \code{dvs},
#' independent variable in \code{ivsl}, and each unique ID in \code{ids}
#' (overridden using \code{ind.mods=FALSE}) controlling for indepented variables in
#' \code{ivs}. For more options submit a \code{PalyticObj}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param ivs See \code{\link{PersonAlytic}}, noting that the variables in \code{ivs}
#' cannot also be in \code{ivsl}.
#' @param ivls Independent variables that are iterated over one at a time. Effects for
#' these variables are labeled as 'target predictor' in the output.
#' @param interactions See \code{\link{PersonAlytic}}.
#' @param time_power See \code{\link{PersonAlytic}}. Ignored if \code{detectAR=TRUE}
#' @param correlation See \code{\link{PersonAlytic}}. Ignored if \code{detectTO=TRUE}
#' @param family See \code{\link{PersonAlytic}}.
#' @param subgroup See \code{\link{PersonAlytic}}.
#' @param standardize See \code{\link{PersonAlytic}}. The default is \code{TRUE}
#' which makes parameter estimate magnitudes comparable across individuals, outcomes in
#' \code{dvs}, and covariates in \code{ivls}.
#' @param package See \code{\link{PersonAlytic}}.
#' @param ind.mods Logical, defaults to \code{TRUE}. Should individual models be
#' fit for each ID?
#' @param grp.mod Logical, defaults to \code{FALSE}. Should a group level model be fit
#' across all IDs? If both \code{ind.mods} and \code{grp.mod} are \code{FALSE},
#' \code{grp.mod} will be changed to \code{TRUE}.
#' @param PalyticObj See \code{\link{PalyticHTP}}. If \code{PalyticObj} is submitted
#' then only \code{dvs}, \code{ivsl}, \code{ind.mods}, and \code{grp.mod} will be
#' used. This allows users access to additional options including generalized linear
#' mixed effects models via the \code{family} option, user specified \code{correlation}
#' structures (in non-\code{NULL} this will override the automated correlation structure
#' search), and user specified models via \code{formula}.
#' @param detectAR Logical, defaults to \code{TRUE}. Should the autoregressive structure
#' be automatically detected? If the \code{time} variable is equally spaced, this is
#' done using the function \code{\link{forecast}}. If the \code{time} variable is not
#' equally spaced, this is done using likelihood ratio tests with mixed effects models
#' using the specified \code{package} under maximum likelihood. This is done separately
#' for each case in \code{ids} if \code{ind.mods=TRUE}.
#' @param detectTO Logical, defaults to \code{TRUE}. Should the \code{time_power} value
#' be automatically detected? Values from 1 to \code{maxOrder} will be tested using
#' likelihood ratio tests with mixed effects models using the specified \code{package}
#' under maximum likelihood. This is done separately for each case in \code{ids}
#' if \code{ind.mods=TRUE}.
#' @param  maxOrder Numeric, defaults to 3. What is the highest order of \code{time}
#' that sholud be tested for both fixed and randmo effects, e.g., \code{time +
#' I(time^2)+...+I(time^maxOrder)}.
#' @param charSub list of paired character strings for character substitution.
#' If the names of the target predictors
#' in \code{ivsl} had to be edited to make valid variable names, this parameter allows
#' users put the illegal characters back in. For example, if the original variable name
#' was "17.00_832.2375m/z", a letter would need to prefix the variable name and the
#' "/" would need to be replaced with another character, e.g., "X17.00_832.2375m.z".
#' To get the row names of the output back to original varibale name, use
#' \code{charSub=list(c("X", ""), c("m.z", "m/z"))}. Note that inputs to charSub
#' must be in double quotes and are case sensitive. All duplicates will be substituted.
#' For example, if the variable name was "X1X23.x" and \code{charSub=list(c("X", ""))},
#' the resulting row label for this variable would be "123.x".
#' @param sigma.formula A formula for the variance under \code{\link{gamlss}}. Static.
#' It will not change dynamically over iterations nor will it be updated by \code{time_power}
#' or \code{detectTO}. If model fitting using this option fails, another attempt will be
#' made after reseting it to its defaul, i.e., \code{~1}.
#'
#' @examples
#' # group model
#' t1 <- PersonAlyticHTP(data=PersonAlyticsLite::OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="TimeSin",
#'                  package='nlme',
#'                  ind.mods=FALSE,
#'                  grp.mod=TRUE)
#' # individual models (using defaults)
#' t1 <- PersonAlyticHTP(data=PersonAlyticsLite::OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="TimeSin",
#'                  package='nlme')
#' hist(t1$Phase.Std.Error, main='Distribution of Phase Effect', xlab='Phase Effect')
#' \dontrun{
#' # if you wish to delete the automatically created csv file, run
#' #NOT IMPLEMENTED YET
#' }

PersonAlyticHTP <- function(file=NULL                ,
                            data=NULL                ,
                            ids                      ,
                            dvs                      ,
                            time                     ,
                            phase=NULL               ,
                            ivs=NULL                 ,
                            ivsl=NULL                ,
                            interactions=NULL        ,
                            time_power=1             ,
                            correlation=NULL         ,
                            family=gamlss.dist::NO() ,
                            subgroup=NULL            ,
                            standardize=TRUE         ,
                            package='gamlss'         ,
                            ind.mods=TRUE            ,
                            grp.mod=FALSE            ,
                            PalyticObj=NULL          ,
                            detectAR=TRUE            ,
                            detectTO=TRUE            ,
                            maxOrder=3               ,
                            charSub=NULL             ,
                            sigma.formula=~1         ,
                            debugforeach = FALSE     ,
                            p.method = "BY"          ,
                            alpha = .05               )
{
  #
  if(is.null(file))
  {
    file <- gsub(":", ".", paste(Sys.time(), 'PersonAlyticHTP_Output.csv'))
  }

  # check that dvs, ivls are lists, if not, force
  if( ! "list" %in% class(dvs) ) dvs <- as.list(dvs)
  if( ! "list" %in% class(ivs) ) ivs <- as.list(ivs)

  # check that inputs conform. this is also done when creating a PalyticHTP
  # object, but we do it early on here to avoid problems after loops start.
  # This is why `clean()` has inputs that apply to PersonAlyticHTP but not to
  # PersonAlytic
  data <- PersonAlyticsLite:::clean(data, ids, dv=NULL, time, phase, ivs,
                                    fixed=NULL, random=NULL, formula=NULL,
                                    correlation,
                                    dvs, ivsl, standardize)

  # subgroup the data and delete the parameter, after this point, it is only
  # used to subgroup to unique ids
  if( is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))
  if(!is.null(data)) data <- data[subgroup,]; rm(subgroup)

  ###########################################################################
  # 20180728 - commented out by Stephen Tueller when debugging metabolomics
  # some DSST output was missing time and phase (but not ids). This code is
  # not yet essential and may be the culprit. PalyticObj is use by utils in
  # PersonAlyticsLite and we may have issues with scope allowing this to be
  # non-null (unlikely though, b/c this is before the loops and only some
  # loops are affected).
  ###########################################################################
  # if a PalyticObj is given, overwrite other objects, avoid if possible,
  # but we want to leave ids, dvs, phase, time required unless PalyticObj is
  # provided
  #if(!is.null(PalyticObj))
  #{
  #  if(! class(PalyticObj) %in% 'PalyticHTP')
  #  {
  #    stop('PalyticObj is not a PalyticHTP object. See ?PalyticHTP')
  #  }
  #  ids=NULL
  #  phase=NULL
  #  time=NULL
  #}

  # check whether any variables in ivs are in ivsl -
  # in the future, split them out automatically
  if(any(ivs %in% ivsl) | any(ivsl %in% ivs))
  {
    stop('ivsl and ivs cannot share any variables.')
  }

  ## if no data are given, use a test data set
  if(is.null(data))
  {
    data   <- PersonAlyticsLite::OvaryICT
    dvs    <- "follicles"
    phase  <- "Phase"
    ids    <- "Mare"
    time   <- "TimeSin"
  }

  # unique ids
  uids <- sort(as.numeric(unique(data[[ids]])))

  # dimensions for loops
  ID <- uids
  IV <- 1:length(ivsl); if(is.null(ivsl)) IV <- 1
  DV <- 1:length(dvs)
  dims <- list(ID=ID, IV=IV, DV=DV)

  #
  if( ind.mods & grp.mod )
  {
    warning('Group and individual models must be run separately, ',
            'the current run will use `grp.mod`=TRUE and `ind.mods=FALSE`. ',
            'To get individual models use `grp.mod=FALSE`` and `ind.mods=TRUE`.')
  }

  #
  if( ind.mods )
  {
    DVout <- htp.foreach(data, dims, dvs, phase, ids, uids, time, ivs, ivsl,
                         interactions, time_power, correlation,
                         family = family, standardize, package,
                         detectAR, detectTO, maxOrder, sigma.formula, debugforeach)
  }
  if( grp.mod )
  {
    grp.dims <- dims
    grp.dims$ID <- "All Cases"

    DVout <- htp.foreach(data, grp.dims, dvs, phase, ids, uids, time, ivs, ivsl,
                         interactions, time_power, correlation,
                         family = family, standardize, package,
                         detectAR, detectTO, maxOrder, sigma.formula, debugforeach)

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

  if(!is.null(p.method) & length(ivsl) > 1) DVout <- psuite(DVout,
                                                            rawdata=data,
                                                            method=p.method,
                                                            alpha=alpha)


  return(DVout)
}
