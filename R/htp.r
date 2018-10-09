#' high throughput
#'
#' @keywords internal

# the current state of using Palytic is to create one object for
# loops across individuals, but overwrite the the Palytic object
# for loops across dvs/ivs
htp <- function(data                       ,
                dims                       ,
                dvs                        ,
                phase                      ,
                ids                        ,
                uids                       ,
                time                       ,
                ivs                        ,
                target_ivs                 ,
                interactions=NULL          ,
                time_power=1               ,
                correlation=NULL           ,
                family = gamlss.dist::NO() ,
                standardize=TRUE           ,
                package='gamlss'           ,
                detectAR = TRUE            ,
                PQ = c(3, 3)               ,
                whichIC = "BIC"            ,
                detectTO = TRUE            ,
                maxOrder=3                 ,
                sigma.formula=~1           ,
                debugforeach = FALSE       ,
                userFormula = list(
                fixed=NULL,
                random=NULL,
                formula=NULL)              )
{
  library(foreach)

  ##############################################################################
  # parralelization: encapsulate ivs and ids in one foreach
  ##############################################################################
  DIM <- expand.grid(ID=dims$ID, IV=dims$IV)
  funcs <- c()
  pkgs  <- c("gamlss", "nlme")
  pb <- txtProgressBar(max = nrow(DIM), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  ##############################################################################
  # outer loop is DV
  ##############################################################################
  DVout <- list()
  for(dv in dims$DV)
  {
    #...........................................................................
    # set up the parent palytic object
    # -- do NOT use $clone() as it has failed in tests (inheritance issues)
    # -- there are inheritance issues, t0 can be modified by its child t0,
    #    hence the save* objects
    #...........................................................................
    t0 <- Palytic$new(data=data,
                      ids=ids,
                      dv=dvs[[dv]],
                      time=time,
                      phase=phase,
                      ivs=ivs, # target_ivs added later
                      interactions=interactions,
                      time_power=time_power,
                      correlation=correlation,
                      family=family,
                      method="ML" # requested method used in final estimation
                     )

    # allow for formula override so that we can test intercept only and
    # slope only models
    if( any(unlist(lapply(userFormula, function(x) !is.null(x)))) )
    {
      if( isNullOrForm(e$userFormula$fixed) ) t0$fixed <- userFormula$fixed
      if( isNullOrForm(e$userFormula$random) ) t0$random <- userFormula$random
      if( isNullOrForm(e$userFormula$formula) ) t0$formula <- userFormula$formula
    }

    # save the inputs as they may be overwritten by $lme, $gamlss, or $arma
    # TODO(Stephen): part of this is due to problems in forms() and other
    # functions, we need to refactor and clean the problem up
    saveMethod  <- t0$method
    saveFormula <- t0$formula
    saveRandom  <- t0$random
    saveIVS     <- t0$ivs
    savePhase   <- t0$phase

    #...........................................................................
    # get the TO and AR
    #...........................................................................
    if(dims$ID[1]!="All Cases")
    {
      if(detectTO) t0$getTime_Power(maxOrder, whichIC[1])
      #t0$time_powers # things like this should be changed to unit tests

      if( detectAR & package != "arma") t0$getAR_order(PQ[1], PQ[2], whichIC[1])
      if(!detectAR)
      {
        t0$corStructs <- data.frame(ids=dims$ID,
                                    arma=rep( ifelse(is.null(t0$correlation),
                                                     "NULL", t0$correlation),
                                                     length(dims$ID)) )
      }
      #t0$corStructs

    }
    if(dims$ID[1]=="All Cases")
    {
      if(detectTO) t0$GroupTime_Power(maxOrder, whichIC[1])
      if(detectAR) t0$GroupAR_order(PQ[1], PQ[2], whichIC[1])
    }

    #...........................................................................
    # start parralelization run
    #...........................................................................

    message('\n\nFitting models of the dependent variable `', dvs[dv],
            '` for ',
            ifelse( length(dims$ID)>0,
                    paste(length(dims$ID), 'cases in `', ids, '`'),
                    dims$ID),
            ifelse( length(target_ivs)>0,
                    paste('\n and for', length(target_ivs),
                          'target indepented variables in `target_ivs`.\n'),
                    ".")
            )
    start <- Sys.time()

    # these must be rerun after each call to stopCluster
    cl <- snow::makeCluster(parallel::detectCores(), type="SOCK")
    snow::clusterExport(cl, funcs)
    doSNOW::registerDoSNOW(cl)

    #IDout <- list()
    IDout <- foreach( id=DIM$ID, iv=DIM$IV,
            .packages = pkgs, .options.snow = opts) %dopar%
    #for(i in 1:nrow(DIM))
    {
      #id=DIM$ID[i]; iv=DIM$IV[i]
      #-------------------------------------------------------------------------
      # for the current id, select rows
      #-------------------------------------------------------------------------
      if(dims$ID[1]!="All Cases")
      {
        useObs <- t0$data[[t0$ids]]==id
        wid    <- which(dims$ID==id)
      }
      if(dims$ID[1]=="All Cases")
      {
        useObs <- rep(TRUE, nrow(t0$data))
        wid    <- 1:nrow(t0$data)
      }

      #-------------------------------------------------------------------------
      # accumulate inputs and errors for the output
      #-------------------------------------------------------------------------
      err_id <- htpErrors(t0, id, dims, useObs)

      #-------------------------------------------------------------------------
      # fit the model w/o target ivs
      # TODO (Stphen): escape if no variance, this may increase speed
      #-------------------------------------------------------------------------
      mod1 <- fitWithTargetIV(t0, package, useObs, dims, PQ=PQ)

      #-------------------------------------------------------------------------
      # for the current id, estimate a full model with the current target IV
      #-------------------------------------------------------------------------
      err_id['target_ivVar'] <- "No target_ivs"
      if( length( target_ivs[iv] ) > 0 )
      {
        # check for 0 variance in the target iv
        tivv <- !all(duplicated(data[[target_ivs[iv]]][data[[ids]]==id])[-1L])
        err_id['target_ivVar'] <- paste('The variance of `',
                                      target_ivs[iv], '` is ',
                                      ifelse(tivv, '> 0 ', '= 0 '),
                                      sep='')

        # add the target IV
        err_id['target_iv'] <- toString( target_ivs[iv] )
        ivs.temp <- unlist(c(ivs, target_ivs[iv]))
        if( is.null(ivs.temp) )
        {
          t0$ivs <- ivs.temp
        }
        if(!is.null(ivs.temp) )
        {
          t0$ivs <- gsub(" ", "", ivs.temp)
        }

        #TODO(Stephen): override correlation search for ARMA?
        fitOutput <- fitWithTargetIV(t0, package, useObs, dims, mod1$modid, PQ)
        err_id <- c(err_id, fitOutput$err_id)
        modid  <- fitOutput$modid
        rm(fitOutput)
      }
      if( length( target_ivs[iv] ) == 0 )
      {
        modid <- mod1$modid
        err_id <- c(err_id, mod1$err_id)
      }

      #-------------------------------------------------------------------------
      # descriptive statistics for target iv
      #-------------------------------------------------------------------------
      err_id <- c(err_id, t0$describe(useObs))

      #-------------------------------------------------------------------------
      # re-fit models with REML (unless arma)
      # TODO(Stephen): this *shouldn't* need the errror accumulation from
      # fitWithTargetIV unless the REML fit leads to a different model than
      # the ML fit
      #-------------------------------------------------------------------------
      Model = "NA"
      if(any(c("gamlss", "lme") %in% class(modid)))
      {
        t0$method <- "REML"
        if("gamlss" %in% class(modid)) Model <- t0$gamlss( useObs )
        if("lme"    %in% class(modid)) Model <- t0$lme( useObs )
      }
      if( any( c("ARIMA", "Arima") %in% class(modid$arima) ) )
      {
        Model <- modid
      }

      #-------------------------------------------------------------------------
      # add final entries to err_id, these may depend on final results
      #-------------------------------------------------------------------------
      err_id["target_iv"]   = nullString(unlist(target_ivs)[iv])
      err_id["fixed"]       = rmSpecChar(t0$fixed)
      err_id["random"]      = rmSpecChar(t0$random)
      err_id["correlation"] = ifelse(all(dims$ID=="All Cases"),
                                     rmSpecChar(t0$correlation),
                                     ifelse(package=="arma",
                                            gerARIMAorder(modid$arima),
                                            rmSpecChar(t0$corStructs[id,2])))
      err_id["formula"]     = rmSpecChar(t0$formula)
      err_id["directory"]   = normalizePath(getwd())
      err_id["date"]        = toString( Sys.time() )

      #-------------------------------------------------------------------------
      # reclaim inputs in case there are inheritance issues due to $gamlss or
      # $lme dropping terms
      #-------------------------------------------------------------------------
      t0$method  <- saveMethod
      t0$formula <- saveFormula
      t0$random  <- saveRandom
      if(length(saveIVS) > 0)  t0$ivs <- saveIVS
      if(length(saveIVS) == 0) t0$ivs <- NULL
      t0$phase <- savePhase

      #-------------------------------------------------------------------------
      # for a clean print after the progress bar
      #-------------------------------------------------------------------------
      cat('\n\n')

      #-------------------------------------------------------------------------
	    # return to foreach
      #-------------------------------------------------------------------------
      # this line stays commented  out except for testing
      #IDout[[i]] <- list( Messages=err_id, Model=Model )
      return( list( Messages=err_id, Model=Model ) )

    } # end of foreach
    # stop the cluster
    parallel::stopCluster(cl)

    message('\n\nModel fitting of the dependent variable `', dvs[dv],
            '` took: ', capture.output(Sys.time() - start), ".\n\n")

    #...........................................................................
    # dis aggregate messages
    #...........................................................................
    IDmsg <- lapply( IDout, function(x) data.frame(x$Messages) )
    IDmsg <- plyr::rbind.fill(IDmsg)

    #...........................................................................
    # post-estimation model extraction
    #...........................................................................
    IDout <- lapply( IDout, function(x) x$Model )
    if(dims$ID[1]!="All Cases") names(IDout) <- paste(ids, uids, sep=".")
    if(package=='gamlss')
    {
      IDoutSum <- lapply( IDout, function(x)
        {
          if(package %in% class(x))
          {
            capture.output( temp <- summary(x), file='NUL' )
          }
          if( all( x == "NA" ) ) temp <- NA
          return(temp)
        } )
    }
    if(package=='nlme') # can't use if(package %in% class(x)) b/c lme!=nlme
    {
      IDoutSum <- lapply( IDout, function(x){ if("lme" %in% class(x)){
        summary(x)$tTable} else NA})
    }
    if(package=="arma")
    {
      IDoutSum <- lapply( IDout, getARMAtbl)
    }

    #...........................................................................
    # turn the numeric output into vectors & rbind them together
    #...........................................................................
    IDoutSumm <- lapply(IDoutSum, rcm, target_ivs=target_ivs)
    IDoutSumm <- plyr::rbind.fill(IDoutSumm)

    DVout[[dvs[[dv]]]] <- data.frame(IDmsg, IDoutSumm)
  } # end of dv loops


  # final post-processing
  outmat <- plyr::rbind.fill(DVout)
  row.names(outmat) <- NULL
  return( outmat )
}

#' htpErrors: accumulate inputs and errors
#' @keywords internal
htpErrors <- function(t0, id, dims, useObs)
{
  err_id <- list()

  # temporary data for checking valid cases
  temp   <- t0$data[useObs, all.vars(t0$formula)]
  nrt    <- nrow(temp)

  # identify rows
  err_id[t0$ids] <- id

  # identify inputs
  err_id['ids']          <- t0$ids
  err_id['dvs']          <- t0$dvs
  err_id['family']       <- t0$family[[1]][2]
  err_id['package']      <- t0$package
  err_id['time']         <- t0$time[1]
  err_id['phase']        <- t0$phase
  err_id['interactions'] <- toString( t0$interactions )
  err_id['time_power']   <- t0$time_power

  # adequate data
  err_id['Nobs'] <- paste('There are', nrow(na.omit(temp)),
                          'observations with complete data.')

  # outcome variance
  err_id['dvVar'] <- paste( 'The variance of `', t0$dv, '` is ',
                            round(var(temp[[t0$dv]], na.rm = TRUE),3),
                            ifelse(all(dims$ID=="All Cases"), ".",
                                   paste(' for ', t0$ids, ' = ', id, '.', sep='')),
                            sep='')
  # time should be monotonically increasing
  err_id['timeVar'] <- paste("`", t0$time, "` is",
                             ifelse(all(t0$monotone$monotonic), "" , "NOT"),
                             "monotonically increasing for `", t0$ids, "` =", id, ".")

  # ivs (!target_ivs yet) and ivs variance
  err_id['ivs']       <- toString( t0$ivs         )
  ivvs <- c(t0$ivs, t0$phase)
  wsplit <- unlist(lapply(ivvs, grepl, pattern=":"))
  if( any(wsplit) )
  {
    ivvs <- unique( c(ivvs[!wsplit], unlist( strsplit(unlist(ivvs[wsplit]), ":") ) ) )
  }
  wsplit <- unlist(lapply(ivvs, grepl, pattern="\\*"))
  if( any(wsplit) )
  {
    ivvs <- unique( c(ivvs[!wsplit],
                      unlist( strsplit(unlist(ivvs[wsplit]), "\\*") ) ) )
  }
  ivvs <- unique( gsub(" ", "", unlist(ivvs)) )
  ivv  <- unlist( lapply(data.frame(temp[,ivvs]),
                         function(x) !all(duplicated(x)[-1L])) )
  err_id['ivVar'] <- paste( paste('The variance of `', ivvs, '` is ',
                                  ifelse(ivv, '> 0 ', '= 0 '),
                                  sep=''), collapse='; ')

  return(err_id)
}

#' rcm: function to force values to data.frame vectors for later stacking
#'
#' @keywords internal
rcm <- function(x, target_ivs)
{
  if(is.matrix(x))
  {
    xo <- data.frame( t(as.vector(t(x)) ) )
    row.names(xo) <- NULL
    colnames(xo)  <- rcn(x, target_ivs)
  }
  if(!is.matrix(x))
  {
    xo <- data.frame(NA)
  }
  return( xo )
}

#' rcn: function to prep row names
#'
#' @keywords internal
rcn <- function(x, target_ivs)
{
  if(is.matrix(x))
  {
    # extract the name(s) of the rows corresponding to the target predictor
    Rns  <- lapply(target_ivs, function(y,x) which( grepl(y, row.names(x)) ), x)
    if(length(Rns)==1)
    {
      Tnms <- row.names(x)[ unlist(Rns[lapply(Rns, length) > 0]) ]

      # if target predictor ! factor, set the name
      if( length(Tnms) == 1 )
      {
        row.names(x)[ which(row.names(x) %in% Tnms) ] <- 'Targ'
      }

      # if target predicor is a factor variable, get the reference category
      # and set the names
      w <- which(lapply(Rns, length) > 0)
      if( length(w) > 0 )
      {
        if( is.factor(data[[target_ivs[[w]]]]) )
        {
          uT     <- levels(data[[target_ivs[[w]]]])
          comCat <- unlist( strsplit(Tnms, target_ivs[[w]]) )
          comCat <- comCat[comCat != ""]
          refCat <- uT[! uT %in% comCat]
          NewNms <- paste('Targ', comCat, 'vs', refCat, sep='_')
          row.names(x)[ which(row.names(x) %in% Tnms) ] <- NewNms
          rm(uT, comCat, refCat, NewNms)
        }
      }

      rm(Tnms)
    }

    xn <- apply(expand.grid(colnames(x), row.names(x)), 1,
                function(x) paste(x[2], x[1], sep=" "))
  }
  if(!is.matrix(x))
  {
    xn <- data.frame(NA)
  }
  return( xn )
}

#' rmSpecChar - strip out special characters from strings
#' @keywords internal
rmSpecChar <- function(x)
{
  gsub("\n|\\s\"|\'", "", paste(format(x), collapse=" "))
}

#' nullString - turn NULL to "NULL"
#' @keywords internal
nullString <- function(x)
{
  if(is.null(x)) return("NULL")
  else return(x)
}

#' gerARIMAorder
#' @keywords internal
gerARIMAorder <- function(x)
{
  if( "Arima" %in% class(x) )
  {
    ao <- forecast::arimaorder(x)
    return( paste('auto.arima::ARMA(p=', ao[1],
                  ', q=', ao[2], ')', sep='') )
  }
  if(! "Arima" %in% class(x) )
  {
    return( NA )
  }
}

#' getARMAtbl: get ARMA tTable
#' @keywords internal
getARMAtbl <- function(x)
{
  if( is.list(x) )
  {
    if("ARIMA" %in% class(x$arima)) return(x$tTable)
  }
  else return(NA)
}

#' isNullOrForm: is an object non-NULL and a formula?
isNullOrForm <- function(x)
{
  if( !is.null(x) )
  {
    return( is.formula(x) )
  }
  else return(FALSE)
}

#TODO(Stephen): add ... to these functions
#' fit models with the target iv & calculate LRT
#' @keywords internal
fitWithTargetIV <- function(t0, package, useObs, dims, mod1=NULL, PQ=c(3,3))
{
  # fit model with targe iv
  if(package=="nlme")   modid <- fitWithTargetIVlme(t0, useObs, dims)
  if(package=="arma")   modid <- fitWithTargetIVarma(t0, useObs, dims, PQ)
  if(package=="gamlss") modid <- fitWithTargetIVgamlss(t0, useObs, dims)

  err_id <- modid$err_id
  modid <- modid$modid

  if(!is.null(mod1))
  {
    # LRT for target iv ####
    lrtp <- NA
    if(any(c("gamlss", "lme") %in% class(mod1) ) &
       any(c("gamlss", "lme") %in% class(modid)) &
       length(dims$IV) > 0
    )
    {
      mod1$call$fixed    <- mod1$PalyticSummary$fixed
      modid$call$fixed   <- modid$PalyticSummary$fixed
      mod1$call$formula  <- mod1$PalyticSummary$formula
      modid$call$formula <- modid$PalyticSummary$formula
      lrt  <- anova(mod1, modid)
      if(nrow(lrt)==2 & "p-value" %in% names(lrt))
      {
        lrtp <- lrt$"p-value"[2]
      }
    }
    if( any(c("ARIMA", "Arima") %in% class(mod1$arima) ) &
        any(c("ARIMA", "Arima") %in% class(modid$arima))  )
    {
      l0 <- logLik(mod1$arima)
      l1 <- logLik(modid$arima)
      df0 <- strsplit( unlist( strsplit(capture.output(l0), "=") )[2] , ")")
      df1 <- strsplit( unlist( strsplit(capture.output(l1), "=") )[2] , ")")
      df0 <- as.numeric( unlist(df0) )
      df1 <- as.numeric( unlist(df1) )
      lrtest <- as.numeric(2*(l1-l0))
      lrtp <- pchisq(lrtest, df1-df0, lower.tail = FALSE)
    }
    err_id['targ_ivs_lrt_pvalue'] <- lrtp
  }

  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVlme
#' @keywords internal
fitWithTargetIVlme <- function(t0, useObs, dims)
{
  err_id <- list()
  modid <- t0$lme( useObs )
  if(! "lme"  %in%  class(modid) )
  {
    err_id['converge']   <- modid
    err_id['estimator']  <- t0$method
    err_id['analyzed_N'] <- NA
    err_id['call'] <- paste(Reduce( paste, deparse( t0$fixed ) ),
                            Reduce( paste, deparse( t0$random ) ),
                            t0$correlation,
                            sep='; ')
    # here is a placeholder for getting error messsages from lme
    # which needs to be updated in the error handling for
    # lme in Palytic
  }
  if(  "lme"  %in%  class(modid) )
  {
    err_id['converge']   <- 'Convergence is `TRUE`'
    err_id['estimator']  <- modid$PalyticSummary$method
    err_id['analyzed_N'] <- paste(modid$dims$N, 'cases were analyzed.')
    err_id['call'] <- paste( Reduce( paste, deparse(modid$PalyticSummary$fixed) ),
                             Reduce( paste, deparse(modid$PalyticSummary$random) ),
                             modid$PalyticSummary$correlation,
                             sep='; ')
  }
  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVarma
#' @keywords internal
fitWithTargetIVarma <- function(t0, useObs, dims, PQ)
{
  err_id <- list()
  modid <- t0$arma( useObs, max.p=PQ[1], max.q=PQ[2] )
  if(! "coeftest"  %in%  class(modid$tTable) )
  {
    err_id['converge']   <- modid$arima
    #t0$method --- ML currently default in arma, REML not an option
    err_id['estimator']  <- "ML"
    err_id['analyzed_N'] <- NA
    err_id['call'] <- paste(Reduce( paste, deparse( t0$fixed ) ),
                            Reduce( paste, deparse( t0$random ) ),
                            t0$correlation,
                            sep='; ')
  }
  if(  "coeftest"  %in%  class(modid$tTable) )
  {
    err_id['converge']   <- 'Convergence is `TRUE`'
    err_id['estimator']  <- "ML" #modid$PalyticSummary$method
    err_id['analyzed_N'] <- paste(modid$arima$nobs, 'cases were analyzed.')
    armaOrder <- arimaorder(modid$arima)
    err_id['call'] <- paste( "auto.arima(y=", t0$dv,
                             ", order=c(", armaOrder[1], 0, armaOrder[3],
                             "), xreg=rxeg) where `xreg` includes ",
                             modid$xregs, ".")
  }
  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVgamlss
#' @keywords internal
fitWithTargetIVgamlss <- function(t0, useObs, dims)
{
  err_id <- list()
  modid <- t0$gamlss( useObs )

  if(! "gamlss"  %in%  class(modid) )
  {
    err_id['converge']         <- "Model failed to converge"
    err_id['re_estimator']     <- NA
    err_id['gamlss_estimator'] <- NA
    err_id['analyzed_N']       <- NA
  }
  if(  "gamlss"  %in%  class(modid) )
  {
    err_id['converge'] = paste('Convergence is `', modid$converged,
                               '`', sep='')
    has_method <- grepl("method", t0$formula)
    if(any(has_method))
    {
      frm_txt <- as.character(t0$formula)[which(has_method)]
      frm_txt <- unlist( strsplit(frm_txt, ',') )
      frm_txt <- frm_txt[which( grepl("method", frm_txt) )]
      err_id['re_estimator'] <- gsub("method|\"| |=", "", frm_txt)
    }
    if(!any(has_method)) err_id['re_estimator'] <- "Cannot be determined"

    err_id['gamlss_estimator'] <- modid$PalyticSummary$method
    err_id['analyzed_N'] <- paste(modid$N, 'cases were analyzed.')
    err_id['call'] <- modid$PalyticSummary$formula
  }
  return( list(err_id = err_id, modid = modid) )
}

