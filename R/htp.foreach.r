#' high throughput
#'
#' @keywords internal

# the current state of using Palytic is to create one object for
# loops across individuals, but overwrite the the Palytic object
# for loops across dvs/ivs
htp.foreach <- function(data                       ,
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
                        debugforeach = FALSE       )
{
  library(foreach)


  ##############################################################################
  # outer loop is DV
  ##############################################################################
  DVout <- list()
  for(dv in dims$DV)
  {
    #...........................................................................
    # set up the parent palytic object
    # -- do NOT use $clone() as it has failed in tests (inheritance issues)
    # -- there are inheritance issues, t0 can be modified by its child t1,
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
                      method="ML"
                     )
    saveMethod  <- t0$method
    saveFormula <- t0$formula

    #...........................................................................
    # get the TO and AR
    #...........................................................................
    if(dims$ID[1]!="All Cases")
    {
      if(detectTO) t0$getTime_Power(maxOrder, whichIC[1])
      #t0$time_powers

      # individual AR orders are pending deprecation -- 20181001
      #if( detectAR) t0$getAR_order(PQ[1], PQ[2], whichIC[1])
      #if(!detectAR) t0$corStructs <- data.frame(ids=dims$ID,
      #                                          arma=rep(ifelse(is.null(t0$correlation),
      #                                                      "NULL", t0$correlation),
      #                                                   length(dims$ID)))
      #t0$corStructs

    }
    if(dims$ID[1]=="All Cases")
    {
      if(detectTO) t0$GroupTime_Power(maxOrder, whichIC[1])
      if(detectAR) t0$GroupAR_order(PQ[1], PQ[2], whichIC[1])
    }


    #...........................................................................
    # parralelization: encapsulate ivs and ids in one foreach
    #...........................................................................
    DIM <- expand.grid(ID=dims$ID, IV=dims$IV)
    funcs <- c()
    cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
    snow::clusterExport(cl, funcs)
    doSNOW::registerDoSNOW(cl)
    pkgs  <- c("gamlss", "nlme")
    pb <- txtProgressBar(max = nrow(DIM), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
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
      # copy the palytic object - this may not be neccessary; the goal is to
      # be able to relaim changes made by fitting functinos, but inheritent may
      # be making changes anyways
      #-------------------------------------------------------------------------
      t1 <- t0

      #-------------------------------------------------------------------------
      # initialize input and error accumulation
      #-------------------------------------------------------------------------
      err_id <- list()

      # temporary data for checking valid cases
      temp   <- t1$data[useObs,all.vars(t1$formula)]
      nrt    <- nrow(temp)

      # identify rows
      err_id[ids] <- id

      # identify inputs
      err_id['ids']          <- t1$ids
      err_id['dvs']          <- t1$dvs
      err_id['family']       <- t1$family[[1]][2]
      err_id['package']      <- t1$package
      err_id['time']         <- t1$time[1]
      err_id['phase']        <- t1$phase
      err_id['interactions'] <- toString( t1$interactions )
      err_id['time_power']   <- t1$time_power

      # identify package
      err_id['package'] <- package

      # adequate data
      err_id['Nobs'] <- paste('There are', nrow(na.omit(temp)),
                              'observations with complete data.')

      # outcome variance
      err_id['dvVar'] <- paste( 'The variance of `', t1$dv, '` is ',
                                round(var(temp[[t1$dv]], na.rm = TRUE),3),
                                ifelse(all(dims$ID=="All Cases"), ".",
                                       paste(' for ', t1$ids, ' = ', id, '.', sep='')),
                                sep='')

      # time should be monotonically increasing
      err_id['timeVar'] <- paste("`", time, "` is",
                                 ifelse(all(t1$monotone$monotonic), "" , "NOT"),
                                 "monotonically increasing for `", ids, "` =", id, ".")

      # iv variance
      ivvs <- c(t1$ivs, t1$phase)
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

      #-------------------------------------------------------------------------
      # fit the model w/o target ivs
      #-------------------------------------------------------------------------
      if(package=="gamlss") mod1 <- t1$gamlss( useObs )
      if(package=="nlme"  ) mod1 <- t1$lme(    useObs )
      if(package=="arma"  ) mod1 <- t1$arma(   useObs, max.p=PQ[1], max.q=PQ[2] )

      #-------------------------------------------------------------------------
      # for the current id, estimate a full model with the current target IV
      #-------------------------------------------------------------------------
      # add the target IV
      ivs.temp <- unlist(c(ivs, target_ivs[iv]))
      if( is.null(ivs.temp) )
      {
        t1$ivs <- ivs.temp
      }
      if(!is.null(ivs.temp) )
      {
        t1$ivs <- gsub(" ", "", ivs.temp)
      }
      err_id['ivs']       <- toString( t1$ivs         )
      err_id['target_iv'] <- toString( target_ivs[iv] )
      # fit models
      if(package=="gamlss")
      {
        modid <- t1$gamlss( useObs )

        if(debugforeach){
          cat('modid for id = ', id, ' is a ', class(modid), '\n\n',
              file = paste('htp', dvs[dv], target_ivs[iv], 'log', sep='.'),
              append = TRUE)
        }

        if(! "gamlss"  %in%  class(modid) )
        {
          err_id['converge'] <- "Model failed to converge"
          err_id['re_estimator'] <- NA
          err_id['gamlss_estimator'] <- NA
          err_id['analyzed_N'] <- NA
        }
        if(  "gamlss"  %in%  class(modid) )
        {
          err_id['converge'] = paste('Convergence is `', modid$converged,
                                     '`', sep='')
          has_method <- grepl("method", t1$formula)
          if(any(has_method))
          {
            frm_txt <- as.character(t1$formula)[which(has_method)]
            frm_txt <- unlist( strsplit(frm_txt, ',') )
            frm_txt <- frm_txt[which( grepl("method", frm_txt) )]
            err_id['re_estimator'] <- gsub("method|\"| |=", "", frm_txt)
          }
          if(!any(has_method)) err_id['re_estimator'] <- "Cannot be determined"

          err_id['gamlss_estimator'] <- modid$PalyticSummary$method
          err_id['analyzed_N'] <- paste(modid$N, 'cases were analyzed.')
          err_id['call'] <- modid$PalyticSummary$formula
        }
      }

      if(package=="nlme")
      {
        modid <- t1$lme( useObs )
        if(! "lme"  %in%  class(modid) )
        {
          err_id['converge']   <- modid
          err_id['estimator']  <- t1$method
          err_id['analyzed_N'] <- NA
          err_id['call'] <- paste(Reduce( paste, deparse( t1$fixed ) ),
                                  Reduce( paste, deparse( t1$random ) ),
                                  t1$correlation,
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
      }

      if(package=="arma")
      {
        modid <- t1$arma( useObs, max.p=PQ[1], max.q=PQ[2] )
        if(! "coeftest"  %in%  class(modid$tTable) )
        {
          err_id['converge']   <- modid
          err_id['estimator']  <- t1$method
          err_id['analyzed_N'] <- NA
          err_id['call'] <- paste(Reduce( paste, deparse( t1$fixed ) ),
                                  Reduce( paste, deparse( t1$random ) ),
                                  t1$correlation,
                                  sep='; ')
        }
        if(  "coeftest"  %in%  class(modid$tTable) )
        {
          err_id['converge']   <- 'Convergence is `TRUE`'
          err_id['estimator']  <- modid$PalyticSummary$method
          err_id['analyzed_N'] <- paste(modid$dims$N, 'cases were analyzed.')
          err_id['call'] <- paste( Reduce( paste, deparse(modid$PalyticSummary$fixed) ),
                                   Reduce( paste, deparse(modid$PalyticSummary$random) ),
                                   modid$PalyticSummary$correlation,
                                   sep='; ')
        }
      }

      #-------------------------------------------------------------------------
      # LRT for target iv
      #-------------------------------------------------------------------------
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
      err_id[['targ_ivs_lrt_pvalue']] <- lrtp

      #-------------------------------------------------------------------------
      # populate IVout - consider writing to file here
      #-------------------------------------------------------------------------
      cln <- function(x)
      {
        gsub("\n|\\s\"|\'", "", paste(format(x), collapse=" "))
      }
      cl2 <- function(x)
      {
        if(is.null(x)) return("NULL")
        else return(x)
      }
      cla <- function(x)
      {
        ao <- forecast::arimaorder(x)
        paste('auto.arima::ARMA(p=', ao[1], ', q=', ao[2], ')', sep='')
      }
      err_id[["target_iv"]]   = cl2(unlist(target_ivs)[iv])
      err_id[["fixed"]]       = cln(t1$fixed)
      err_id[["random"]]      = cln(t1$random)
      err_id[["correlation"]] = ifelse(all(dims$ID=="All Cases"),
                                  cln(t1$correlation),
                                  ifelse(package=="arma",
                                         cla(modid$arima),
                                         cln(t1$corStructs[id,2])))
      err_id[["formula"]]     = cln(t1$formula)
      err_id[["directory"]]   = normalizePath(getwd())
      err_id[["date"]]        = toString( Sys.time() )

      #-------------------------------------------------------------------------
      # return messages and re-run models with REML (unless arma)
      #-------------------------------------------------------------------------
      Model = "NA"
      if(any(c("gamlss", "lme") %in% class(modid)))
      {
        t1$method <- "REML"
        if("gamlss" %in% class(modid)) Model <- t1$gamlss( useObs )
        if("lme"    %in% class(modid)) Model <- t1$lme( useObs )
      }
      if( any( c("ARIMA", "Arima") %in% class(modid$arima) ) )
      {
        Model <- modid
      }

      return( list( Messages=err_id, Model=Model ) )

      # reclaim inputs in case there are inheritance issues due to $gamlss or
      # $lme dropping terms
      t0$method  <- saveMethod
      t0$formula <- saveFormula

      # clean up
      rm(t1)
    } # end of foreach
    # stop the cluster
    parallel::stopCluster(cl)

    #...........................................................................
    # dis aggregate messages
    #...........................................................................
    IDmsg <- lapply( IDout, function(x) x$Messages )
    IDmsg <- data.frame( do.call(rbind, IDmsg) )

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
      IDoutSum <- lapply( IDout, function(x){ if("ARIMA" %in% class(x$arima)){
        x$tTable} else NA})
    }

    #...........................................................................
    # function to prep row names
    # this function is fragile, we need it to adapt to any number of fixed ivs
    # and use a common label for each target_ivs
    #...........................................................................
    rcn <- function(x)
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

    #...........................................................................
    # function to force values to data.frame vectors for later stacking
    #...........................................................................
    rcm <- function(x)
    {
      if(is.matrix(x))
      {
        xo <- data.frame( t(as.vector(t(x)) ) )
        row.names(xo) <- NULL
        colnames(xo)  <- rcn(x)
      }
      if(!is.matrix(x))
      {
        xo <- data.frame(NA)
      }
      return( xo )
    }

    #...........................................................................
    # turn the numeric output into vectors & rbind them together
    #...........................................................................
    IDoutSumm <- lapply(IDoutSum, rcm)
    IDoutSumm <- plyr::rbind.fill(IDoutSumm)

    DVout[[dvs[[dv]]]] <- data.frame(IDmsg, IDoutSumm)
  } # end of dv loops


  # final post-processing
  outmat <- plyr::rbind.fill(DVout)
  row.names(outmat) <- NULL
  return( outmat )
}

