#' high throughput
#'
#' @keywords internal

# the current state of using PalyticHTP is to create one object for
# loops across individuals, but overwrite the the PalyticHTP object
# for loops across dvs/ivs
htp <- function(data, dims, dvs, phase, ids, uids, time, ivs, ivsl,
                interactions=NULL, time_power=1, correlation=NULL,
                family = gamlss.dist::NO(), standardize=TRUE, package='gamlss',
                detectAR = TRUE, detectTO = TRUE, sigma.formula=~1)
{

  # function to preclean the data to allow for escapes in loops
  # may not be needed, escape via 'model didn't converge'
  preclean <- function()
  {

  }

  # outer loop for dependent variables
  DVout <- list()
  for(dv in dims$DV)
  {
    t1 <- PalyticHTP$new(data=data,
                      ids=ids,
                      dv=dvs[[dv]],
                      time=time,
                      phase=phase,
                      ivs=NULL, # NULL for now, populate later
                      interactions=interactions,
                      time_power=time_power,
                      correlation=correlation,
                      family=family,
                      method="REML"
                      )
    if(detectAR) t1$getAR_order(dvs[[dv]])
    if(!detectAR) t1$corStructs <- data.frame(ids=dims$ID,
                                              arma=rep("NULL", length(dims$ID)))
    #t1$corStructs

    if(detectTO) t1$getTime_Power()
    #t1$time_powers

    # middle loop for covariates
    IVout <- list()
    for(iv in dims$IV)
    {
      # note that ivs should NOT be t1$ivs, pull from the input, not the t1
      # object, otherwise you'll get accumulation
      ivs.temp <- unlist(c(ivs, ivsl[iv]))
      t1$ivs <- gsub(" ", "", ivs.temp)

      # inner loop for participants
      IDout <- list()
      for(id in dims$ID)
      {
        cat('starting id = ', id, '\n\n',
            file = paste('htp', dvs[dv], ivsl[iv], 'log', sep='.'),
            append = TRUE)

        # select rows
        useObs <- t1$data[[t1$ids]]==id
        wid    <- which(dims$ID==id)
        ### error handling
        err_id <- list()
        temp   <- t1$data[useObs,all.vars(t1$formula)]
        nrt    <- nrow(temp)
        # identify rows
        err_id[ids] <- id
        # adequate data
        err_id['Nobs'] <- paste('There are', nrow(na.omit(temp)),
                                'time points with complete data.')
        # outcome variance
        err_id['dvVar'] <- paste( 'The variance of `', t1$dv, '` is ',
                                  round(var(temp[[t1$dv]], na.rm = TRUE),3),
                                  ' for ', t1$ids, ' = ', id, '.', sep='')
        # time should be monotonically increasing
        err_id['timeVar']  <- paste( '`', t1$time, '` is monotonically increasing.')
        if(! all( diff(temp[[t1$time]][temp[[t1$ids]]==id]) >= 0 ) )
        {
          err_id['timeVar']  <- paste( '`', t1$time, '` is NOT monotonically increasing.')
        }
        # iv variance
        ivvs <- c(t1$ivs, t1$phase)
        wsplit <- unlist(lapply(ivvs, grepl, pattern=":"))
        ivvs <- unique( c(ivvs[!wsplit], unlist( strsplit(ivvs[wsplit], ":") ) ) )
        wsplit <- unlist(lapply(ivvs, grepl, pattern="\\*"))
        ivvs <- unique( c(ivvs[!wsplit], unlist( strsplit(ivvs[wsplit], "\\*") ) ) )
        ivv  <- unlist( lapply(data.frame(temp[,ivvs]),
                               function(x) !all(duplicated(x)[-1L])) )
        err_id['ivVar'] <- paste( paste('The variance of `', ivvs, '` is ',
                                        ifelse(ivv, '> 0 ', '= 0 '),
                                        sep=''), collapse='; ')
        # update correlation and fit model
        if(is.null(t1$corStructs[wid,2]) |
           length(t1$corStructs[wid,2]) == 0 ) t1$correlation <- "NULL"
        if(  is.na(t1$corStructs[wid,2]) )     t1$correlation <- "NULL"
        if(! is.null(t1$corStructs[wid,2]) )   t1$correlation <-
          as.character( t1$corStructs[wid,2] )

        cat('correllation for id = ', id, ' is ', t1$correlation, '\n\n',
            file = paste('htp', dvs[dv], ivsl[iv], 'log', sep='.'),
            append = TRUE)

        if(package=="gamlss")
        {
          # Not sure why this fails, it prints the correct deviance, and yet
          # the model did not converge error is returned, even when running
          # the subsequent line of code works
          modid <- t1$gamlss( useObs )
          #modid <- try( gamlss::gamlss(t1$formula,
          #                        data=na.omit(t1$data[t1$data[[t1$ids]]==id,
          #                                             all.vars(t1$formula)])),
          #              TRUE )

          cat('modid for id = ', id, ' is a ', class(modid), '\n\n',
              file = paste('htp', dvs[dv], ivsl[iv], 'log', sep='.'),
              append = TRUE)

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
            # The following only makes sense if we use palytic$gamlss and
            # REML happens to be changed to ML; if this is restored, t1 in the
            # if()... needs to be changet to modid$call$formula
            #has_method <- grepl("method", modid$call$formula)
            has_method <- grepl("method", t1$formula)
            if(any(has_method))
            {
              frm_txt <- as.character(t1$formula)[which(has_method)]
              frm_txt <- unlist( strsplit(frm_txt, ',') )
              frm_txt <- frm_txt[which( grepl("method", frm_txt) )]
              err_id['re_estimator'] <- gsub("method|\"| |=", "", frm_txt)
            }
            if(!any(has_method)) err_id['re_estimator'] <- "Cannot be determined"
            err_id['gamlss_estimator'] <- as.character( modid$method )
            err_id['analyzed_N'] <- paste(modid$N, 'cases were analyzed.')
          }
        }
        if(package=="nlme")
        {
          modid <- t1$lme( useObs )
          if(! "lme"  %in%  class(modid) )
          {
            err_id['converge'] <- modid
            err_id['estimator'] <- NA
            err_id['analyzed_N'] <- NA
            # this is a placeholder for getting error messsages from lme
            # which needs to be updated in the error handling for
            # lme in PalyticHTP
          }
          if(  "lme"  %in%  class(modid) )
          {
            err_id['converge'] <- paste('Convergence is `', 'TRUE`', sep='')
            err_id['estimator'] <- modid$method
            err_id['analyzed_N'] <- paste(modid$dims$N, 'cases were analyzed.')
          }
        }
        # return messages and model
        if(any(c("gamlss", "lme") %in% class(modid))) Model = modid
        else Model = "NA"
        IDout[[paste(ids,id,sep='')]] <- list(Messages=err_id, Model=Model )
      }

      #parallel::stopCluster(cl)

      # dis aggregate messages and models
      IDmsg <- lapply( IDout, function(x) x$Messages )
      IDout <- lapply( IDout, function(x) x$Model )

      # restructure messages into data.frame
      IDmsg <- plyr::ldply(lapply(IDmsg, unlist), rbind)
      #IDmsg <- do.call(rbind, IDmsg)

      # post-estimation aggregation
      names(IDout) <- paste(ids, uids, sep=".")
      # find a way to supress console printing for gamlss
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

      rcm <- function(x)
      {
        if(is.matrix(x))
        {
          xo <- data.frame( t(as.vector(t(x)) ) )
          row.names(xo) <- NULL
        }
        if(!is.matrix(x))
        {
          xo <- data.frame(NA)
        }
        return( xo )
      }

      # this function is fragile, we need it to adapt to any number of fixed ivs
      # and use a common label for each ivsl, potentially as a factor!
      rcn <- function(x)
      {
        if(is.matrix(x))
        {
          row.names(x)[which( row.names(x) == ivsl[[iv]] )] <- 'TargetPredictor'
          xn <- apply(expand.grid(colnames(x), row.names(x)), 1,
                      function(x) paste(x[2], x[1], sep=" "))
        }
        if(!is.matrix(x))
        {
          xn <- data.frame(NA)
        }
        return( xn )
      }

      IDoutSumm <- lapply(IDoutSum, rcm)
      IDoutSumn <- lapply(IDoutSum, rcn)
      IDoutSumm <- plyr::rbind.fill(IDoutSumm)
      firstNonMissResult <- which(!unlist(lapply(lapply(IDoutSumn, is.na), any)))[1]
      # need an escape here for all missing
      if(length(firstNonMissResult)<1) firstNonMissResult <- 1
      names(IDoutSumm) <- unlist( IDoutSumn[firstNonMissResult] )

      # populate IVout - consider writing to file here
      cln <- function(x)
      {
        gsub("\n|\\s\"|\'", "", paste(format(x), collapse=" "))
      }
      cl2 <- function(x)
      {
        if(is.null(x)) return("NULL")
        else return(x)
      }
      IDoutSumm <- data.frame(#ids        = row.names(IDmsg),
                              "ivs"       = cl2(paste(t1$ivs, collapse=', ')),
                             "target_iv" = cl2(unlist(ivsl)[iv]),
                             fixed       = cln(t1$fixed),
                             random      = cln(t1$random),
                             correlation = cln(t1$corStructs[iv,2]),
                             formula     = cln(t1$formula),
                             directory   = normalizePath(getwd()),
                             date        = Sys.time(),
                             IDmsg,
                             IDoutSumm)
      #names(IDoutSumm)[1] <- ids
      IVout[[iv]] <- IDoutSumm
    }
    ncols <- unlist( lapply(IVout, function(x) dim(x)[2]) )
    DVout[[dvs[[dv]]]] <- data.frame(dv=dvs[[dv]], plyr::rbind.fill(IVout))
  }
  # final post-processing
  outmat <- plyr::rbind.fill(DVout)
  row.names(outmat) <- NULL
  return( outmat )

}

