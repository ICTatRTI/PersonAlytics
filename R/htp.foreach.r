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

  ######################################
  # outer loop for dependent variables
  ######################################
  DVout <- list()
  for(dv in dims$DV)
  {
    message("\n\nPersonAlytics: Starting analyses for ", dvs[[dv]], " \n\n")

    # initialize the null object to be copied (do NOT use $clone() as it
    # has failed in tests) ---- note that there are inheritance issues and
    # t0 can be modified by its child t1, hence the save* objects
    t0 <- Palytic$new(data=data,
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
    saveMethod  <- t0$method
    saveFormula <- t0$formula

    if(dims$ID[1]!="All Cases")
    {
      if(detectTO) t0$getTime_Power(maxOrder, whichIC[1])
      #t0$time_powers

      if( detectAR) t0$getAR_order(PQ[1], PQ[2], whichIC[1])
      if(!detectAR) t0$corStructs <- data.frame(ids=dims$ID,
                                     arma=rep(t0$correlation, length(dims$ID)))
      #t0$corStructs

    }
    if(dims$ID[1]=="All Cases")
    {
      if(detectTO) t0$GroupTime_Power(maxOrder, whichIC[1])
      if(detectAR) t0$GroupAR_order(PQ[1], PQ[2], whichIC[1])
    }

    # parralelization
    funcs <- c("mean") # c(".eds") -- not importing from PersonAlytic correctly
    cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
    snow::clusterExport(cl, funcs)
    doSNOW::registerDoSNOW(cl)
    pkgs  <- c("gamlss", "nlme")
    pb <- txtProgressBar(max = length(dims$IV), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    ##############################
    # middle loop for covariates
    ##############################
    message("\n\nPersonAlytics: Starting target_ivs loops \n\n")
    IVout <- foreach(iv = dims$IV, .packages = pkgs,
                     .options.snow = opts) %dopar%
    {
      t1 <- t0
      t1$method <- "REML"
      ivs.temp <- unlist(c(ivs, target_ivs[iv]))
      if( is.null(ivs.temp) )
      {
        t1$ivs <- ivs.temp
      }
      if(!is.null(ivs.temp) )
      {
        t1$ivs <- gsub(" ", "", ivs.temp)
      }

      ##############################
      # inner loop for participants
      ##############################
      IDout <- list()
      for(id in dims$ID)
      {
        print(id)
        if(debugforeach){
          cat('starting id = ', id, '\n\n',
              file = paste('htp', dvs[dv], target_ivs[iv], 'log', sep='.'),
              append = TRUE)
        }

        # select rows
        if(dims$ID[1]!="All Cases")
        {
          useObs <- t1$data[[t1$ids]]==id
          wid    <- which(dims$ID==id)
        }
        if(dims$ID[1]=="All Cases")
        {
          useObs <- rep(TRUE, nrow(t1$data))
          wid    <- 1:nrow(t1$data)
        }

        ### error handling
        err_id <- list()

        # temporary data for checking valid cases
        temp   <- t1$data[useObs,all.vars(t1$formula)]
        nrt    <- nrow(temp)

        # identify rows
        err_id[ids] <- id

        # identify inputs
        err_id['ids']          <- ids
        err_id['ids']          <- t1$ids
        err_id['dvs']          <- t1$dvs
        err_id['family']       <- t1$family[[1]][2]
        err_id['package']      <- t1$package
        err_id['time']         <- t1$time
        err_id['phase']        <- t1$phase
        err_id['ivs']          <- toString( t1$ivs          )
        #err_id['target_ivs']         <- toString( t1$target_ivs         )
        err_id['interactions'] <- toString( t1$interactions )
        err_id['time_power']   <- t1$time_power

        # identify package
        err_id['package'] <- t1$package

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

        # some clean up
        rm(temp)

        # update correlation for individuals if individual models are being fit
        if(dims$ID[1]!="All Cases")
        {
          if(is.null(t1$corStructs[wid,2]) | length(t1$corStructs[wid,2]) == 0 |
             t1$corStructs[wid,2] == "NULL" | is.na(t1$corStructs[wid,2]) )
          {
            t1$correlation <- "NULL"
          }
          if(! is.null(t1$corStructs[wid,2]) )
          {
            t1$correlation <- as.character( t1$corStructs[wid,2] )
          }

          if(debugforeach)
          {
            cat('correllation for id = ', id, ' is ', t1$correlation, '\n\n',
              file = paste('htp', dvs[dv], target_ivs[iv], 'log', sep='.'),
              append = TRUE)
          }
        }


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
            err_id['converge'] <- modid
            err_id['estimator'] <- t1$method
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

        #############################
        # return messages and model
        #############################
        Model = "NA"
        if(any(c("gamlss", "lme") %in% class(modid))) Model = modid
        IDout[[paste(ids, id, sep='.')]] <- list( Messages=err_id, Model=Model )

        # reclaim inputs in case there are inheritance issues due to $gamlss or
        # $lme dropping terms
        t0$method  <- saveMethod
        t0$formula <- saveFormula

      } # end of for IDout

      # dis aggregate messages and models
      IDmsg <- lapply( IDout, function(x) x$Messages )
      IDout <- lapply( IDout, function(x) x$Model )

      # restructure messages into data.frame
      #IDmsg <- plyr::ldply(lapply(IDmsg, unlist), rbind)
      IDmsg <- do.call(rbind, IDmsg)

      # post-estimation aggregation
      if(dims$ID[1]!="All Cases") names(IDout) <- paste(ids, uids, sep=".")
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

      # function to prep row names
      # this function is fragile, we need it to adapt to any number of fixed ivs
      # and use a common label for each target_ivs
      rcn <- function(x)
      {
        if(is.matrix(x))
        {
          # extract the name(s) of the rows corresponding to the targe predictor
          Tnms <- row.names(x)[ which( grepl(target_ivs[[iv]], row.names(x)) ) ]

          # if target predictor ! factor, set the name
          if( length(Tnms) == 1 )
          {
            row.names(x)[ which(row.names(x) %in% Tnms) ] <- 'Targ'
          }

          # if target predicor is a factor variable, get the reference category
          # and set the names
          if( is.factor(t1$data[[target_ivs[[iv]]]]) )
          {
            uT     <- levels(t1$data[[target_ivs[[iv]]]])
            comCat <- unlist( strsplit(Tnms, target_ivs[[iv]]) )
            comCat <- comCat[comCat != ""]
            refCat <- uT[! uT %in% comCat]
            NewNms <- paste('Targ', comCat, 'vs', refCat, sep='_')
            row.names(x)[ which(row.names(x) %in% Tnms) ] <- NewNms
            rm(uT, comCat, refCat, NewNms)
          }
          rm(Tnms)

          xn <- apply(expand.grid(colnames(x), row.names(x)), 1,
                      function(x) paste(x[2], x[1], sep=" "))
        }
        if(!is.matrix(x))
        {
          xn <- data.frame(NA)
        }
        return( xn )
      }

      # function to force values to data.frame vectors for later stacking
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

      # turn the numeric output into vectors & rbind them together
      IDoutSumm <- lapply(IDoutSum, rcm)
      IDoutSumm <- plyr::rbind.fill(IDoutSumm)

      # functions to clean inputs, can cln be replaced with Reduce(paste, x)?
      cln <- function(x)
      {
        gsub("\n|\\s\"|\'", "", paste(format(x), collapse=" "))
      }
      cl2 <- function(x)
      {
        if(is.null(x)) return("NULL")
        else return(x)
      }

      # populate IVout - consider writing to file here
      IDoutSumm <- data.frame("target_iv" = cl2(unlist(target_ivs)[iv]),
                              fixed       = cln(t1$fixed),
                              random      = cln(t1$random),
                              correlation = ifelse(all(dims$ID=="All Cases"),
                                                   cln(t1$correlation),
                                                   cln(t1$corStructs[id,2])),
                              formula     = cln(t1$formula),
                              directory   = normalizePath(getwd()),
                              date        = Sys.time(),
                              IDmsg,
                              IDoutSumm)
      ###
      return( IDoutSumm )
      rm( t1 )
    } # end of for each IVout

    # populate DVout
    ncols <- unlist( lapply(IVout, function(x) dim(x)[2]) )
    DVout[[dvs[[dv]]]] <- data.frame(dv=dvs[[dv]], plyr::rbind.fill(IVout))
  }
  # stop the cluster
  parallel::stopCluster(cl)

  # final post-processing
  outmat <- plyr::rbind.fill(DVout)
  row.names(outmat) <- NULL
  return( outmat )
}

