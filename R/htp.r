# high throughput functions ####

#' high throughput
#'
#' @keywords internal

# the current state of using Palytic is to create one object for
# loops across individuals, but overwrite the the Palytic object
# for loops across dvs/ivs
htp <- function(data                                                   ,
                dims                                                   ,
                dvs                                                    ,
                phase                                                  ,
                ids                                                    ,
                uids                                                   ,
                time                                                   ,
                ivs                                                    ,
                target_ivs                                             ,
                interactions=NULL                                      ,
                time_power=1                                           ,
                alignPhase='none'                                      ,
                correlation=NULL                                       ,
                family = gamlss.dist::NO()                             ,
                standardize = list(dvs=FALSE,ivs=FALSE,byids=FALSE)    ,
                fpc = FALSE                                            ,
                popsize2 = 0                                           ,
                package='gamlss'                                       ,
                autoSelect = list(AR=list(P=3, Q=3)    ,
                               TO=list(polyMax=3)      ,
                               DIST=list())            ,
                whichIC = "BIC"                                        ,
                sigma.formula=~1                                       ,
                debugforeach = FALSE                                   ,
                userFormula = list(
                fixed=NULL,
                random=NULL,
                formula=NULL)                                          ,
                cores=parallel::detectCores()-1                        )
{

  PQ <- c(autoSelect$AR$P, autoSelect$AR$Q)

  #----------------------------------------------------------------------------#
  # log files - create a log directory and overwrite all logs
  #----------------------------------------------------------------------------#
  if(debugforeach)
  {
    if(!file.exists('PAlogs')) dir.create('PAlogs')
    cat( correlation, '\n\n', file = './PAlogs/getARnEQ1run.log', append=FALSE)
    cat( 'Start $lme.log\n\n', file = "./PAlogs/$lme.log", append=FALSE)
    cat( 'Start formula.log\n\n', file = "./PAlogs/formula.log", append=FALSE)
  }

  #----------------------------------------------------------------------------#
  # determine which loop type to use, if dvLoop==TRUE, a non-parallelized
  # loop for the DVs is used
  #----------------------------------------------------------------------------#
  dvLoop  <- TRUE
  if( length(dvs) > 1 & length(target_ivs) <= 1 & dims$ID[1]=="All Cases" )
  {
    dvLoop <- FALSE
  }

  #----------------------------------------------------------------------------#
  # autoSelect$TO not implemented for alignPhase == "piecewise"
  #----------------------------------------------------------------------------#
  if(alignPhase == "piecewise" & !is.null(autoSelect$TO))
  {
    autoSelect$TO <- NULL
    message("\n`alignPhase=='piecewise'`, hence `autoSelect$TO` has been set to NULL",
            "\nbecause automatic decection of the time order is not implemented",
            "within phases.")
  }

  #----------------------------------------------------------------------------#
  # functions to export
  #----------------------------------------------------------------------------#
  exports <- c()

  # parralelization options could be implemented as methods for a generic,
  # would that be faster? probably not, and since the user con't touch them,
  # the generic doesn't help anyone
  #----------------------------------------------------------------------------#
  # parralelization option 1: iterate over dvs only ####
  #----------------------------------------------------------------------------#
  if(!dvLoop )
  {
    DIM <- expand.grid(ID=dims$ID, IV=dims$IV)
    if(is.factor(DIM$ID)) DIM$ID <- as.character(DIM$ID)
    capture.output( pb <- txtProgressBar(max = length(dvs), style = 3),
                    file='NUL')
    progress <- function(n) setTxtProgressBar(pb, n)
    opts     <- list(progress = progress)

    # parralelization setup
    pkgs  <- c("gamlss", "nlme", "foreach")
    cl <- snow::makeCluster(cores, type="SOCK", outfile="")
    snow::clusterExport(cl, list())
    doSNOW::registerDoSNOW(cl)

    start <- messenger(dvLoop)
    DVout <- list()
    IDout <- foreach( dv=seq_along(dvs), .packages = pkgs,
                      .options.snow = opts) %dopar%
    {
      # pre-clean loop
      if(exists("t0")) rm(t0)
      # set up the  palytic object
      t0 <- PersonAlytics::Palytic$new(data=data  ,
                        ids=ids                   ,
                        dv=dvs[[dv]]              ,
                        time=time$raw             ,
                        phase=phase               ,
                        ivs=ivs                   ,
                        interactions=interactions ,
                        standardize=standardize   ,
                        autoSelect=autoSelect     ,
                        time_power=time_power     ,
                        alignPhase=alignPhase     ,
                        correlation=correlation   ,
                        family=family             ,
                        method="ML"

      )

      # autoselect
      t0$detect(subgroup    = NULL        ,
                model       = NULL        ,
                parallel    = "no"        ,
                plot        = FALSE       ,
                userFormula = userFormula ,
                dims        = dims        ,
                package     = package     )
      .htp(t0, id=1, iv=1, dv, dvs, ivs,
           dims, package, target_ivs, PQ, family, fpc, popsize2, debugforeach)
    }# end of foreach
    # stop the cluster
    parallel::stopCluster(cl)

    # the next several lines are repeated in if(dvLoop), could be moved to
    # function, but unpacking the results will take as much code a repeating
    # the lines
    # disaggregate messages
    IDmsg <- lapply( IDout, function(x) data.frame(x$Messages) )
    IDmsg <- plyr::rbind.fill(IDmsg)

    # disaggregate depcriptive statistics
    IDdesc <- lapply( IDout, function(x) data.frame(x$Describe))
    IDdesc <- plyr::rbind.fill(IDdesc)

    # parameter estimates
    if(dims$ID[1]!="All Cases") names(IDout) <- paste(ids, uids, sep=".")
    IDoutSum <- lapply( IDout, function(x) data.frame(x$IDoutSum) )
    IDoutSumm <- plyr::rbind.fill(IDoutSum)

    # fpc
    if(fpc & package == "nlme")
    {
      IDoutSumFPC <- lapply( IDout, function(x) data.frame(x$IDoutSumFPC) )
      IFoutSummFPC <- plyr::rbind.fill(IDoutSumFPC)
    }

    if(!fpc) DVout[[1]] <- list(IDmsg=IDmsg, IDdesc=IDdesc, IDoutSumm=IDoutSumm)
    if(fpc & package == "nlme")
    {
      DVout[[1]] <- list(IDmsg=IDmsg, IDdesc=IDdesc, IDoutSumm=IDoutSumm,
                         IDoutSummFPC=IFoutSummFPC)
    }

    message('\n\nModel fitting took:\n ', capture.output(Sys.time() - start), ".\n\n")

    rm(IDoutSumm, IDmsg, IDdesc )
  }

  #----------------------------------------------------------------------------#
  # parralelization option 2: outer loop is DV ####
  # so that correlation and time order searches only need to happen 1x/DV
  #----------------------------------------------------------------------------#
  if( dvLoop )
  {
    DIM <- expand.grid(ID=dims$ID, IV=dims$IV)
    if(is.factor(DIM$ID)) DIM$ID <- as.character(DIM$ID)
    capture.output( pb <- txtProgressBar(max = nrow(DIM), style = 3), file='NUL')
    progress <- function(n) setTxtProgressBar(pb, n)
    opts     <- list(progress = progress)

    DVout <- list()
    for(dv in dims$DV)
    {
      # set up the parent palytic object
      t0 <- Palytic$new(data         = data         ,
                        ids          = ids          ,
                        dv           = dvs[[dv]]    ,
                        time         = time$raw     ,
                        phase        = phase        ,
                        ivs          = ivs          , # target_ivs added later
                        interactions = interactions ,
                        standardize  = standardize  ,
                        autoSelect   = autoSelect   ,
                        time_power   = time_power   ,
                        alignPhase   = alignPhase   ,
                        correlation  = correlation  ,
                        family       = family       ,
                        method       = "ML"         # requested method used later
      )

      # autoselect
      t0$detect(subgroup    = NULL        ,
                model       = NULL        ,
                parallel    = "snow"      ,
                plot        = FALSE       ,
                userFormula = userFormula ,
                dims        = dims        ,
                package     = package     )

      # parralelization setup -- must reoccur for each dv in dims$DV
      pkgs  <- c("gamlss", "nlme", "foreach")
      cl <- snow::makeCluster(cores, type="SOCK", outfile="")
      snow::clusterExport(cl, list())
      doSNOW::registerDoSNOW(cl)

      start <- messenger(dvLoop, dvs, dv, dims, ids, target_ivs)
      IDout <- foreach( id=DIM$ID, iv=DIM$IV, .export = exports,
                        .packages = pkgs, .options.snow = opts) %dopar%
      {
        #i=1; id<-DIM$ID[i]; iv<-DIM$IV[i]
        t1 <- t0$clone(deep=TRUE)
        .htp(t1, id, iv, dv, dvs, ivs,
             dims, package, target_ivs, PQ, family, fpc, popsize2, debugforeach)

      } # end of foreach
      # stop the cluster
      parallel::stopCluster(cl)

      message('\n\nModel fitting for the dependent variable `', dvs[dv],
              '` took:\n ', capture.output(Sys.time() - start), ".\n\n")

      cat('\n\n') #TODO is this needed?
    }

    # disaggregate messages
    IDmsg <- lapply( IDout, function(x) data.frame(x$Messages) )
    IDmsg <- plyr::rbind.fill(IDmsg)

    # disaggregate depcriptive statistics
    IDdesc <- lapply( IDout, function(x) data.frame(x$Describe))
    IDdesc <- plyr::rbind.fill(IDdesc)

    # parameter estimates
    if(dims$ID[1]!="All Cases") names(IDout) <- paste(ids, uids, sep=".")
    IDoutSum <- lapply( IDout, function(x) data.frame(x$IDoutSum))
    IDoutSumm <- plyr::rbind.fill(IDoutSum)

    # fpc
    if(fpc & package == "nlme")
    {
      IDoutSumFPC <- lapply( IDout, function(x) data.frame(x$IDoutSumFPC) )
      IFoutSummFPC <- plyr::rbind.fill(IDoutSumFPC)
    }

    # put outputs in a list, to be aggregated outside of the dv loop
    if(!fpc) DVout[[dvs[[dv]]]] <- list(IDmsg=IDmsg, IDdesc=IDdesc,
                                        IDoutSumm=IDoutSumm)
    if(fpc & package == "nlme")
    {
      DVout[[dvs[[dv]]]] <- list(IDmsg=IDmsg, IDdesc=IDdesc,
                                 IDoutSumm=IDoutSumm,
                                 IDoutSummFPC=IFoutSummFPC)
    }

    rm(t0, IDoutSumm, IDmsg, IDdesc )

  } # end of dv loops

  # final post-processing
  IDmsg     <- plyr::rbind.fill( lapply(DVout, function(x) x$IDmsg) )
  IDdesc    <- plyr::rbind.fill( lapply(DVout, function(x) x$IDdesc) )
  IDoutSumm <- plyr::rbind.fill( lapply(DVout, function(x) x$IDoutSumm) )

  if(fpc & package == "nlme")
  {
    IDoutSummFPC <- plyr::rbind.fill( lapply(DVout, function(x) x$IDoutSummFPC) )
    names(IDoutSummFPC) <- paste(names(IDoutSummFPC), 'FPC', sep='.')
    outmat <- cbind(IDmsg, IDdesc, IDoutSumm, IDoutSummFPC)
  }
  if(!fpc) outmat <- cbind(IDmsg, IDdesc, IDoutSumm)

  row.names(outmat) <- NULL

  # drop redundant columns that may crop up under some combinations of settings
  outmat <- outmat[,! names(outmat) %in% c('converge.1',	'estimator.1',
              'analyzed_N.1',	'call.1',	'targ_ivs_lrt_pvalue.1')]
  stats <- which( grepl('statName', names(outmat)) | grepl('statValue', names(outmat)) )
  nostats <- stats[ unlist(lapply(outmat[,stats], function(x) all(is.na(x)))) ]
  if(length(nostats)>0) outmat <- outmat[,-nostats]

  return( outmat )
}

#' messenger
#' @keywords internal
messenger <- function(dvLoop, dvs=NULL, dv=NULL,
                      dims=NULL, ids=NULL, target_ivs=NULL)
{
  if(!dvLoop)
  {
    message('\n\nModel fitting starting...\n\n')
  }
  if(dvLoop)
  {
    message('\n\nFitting models of the dependent variable `', dvs[dv],
            '` for ',
            ifelse( dims$ID[1]=="All Cases",
                    paste(length(dims$indID), ' cases in `', ids, '`', sep=''),
                    paste("All indidivuals in `", ids, '`', sep='')
            ),
            ifelse( length(target_ivs)>0,
                    paste('\n and for', length(target_ivs),
                          'target independent variables in `target_ivs`.\n'),
                    ".")
    )
  }
  start <- Sys.time()
  return(start)
}

#' .htp
#' @keywords internal
.htp <- function(t1, id, iv, dv, dvs, ivs, dims,
            package, target_ivs, PQ, family, fpc, popsize2, debugforeach)
{
  #-------------------------------------------------------------------------
  # save information needed to help debug
  #-------------------------------------------------------------------------
  if(debugforeach)
  {
    cat(Reduce(paste, deparse(t1$time)),
        '\n\n', file = "./PAlogs/formula.log", append = TRUE)

    cat(Reduce(paste, deparse(t1$formula)),
        '\n\n', file = "./PAlogs/formula.log", append = TRUE)
  }

  #-------------------------------------------------------------------------
  # initialize the model as NA
  #-------------------------------------------------------------------------
  modid <- NA

  #-------------------------------------------------------------------------
  # for the current id, select potential rows, useObs will be updated
  # based on missing data determined by htpErrors()
  #-------------------------------------------------------------------------
  if(dims$ID[1]!="All Cases")
  {
    useObs <- t1$data[[t1$ids]]==id
    wid    <- which(dims$ID==id)
    errid  <- wid
  }
  if(dims$ID[1]=="All Cases")
  {
    useObs <- rep(TRUE, nrow(t1$data))
    wid    <- 1:nrow(t1$data)
    errid  <- dims$ID[1]
  }

  #-------------------------------------------------------------------------
  # accumulate inputs and errors for the output, results are used in
  # row selection `useObs`
  #-------------------------------------------------------------------------
  htpErr <- htpErrors(t1=t1, id=errid, dv=dvs[[dv]], dims=dims,
                      package=package, useObs=useObs,
                      target_iv=target_ivs[[iv]])
  tivv   <- htpErr$tivv
  dvVar  <- htpErr$dvVar
  err_id <- htpErr$err_id
  rm(htpErr)

  #-------------------------------------------------------------------------
  # populate time_power
  #-------------------------------------------------------------------------
  if(dims$ID[1]!="All Cases")
  {
    t1$time_power     <- as.numeric( t1$time_powers[wid,2] )
    err_id$time_power <- t1$time_powers[wid,2]
  }
  if(dims$ID[1]=="All Cases")
  {
    err_id$time_power <- t1$time_power
  }

  #-------------------------------------------------------------------------
  # for the current id, estimate a full model with the current target IV
  #-------------------------------------------------------------------------
  if( length( target_ivs[[iv]] ) > 0 & !is.na(tivv) & tivv & dvVar>0 )
  {
    # add the target IV
    ivs.temp <- unlist(c(ivs, target_ivs[[iv]]))
    if( is.null(ivs.temp) )
    {
      t1$ivs <- ivs.temp
    }
    if(!is.null(ivs.temp) )
    {
      t1$ivs <- gsub(" ", "", ivs.temp)
    }

    #TODO(Stephen): override correlation search for ARMA?
    fitOutput <- fitWithTargetIV(t1, package, useObs, dims,
                                 dropVars=target_ivs[[iv]], PQ,
                                 fpc=fpc, popsize2=popsize2)
    err_id <- c(err_id, fitOutput$err_id)
    modid  <- fitOutput$modid
    rm(fitOutput)

    #cat(modid, "\n\n", file='FitWithTargetIV.txt', append=TRUE)

    #test <- t1$lme(useObs, target_ivs[[iv]], PQ)
    #cat(class(test), '\nuseObs: ', toString(table(useObs)),
    #    '\ntarget_iv: ', target_ivs[[iv]], '\nPQ: ', toString(PQ),
    #    "\n\n", file='FitWithTargetIVtest.txt', append=TRUE)
  }
  # the target iv variance was 0
  if( is.na(tivv) ) tivv <- FALSE
  if( length( unlist(target_ivs[iv]) ) != 0 & !tivv )
  {
    #modid <- NA
    err_id$converge    <- paste('No variance in `', target_ivs[iv], '`.')
    err_id$estimator   <- toString( NA )
    err_id$analyzed_N  <- "0 cases were analyzed."
    err_id$call        <- toString( NA )
    err_id$targ_ivs_lrt_pvalue <- as.numeric( NA )
  }

  #-------------------------------------------------------------------------
  # if there was no target IV, fit the null model
  #-------------------------------------------------------------------------
  if( length( target_ivs[[iv]] ) == 0 | !tivv )
  {
    mod1   <- fitWithTargetIV(t1, package, useObs, dims,
                              dropVars=NULL, PQ=PQ,
                              fpc=fpc, popsize2=popsize2)
    modid  <- mod1$modid
    err_id <- c(err_id, mod1$err_id)
    rm(mod1)
  }

  #-------------------------------------------------------------------------
  # if the dv variance was 0
  #-------------------------------------------------------------------------
  if( is.na(dvVar) ) dvVar <- 0
  if( dvVar==0 )
  {
    #modid <- NA
    err_id$converge    <- paste('No variance in `', dvs[dv], '`.')
    err_id$estimator   <- toString( NA )
    err_id$analyzed_N  <- "0 cases were analyzed."
    err_id$call        <- toString( NA )
    err_id$targ_ivs_lrt_pvalue <- as.numeric( NA )
  }

  #-------------------------------------------------------------------------
  # descriptive statistics
  #-------------------------------------------------------------------------
  descr_id <- t1$describe(useObs)

  #-------------------------------------------------------------------------
  # re-fit models with REML (unless arma)
  # TODO(Stephen): this *shouldn't* need the errror accumulation from
  # fitWithTargetIV unless the REML fit leads to a different model than
  # the ML fit
  #-------------------------------------------------------------------------
  Model <- data.frame(NA)
  if(any(c("gamlss", "lme") %in% class(modid)))
  {
    try( cat("For id ", id, " and iv ", iv, " htp tried to REML lme\n\n",
             file='REMLlme.txt', append=TRUE), silent = TRUE )
    t1$method <- "REML"
    t1$family <- family #TODO(Stephen): prior line drops family, why??
    if("gamlss" %in% class(modid)) Model <- t1$gamlss( useObs, sigma.formula=sigma.formula )
    if("lme"    %in% class(modid)) Model <- t1$lme( useObs,
                                      fpc=fpc, popsize2 = popsize2)
    if( any(c("gamlss", "lme") %in% class(Model)) )
    {
      err_id$method <- err_id$estimator <- "REML"
    }
    if( any("Model did not converge" %in% Model) )
    {
      Model <- modid
    }
  }
  if( !any(c("gamlss", "lme") %in% class(modid)) &
      !any("Model did not converge" %in% modid)  &
      !any(is.na(modid)) )
  {
    if( any( c("ARIMA", "Arima") %in% class(modid$arima) ) )
    {
      try( cat("For id ", id, " and iv ", iv, "I made Model<-modid\n\n",
               file='NotREMLarma.txt', append=TRUE), silent = TRUE )
      Model <- modid
    }
  }
  rm(modid)

  #-------------------------------------------------------------------------
  # add final entries to err_id, these may depend on final results
  #-------------------------------------------------------------------------
  err_id <- htpForms(err_id, t1, dims, id, package, modid=Model)

  #-------------------------------------------------------------------------
  # reduce the size of Model
  #-------------------------------------------------------------------------
  IDoutSum <- getParameters(Model, package, target_ivs[[iv]], t1$datac,
                            fpc=FALSE)
  if(fpc & package == "nlme")
  {
    IDoutSumFPC <- getParameters(Model, package, target_ivs[[iv]], t1$datac,
                                fpc=fpc)
  }
  rm(Model)

  #-------------------------------------------------------------------------
  # return to foreach
  #-------------------------------------------------------------------------
  # this line stays commented  out except for testing
  #IDout[[i]] <- list( Messages=err_id, IDoutSum=IDoutSum, Describe=descr_id)
  if( !fpc | package != "nlme" )
  {
    return( list( Messages=err_id, IDoutSum=IDoutSum, Describe=descr_id ) )
  }
  if(fpc & package == "nlme")
  {
    return( list( Messages=err_id, IDoutSum=IDoutSum, Describe=descr_id,
                  IDoutSumFPC=IDoutSumFPC) )
  }
}

#' getParameters: extract parameter table, the `Model` variable is too big
#'
#' @param Model Statistical model output from \code{\link{Palytic}} methods
#' \code{$lme}, \code{$gamlss}, or \code{$arma}
#' @param package See \code{package} in \code{\link{PersonAlytic}}
#' @param target_ivs See \code{targe_ivs} in \code{\link{PersonAlytic}}
#' @param data A data.frame, for now only pass t1$datac in \code{\link{htp}}
#'
#' @keywords internal
getParameters <- function(Model, package, target_iv, data, fpc)
{
  IDoutSum <- data.frame(NA)
  if(!all(is.na(Model)))
  {
    if(package=='gamlss')
    {
      if(! "gamlss" %in% class(Model)) IDoutSum <- NA
      if(  "gamlss" %in% class(Model))
      {
        IDoutSum <- rcm( summary(Model),
                         target_iv, data)
      }
    }

    if(package=='nlme')
    {
	    if(! "lme" %in% class(Model)) IDoutSum <- NA
      if(  "lme" %in% class(Model))
      {
        if(!fpc)
        {
          IDoutSum <- rcm(Model$tTable, target_iv, data)
        }

        if(fpc)
        {
          IDoutSum <- rcm(Model$FPCtTable, target_iv, data)
        }
      }
    }

    if(package=="arma")
    {
      if(! "coeftest"  %in%  class(Model$tTable)) IDoutSum <- NA
      if(  "coeftest"  %in%  class(Model$tTable))
      {
        IDoutSum <- rcm(getARMAtbl(Model), target_iv, data)
      }
    }
  }

  return(IDoutSum)
}


#' htpForms: accumulate formula information
#' @keywords internal
htpForms <- function(err_id, t1, dims, id, package, modid)
{

  err_id["fixed"]        <- toString( NA )
  err_id["random"]       <- toString( NA )
  err_id["formula"]      <- toString( NA )
  err_id["correlation0"] <- toString( NA )

  correlation0 <- t1$correlation0
  if(is.null(correlation0)) correlation0 <- "NULL"

  if(!any(is.na(modid)) & !any(modid == "Model did not converge"))
  {
    err_id["fixed"]       <- rmSpecChar(modid$PalyticSummary$fixed)
    err_id["random"]      <- rmSpecChar(modid$PalyticSummary$random)
    err_id["correlation0"] <- ifelse(all(dims$ID=="All Cases"),
                                   #rmSpecChar(modid$PalyticSummary$correlation),
                                   correlation0,
                                   ifelse(package=="arma",
                                          #TODO() this doesn't get the ride order
                                          #gerARIMAorder(modid$arima),
                                          err_id$correlation, # placeholder
                                          rmSpecChar(t1$corStructs[id,2])))
    err_id["formula"]     <- rmSpecChar(modid$PalyticSummary$formula)
  }

  err_id["directory"]   <- normalizePath(getwd())
  err_id["date"]        <- toString( Sys.time() )

  return(err_id)
}

#' htpErrors: accumulate inputs and errors
#' @keywords internal
htpErrors <- function(t1, id, dv, dims, package, useObs, target_iv)
{
  err_id <- list()

  # temporary data for checking valid cases
  temp   <- na.omit( t1$datac[useObs, c(all.vars(t1$formula), target_iv)] )

  # identify rows
  err_id[[t1$ids]] <- id

  # identify inputs
  err_id['ids']          <- t1$ids
  err_id['dv']           <- dv #TODO(Stephen) column not appearing in output??
  err_id['time']         <- t1$time[1]
  err_id['phase']        <- t1$phase
  err_id['ivs']          <- toString( t1$ivs )
  err_id['target_iv']    <- toString( target_iv )
  err_id['interactions'] <- toString( t1$interactions )
  err_id['time_power']   <- t1$time_power # updated later
  err_id['alignPhase']   <- toString( t1$alignPhase )
  if(package=="arma") err_id['correlation']  <- "See 'call' column'"
  if(package!="arma") err_id['correlation']  <- t1$correlation
  err_id['family']        <- t1$family[[1]][2]
  err_id['gamlss.family'] <- t1$family[[1]][1]
  err_id['standardize']   <- paste(paste(names(t1$standardize),
                                        t1$standardize, sep='='), collapse=', ')
  err_id['method']       <- t1$method
  err_id['package']      <- package
  err_id['Personalytics']<- paste('Version', packageVersion("PersonAlytics"))
  err_id['Date_Time']    <- format(Sys.time(), format='%Y%m%d_%H.%M%p')

  # id and time series lengths
  nid <- length(unique(temp[[t1$ids]]))
  err_id['N_participants'] <- nid
  ntp <- length(unique(temp[[t1$time$raw]]))
  err_id['N_time_points'] <- ntp

  # check for adequate data
  nrt            <- length(unique(temp[[t1$time$analysis[1]]]))
  err_id['N_time_points_complete'] <- paste('There are', nrt, 'time points w/ complete data')

  # outcome variance
  dvVar <- round(var(temp[[t1$dv]], na.rm = TRUE),3)
  err_id['dvVariance'] <- paste( 'The variance of `', t1$dv, '` is ', dvVar,
                            ifelse(all(dims$ID=="All Cases"), ".",
                                   paste(' for ', t1$ids, ' = ', id, '.', sep='')),
                            ' If unexpected, check your settings in `standardize`.',
                            sep='')

  # time should be monotonically increasing
  err_id['timeVariance'] <- paste("`", t1$time[1], "` is",
                             ifelse(t1$monotone$monotonic[id], "" , "NOT"),
                             " monotonically increasing for `",
                             t1$ids, "` =", id, ".")

  # ivs and ivs variance
  #TODO() add value pass for iv variances and add to output
  err_id['ivVariance'] <- "There are no variables in `ivs`."
  if( length(t1$ivs) > 0)
  {
    err_id['ivs']       <- toString( t1$ivs )
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
    wsplit <- unlist(lapply(ivvs, grepl, pattern="\\+"))
    if( any(wsplit) )
    {
      ivvs <- unique( c(ivvs[!wsplit],
                        unlist( strsplit(unlist(ivvs[wsplit]), "\\+") ) ) )
    }
    ivvs <- unique( gsub(" ", "", unlist(ivvs)) )
    ivvs <- ivvs[!grepl("\\(", ivvs)]
    ivvs <- ivvs[!grepl("\\^", ivvs)]
    ivv  <- unlist( lapply(data.frame(temp[,ivvs]),
                           function(x) !all(duplicated(x)[-1L])) )

    err_id['ivVarariance'] <- paste( paste('The variance of `', ivvs, '` is ',
                                    ifelse(ivv, '> 0 ', '= 0 '),
                                    sep=''), collapse='; ')
  }

  # target iv and target iv variance
  err_id['target_ivVariance'] <- "No target_ivs"
  tivv <- NA

  # check for 0 variance in the target iv
  if(!is.null(target_iv))
  {
    if(is.character(temp[[target_iv]]) | is.factor(temp[[target_iv]]))
    {
      tivv <- !all(duplicated(temp[[target_iv]])[-1L])
      err_id['target_ivVariance'] <- paste('The variance of `',
                                    target_iv, '` is ',
                                    ifelse(tivv, '> 0 ', '= 0 '),
                                    sep='')
    }
    if(is.numeric(temp[[target_iv]]))
    {
      tivv <- round(var(temp[[target_iv]], na.rm=TRUE), 2)
      err_id['target_ivVariance'] <- paste('The variance of `',
                                    target_iv, '` is ', tivv)
    }

  }

  return( list(err_id=err_id, tivv=tivv, dvVar=dvVar) )
}

#' rcm: function to force values to data.frame vectors for later stacking
#'
#' @param x is a matrix-like ANOVA-like table from lme, gamlss, or coeftest in
#' arima
#' @param target_ivs See \code{targe_ivs} in \code{\link{PersonAlytic}}
#' @param data A data.frame, for now only pass t1$datac in \code{\link{htp}}
#'
#' @keywords internal
rcm <- function(x, target_ivs, data)
{
  if(is.matrix(x))
  {
    xo <- data.frame( t(as.vector(t(x)) ) )
    row.names(xo) <- NULL
    colnames(xo)  <- rcn(x, target_ivs, data)
  }
  if(!is.matrix(x))
  {
    xo <- data.frame(NA)
  }
  return( xo )
}

#' rcn: function to prep row names
#'
#' @param x See `x` for \code{\link{rcm}}
#' @param target_ivs See \code{targe_ivs} in \code{\link{PersonAlytic}}
#' @param data A data.frame, for now only pass t1$datac in \code{\link{htp}}
#'
#' @keywords internal
rcn <- function(x, target_ivs, data)
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

      # if target predictor is a factor variable, get the reference category
      # and set the names
      w <- which(lapply(Rns, length) > 0)
      if( length(w) > 0 )
      {
        #cat(toString(target_ivs), file = "htpLine431.txt")
        #cat(capture.output(target_ivs), '\n',
        #    capture.output(target_ivs[[w]]), '\n',
        #    capture.output(is.factor(data[[target_ivs[[w]]]])),
        #    file = "htpLine431.txt")
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
    return( paste('arima(p=', ao[1], ", d=0",
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
  #TOTO() what kind of NA should it return? data.frame? character?
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

#TODO(Stephen): add ... to these functions to pass to Palytic methods
#' fit models with the target iv & calculate LRT
#' @keywords internal
fitWithTargetIV <- function(t1, package, useObs, dims, dropVars=NULL, PQ=c(3,3),
                            fpc, popsize2)
{
  # fit model with targe iv
  if(package=="nlme")
  {
    modid <- fitWithTargetIVlme(t1, useObs, dims, dropVars, PQ,
                                fpc, popsize2)
  }
  if(package=="arma")
  {
    modid <- fitWithTargetIVarma(t1, useObs, dims, dropVars, PQ)
  }
  if(package=="gamlss")
  {
    modid <- fitWithTargetIVgamlss(t1, useObs, dims, dropVars)
  }

  err_id <- modid$err_id
  modid  <- modid$modid

  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVlme
#' @keywords internal
fitWithTargetIVlme <- function(t1, useObs, dims, dropVars, PQ, fpc, popsize2)
{
  err_id <- list()
  modid  <- t1$lme( useObs, dropVars, PQ,
                    fpc=fpc, popsize2 = popsize2 )
  if(! "lme"  %in%  class(modid) )
  {
    err_id['converge']   <- modid
    err_id['estimator']  <- toString( t1$method )
    err_id['analyzed_N'] <- toString( NA )
    err_id['call']       <- toString( NA )
    err_id['wasLRTrun']  <- FALSE
    err_id['targ_ivs_lrt_pvalue'] <- as.numeric( NA )
  }
  if(  "lme"  %in%  class(modid) )
  {
    err_id['converge']   <- 'Convergence is `TRUE`'
    err_id['estimator']  <- toString( modid$PalyticSummary$method )
    err_id['analyzed_N'] <- paste(modid$dims$N, 'observations were analyzed.')
    err_id['call'] <- paste(
      paste( Reduce( paste, deparse(modid$PalyticSummary$fixed) ),
             Reduce( paste, deparse(modid$PalyticSummary$random) ),
             modid$PalyticSummary$correlation,
           sep='; '), collapse = '')
    err_id['wasLRTrun']  <- modid$lrt$wasLRTrun
    err_id['targ_ivs_lrt_pvalue'] <- modid$lrt$lrtp
  }
  # here is a placeholder for getting error messsages from lme
  # which needs to be updated in the error handling for
  # lme in Palytic


  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVarma
#' @keywords internal
fitWithTargetIVarma <- function(t1, useObs, dims, dropVars, PQ)
{
  err_id <- list()
  modid <- t1$arma( useObs, max.p=PQ[1], max.q=PQ[2], dropVars )
  if(! "coeftest"  %in%  class(modid$tTable) )
  {
    err_id['converge']   <- modid$arima
    #t1$method --- ML currently default in arma, REML not an option
    err_id['estimator']  <- toString( "ML" )
    err_id['analyzed_N'] <- toString( NA )
    err_id['call']       <- toString( NA )
  }
  if(  "coeftest"  %in%  class(modid$tTable) )
  {
    err_id['converge']   <- 'Convergence is `TRUE`'
    err_id['estimator']  <- toString( "ML" )
    err_id['analyzed_N'] <- paste(modid$arima$nobs, 'observations were analyzed.')
    armaOrder <- forecast::arimaorder(modid$arima)
    err_id['call'] <- paste( "arima(y=", t1$dv,
                             ", order=c(",
                             armaOrder[1], ",",
                             armaOrder[2], ",",
                             armaOrder[3], "), ",
                             "xreg=rxeg) where `xreg` includes: ",
                             paste(modid$xregs, collapse=", "),
                             ".")
  }
  err_id['wasLRTrun']  <- modid$lrt$wasLRTrun
	err_id['targ_ivs_lrt_pvalue'] <- modid$lrt$lrtp
  return( list(err_id = err_id, modid = modid) )
}

#' fitWithTargetIVgamlss
#' @keywords internal
fitWithTargetIVgamlss <- function(t1, useObs, dims, dropVars)
{
  err_id <- list()
  modid  <- t1$gamlss( useObs, sigma.formula=sigma.formula, dropVars=dropVars )

  if(! "gamlss"  %in%  class(modid) )
  {
    err_id['converge']         <- "Model failed to converge"
    err_id['re_estimator']     <- toString( NA )
    err_id['gamlss_estimator'] <- toString( NA )
    err_id['analyzed_N']       <- toString( NA )
    err_id['call']             <- toString( NA )
    err_id['wasLRTrun']           <- FALSE
    err_id['targ_ivs_lrt_pvalue'] <- as.numeric( NA )
  }
  if(  "gamlss"  %in%  class(modid) )
  {
    err_id['converge'] <- paste('Convergence is `', modid$converged,
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

    err_id['gamlss_estimator']    <- toString( modid$PalyticSummary$method )
    err_id['analyzed_N']          <- paste(modid$N, 'observations were analyzed.')
    err_id['call']                <- rmSpecChar(modid$PalyticSummary$formula)
    err_id['wasLRTrun']           <- modid$lrt$wasLRTrun
    err_id['targ_ivs_lrt_pvalue'] <- modid$lrt$lrtp
  }

  return( list(err_id = err_id, modid = modid) )
}

