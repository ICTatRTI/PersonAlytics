#' Hight throughput methods for the \code{PalyticHTP} class generator,
#' which inherits from \code{\link{Palytic}}.
#'
#' @name PalyticHTP
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @import PersonAlyticsLite
#' @importFrom forecast auto.arima
#'
#' @export
#' @format A \code{PalyticHTP} generator
#' @keywords data
#'
#' @usage PalyticHTP(data, ids, y, time)
#'
#' @details
#' High throughput methods for the \code{\link{Palytic}} class generator in the
#' \code{\link{PersonAlyticsLite}} package and the documented here.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{getAR_order(dv, maxAR=3, maxMA=3, crit="BIC", lrt=TRUE, alpha=.05)}}{
#'   This method automates the task of determining the correlation structure for each case in
#'   \code{ids} (see \code{\link{PersonAlytic} or \code{\link{PersonAlytic}}}). \code{maxAR} and
#'   \code{maxMA} set the highest autoregressive and moving average parameters to be tested.
#'   If the time variable is approximatetly equally spaced, \code{crit} is the criterion used for
#'   determining the correlation structure for each \code{ids} using the \code{\link{auto.arima}}
#'   function. If the time variable is unequally spaced, \code{crit} as also the criterion for
#'   model selection via mixed effects models using \code{\link{lme}} if \code{lrt=FALSE}.
#'   If \code{lrt=TRUE} (the default), likelihood ratios are used via the \code{\link{anova}}
#'   method for \code{\link{lme}} objects. Calling \code{getAR_order} populates the
#'   \code{corStructs} field of a \code{PalyticHTP} object. For usage, see the examples.}
#'   \item{\code{getTime_Power(subset, maxOrder)}}{This method automates the task of determining
#'   \code{time_power} for each case in \code{ids} (see \code{\link{PersonAlytic} or
#'   \code{\link{PersonAlytic}}}). For example,
#'   if \code{getTime_Power} returns \code{time_power=3}, then \code{time + time^2 + time^3}
#'   will be added to the fixed effects of the model. Calling \code{getTime_Power} populates the
#'   \code{time_powers} field of a \code{PalyticHTP} object. For usage, see the examples.}
#' }
#'
#' @examples
#' t1 <- PalyticHTP$new(data = PersonAlyticsLite::OvaryICT, ids='Mare',
#'                   dv='follicles', time='TimeSin', phase='Phase')
#' t1$getTime_Power()
#' t1$time_powesr
#' # getAR_order works on one case at a time
#' t1$getAR_order(t1$dv)
#' t1$corStructs

PalyticHTP <- R6::R6Class("PalyticHTP",
                       inherit = PersonAlyticsLite::Palytic)

# this will only be applied to one participant at a time
# crit can take on AIC or BIC
PalyticHTP$set("public", "getAR_order",
            function(dV         ,
                     maxAR=3    ,
                     maxMA=3    ,
                     crit="BIC" ,
                     lrt=TRUE   ,
                     alpha=.05  )
            {
              # check whether time is (approximately) equally spaced, use the
              # first participant for now; the method is the standard deviation
              # of the 1st order difference score
              uid <- sort( as.numeric( unique(self$data[[self$ids]]) ) )
              tt  <- self$data[[self$time]][self$data[[self$ids]]==uid[1]]
              ttd <- diff( tt )
              eqSpace <- sd(ttd) < .1
              #
              if(eqSpace)
              {
                if(crit=="BIC") ic = "bic"
                if(crit=="AIC") ic = "aic"
                #frmToChar(self$fixed) # 20180816 not sure why this is here, depricate
                AR_orders <- by(self$data[[dV]],
                                self$data[[self$ids]],
                                FUN = function(x) forecast::auto.arima(x,
                                                  ic=ic)$arma[c(1,3)])
                AR_orders <- lapply(AR_orders, function(x)as.data.frame(t(x)))
                AR_orders <- plyr::rbind.fill(AR_orders)
                AR_orders <- data.frame(ids=as.numeric(row.names(AR_orders)),
                                        p=AR_orders[,1],
                                        q=AR_orders[,2])
                AR_orders <- AR_orders[order(AR_orders$ids),]
                row.names(AR_orders) <- NULL

                AR_orders <- data.frame(ids=AR_orders$ids,
                                        arma=paste("corARMA(p=", AR_orders$p, ",q=",
                                        AR_orders$q, ")", sep=""))
                AR_orders$arma <- as.character( AR_orders$arma )
                AR_orders[AR_orders=="corARMA(p=0,q=0)"] <- "NULL"

                AR_orders[,2]   <- as.character( AR_orders[,2] )
                self$corStructs <- AR_orders
              }
              if(!eqSpace)
              {
                temp    <- self
                temp$dv <- dV
                bestCors <- list()

                # this is supposedly taboo but I've been unable to work around it b/c we
                # cannot :: the %dopar% operator
                require(foreach)
                # parralelization
                funcs <- c("mean") # c(".eds") -- not importing from PersonAlytic correctly
                cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
                snow::clusterExport(cl, funcs)
                doSNOW::registerDoSNOW(cl)
                pkgs  <- c("gamlss", "nlme")

                bestCors <- foreach(id=uid, .packages = pkgs)  %dopar%
                {
                  corModsid <- list(); cc = 1
                  temp$correlation <- "NULL"
                  nullMod <- temp$lme(temp$data[[temp$ids]]==id)
                  if( "lme" %in% class(nullMod) )
                  {
                    for(p in 1:maxAR)
                    {
                      for(q in 1:maxMA)
                      {
                        # will this automatically update the fixed effects? it doesn't
                        # need to for nlme which takes `correlation` directly, but would
                        # need to be updated for gamlss; hence, lme for now (faster too)
                        cortemp <- paste("nlme::corARMA(p=", p, ",
                                                  q=", q, ")", sep="")
                        cortemp <- gsub('\n', '', cortemp)
                        cortemp <- gsub(' ', '', cortemp)
                        temp$correlation <- cortemp
                        corModsid[[cc]]  <- temp$lme(temp$data[[temp$ids]]==id)
                        if( any(corModsid[[cc]]=="Model did not converge") )
                        {
                          corModsid[[cc]] <- NULL
                        }
                        else cc = cc + 1
                      }
                    }
                    names(corModsid) <- unlist( lapply(corModsid,
                                                       function(x) x$call$correlation) )

                    if(lrt)
                    {
                      lrts <- lapply( lapply( corModsid,
                                              function(x) anova(x, nullMod)),
                                      function(x) x[2,])
                      lrts <- data.frame(cor=names(corModsid), plyr::rbind.fill(lrts))
                      wlrt <- which(lrts$`p.value`<=alpha)
                      if( length(wlrt) > 1 )
                      {
                        # we may need to do this recursively until nothing is
                        # sig, e.g., while loop
                        newnullmod <- corModsid[wlrt[1]]
                        nmnnm      <- names(newnullmod)
                        newnullmod <- newnullmod[[1]]
                        compmods   <- corModsid[wlrt[2:length(wlrt)]]
                        newlrts    <- lapply( lapply( compmods, function(x)
                                              anova(x, newnullmod)),
                                              function(x) x[2,])
                        wnl <- unlist(lapply(newlrts, length))
                        wnl <- which(wnl < max(wnl, na.rm=TRUE))
                        newlrts[wnl] <- NULL
                        compmods[wnl] <- NULL
                        newlrts    <- data.frame(cor=names(compmods),
                                                 plyr::rbind.fill(newlrts))
                        if( all(newlrts$p.value<.05) )
                        {
                          # minimum p-value is a bad criterion,
                          # serving only as a placeholder for
                          # a recursive search
                          return( as.character( names(compmods)[which.min(newlrts$p.value)] ) )
                        }
                        else return( nmnnm )
                      }
                      if( length(wlrt)==1 )
                      {
                        return( as.character( lrts$cor[wlrt] ) )
                      }
                      if( length(wlrt)==0 )
                      {
                        return( "NULL" )
                      }
                      #else stop('Failure in getAR_order, likelihood ratio test (lrt)')
                    }

                    if(!lrt)
                    {
                      if(crit=="AIC") ICs <- data.frame( unlist( lapply(corModsid, AIC) ) )
                      if(crit=="BIC") ICs <- data.frame( unlist( lapply(corModsid, BIC) ) )
                      else( stop( paste(crit, "is not a valid value for `crit`,",
                                        "use AIC or BIC.")))
                      ICs <- rbind(NULL = ifelse(crit=="AIC", AIC(nullMod),
                                                 BIC(nullMod)), ICs)
                      return( row.names(ICs)[which.min(ICs[,1])] )
                    }
                  }
                  if(! "nlme" %in% class(nullMod) )
                  {
                    return( "NULL" )
                  }
                }
                #parallel::stopCluster(cl)

                self$corStructs <- data.frame(ids=uid,
                                              arma=as.character( unlist(bestCors)) )
              } #oef !eqspace
              #self$corStructs$arma[self$corStructs$arma=="NULL"] <- NULL
            },
            overwrite = TRUE)

# crit can take on AIC or BIC
PalyticHTP$set("public", "GroupAR_order",
               function(dV, maxAR=3, maxMA=3, crit="BIC", lrt=TRUE, alpha=.05,
                        subgroup=NULL)
               {

                 corMods <- list(); cc <- 1
                 self$correlation <- "NULL"
                 nullMod <- self$lme(subgroup)
                 if( "lme" %in% class(nullMod) )
                 {
                   for(p in 1:maxAR)
                   {
                     for(q in 1:maxMA)
                     {
                       # will this automatically update the fixed effects? it doesn't
                       # need to for nlme which takes `correlation` directly, but would
                       # need to be updated for gamlss; hence, lme for now (faster too)
                       cortemp <- paste("nlme::corARMA(p=", p, ",
                                        q=", q, ")", sep="")
                       cortemp <- gsub('\n', '', cortemp)
                       cortemp <- gsub(' ', '', cortemp)
                       self$correlation <- cortemp
                       corMods[[cc]]  <- self$lme(subgroup)
                       if( any(corMods[[cc]]=="Model did not converge") )
                       {
                         corMods[[cc]] <- NULL
                       }
                       else cc = cc + 1
                     }
                   }
                   names(corMods) <- unlist( lapply(corMods,
                                                    function(x) x$call$correlation) )

                   if(lrt)
                   {
                     lrts <- lapply( lapply( corMods,
                                             function(x) anova(x, nullMod)),
                                     function(x) x[2,])
                     #print(lrts)
                     lrts <- data.frame(cor=names(corMods), plyr::rbind.fill(lrts))
                     wlrt <- which(lrts$`p.value`<=alpha)
                     if( length(wlrt) > 1 )
                     {
                       # we may need to do this recursively until nothing is
                       # sig, e.g., while loop
                       newnullmod <- corMods[wlrt[1]]
                       nmnnm      <- names(newnullmod)
                       newnullmod <- newnullmod[[1]]
                       compmods   <- corMods[wlrt[2:length(wlrt)]]
                       newlrts    <- lapply( lapply( compmods, function(x)
                         anova(x, newnullmod)),
                         function(x) x[2,])
                       wnl <- unlist(lapply(newlrts, length))
                       wnl <- which(wnl < max(wnl, na.rm=TRUE))
                       newlrts[wnl] <- NULL
                       compmods[wnl] <- NULL
                       newlrts    <- data.frame(cor=names(compmods),
                                                plyr::rbind.fill(newlrts))
                       if( all(newlrts$p.value<.05) & nrow(newlrts) > 1 )
                       {
                         # minimum p-value is a bad criterion,
                         # serving only as a placeholder for
                         # a recursive search
                         bestCor <- as.character( names(compmods)[which.min(newlrts$p.value)] )
                       }
                       else bestCor <- nmnnm
                     }
                     if( length(wlrt)==1 )
                     {
                       bestCor <- as.character( lrts$cor[wlrt] )
                     }
                     if( length(wlrt)==0 )
                     {
                       bestCor <- "NULL"
                     }
                   }

                   if(!lrt)
                   {
                     if(crit=="AIC") ICs <- data.frame( unlist( lapply(corMods, AIC) ) )
                     if(crit=="BIC") ICs <- data.frame( unlist( lapply(corMods, BIC) ) )
                     else( stop( paste(crit, "is not a valid value for `crit`,",
                                       "use AIC or BIC.")))
                     ICs <- rbind("NULL" = ifelse(crit=="AIC", AIC(nullMod),
                                                  BIC(nullMod)), ICs)
                     bestCor <- c("NULL", names(corMods))[which.min(ICs[,1])]
                   }
                 }
                 if(! "nlme" %in% class(nullMod) )
                 {
                   bestCor <- "NULL"
                 }

                 self$corStructs <- bestCor
               },
               overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
PalyticHTP$set("public", "getTime_Power",
               function(maxOrder=3)
               {
                 uid <- sort( as.numeric( unique(self$data[[self$ids]]) ) )
                 time_powers <- list()
                 temp <- PalyticHTP$new(self$data, self$ids, self$dv,
                                        self$time)
                 for(id in uid)
                 {
                   aics <- list()
                   for(i in 1:maxOrder)
                   {
                     temp$time_power <- i
                     mod0 <- temp$lme(self$data[[self$ids]]==id)
                     if("lme" %in% class(mod0))
                     {
                       aics[[i]] <- AIC( mod0 )
                     }
                     else aics[[i]] <- NA

                   }
                   aics <- unlist(aics)
                   time_powers[[id]] <- ifelse( all( is.na(aics) ),
                                                1,
                                                which.min( aics ) )
                 }
                 self$time_powers <- data.frame(ids=uid, time_power=unlist(time_powers))
               },
               overwrite = TRUE)

# hard coded lme at this point, option for gamlss later
PalyticHTP$set("public", "GroupTime_Power",
               function(maxOrder=3)
               {
                 uid <- sort( as.numeric( unique(self$data[[self$ids]]) ) )
                 time_powers <- list()
                 temp <- PalyticHTP$new(self$data, self$ids, self$dv,
                                        self$time)
                 aics <- list()
                 for(i in 1:maxOrder)
                 {
                   temp$time_power <- i
                   mod0 <- temp$lme()
                   if("lme" %in% class(mod0))
                   {
                     aics[[i]] <- AIC( mod0 )
                   }
                   else aics[[i]] <- NA

                 }
                 aics <- unlist(aics)
                 self$time_power <- which.min( aics )
               },
               overwrite = TRUE)


