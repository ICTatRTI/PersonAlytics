#' psuite - apply p.adjust and add output and graphics
#'
#' @name psuite
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#' @import gridExtra
#' @import plyr
#'
#' @description
#' Apply p.adjust and add output and graphics
#'
#' @param DVout A data frame created by \code{\link{PersonAlytic}}.
#' @param ids The variable in \code{DVout} across which to adjust.
#' @param output A character string for labeling the output
#' @param method One of \code{\link{p.adjust.methods}}.
#' @param nbest The \code{nbest} largest parameters for the target predictor with
#' \code{p<alpha} to be print to output.
#' @param alpha The type I error rate.
#' @param rawdata The raw data for plotting (not implemented).
#'

psuite <- function(DVout, ids, output, method="BY", nbest=NULL, alpha=.05,
                   rawdata=NULL)
{

  pav <- paste("-PAv", packageVersion("PersonAlytics"), "-", sep='')

  if(! method %in% p.adjust.methods)
  {
    stop("`method` must be one of: ", paste(p.adjust.methods, collapse = ', '))
  }

  # find the lrt p-value columns
  wc <- which( names(DVout) == 'targ_ivs_lrt_pvalue')

  # apply the adjustments by ids and dvs, if
  # - group based, byVariable will be dvs only
  # - only 1 dv, we still get ids within dv
  byVariable <- paste(DVout$dv, DVout[[ids]], sep="_")
  DVoutadj   <- plyr::rbind.fill(as.list(by(data = DVout[, c(ids, 'dv', 'target_iv',
                                                        'targ_ivs_lrt_pvalue')],
                                            INDICES = byVariable,
                                            FUN = adjuster, method=method)))
  DVoutadj   <- merge(DVout, DVoutadj)
  #all(DVoutadj$targ_ivs_lrt_pvalue== DVoutadj$targ_ivs_lrt_pvalue_raw,na.rm=TRUE)
  rm(DVout)

  if(!is.null(nbest))
  {
    #st <- unlist( strsplit(as.character( Sys.time()), " ") )
    #st[2] <- gsub(":", "-", st[2])
    #st <- paste(st, collapse=".")

    if(!exists("output")) stop("Please provide a name for your output.")

    dn <- paste('PersonAlytics output for', output)
    if(!dir.exists(dn))  dir.create(dn)

    fn <- paste(paste("./", dn, "/", sep=""),
                paste("Best", nbest, 'pvalues', sep="_"), sep="")

    rawColumn   <- which( names(DVoutadj) == "targ_ivs_lrt_pvalue_raw"  )
    targColumns <- which( names(DVoutadj) %in%
                            paste("targ_ivs_lrt_pvalue", method, "FDR", sep="_") )
    wc <- c(rawColumn, targColumns)

    # only works for continuous predictors
    targEffecSize <- which( names(DVoutadj) %in% c("Targ.Estimate", "Targ.Value") )

    #TODO(Stephen) this worked for the metabolomics data b/c we had a large number
    # of continuous covariates. For curelator, this is more challegenging b/c
    # some covariates are factors and some are continuous. This is on hold for
    # now b/c curelator doesn't need it. When yo uhave time to fix it, turn it
    # into a separate function.

    sink(file = paste(fn, pav, "txt", sep=".") )
    ln <- paste("\n", paste(rep("-", 80), collapse=""), "\n\n", sep = "")
    for(d in levels(DVoutadj$dv))
    {
      temp <- DVoutadj[DVoutadj$dv==d,]
      for(i in seq_along(wc))
      {
        # method
        if(i==1) mthd <- "Unadjusted"
        if(i >1) mthd <- method[i-1]

        fn <- paste(paste("./", dn, "/", sep=""),
                    paste("Best", nbest, mthd, 'pvalues_for', d, sep="_"), sep="")

        # escape if all p-values are missing
        if(all(is.na(temp[,wc[i]])))
        {
          cat("All p-values are NA\n", ln)
        }
        if(!all(is.na(temp[,wc[i]])))
        {
          # summary stats
          ss <- mean(temp[,wc[i]]<alpha, na.rm=TRUE)

          # section title
          cat(ln, nbest, "largest effects for `", d, "` with the", mthd, "p-value.\n\n",
              'Proportion p < alpha = ', alpha, ": ", ss, '\n',
              ln)
          # select rows with p < alpha
          tivloc <- which(names(temp)=='target_iv')
          best <- temp[temp[,wc[i]]<alpha,c(tivloc, targEffecSize, wc[i])]
          # select the nbest best
          nbestr <- nbest
          if(nbest > nrow(best)) nbestr <- nrow(best)
          if(nbestr > 0)
          {
            best <- best[order(abs(best[[names(temp)[targEffecSize]]]),
                               decreasing = TRUE),][seq(1,nbestr),]
          }
          best <- data.frame(best)
          row.names(best) <- NULL
          if(!all(is.na(best)))
          {
            print( best, digits = 3 )

            ## graphics
            #if(!is.null(rawdata))
            #{
            #  # add foreach %dopar% here
            #  bestplots <- list()
            #  for(p in seq_along(best))
            #  {
            #    w  <- which(temp$target_iv==best$target_iv[p])
            #    iv <- unlist( strsplit( toString(temp$ivs[w]), ", " ) )
            #    g  <- try(trajplot(data=rawdata                                         ,
            #                       ids=toString(temp$ids[w])                            ,
            #                       dv=toString(temp$dv[w])                              ,
            #                       time=toString(temp$time[w])                          ,
            #                       phase=toString(temp$phase[w])                        ,
            #                       ivs=ifelse(length(iv)>1, iv[seq_along(iv)[-1]], NULL) ,
            #                       target_iv=iv[length(iv)]                             ,
            #                       target_nm=toString(temp$target_iv[w])
            #    ), TRUE)
            #    if(! "try-error" %in% class(g) ) bestplots[[p]] <- g

            #  }
            #  pdf(paste(fn, pav, "pdf", sep="."), onefile = TRUE)
            #  lapply(bestplots, gridExtra::grid.arrange)
            #  dev.off()
            #}
          }
          if(all(is.na(best)))
          {
            cat("All p-values are >", alpha, "\n\n\n")
          }

          write.csv(best,  paste(fn, pav, "csv", sep="."))
        }
      }
    }
    while( sink.number() > 0 ) sink()

  }

  return( DVoutadj )
}

#' adjuster - a function to adjust p-values for a given column `wc` in the data
#' `x`.
#'
#' note that with missing data, the parameter `n` in p.adjust will be
#' length(temp), not length(na.omit(temp)). If the latter is used, results
#' will be more conservative, i.e., adjusting for the number of converged
#' analyses rather than the total number of analyses
adjuster <- function(x, method="BY")
{
  adj.ps <- p.adjust(x$targ_ivs_lrt_pvalue, method=method)
  output <- data.frame(x, adj.ps)
  names(output) <- c(names(x)[seq(1,3)],
                     paste(names(x)[4], 'raw', sep='_'),
                     paste(names(x)[4], method, 'FDR', sep='_'))
  return(output)
}

### metabolomics version, updates needed include
# - allow user to specify different time variables for analysis and visualization

#' trajplot - plot trajectories
#'
#' @name trajplot
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#' @import ggplot2
#'
#' @description
#' trajectoy plots
#'
#' @param data Inputs from psuite()
#' @param ids See \code{\link{PersonAlytics}}
#' @param dv See \code{\link{PersonAlytics}}
#' @param time See \code{\link{PersonAlytics}}
#' @param phase See \code{\link{PersonAlytics}}
#' @param ivs  See \code{\link{PersonAlytics}}
#' @param target_iv See \code{\link{PersonAlytics}}
#' @param target_nm The name of the target_iv
#'
trajplot <- function(data, ids, dv, time, phase, ivs, target_iv, target_nm)
{
  theColumns <- c(ids, dv, time, phase, ivs, target_iv)
  theColumns <- theColumns[! theColumns %in% c(1, "NULL", NA)]
  tempdata <- data[,theColumns]
  names(tempdata) <- c('id', 'dv', 'time', 'phase', 'ivs', 'target_iv')

  tempdata$dvr <- tempdata$dv
  tempdata$dv <- scale(tempdata$dv)
  tempdata$target_ivr <- tempdata$target_iv
  tempdata$target_iv <- scale(tempdata$target_iv)

  t1 <- Palytic$new(data = tempdata,
                    ids = "id",
                    dv  = "dv",
                    time = "time",
                    phase = "phase",
                    ivs = list("target_iv", "ivs"),
                    standardize = TRUE)

  t1$GroupTO()
  t1$GroupAR()

  # predict is looking for self, for now est directly
  cor <- eval(parse(text = ifelse(!is.null(t1$correlation),
                                  t1$correlation,
                                  'NULL')))
  t1.lme <- nlme::lme(t1$fixed, data = tempdata,
                      random = t1$random,
                      correlation = cor,
                      na.action = na.omit)
  t1.lme$call$fixed  <- t1$fixed
  t1.lme$call$random <- t1$random
  t1.lme$call$correlation <- cor

  # level = 0 is fixed effects only, level = 1 includes random effects
  preddat <- predict(t1.lme, tempdata, level = 1)

  # structure plotting data
  dvdat <- ivdat <- prdat <- tempdata

  dvdat$g <- t1$dv
  ivdat$g <- t1$ivs[[1]]
  prdat$g <- 'Predicted'

  dvdat$y <- t1$data[[t1$dv]]
  ivdat$y <- t1$data[[t1$ivs[[1]]]] # scaled above, do so here too?
  prdat$y <- preddat

  pdat <- rbind(dvdat, ivdat, prdat)
  pdat$id <- paste(pdat$id, pdat$g)


  g <- ggplot(pdat, aes(x=time, y=y, group=id, col=g)) +
    geom_line(position=position_jitter(height=.05), alpha=.25) +
    ylab("Standardized Scale (mean 0, unit variance)") +
    ggtitle(target_nm) +
    stat_summary(aes(group=g), fun.y=mean, geom="line", size=2) +
    facet_grid(.~phase, scales = "free") +
    guides(col=guide_legend(title=NULL)) +
    scale_colour_manual(values = c('blue', 'grey', 'red'),
                        labels=c(paste(dv, "Observed Values"),
                                 paste(dv, "Predicted Values"),
                                 target_nm))# +
  #scale_linetype_manual(values = c('dotted', 'dashed', 'solid'))


  return( g )

  #pdat$y[pdat$g==t1$ivs[[1]]] <- 2^( tempdata$X19.43_887.1030m.zr )
  #
  #ggplot(pdat, aes(x=Time, y=y, group=id, col=g)) +
  #  geom_line(position=position_jitter(height=.05), alpha=.25) +
  #  ylab("Standardized Scale (mean 0, unit variance)") +
  #  ggtitle(paste(as.character(t1$fixed)[c(2,1,3)], collapse = "")) +
  #  stat_summary(aes(group=g), funy.y=mean, geom="line", size=2.5) +
  #  facet_grid(.~Tx, scales = "free") +
  #  guides(col=guide_legend(title=NULL)) +
  #  scale_colour_discrete(name="", labels=c("PASAT Predicted Values", t1$ivs[[1]]))

}


if(1==2)
{
  # compound only
  ldat <- tempdata
  ldat$target_iv <- 2^ldat$target_iv
  ldat$transform=3
  tempdata$transform=2
  tempdata <- rbind(tempdata,ldat)

  ldat <- tempdata
  ldat$transform <- 1
  ldat$target_iv <- ldat$PASAT
  tempdata <- rbind(ldat, tempdata)
  tempdata$transform <- factor(tempdata$transform, labels=c('PASAT', 'Transformed Compound',
                                                            'Raw Compound'))

  ggplot(tempdata, aes(x=Time, y=target_iv, group=interaction(ID, Tx),
                       color=Tx)) + ylab('') +
    geom_line(position=position_jitter(height=.1), alpha=.5) +
    stat_summary(aes(group=Tx), fun.y=median, geom="line", size=3) +
    facet_grid(transform~., scales="free")
}
