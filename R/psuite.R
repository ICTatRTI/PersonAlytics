#' psuite - apply p.adjust and add output and graphics
#'
#' @name psuite
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @usage Palytic(DVout, method)
#'
#' @details
#' @description apply p.adjust and add output and graphics
#'
#' @param DVout A data frame created by \code{\link{PersonAlyticHTP}}.
#' @param method One of \code{\link{p.adjust.methods}}.
#' @param nbest The \code{nbest} largest pameters for the target predictor with
#' \code{p<alpha} to be print to output.
#' @param alpha The type I error rate.
#' @param rawdata The raw data for plotting.
#'
#' @examples

psuite <- function(DVout, method="BY", nbest=25, alpha=.05,
                   rawdata=NULL)
{
  library(gridExtra)
  library(ggplot2)

  if(! method %in% p.adjust.methods)
  {
    stop("`method` must be one of", p.adjust.methods)
  }

  # find all the p-value columns
  wc <- which(! lapply( strsplit(names(DVout), '.p.value'),
                  function(DVout) paste(DVout, collapse='')) == names(DVout))

  adjuster <- function(x, wc, method)
  {
    temp <- x[,wc]
    adj.ps <- lapply(temp, function(x) p.adjust(x, method=method))
    names(adj.ps) <- paste(names(x)[wc], method, sep='_')
    return(data.frame(x, adj.ps))
  }
  DVoutadj <- do.call(rbind, as.list(by(DVout, DVout$dv,
                              adjuster, wc=wc, method=method)))

  tn <- which(names(DVoutadj)=='target_iv')
  te <- which(grepl("TargetPredictor.Value", names(DVoutadj)))
  tp <- which(grepl(glob2rx("TargetPredictor*p.value*"), names(DVoutadj)))

  st <- unlist( strsplit(as.character( Sys.time()), " ") )
  st[2] <- gsub(":", "-", st[2])
  st <- paste(st, collapse="_")

  dn <- paste('PersonAlyticsPro p-value report', st)
  if(!dir.exists((dn)))  dir.create(dn)

  fn <- paste(paste("./", dn, "/", sep=""),
              paste("Best", nbest, 'pvalues', sep="_"), sep="")
  print(fn)

  sink(file = paste(fn, "txt", sep=".") )
  ln <- paste("\n", paste(rep("-", 80), collapse=""), "\n\n", sep = "")
  for(d in unique(DVoutadj$dv))
  {
    temp <- DVoutadj[DVoutadj$dv==d,]
    for(i in 1:length(tp))
    {
      # method
      if(i==1) mthd <- "Unadjusted"
      if(i >1) mthd <- method[i-1]

      fn <- paste(paste("./", dn, "/", sep=""),
                  paste("Best", nbest, mthd, 'pvalues_for', d, sep="_"), sep="")

      # section stats
      ss <- mean(temp[,tp[i]]<alpha, na.rm=TRUE)

      # section title
      cat(ln, nbest, "largest effects for `", d, "` with the", mthd, "p-value.\n\n",
          'Proportion p < alpha = ', alpha, ": ", ss, '\n',
          ln)
      # select p < alpha
      best <- temp[temp[,tp[i]]<alpha,c(tn,te,tp[i])]
      # select the nbest best
      nbestr <- nbest; if(nbest > nrow(best)) nbestr <- nrow(best)
      if(nbestr > 0)
      {
        best <- best[order(abs(best[,(1+i)]), decreasing = TRUE),][1:nbestr,]
      }
      best <- data.frame(best)
      row.names(best) <- NULL
      if(!all(is.na(best)))
      {
        print( best )

        # graphics
        if(!is.null(rawdata))
        {
          # add foreach %dopar% here
          bestplots <- list()
          for(p in 1:nrow(best))
          {
            w  <- which(temp$target_iv==best$target_iv[p])
            iv <- unlist( strsplit( toString(temp$ivs[w]), ", " ) )
            g  <- try(trajplot(data=rawdata                                         ,
                               ids=toString(temp$ids[w])                            ,
                               dv=toString(temp$dv[w])                              ,
                               time=toString(temp$time[w])                          ,
                               phase=toString(temp$phase[w])                        ,
                               ivs=ifelse(length(iv)>1, iv[1:(length(iv)-1)], NULL) ,
                               target_iv=iv[length(iv)]                             ,
                               target_nm=toString(temp$target_iv[w])
            ), TRUE)
            if(! "try-error" %in% class(g) ) bestplots[[p]] <- g

          }
          pdf(paste(fn, "pdf", sep="."), onefile = TRUE)
          lapply(bestplots, grid.arrange)
          dev.off()
        }
      }
      if(all(is.na(best)))
      {
        cat("All p-values are >", alpha, "\n\n\n")
      }

      write.csv(best,  paste(fn, "csv", sep="."))
    }
  }
  sink()


  return( DVoutadj )
}

### metabolomics version, updates needed include
# - allow user to specify different time variables for analysis and visualization

trajplot <- function(data, ids, dv, time, phase, ivs, target_iv, target_nm)
{
  library(ggplot2)

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
                   ivs = list("target_iv", "ivs"))

  t1$GroupAR_order("dv"); t1$correlation
  t1$GroupTime_Power(); t1$time_power
  t1.lme <- t1$lme()

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
    geom_line(position=position_jitter(h=.05), alpha=.25) +
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
  #  geom_line(position=position_jitter(h=.05), alpha=.25) +
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
    geom_line(position=position_jitter(h=.1), alpha=.5) +
    stat_summary(aes(group=Tx), fun.y=median, geom="line", size=3) +
    facet_grid(transform~., scales="free")
}