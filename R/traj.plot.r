traj.plot <- function(plotName	      = paste("_PACT_trajectory_plot.jpg", sep=""),
                      Individual	    = c(TRUE , "black", "solid", .5),
                      Avg.Observed	  = c(FALSE, "red", "solid", 2),
                      IC.Observed	    = c(FALSE, "pink", "solid", 1),
                      Avg.Fitted	    = c(FALSE, "yellow", "solid", 2),
                      Avg.IC		      = c(FALSE, "blue", "solid", 1),
                      Avg.PI		      = c(FALSE, "green", "solid", 1),
                      Fit.Ind         = c(FALSE, "black", "solid", .5),
                      Ind.Smooth		  = c(FALSE, "red", "solid", 2),
                      Fit.Smooth		  = c(TRUE, "yellow", "solid", 2),
                      Ind.Smooth.IC	  = c(FALSE, "pink", "solid", 2),
                      Fit.Smooth.IC	  = c(TRUE, "blue", "solid", 2),
                      time.power.plot = 1,
                      xlabel		      = "",
                      ylabel		      = "",
                      boxlabel	      = "",
                      phaselable	    = "Phase",
                      w	              = 11,
                      h		            = 8.5,
                      IC.transparency = .50,
                      conf.int	      = .95,
                      boxCols			    = NULL,
                      outlier.color	  = NULL,
                      subgroup.var    = NULL,
                      subgroup        = NULL)
{
  #load("first_run_plotting_inputs.Rdata") #testing defaults only
  load("plotData.Rdata")

  ### subset the data if requested
  if(!is.null(subgroup))
  {
    w <- dat[,subgroup.var]==subgroup
    dat <- dat[w,]
    plotdat$`Individual Observed` <- plotdat$`Individual Observed`[w,]
    plotdat$`Individual Model Based`   <- plotdat$`Individual Model Based`[w,]
  }

  ### augment data with plotting inputs
  #   get fitted data and standard errors
  cmult <- qnorm(conf.int+.5*(1-conf.int))

  #   IC for each person
  plotdat$"Model Based Avg. IC" <- data.frame(group=dat[,id],
                                              x=dat[,time],
                                              ymin=plotdat$"Individual Model Based"$y-dat$SE*cmult,
                                              ymax=plotdat$"Individual Model Based"$y+dat$SE*cmult)

  #   PI for each person
  plotdat$"Model Based Avg.  PI" <- data.frame(group=dat[,id],
                                               x=dat[,time],
                                               ymin=plotdat$"Individual Model Based"$y-dat$SE2*cmult,
                                               ymax=plotdat$"Individual Model Based"$y+dat$SE2*cmult)

  #   IC mean at each time point
  plotdat$"Model Based Avg. IC" <- aggregate(plotdat$"Model Based Avg. IC"[,3:4],
                                             by=list(plotdat$"Model Based Avg. IC"$x), FUN=mean)
  names(plotdat$"Model Based Avg. IC")[1] <- "x"
  plotdat$"Model Based Avg. IC" <- data.frame(plotdat$"Model Based Avg. IC")

  #   PI mean at each time point
  plotdat$"Model Based Avg.  PI" <- aggregate(plotdat$"Model Based Avg.  PI"[,3:4],
                                              by=list(plotdat$"Model Based Avg.  PI"$x), FUN=mean)
  names(plotdat$"Model Based Avg.  PI")[1] <- "x"
  plotdat$"Model Based Avg.  PI" <- data.frame(plotdat$"Model Based Avg.  PI")

  #   Observed mean at each time point
  plotdat$"Observed Avg." <- aggregate(plotdat$"Individual Observed"$y,
                                       list(plotdat$"Individual Observed"$x), FUN=mean)
  plotdat$"Observed Avg." <- data.frame(x=plotdat$"Observed Avg."$Group.1,
                                        y=plotdat$"Observed Avg."$x)

  #   Observed mean IC at each time point
  omeansd <- aggregate(plotdat$`Individual Observed`$y,
                       list(plotdat$`Individual Observed`$x), FUN=sd, na.rm=TRUE)
  names(omeansd) <- c("x", "sd")
  omeann  <- aggregate(plotdat$`Individual Observed`$y,
                       list(plotdat$`Individual Observed`$x), FUN=length)
  names(omeann) <- c("x", "n")
  omeanse <- merge(plotdat$"Observed Avg.", merge(omeansd, omeann))
  omeanse$se <- omeanse$sd/sqrt(omeanse$n)
  plotdat$`Observed Avg. IC` <- data.frame(x=omeanse$x,
                                           y=omeanse$y,
                                           ymin=omeanse$y - omeanse$se*cmult,
                                           ymax=omeanse$y + omeanse$se*cmult)


  #   Fitted mean at each time point
  plotdat$"Model Based Avg." <- aggregate(plotdat$"Individual Model Based"[,1:2],
                                          list(plotdat$"Individual Model Based"$x), FUN=mean)
  plotdat$"Model Based Avg." <- data.frame(x=plotdat$"Model Based Avg."$x,
                                           y=plotdat$"Model Based Avg."$y)

  #   Smoothed observed data & IC
  obs.gam <- gam(y ~ s(x, bs = "cs"), data=plotdat$`Individual Observed`)
  obs.gam <- predict(obs.gam, se.fit=TRUE)
  obs.gam.y  <- aggregate(obs.gam$fit, list(dat[,time]), mean, na.rm=TRUE)
  obs.gam.se <- aggregate(obs.gam$se.fit, list(dat[,time]), mean, na.rm=TRUE)
  xna <- rep(NA, nrow(obs.gam.y))
  xob <- rep("Smooth Observed Avg.", nrow(obs.gam.y))
  plotdat$"Smooth Observed Avg." <- data.frame(x=obs.gam.y$Group.1,
                                               y=obs.gam.y$x,
                                               ymin=obs.gam.y$x - obs.gam.se$x*cmult,
                                               ymax=obs.gam.y$x + obs.gam.se$x*cmult)

  #   Smoothed Fitted data & IC (use observed se for now)
  fit.gam <- gam(y ~ s(x, bs = "cs"), data=plotdat$`Individual Model Based`)
  fit.gam <- predict(fit.gam, se.fit=TRUE)
  fit.gam.y  <- aggregate(fit.gam$fit, list(dat[,time]), mean, na.rm=TRUE)
  fit.gam.se <- aggregate(fit.gam$se.fit, list(dat[,time]), mean, na.rm=TRUE)
  xna <- rep(NA, nrow(fit.gam.y))
  xob <- rep("Smooth Model Based", nrow(fit.gam.y))
  plotdat$"Smooth Model Based" <- data.frame(x=fit.gam.y$Group.1,
                                             y=fit.gam.y$x,
                                             ymin=fit.gam.y$x - obs.gam.se$x*cmult,
                                             ymax=fit.gam.y$x + obs.gam.se$x*cmult)

  # get average phase times
  phase.locs <- phase.loc(dat[,id], dat[,phase], dat[,time])
  phase_loc  <- phase.locs$phase_loc
  phase_tic  <- phase.locs$phase_tic
  phase_nms  <- phase.locs$phase_nms

  ### general plotting inputs
  if( xlabel=="" ) xlabel <- time
  if( ylabel=="" ) ylabel <- dv
  ylim <- range(c(plotdat$`Individual Observed`$y, plotdat$`Individual Model Based`$y), na.rm=TRUE)
  xlim <- range(c(plotdat$`Individual Observed`$x, plotdat$`Individual Model Based`$x), na.rm=TRUE)
  ylim <- range(pretty(ylim))
  xlim <- range(pretty(xlim))

  ### adjust color transparency
  Avg.IC[2] <- adjustcolor(Avg.IC[2], alpha.f=IC.transparency)
  Avg.PI[2] <- adjustcolor(Avg.PI[2], alpha.f=IC.transparency)
  Ind.Smooth.IC[2] <- adjustcolor(Ind.Smooth.IC[2], alpha.f=IC.transparency)
  Fit.Smooth.IC[2] <- adjustcolor(Fit.Smooth.IC[2], alpha.f=IC.transparency)

  ### set up blank plot
  jpeg(file = plotName, width = w, height = h, units = "in", res=500)
  par(mfrow=c(1,1), bg="white", mar=c(5.1, 4.1, 4.1, 16.1), xpd=TRUE)
  plot(plotdat$`Individual Observed`$x, plotdat$`Individual Observed`$y,
       col="transparent", axes=FALSE, xlab=xlabel, ylab=ylabel,
       xlim=xlim, ylim=ylim)
  alim <- par("usr")
  rect(alim[1], alim[3], alim[2], alim[4], col="lightgray", border="NA")

  ### set up axes
  x1 <- plotdat$`Individual Observed`$x
  axis(1, at = pretty(xlim))
  axis(2, at = pretty(ylim))
  axis(3, at=phase_loc, labels=paste(phase, "=", phase_nms), tick = FALSE)
  axis(3, at=phase_tic, labels=FALSE, tick = TRUE)
  #mtext(phaselable, 3, line = 2)
  par(xpd=FALSE)
  grid(col="white", lty=1)

  ### pre-allocate memory for legend inputs
  cols.i <- ltys.i <- lwds.i <- vector("list", 0)
  cols.a <- ltys.a <- lwds.a <- vector("list", 0)

  ### plot trajectories
  if(Individual[1])
  {
    u.i <- unique(plotdat$`Individual Observed`$group)
    for(i in seq_along(u.i))
    {
      temp.dat <- plotdat$`Individual Observed`[plotdat$`Individual Observed`$group==u.i[i],]
      lines(temp.dat$x, temp.dat$y, col=Individual[2], lty=Individual[3], lwd=Individual[4])
    }
    cols.i$`Observed` <- Individual[2]
    ltys.i$`Observed` <- Individual[3]
    lwds.i$`Observed` <- Individual[4]
  }
  if(Fit.Ind[1]=="TRUE" & Individual[1]!="TRUE")
  {
    u.i <- unique(plotdat$`Individual Model Based`$group)
    for(i in seq_along(u.i))
    {
      temp.dat <- plotdat$`Individual Model Based`[plotdat$`Individual Model Based`$group==u.i[i],]
      lines(temp.dat$x, temp.dat$y, col=Fit.Ind[2], lty=Fit.Ind[3], lwd=Fit.Ind[4])
    }
    cols.i$`Model Based` <- Fit.Ind[2]
    ltys.i$`Model Based` <- Fit.Ind[3]
    lwds.i$`Model Based` <- Fit.Ind[4]
  }
  if(Avg.Observed[1])
  {
    lines(plotdat$`Observed Avg.`$x, plotdat$`Observed Avg.`$y,
          col=Avg.Observed[2], lty=Avg.Observed[3], lwd=Avg.Observed[4])
    cols.a$`Observed` <- Avg.Observed[2]
    ltys.a$`Observed` <- Avg.Observed[3]
    lwds.a$`Observed` <- Avg.Observed[4]
  }
  if(IC.Observed[1])
  {
    polygon(c(plotdat$`Observed Avg. IC`$x, plotdat$`Observed Avg. IC`$x),
            c(plotdat$`Observed Avg. IC`$ymin, plotdat$`Observed Avg. IC`$ymax),
            col=IC.Observed[2], lty=IC.Observed[3], lwd=IC.Observed[4], border=NA)
    cols.a$`Observed Avg. IC` <- IC.Observed[2]
    ltys.a$`Observed Avg. IC` <- IC.Observed[3]
    lwds.a$`Observed Avg. IC` <- IC.Observed[4]
  }
  if(Avg.PI[1])
  {
    polygon(c(plotdat$`Model Based Avg.  PI`$x, plotdat$`Model Based Avg.  PI`$x[nrow(plotdat$`Model Based Avg.  PI`):1]),
            c(plotdat$`Model Based Avg.  PI`$ymax, plotdat$`Model Based Avg.  PI`$ymin[nrow(plotdat$`Model Based Avg.  PI`):1]),
            col=Avg.PI[2], lty=Avg.PI[3], lwd=Avg.PI[4], border=NA)
    cols.a$`PI Model Based` <- Avg.PI[2]
    ltys.a$`PI Model Based` <- Avg.PI[3]
    lwds.a$`PI Model Based` <- Avg.PI[4]
  }
  if(Avg.IC[1])
  {
    polygon(c(plotdat$`Model Based Avg. IC`$x, plotdat$`Model Based Avg. IC`$x[nrow(plotdat$`Model Based Avg. IC`):1]),
            c(plotdat$`Model Based Avg. IC`$ymax, plotdat$`Model Based Avg. IC`$ymin[nrow(plotdat$`Model Based Avg. IC`):1]),
            col=Avg.IC[2], lty=Avg.IC[3], lwd=Avg.IC[4], border=NA)
    cols.a$`IC Model Based` <- Avg.IC[2]
    ltys.a$`IC Model Based` <- Avg.IC[3]
    lwds.a$`IC Model Based` <- Avg.IC[4]
  }
  if(Avg.Fitted[1])
  {
    lines(plotdat$`Model Based Avg.`$x, plotdat$`Model Based Avg.`$y,
          col=Avg.Fitted[2], lty=Avg.Fitted[3], lwd=Avg.Fitted[4])
    cols.a$`Model Based` <- Avg.Fitted[2]
    ltys.a$`Model Based` <- Avg.Fitted[3]
    lwds.a$`Model Based` <- Avg.Fitted[4]
  }
  if(Ind.Smooth.IC[1])
  {
    polygon(c(plotdat$`Smooth Observed`$x, plotdat$`Smooth Observed`$x[nrow(plotdat$`Smooth Observed`):1]),
            c(plotdat$`Smooth Observed`$ymax, plotdat$`Smooth Observed`$ymin[nrow(plotdat$`Smooth Observed`):1]),
            col=Ind.Smooth.IC[2], lty=Ind.Smooth.IC[3], lwd=Ind.Smooth.IC[4], border=NA)
    cols.a$`IC Smooth Observed` <- Ind.Smooth.IC[2]
    ltys.a$`IC Smooth Observed` <- Ind.Smooth.IC[3]
    lwds.a$`IC Smooth Observed` <- Ind.Smooth.IC[4]
  }
  if(Ind.Smooth[1])
  {
    lines(plotdat$`Smooth Observed`$x, plotdat$`Smooth Observed`$y,
          col=Ind.Smooth[2], lty=Ind.Smooth[3], lwd=Ind.Smooth[4])
    cols.a$`Smooth Observed` <- Ind.Smooth[2]
    ltys.a$`Smooth Observed` <- Ind.Smooth[3]
    lwds.a$`Smooth Observed` <- Ind.Smooth[4]
  }
  if(Fit.Smooth.IC[1])
  {
    polygon(c(plotdat$`Smooth Model Based`$x, plotdat$`Smooth Model Based`$x[nrow(plotdat$`Smooth Model Based`):1]),
            c(plotdat$`Smooth Model Based`$ymax, plotdat$`Smooth Model Based`$ymin[nrow(plotdat$`Smooth Model Based`):1]),
            col=Fit.Smooth.IC[2], lty=Fit.Smooth.IC[3], lwd=Fit.Smooth.IC[4], border=NA)
    cols.a$`IC Smooth Model Based` <- Fit.Smooth.IC[2]
    ltys.a$`IC Smooth Model Based` <- Fit.Smooth.IC[3]
    lwds.a$`IC Smooth Model Based` <- Fit.Smooth.IC[4]
  }
  if(Fit.Smooth[1])
  {
    lines(plotdat$`Smooth Model Based`$x, plotdat$`Smooth Model Based`$y,
          col=Fit.Smooth[2], lty=Fit.Smooth[3], lwd=Fit.Smooth[4])
    cols.a$`Smooth Model Based` <- Fit.Smooth[2]
    ltys.a$`Smooth Model Based` <- Fit.Smooth[3]
    lwds.a$`Smooth Model Based` <- Fit.Smooth[4]
  }

  ### legend
  cols.i <- unlist(cols.i)
  ltys.i <- unlist(ltys.i)
  lwds.i <- unlist(lwds.i)
  cols.a <- unlist(cols.a)
  ltys.a <- unlist(ltys.a)
  lwds.a <- unlist(lwds.a)
  par(xpd=TRUE)

  if(Individual[1]=="TRUE" | Fit.Ind[1]=="TRUE")
  {
    legend(xlim[2] + .05*(xlim[2]-xlim[1]), max(ylim), legend=names(cols.i),
           fill=cols.i, lty=ltys.i,# merge=TRUE, lwd=lwds,
           bg="lightgrey", title = "Individual Trajectories")
  }
  legend(xlim[2] + .05*(xlim[2]-xlim[1]), mean(ylim), legend=names(cols.a),
         fill=cols.a, lty=ltys.a, #merge=TRUE, lwd=lwds,
         bg="lightgrey", title = "Average Trajectory")

  legend(x=xlim[2] + .05*(xlim[2]-xlim[1]), y=ylim[1],
         legend="Note: Phases are deliminated at\n average start, transition,\nand end points",
         bty="n")

  par(xpd=FALSE)
  ### turn off jpeg device
  dev.off()


  ### box plot parameters
  if(is.null(boxCols))
  {
    n.cols <- n.phases <- length(unique(dat[,phase]))
    if(n.phases<3) n.cols <- 3
    boxCols  <- brewer.pal(n.cols, "Accent")
    if(n.phases!=length(boxCols)) boxCols <- boxCols[seq(1,n.phases)]
  }
  if(is.null(outlier.color))
  {
    outlier.color <- "red"
  }
  if(boxlabel=="")
  {
    boxlabel <- phase
  }

  ### box plot at each phase
  b1nm <- "PACT_boxplot.jpg"
  if(!file.exists(b1nm))
  {
    print(
      ggplot(mapping = aes(y=y, x=as.factor(phase)), data=plotdat$"Individual Observed") +
        stat_boxplot(geom ="errorbar", width = 0.5) +
        geom_boxplot(aes(fill = as.factor(phase)), outlier.color=outlier.color) +
        scale_fill_manual(values = boxCols, name = boxlabel) + xlab("") + ylab(ylabel) +
        ggtitle(paste("Observed Mean", dv, "by Phase")) +
        labs(caption = paste(
          "The line in the middle of the box is the median observed value.\n",
          "The edges of the box are the 1st and 3rd quantiles (25th and 75th percentiles).\n",
          "The whiskers extend from the edge of the box to the observed data point that is, at most, 1.5 times the interquartile range.")) +
        theme(plot.caption=element_text(hjust=0))
    )
    ggsave(file = b1nm, device="jpeg", width = w, height = h, units = "in", dpi = 300)
  }

  ### Ty"s custom boxplot
  # Having another plot that summarizes results for efficacy estimates would be especially
  # appealing to our users…  So, the features described above would be retained, but the
  # centering line would representing the mean (based on the model’s predicted outcomes
  # and the median appearing as a dot), the upper and lower ends of the box representing
  # 95% IC of the mean, and the whiskers presenting standard deviations.  I’m hoping that
  # for this second plot, you could use existing R code and simply replace the typical
  # box-and-whiskers values with these alternatives?
  TyStats <- function(x, conf.int)
  {
    x.sd <- sd(x)
    x.IC <- as.data.frame( t(IC(x, conf.int)) )
    r <- c(x.IC$mean-x.sd, x.IC$lower, x.IC$mean, x.IC$upper, x.IC$mean + x.sd)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  b2nm <- "PACT_boxplot2.jpg"
  if(!file.exists(b2nm))
  {
    print(
      ggplot(mapping = aes(x=as.factor(phase), y=y, fill=as.factor(phase)),
             data=plotdat$"Individual Model Based") +
        stat_summary(fun.data = TyStats, geom = "boxplot", fun.args = list(conf.int=conf.int)) +
        #stat_boxplot(geom ="errorbar", width = 0.25) +
        scale_fill_manual(values = boxCols, name = boxlabel) + xlab("") + ylab(ylabel) +
        ggtitle(paste("Modelled Average", dv, "by Phase")) +
        labs(caption = paste("The line in the middle of the box is the mean of the fitted values.\n",
                             "The edges of the box are at the ends of the 95% IC for the mean of the fitted values.\n",
                             "The whiskers extend from the line in the middle of the box to +/- 1 SD of the fitted values."))+
        theme(plot.caption=element_text(hjust=0))
    )
    ggsave(file = b2nm, device="jpeg", width = w, height = h, units = "in", dpi = 300)
  }
}
