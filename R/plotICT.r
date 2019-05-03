# the reason this isn't in the Palytic method $plot() is because it has multiple
# calls in that method

#' plotICT - function to plot ICT data.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param out An \code{lme} object
#'
#' @export
#' @import ggplot2
#' @import gridExtra
#'
#' @examples
#'
#' t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
#'         time='Time', phase='Phase')
#' t1$plot()
#' t1$plot(groupvar = 'Target1')

# Note that this is implemented as a
# function rather than a method for Palytic so that other packages can use
# it's functionality. Used in Palytic$plot()...wait...it only takes a Palytic
# object 'self'...maybe move this back to the method?? This is still a utility
# function that gets multiple calls
plotICT <- function(self, data, legendName=NULL,
                    type=c('density', 'histogram', 'freqpoly'),
                    ylim=NULL)
{
  # check whetehr the time variable is continuous, if yes it
  # needs to be aggregated so that we don't have time points with only
  # one observation
  time <- data[[self$time]]
  if(!all(round(time) == time, na.rm = TRUE))
  {
    breaks <- hist(time, plot = FALSE)$breaks
    r      <- nchar(strsplit(as.character(breaks), "\\.")[[1]][2])
    data[[self$time]] <- round(time, r)
  }
  rm(time)

  # clean out observations with missing time variables
  data <- data[!is.na(data[[self$time]]),]

  # summarize the data for plotting
  summData <- summarySE(data=data, measurevar=self$dv, groupvars=self$time,
                        phase=self$phase, na.rm=TRUE)
  # set up grouping variable(s)
  group <- 1
  #*#if(length(groupvars)>1) group <- groupvars[2] # TODO general to 2+ groupvars

  # calculate dodging & range
  xrange <- range(data[[self$time]], na.rm=TRUE)
  pd <- position_dodge(0.1*(xrange[2]-xrange[1]))
  #ylim <- range(c(summData$sdlo, summData$sdhi), na.rm=TRUE) # this results in truncated tails
  if(is.null(ylim)) ylim <- range(self$datac[[self$dv]], na.rm=TRUE)

  # set up phase colors / borrowed from ICTviz() in PersonAlyticsPower
  # see'rects' in PersonAlyticsPower setup.r ICTviz()
  #TODO this won't work for non-unique phase, eg ABA
  if(!is.null(self$phase))
  {
    phaseDup  <- !duplicated(summData[[self$phase]])
    rect      <- summData[phaseDup,]
    phaseEnd  <- which(phaseDup) - 1
    rect$xmax <- summData[[self$time]][c(phaseEnd[2:length(phaseEnd)]+1,
                                         length(phaseDup))]
    # + 1 b/c requires 3+ colors, this allows for 1 & 2 phase studies
    rect$cols <- RColorBrewer::brewer.pal(nrow(rect) + 1, 'Accent')[1:nrow(rect)]

    # phase colors will be incorrect if the phase levels are not in
    # the data order
    rect[[self$phase]] <- factor(as.character(rect[[self$phase]]),
                                 levels=as.character(rect[[self$phase]]))
  }

  # legend name
  if(!is.null(legendName)) legendName <- paste(self$phase, legendName, sep=': ')
  if( is.null(legendName)) legendName <- self$phase

  # theme
  plotTheme <- theme(panel.background = element_rect(fill="white"),
                     panel.grid.major = element_line(color="lightgrey"))
  #panel.grid.minor = element_line(color="lightgrey")


  # trajectory plot
  s <- ggplot(summData, aes_string(x=self$time, y=self$dv,
                                   group=group, col=group)) +
    geom_errorbar(aes(ymin=sdlo, ymax=sdhi),
                  width = .1, position = pd) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    ylab('') + ylim(ylim) + plotTheme

  #
  if(!is.null(self$phase))
  {
    s <- s + geom_rect(data=rect,
                       aes_string(xmin=self$time, xmax='xmax',
                                  ymin=-Inf, ymax=Inf,
                                  fill=self$phase),
                       alpha=0.4, color=NA, inherit.aes=F) +
         scale_fill_manual(values=rect$cols, name=legendName,
                        labels=rect[[self$phase]])
  }

  # turn off line legend
  if(group==1) s <- s + scale_color_continuous(guide = FALSE)

  # density with no phase
  if( is.null(self$phase))
  {
    densdat <- self$datac
    d <- ggplot(data=densdat,
                aes_string(x=self$dv)) +
      theme(legend.position="none") +
      xlim(ylim) + plotTheme
    d <- d + coord_flip() + scale_y_reverse()
    if(type=='density')   d <- d + geom_density(alpha=.5, col='darkblue', fill='lightblue')
    if(type=='histogram') d <- d + geom_histogram(alpha=.5, col='darkblue', fill='lightblue')
    if(type=='freqpoly')  d <- d + geom_freqpoly(alpha=.5, col='darkblue', fill='lightblue')
  }

  # denisty by phase
  if(!is.null(self$phase))
  {
    densdat <- self$datac
    densdat[[self$phase]] <- factor(densdat[[self$phase]])
    d <- ggplot(data=densdat,
                aes_string(x=self$dv,
                           colour=self$phase,
                           fill=self$phase)) +
      scale_fill_manual(values = rect$cols) +
      scale_color_manual(values = rect$cols) +
      theme(legend.position="none") +
      geom_density(alpha=.5) +
      xlim(ylim) + plotTheme
    d <- d + coord_flip() + scale_y_reverse()
    if(type=='density')   d <- d + geom_density(alpha=.5)
    if(type=='histogram') d <- d + geom_histogram(alpha=.5)
    if(type=='freqpoly')  d <- d + geom_freqpoly(alpha=.5)
  }

  return( list(d=d, s=s) )
}
