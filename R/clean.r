#' clean, a function to check data and inputs.
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#' @param dv See \code{\link{PersonAlytic}}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param ivs Used by \code{PersonAlyticPro}.
#' @param dvs Used by \code{PersonAlyticPro}.
#' @param target_ivs Used by \code{PersonAlyticPro}.
#' @param fixed See \code{\link{Palytic}}.
#' @param random See \code{\link{Palytic}}.
#' @param formula See \code{\link{Palytic}}.
#' @param standardize Logical. Should all variables be standardized? Only applies
#' to \code{dv}, \code{ivs}, and \code{target_ivs}.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @examples
#' \dontrun{
#' data <- clean(data=OvaryICT, ids='Mare', dv='follicles',
#'               time='Time', phase='Phase')
#' # illustrate dropping cases with < 2 observations
#' data <- clean(data=OvaryICT[29:nrow(OvaryICT),], ids='Mare', dv='follicles',
#'               time='Time', phase='Phase')
#' }
#'
#' @keywords internal

clean <- function(data		 	            ,
                  ids		 	              ,
                  dv		 	              ,
                  time		 	            ,
                  phase	        =	NULL	,
                  ivs	          =	NULL	,
                  fixed       	=	NULL	,
                  random      	=	NULL	,
                  formula	      =	NULL	,
                  correlation 	=	NULL	,
                  family        = NULL  ,
                  dvs	          =	NULL	,
                  target_ivs    =	NULL	,
                  standardize 	=	FALSE	,
                  sortData	    =	TRUE	,
                  alignPhase 	  =	TRUE  )
{
  # check that variables are in the data set
  vars <- unique( c(ids, dv, time, phase, unlist(ivs), unlist(dvs), unlist(target_ivs),
                    all.vars(fixed), all.vars(random), all.vars(formula)) )
  vars <- vars[which(vars!='1')]
  wvars <- which( ! vars %in% names(data) )
  if( length(wvars) > 0 )
  {
    stop( paste('\n`', vars[wvars], '` is not in the data\n'))
  }

  # if you ever decide to make a numeric id, do so here, for now enforce numeric
  if(!is.numeric(data[[ids]]))
  {
    stop( paste('\n`', ids, '` must be numeric\n'))
  }

  # check correlation structure
  invisible( iscorStruct(correlation) )

  # see issue #12 on github
  # check that time is monotonically increasing - this fails for
  # time that is a trigonometric function of raw time
  #by(data[[time]], data[[ids]], function(x) all( diff(x) > 0) )

  # check for missing data and perform missing data handling.
  # I'm hesitant to do it here b/c it changes the data, do it on the fly
  # in the fitting methods

  # check variance - note that this also must be done for each ids in loops,
  # but follow curelator analyses and return per-person errors in output rather
  # than stopping analyses
  novar <- lapply(data[,unlist(vars[2:length(vars)])], function(x) all(duplicated(x)[-1L]) )
  novar <- unlist(novar)
  if( any(novar) )
  {
    stop( paste('\n`', names(novar)[novar], '` has zero variance' ))
  }

  # standardize the data
  if(standardize)
  {
    # only standardize dv if it is normal
    if( !is.null(dv) &
        (is.null(family) | as.character(family)[1]=="c(\"NO\", \"Normal\")") )
    {
      data[[dv]] <- scale(data[[dv]])
    }
    if( !is.null(dvs) & length(dvs) > 0 )
    {
      for(i in 1:length(dvs)) data[[dvs[[i]]]] <- scale( data[[dvs[[i]]]] )
    }
    if( !is.null(ivs) & length(ivs) > 0 )
    {
      for(i in 1:length(ivs))
      {
        if(!is.factor(data[[ivs[[i]]]])) data[[ivs[[i]]]] <- scale( data[[ivs[[i]]]] )
      }
    }
    if( !is.null(target_ivs) & length(target_ivs) > 0 )
    {
      for(i in 1:length(target_ivs))
      {
        if(!is.factor(data[[target_ivs[[i]]]])) data[[target_ivs[[i]]]] <- scale( data[[target_ivs[[i]]]] )
      }
    }
  }

  # redundant with monotonic(), clear this out
  # check whether any time points are duplicated
  ####dupTime <- lapply(by(data[[time]], data[[ids]], duplicated), any)
  ####dupTimem <- list()
  ####dupTimem[[ids]] <- names(dupTime)
  ####dupTimem[['isDuplicated']] <- unlist( dupTime )
  ####if(any(dupTimem$isDuplicated))
  ####{
  ####  warning(paste('\nThe following have duplicated values for', time, ':\n\n'),
  ####          paste(paste(ids, dupTimem[[ids]][dupTimem$isDuplicated]), collapse='\n'),
  ####          '\n\n1. Check that this is intentional. Duplicated time points\n',
  ####          'should only occur when ',
  ####          paste(time, 'is a trigonometric function of actual time.'),
  ####          '\n2. Data cannot be sorted by PersonAlytic and the user should ensure that\n',
  ####          paste('The data are sorted by', ids, 'then by', time, '. '),
  ####          'Failure to do this can cause \nestimation problems and invalidate ',
  ####          'graphical output.')
  ####  sortData <- FALSE
  ####}

  # sort the data
  # -- this is why ids must be numeric and
  # -- time must be monotonically increasing
  # this doesn't work for trigonometric functions of time
  if(sortData)
  {
    data <- data[order(data[[ids]], data[[time]]),]
  }
  # data <- data[order(data[[ids]]), ] # see issue #12 on github

  # if any cases have only 1 observation, drop them
  tuid <- table(data[[ids]])
  tofew <- names(tuid)[ tuid<2 ]
  if(length(tofew)>0)
  {
    warning(paste('\nThe following have less than 2 observations and will be dropped:\n\n'),
            paste(paste(ids, tofew), collapse='\n'))
    data <- data[!data[[ids]]%in%as.numeric(tofew),]
  }

  # align the data
  if(alignPhase & !is.null(phase))
  {
    data <- alignPhases(dat = data, id = ids, phase = phase, time = time)
  }

  return(data)
}


#' align the data at the transition between the first and second phase
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param id See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param do.plot Logical. Should the resulting data be plotted?
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

alignPhases <- function(dat, id, phase, time, do.plot=FALSE)
{
  time.old <- dat[,time]

  ### only align if phase order is the same for everyone
  phase.order <- aggregate(dat[,time], by=list(dat[,id], dat[,phase]), min, na.rm=TRUE)
  phase.order <- aggregate(phase.order$x, by=list(phase.order$Group.1), order)
  phase.is.ordered <- all(phase.order$x[1,1]==phase.order$x[,1])

  u.id <- unique(dat[,id])

  ### order phases if they are unordered
  if(!phase.is.ordered)
  {
    phase.Levels <- levels(as.factor(dat[,phase]))
    ### this automatically takes the levels alphabetically except for the
    #   selected phase.Level; user input will be needed in the future
    phase.Levels <- c(phase.Level, phase.Levels[phase.Levels!=phase.Level])
    min.times <- max.times <- vector("list", 0)

    ### rescale in.traphase times to start at 0
    for(i in 1:length(phase.Levels))
    {
      temp      <- dat[dat[,phase]==phase.Levels[i], c(id,time)]
      min.times[[i]] <- aggregate(temp[,time], list(temp[,id]), min, na.rm=TRUE)
      max.times[[i]] <- aggregate(temp[,time], list(temp[,id]), max, na.rm=TRUE)
      for(j in 1:length(u.id))
      {
        temp.time <- temp[temp[,id]==u.id[j],time]
        temp.time <- temp.time - min.times[[i]]$x[min.times[[i]]$Group.1==u.id[j]]
        temp[temp[,id]==u.id[j], time] <- temp.time
        rm(temp.time)
      }
      dat[dat[,phase]==phase.Levels[i], time] <- temp[,time]
    }

    ### reaccumulate time for the phase order in phase.Levels
    for(i in 2:length(phase.Levels))
    {
      ### get the interphase gap times
      max.pmi <- cbind(min.times[[i-1]]$Group.1, apply(cbind(min.times[[i-1]]$x,
                                                             min.times[[i]]$x),1,max))
      min.pmi <- cbind(max.times[[i-1]]$Group.1, apply(cbind(max.times[[i-1]]$x,
                                                             max.times[[i]]$x),1,min))
      dif.pmi <- data.frame(id=min.times[[i-1]]$Group.1, diff=max.pmi[,2] - min.pmi[,2])

      ### accumulate
      temp_im1  <- dat[dat[,phase]==phase.Levels[i-1], c(id,time)]
      max.times <- aggregate(temp_im1[,time], list(temp_im1[,id]), max, na.rm=TRUE)
      temp      <- dat[dat[,phase]==phase.Levels[i], c(id,time)]
      min.times <- aggregate(temp[,time], list(temp[,id]), min, na.rm=TRUE)
      for(j in 1:length(u.id))
      {
        temp.time <- temp[temp[,id]==u.id[j],]
        temp.time <- temp.time + (max.times$x[max.times$Group.1==u.id[j]] +
                                    dif.pmi$diff[dif.pmi$id==u.id[j]])
        temp[temp[,id]==u.id[j], time] <- temp.time[,time]
        rm(temp.time)
      }
      dat[dat[,phase]==phase.Levels[i], time] <- temp[,time]
    }
    ### re-sort the data
    dat <- dat[order(dat[,id], dat[,time]),]
  }

  ### alignment to the first phase transition
  if(do.plot)
  {
    par(mfrow=c(2,2))
    hist(dat[,time])
    time.old <- dat[,time]
  }
  for(i in u.id)
  {
    w       <- dat[,id]==i
    phase.i <- factor(dat[w,phase])
    if(is.numeric(phase.i)) phase.t <- table(phase.i!=min(phase.i))
    if(is.factor(phase.i))  phase.t <- table(phase.i)
    time.i  <- dat[w, time]
    time.p  <- time.i - time.i[ phase.t[1]+1 ]
    #time.p2  <- -(phase.t[1]):(phase.t[2]-1)
    dat[w,time] <- time.p
    #cbind(time.i, time.p, time.p2, phase.i)
  }
  if(do.plot)
  {
    plot(time.old, dat[,time])
    hist(dat[,time])
    par(mfrow=c(1,1))
  }


  return(dat)
}
