#' clean, a function to check data and inputs.
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#' @param dv See \code{\link{Palytic}}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param ivs See \code{PersonAlytic}.
#' @param fixed See \code{\link{Palytic}}.
#' @param random See \code{\link{Palytic}}.
#' @param formula See \code{\link{Palytic}}.
#' @param correlation See \code{\link{PersonAlytic}}.
#' @param family See \code{\link{PersonAlytic}}.
#' @param dvs See \code{\link{PersonAlytic}}.
#' @param target_ivs See \code{PersonAlytic}.
#' @param standardize Logical. Should all variables be standardized? Only applies
#' to \code{dvs}, \code{ivs}, and \code{target_ivs}.
#' @param alignPhase See \code{\link{PersonAlytic}}.
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

clean <- function(data		 	                                            ,
                  ids		 	                                              ,
                  dv		 	                                              ,
                  time		 	                                            ,
                  phase	        =	NULL	                                ,
                  ivs	          =	NULL	                                ,
                  fixed       	=	NULL	                                ,
                  random      	=	NULL	                                ,
                  formula	      =	NULL	                                ,
                  correlation 	=	NULL	                                ,
                  family        = NO()                                  ,
                  dvs	          =	NULL	                                ,
                  target_ivs    =	NULL	                                ,
                  standardize 	=	list(dv=FALSE, iv=FALSE, byids=FALSE)	,
                  sortData	    =	TRUE	                                ,
                  alignPhase 	  =	"none"                                ,
                  debugforeach  = debugforeach                          )
{
  # check that the family is a gamlss.family object
  if(! "gamlss.family" %in% class(family))
  {
    stop("\n`family`=", family, " is not a 'gamlss.family' family.")
  }

  # check that variables are in the data set
  vars <- unique( c(ids, dv, time$raw, phase, unlist(ivs), unlist(dvs),
                    unlist(target_ivs), all.vars(fixed), all.vars(random),
                    all.vars(formula)) )
  vars <- vars[which(vars!='1' & vars!='0')]
  #vars <- gsub(" ", "", vars) # this will be a problem if var names have spaces...
  wvars <- which( ! vars %in% names(data) )
  if( length(wvars) > 0 )
  {
    stop( paste('\n`', vars[wvars], '` is not in the data\n'))
  }

  # check id variable
  data[[ids]] <- checkID(data[[ids]], ids)

  # check phase variable
  if(!is.null(phase))
  {
    if(any(is.na(data[[phase]]))) data[[phase]] <- fixPhaseNA(data[[phase]])
  }

  # check correlation structure
  invisible( iscorStruct(correlation) )

  # see issue #12 on github
  # check that time is monotonically increasing - this fails for
  # time that is a trigonometric function of raw time
  #by(data[[time$raw]], data[[ids]], function(x) all( diff(x) > 0) )

  # check for missing data and perform missing data handling.
  # I'm hesitant to do it here b/c it changes the data, do it on the fly
  # in the fitting methods

  # check variance - note that this also must be done for each ids in loops,
  # but follow curelator analyses and return per-person errors in output rather
  # than stopping analyses
  novar <- lapply(data[,unlist(vars[2:length(vars)])],
                  function(x) all(duplicated(x)[-1L]) )
  novar <- unlist(novar)
  if( any(novar) )
  {
    stop( paste('\n`', names(novar)[novar], '` has zero variance' ))
  }

  # Standardize the data. Note that dvs and target_ivs (as part of ivs) are
  # resubmitted to clean() for each Palytic object created in htp, so we don't
  # need to standardize them here
  if(!is.list(standardize))
  {
    stop("`standardize=", standardize, "`,\n",
      "but `standardize` should be a named logical list with at least one of ",
      "`dvs`, `ivs`, or `byids`. For example, `standardize=list(dvs=FALSE,",
      "ivs=TRUE,byids=TRUE)`.")
  }
  data <- pstand(data, standardize, dv, ivs, family, ids)

  # sort the data
  # -- this is why ids must be numeric and
  # -- time must be monotonically increasing
  # this doesn't work for trigonometric functions of time
  if(sortData)
  {
    if(debugforeach)
    {
      print(time)
    }
    data <- data[order(data[[ids]], data[[time$raw]]),]
  }
  # data <- data[order(data[[ids]]), ] # see issue #12 on github

  # if any cases have only 1 observation, drop them
  tuid <- table(data[[ids]])
  toofew <- names(tuid)[ tuid<2 ]
  if(length(toofew)>0)
  {
    warning(paste('\nThe following have less than 2 observations and will be dropped:\n\n'),
            paste(paste(ids, toofew), collapse='\n'))
    data <- data[!data[[ids]]%in%as.numeric(toofew),]
  }

  # align the data
  if(!is.null(phase))
  {
    if(alignPhase == 'align') data <- alignPhases(dat = data, id = ids,
                                      phase = phase, time = time$raw)
    if(alignPhase == 'piecewise')
    {
      data <- data.frame(data, pwtime(time = data[[time$raw]],
                                      phase = data[[phase]])
      )
    }
  }

  return(data)
}

#' fixPhaseNA - function to carry last value forward if phase is missing
#' (this is common in example data sets Ty has worked with)
#'
#' @param phase The \code{phase} variable
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
fixPhaseNA <- function(phase)
{
  for(i in 2:length(phase))
  {
    if(is.na(phase[i])) phase[i] <- phase[i-1]
  }
  phase
}

#' Function to check ID variable and if possible, coerce to numeric
#'
#' @param x The ID variable values as a vector.
#' @param ids The ID variable name. See \code{\link{PersonAlytic}}.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
checkID <- function(x, ids)
{
  if(!is.numeric(x))
  {
    if(is.character(x))
    {
      # length unique x
      lux  <- length(unique(x))
      # length unique numeric x
      luxn <- length(unique(as.numeric(x)))
      if(lux==lunx)
      {
        warning('`', ids, '` is character but will be forced to numeric')
        return(as.numeric(x))
      }
      if(lux==lunx)
      {
        stop('`', ids, '` is character but cannot be forced to numeric.',
             ' For example, id="1.0" and id="1.00" are not unique after',
             ' conversion to numeric.')
      }
    }
    if(!is.character(x))
    {
      stop( paste('\n`', ids, '` must be numeric\n', sep='') )
    }
  }
  if(is.numeric(x)) return(x)
}

#' Function to standardize data with options from the \code{standardize} parameter.
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param standardize See \code{\link{PersonAlytic}}.
#' @param dv See \code{\link{PersonAlytic}}.
#' @param ivs See \code{\link{PersonAlytic}}.
#' @param family See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
pstand <- function(data, standardize, dv, ivs, family, ids)
{
  # determine which variables to standardize
  dostand <- list()

  # only standardize dv if requested and if dv is normal
  if( !is.null(standardize$dv) &
      (is.null(family) | as.character(family)[1]=="c(\"NO\", \"Normal\")") )
  {
    if(standardize$dv) dostand$dv=dv
  }

  # ivs
  if( !is.null(standardize$iv) & length(ivs) > 0 )
  {
    ivs_stand <- list()
    for(i in seq_along(ivs))
    {
      if(is.numeric(data[[ivs[[i]]]])) ivs_stand[[i]] <- ivs[[i]]
    }
    if(standardize$iv) dostand$iv=ivs_stand
  }

  # turn dostand into a character vector
  dostand <- unlist(dostand)

  # group standardization
  dogroup <- FALSE
  if(!is.null(standardize$byids))
  {
    if(standardize$byids) dogroup <- TRUE
  }
  for(i in dostand)
  {
    if(!dogroup) data[[i]] <- scale(data[[i]])
    if( dogroup)
    {
      data[[i]] <- unlist( by(as.numeric( data[[i]] ), data[[ids]], scale))
    }
  }

  return( data )

}


#' align the time variable to be zero at the transition between
#' the first and second phase
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param id See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param do.plot Logical. Should the resulting data be plotted?
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export

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
    for(i in seq_along(phase.Levels))
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
    for(i in seq_along(phase.Levels)[-1])
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
