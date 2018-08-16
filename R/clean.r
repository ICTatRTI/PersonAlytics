#' clean, a function to check data and inputs.
#'
#' @param data See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#' @param dv See \code{\link{PersonAlytic}}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param ivs Used by \code{PersonAlyticPro}.
#' @param dvs Used by \code{PersonAlyticPro}.
#' @param ivsl Used by \code{PersonAlyticPro}.
#' @param fixed See \code{\link{Palytic}}.
#' @param random See \code{\link{Palytic}}.
#' @param formula See \code{\link{Palytic}}.
#' @param standardize Logical. Should all variables be standardized? Only applies
#' to \code{dv}, \code{ivs}, and \code{ivsl}.
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


# consider moving to .utils unless it gets too big
# consider adding this to all active bindings, i.e., recheck the data if a change is made, or just leave it to the active binding itself
clean <- function(data, ids, dv, time, phase=NULL, ivs=NULL,
                  fixed=NULL, random=NULL, formula=NULL, correlation=NULL,
                  dvs=NULL, ivsl=NULL, standardize=FALSE, sortData=TRUE)
{
  # check that variables are in the data set
  vars <- unique( c(ids, dv, time, phase, unlist(ivs), unlist(dvs), unlist(ivsl),
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
    if( !is.null(dv)  ) data[[dv]] <- scale(data[[dv]])
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
    if( !is.null(ivsl) & length(ivsl) > 0 )
    {
      for(i in 1:length(ivsl))
      {
        if(!is.factor(data[[ivsl[[i]]]])) data[[ivsl[[i]]]] <- scale( data[[ivsl[[i]]]] )
      }
    }
  }

  # redundant with monotonic(), clear this out
  # check whether any time points are duplicated
  dupTime <- lapply(by(data[[time]], data[[ids]], duplicated), any)
  dupTimem <- list()
  dupTimem[[ids]] <- names(dupTime)
  dupTimem[['isDuplicated']] <- unlist( dupTime )
  if(any(dupTimem$isDuplicated))
  {
    warning(paste('\nThe following have duplicated values for', time, ':\n\n'),
            paste(paste(ids, dupTimem[[ids]][dupTimem$isDuplicated]), collapse='\n'),
            '\n\n1. Check that this is intentional. Duplicated time points\n',
            'should only occur when ',
            paste(time, 'is a trigonometric function of actual time.'),
            '\n2. Data cannot be sorted by PersonAlytic and the user should ensure that\n',
            paste('The data are sorted by', ids, 'then by', time, '. '),
            'Failure to do this can cause \nestimation problems and invalidate ',
            'graphical output.')
    sortData <- FALSE
  }

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

  return(data)
}
