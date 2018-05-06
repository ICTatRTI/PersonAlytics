#' function for hight throughput (htp)
#'


.htp <- function(data, dvs, phase, ids, uid, time, ivs, ivsl)
{
  # parralelization
  funcs <- c('eds')
  .foreach(funcs)
  pkgs  <- c('gamlss', 'nlme')

  #

  # outer loop for dependent variables
  DVout <- list()
  for(dv in dims$DV)
  {
    # middle loop for covariates
    IVout <- list()
    for(iv in dims$IV)
    {
      # set up Palytic object//needs editing for covariates & situation
      # with no phase
      t1 <- Palytic$new(data=data,
                        fixed  = formula(paste(
                                         paste(dvs[dv], '~', time, '*', phase),
                                         paste(ivs, collapse = '+'),
                                         ivsl[[iv]],
                                         collapse = '+')),
                        random = formula(paste('~', time, '|', ids)) )

      # inner loop for participants
      IDout <- foreach::foreach(id=dims$ID, .packages = pkgs) %dopar%
      {
        t1$gamlss( which(data[[ids]]==id) )
        # post processing//do you want it here? as a method for Palytic objects?
        # note that a method doesn't work (yet) b/c we're not saving the gamlss
        # output the the object t1. We'll need to do this if we want to post-
        # process on results from the gamlss object.
      }
      # post-estimation aggregation
      names(IDout) <- paste(ids, uid, sep='.')
      IDout1 <- lapply(IDout, function(x) summary(x))
      #IDout1 <- lapply(IDout1, functionx(x) x$id=) # fix names
      IDout1 <- do.call(rbind, IDout1)

      IVout[[ifelse(is.null(ivs[iv]), 1, ivs[iv])]] <- IDout1 # aggregated IDout
    }
    DVout[[dvs[dv]]] <- IVout
  }

  # final post-processing
  return(DVout)

}
