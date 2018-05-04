#' function for hight throughput (htp)
#'


.htp <- function(data, dvs, phase, ids, time, ivs)
{
  # parralelization
  funcs <- c('eds')
  .foreach(funcs)
  pkgs  <- c('gamlss', 'nlme')

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
                        fixed  = formula(paste(dvs[dv], '~', time, '*', phase)),
                        random = formula(paste('~', time, '|', ids)) )

      # inner loop for participants
      IDout <- foreach::foreach(id=dims$ID, .packages = pkgs) %dopar%
      {
        t1$gamlss( which(data[[ids]]==id) )
        # post processing//do you want it hear? as a method for Palytic objects?
      }
      # post-estimation aggregation

      IVout[[ivs[iv]]] <- 0 # aggregated IDout
    }
    DVout[[dvs[dv]]] <- IVout
  }
  # final post-processing

}
