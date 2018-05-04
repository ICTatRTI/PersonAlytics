#' .foreach Utility to parallelize loops
#'
#'
#'
#'

.foreach <- function(funcs)
{
  cl     <- snow::makeCluster(parallel::detectCores(), type="SOCK")
  snow::clusterExport(cl, funcs)
  doSNOW::registerDoSNOW(cl)
}

