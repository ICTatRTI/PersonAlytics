#' Citations for PersonAlytics
#'
#' @param load Logical. Should packages be loaded? If \code{FALSE} (the default) citations are returned as a list. Called internally by [link to function].
#'
#' @examples
#' getCitations()
#'

getCitations <- function(load=FALSE)
{
  packages <- list(
    pbkrtest     = "pbkrtest"     ,
    lme4         = "lme4"         ,
    nlme         = "nlme"         ,
    ggplot2      = "ggplot2"      ,
    gridExtra    = "gridExtra"    ,
    RColorBrewer = "RColorBrewer" ,
    Rmisc        = "Rmisc"        ,
    Hmisc        = "Hmisc"        ,
    mgcv         = "mgcv"
    )
  if(load)
  {
    new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(unlist(new.packages))
    lapply(packages, require, character.only = TRUE)
  }
  if(!load)
  {
    return(lapply(packages, citation))
  }
}
