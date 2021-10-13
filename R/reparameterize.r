#' Function to reparameterize ICT data with 2+ groups and 2+ phase
#'
#' @export
#'
#' @param data A \code{data.frame} with a time and phase variables
#' (see \link{\code{PersonAlytic}}), and a grouping variable (which would
#' normally be passed to \code{PersonAlytic} via the \code{ivs} parameter).
#'
#' @param time Character. The name of the time variable.
#'
#' @param phase Character. The name of the phase variable.
#'
#' @param group Character. The name of the group variable.
#'
#' @examples
#'
byPhasebyGroup <- function(data, time, phase, group)
{
  g <- unique(data[[group]])
  p <- unique(data[[phase]])

  dummies <- list()
  for(i in seq_along(g))
  {
    for(j in seq_along(p))
    {
      term <- paste("g", g[i], "p", p[j], sep="")

      wg <- data[[group]] == g[i]
      wp <- data[[phase]] == p[j]

      # slope first so it is easy to drop the first one as referent
      slp  <- paste(term, "slope", sep="_")
      dummies[[slp]] <- as.numeric(wg & wp) * data[[time]]

      # intercepts next
      int  <- paste(term, "int", sep="_")
      dummies[[int]] <- as.numeric(wg & wp)
    }
  }

  newdata <- data.frame(data, do.call(data.frame, dummies))

  list(data=newdata, dummyNames = names(dummies))

}






