#' Function to create dummy codes for phase and group specific effects.
#'
#' In ICT data with 2+ groups and 2+ phases, dummy coding can be used to
#' estimate group and phase specific intercepts and slopes.
#'
#' @export
#'
#' @param data A \code{data.frame} with a time and phase variables
#' (see \code{\link{PersonAlytic}}), and a grouping variable (which would
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
#' # create group and phase specific dummy indicators after creating an
#' # artificial "Group" variable in the Ovary data.
#' OvaryICT2 <- OvaryICT
#' OvaryICT2$Group <- as.numeric(OvaryICT2$Mare<7)
#' OvaryICT2 <- byPhasebyGroup(OvaryICT2, "Time", "Phase", "Group")
#'
#' # view the new variables, group and phase specific intercepts and slopes
#' OvaryICT2$dummyNames
#'
#' # view the right hand side formula of the fixed effects produced by
#' # `byPhasebyGroup`, noting the `-1` elimination of an overall intercept
#' OvaryICT2$fixed
#'
#' # fit a model using the new variables with
#' # - no intercept
#' # - no overall time variable
#' userFormula <- list()
#' userFormula$fixed <- formula(paste("follicles ~ ", OvaryICT2$fixed)) # add dependent variable
#' userFormula$random <- ~ Time | Mare
#'
#' t0 <- PersonAlytic(output      = 'byPhasebyGroup'     ,
#'                    data        = OvaryICT2$data       ,
#'                    ids         = "Mare"               ,
#'                    dvs         = "follicles"          ,
#'                    time        = "Time"               ,
#'                    ivs         = OvaryICT2$dummyNames ,
#'                    package     = "nlme"               ,
#'                    userFormula = userFormula          ,
#'                    correlation = "corARMA(p=1)"       ,
#'                    autoSelect  = list()               ,
#'                    method      = "ML"                 )
#' summary(t0)
#'
#' # note that the following moves the `-1` to the end of the formula, test
#' # that this still works
#' userFormula$fixed <- update(formula(paste(" ~ ", OvaryICT2$fixed)), follicles ~ .)
#' t1 <- PersonAlytic(output      = 'byPhasebyGroup'     ,
#'                    data        = OvaryICT2$data       ,
#'                    ids         = "Mare"               ,
#'                    dvs         = "follicles"          ,
#'                    time        = "Time"               ,
#'                    ivs         = OvaryICT2$dummyNames ,
#'                    package     = "nlme"               ,
#'                    userFormula = userFormula          ,
#'                    correlation = "corARMA(p=1)"       ,
#'                    autoSelect  = list()               ,
#'                    method      = "ML"                 )
#' summary(t1)
#'
#' # if nothing returns, the test passes
#' testthat::expect_equal(summary(t0)$tTable, summary(t1)$tTable)
#'
#' \dontrun{
#'
#' # multiple DVs with no target_ivs and group models
#' # this fails to implement the userFormula, see https://github.com/ICTatRTI/PersonAlytics/issues/25
#' t2 <- PersonAlytic(output          = 'MultiDVnoIDnoIV'          ,
#'                    data            = OvaryICT2$data             ,
#'                    ids             = "Mare"                     ,
#'                    dvs             = names(OvaryICT)[c(3,9:11)] ,
#'                    ivs             = OvaryICT2$dummyNames       ,
#'                    time            = "Time"                     ,
#'                    package         = "nlme"                     ,
#'                    individual_mods = FALSE                      ,
#'                    target_ivs      = NULL                       ,
#'                    autoSelect      = list()                     ,
#'                    userFormula     = userFormula                )
#'
#'
#' }
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

  # create a data frame with the new variables
  newdata <- data.frame(data, do.call(data.frame, dummies))

  # create a fixed effects formula with no intercept
  fixed <- paste("- 1 + ", paste(names(dummies), collapse = "+"))


  list(data=newdata, dummyNames = names(dummies), fixed = fixed)

}






