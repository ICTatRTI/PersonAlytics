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
#' Alternatively, \code{data=NULL} can be used in conjunction with \code{time},
#' \code{phase}, and \code{group} to produce generic variable names and a
#' design matrix.
#'
#' @param time Character. The name of the time variable, or if \code{data=NULL},
#' \code{time} must be a numeric vector with all of time points in a
#' hypothetical study.
#'
#' @param phase Character. The name of the phase variable, or if \code{data=NULL},
#' \code{phase} must be either a numeric integer specifying the number of phases
#' or a character vector with the phase names.
#'
#' @param group Character. The name of the group variable., or if \code{data=NULL},
#' \code{group} must be either a numeric integer specifying the number of groups
#' or a character vector with the group names.
#'
#' @param retime Logical. Default is \code{TRUE}. Should the time (slope)
#' variables be rescaled within each phase and group to go from 0 to the
#' number of time points minus 1.
#'
#' @examples
#'
#' # produce a design matrix for a hypothetical study with 20 time points,
#' # 2 phases, and 3 groups
#' studyDesign <- byPhasebyGroup(NULL, -9:10, 2, 3)
#' head(studyDesign$data)
#' studyDesign$dummyNames
#'
#' # an alternative specification
#' studyDesign <- byPhasebyGroup(NULL, -9:10, c("baseline", "intervention"),
#'                  c("treatment", "control"))
#' head(studyDesign$data)
#' studyDesign$dummyNames
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

byPhasebyGroup <- function(data, time, phase, group, retime = TRUE)
{
  if(is.null(data))
  {
    if(length(phase)==1 & length(group)==1)
    {
      if(!is.numeric(phase) | !is.numeric(group))
      {
        stop("\n`phase` and `group` must be numeric integers, currently:\n\n",
             paste(paste(c("phase", "group"), c(phase, group), sep=" = "),
                   collapse = "\n"))
      }
      phase <- 1:phase
      group <- 1:group
    }
    data <- expand.grid(time=time, phase=phase, group=group)
    time <- "time"
    phase <- "phase"
    group <- "group"
  }

  g <- unique(data[[group]])
  p <- unique(data[[phase]])

  gpref <- ""
  ppref <- ""

  if(is.numeric(g)) gpref <- "g"
  if(is.numeric(p)) ppref <- "p"

  dummies <- list()
  for(i in seq_along(g))
  {
    for(j in seq_along(p))
    {
      term <- paste(gpref, g[i], "_", ppref, p[j], sep="")

      wg <- data[[group]] == g[i]
      wp <- data[[phase]] == p[j]

      # intercepts
      int  <- paste(term, "int", sep="_")
      dummies[[int]] <- as.numeric(wg & wp)

      # slopes
      slp  <- paste(term, "slope", sep="_")
      dummies[[slp]] <- as.numeric(wg & wp) * data[[time]]
    }
  }

  # rescale the slope variables to have min(x)==0
  if(retime) dummies <- retimes(dummies)

  # create a data frame with the new variables
  newdata <- data.frame(data, do.call(data.frame, dummies))

  # create a fixed effects formula with no intercept
  fixed <- paste("- 1 + ", paste(names(dummies), collapse = "+"))


  list(data=newdata, dummyNames = names(dummies), fixed = fixed)

}

#' retime
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @keywords internal
retime <- function(x=30:45)
{
  xnot0 <- which(x!=0)
  xno0  <- x[xnot0] - min(x[xnot0], na.rm = T)
  x[xnot0] <- xno0
  x
}

#' retimes
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @keywords internal
retimes <- function(dummies)
{
  wslopes <- names(dummies)[grepl("_slope", names(dummies))]
  for(i in wslopes)
  {
    dummies[[i]] <- retime(dummies[[i]])
  }
  dummies
}


#' Function to create a formula from a list of group and phase specific
#' intercepts and slopes
#'
#' used by the PersonAlyticsPower GUI
#'
#' @export
#'
#' @param terms Character vector or list.
#'
#' @param random A random effects formula. The default is \code{~ Time | id} as
#' used by PersonAlyticsPower
#'
#' @examples
#'
#' terms <- c("group1_phase1_int",
#'            "group1_phase1_slope",
#'            "group1_phase2_int",
#'            "group1_phase2_slope",
#'            "group1_phase3_int",
#'            "group1_phase3_slope",
#'            "group2_phase1_int",
#'            "group2_phase1_slope",
#'            "group2_phase2_int",
#'            "group2_phase2_slope",
#'            "group2_phase3_int",
#'            "group2_phase3_slope")
#' termsToFormula(terms)
#'
termsToFormula <- function(terms, random = ~Time | id)
{
  if(length(terms) > 0)
  {
    frm <- paste(unlist(terms), collapse = "+")
    frm <- paste("y1 ~ -1", frm, sep="+")
    return(list(fixed = formula(frm),
         random = random))
  }
  if(length(terms)==0)
  {
    return(list())
  }
}






