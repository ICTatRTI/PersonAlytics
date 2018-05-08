#' Personalytic
#'
#' @param data
#' @param ids Character. Name of the ID variable. If left as NULL, an example will run using the Ovary data from the nlme package.
#' @param dvs Character list of dependent variable name(s).
#' @param phase Charcter. Name of the phase variable.
#' @param time Character. Name of the time variable.
#' @param ivs Character list of independent variables. These will all be used in every analysis.
#' @param ivsl Character list of independent variables. These will be looped through one at a time.
#' @param ind.mods Logical. Should individual models be fit for each ID?
#' @param grp.mod Logical. Should a group level model be fit across all IDs?

PersonAlytic <- function(data=NULL,
                         ids,
                         dvs,
                         phase,
                         time,
                         ivs=NULL,
                         ivsl=NULL,
                         ind.mods=TRUE,
                         grp.mod =FALSE)
{
  # if no data are given, use a test data set
  if(is.null(data))
  {
    data <- as.data.frame(nlme::Ovary)
    data$Mare <- factor(data$Mare, ordered = FALSE)
    data$Phase <- as.numeric(data$Time > .5)
    data$TimeSin <- sin(2*pi*data$Time)
    dvs    <- "follicles"
    phase  <- "Phase"
    ids    <- "Mare"
    time   <- "TimeSin"
  }

  # unique ids
  uids <- unique(data[[ids]])

  # dimensions for loops
  ID <- 1:length(uid)
  IV <- 1:length(ivsl); if(is.null(ivsl)) IV <- 1
  DV <- 1:length(dvs)
  dims <- list(ID=ID, IV=IV, DV=DV)

  if( !any( lapply(dims, length) > 1) | grp.mod )
  {
    t1 <- Palytic$new(data=data,
                      fixed  = formula(paste(
                        paste(dvs, "~", time, "*", phase),
                        paste(ivs, collapse = "+"),
                        collapse = "+")),
                      random = formula(paste("~", time, "|", ids)) )
    Grp.out <- t1$gamlss()
  }
  else DVout <- .htp(data, dvs, phase, ids, uid, time, ivs, ivsl)



}
