#'
#'
#'

PersonAlytic <- function(data=NULL,
                         ids,
                         dvs,
                         phase,
                         time,
                         ivs=NULL,
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
    dvs    <- 'follicles'
    phase  <- 'Phase'
    ids     <- 'Mare'
    time   <- 'TimeSin'
  }

  # unique idss
  uids <- unique(data[[ids]])

  # dimensions for loops
  ID <- 1:length(uid)
  IV <- 1:length(ivs); if(is.null(ivs)) IV <- 1
  DV <- 1:length(dvs)

  dims <- list(ID=ID, IV=IV, DV=DV)

  if( any( lapply(dims, length) > 1) )
  {
    .htp()
  }

}
