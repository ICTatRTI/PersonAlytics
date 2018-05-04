#'
#'
#'

PersonAlytic <- function(data=NULL,
                         id,
                         dvs,
                         phase,
                         time,
                         cvs=NULL,
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
    id     <- 'Mare'
    time   <- 'TimeSin'
  }

  # unique ids
  uid <- unique(data[[id]])

  # dimensions for loops
  P <- 1:length(uid)
  C <- 1:length(cvs); if(is.null(cvs)) C <- cvs
  D <- 1:length(dvs)



}
