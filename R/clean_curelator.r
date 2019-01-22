#' clean_curelator, a function to clean Curelator input dataa.
#'
#' @export
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @description To be moved to package with Curelator specific tools (this may be
#' the only thing in it).
#'
#' @param data A \code{\link{data.frame}} with ---specs here---.
#'
#'
clean_curelator <- function(data)
{
  if(class(data)!="data.frame") stop("clean_curelator: `data` must be a data.frame.")
  names(data) <- tolower(names(data))
  # atmospheric variables
  data$windgust      [data$windgust      >= 100] <- NA
  data$windspeedmax  [data$windspeedmax  >= 100] <- NA
  data$windspeedmean [data$windspeedmean >= 100] <- NA

  data$humiditychange    <- data$humiditymax - data$humiditymin
  data$pressurechange    <- data$pressuremax - data$pressuremin
  data$temperaturechange <- data$temperaturemax - data$temperaturemin
  data$windspeedchange   <- data$windspeedmax - data$windspeedmin

  ### fix fauxlogicals
  for(i in seq_along(data))
  {
    if(all(unique(data[,i])%in%c('False', 'True', NA))) data[,i] <- as.logical(data[,i])
  }

  ### create composites using mean imputation
  data$sensitivities <- apply(data[,c('brightlights', 'loudnoise', 'odours',
                                    'skinsensitivity', 'eyestrain')], 1, mean, na.rm=TRUE)*5

  data$totalactivity <- apply(data[,c('activityintense',
                                    'activitymoderate')], 1, mean, na.rm=TRUE)*2

  data$totalmealsmissed <- apply(data[,c('missedbreakfast', 'missedlunch',
                                       'misseddinner', 'missedothermeals')],
                                1, mean, na.rm=TRUE)*4

  data$totalalcohol <- apply(data[,c('beer', 'redwine', 'whitewine',
                                   'sparklingwine', 'spirits')], 1, mean, na.rm=TRUE)*5

  data$totalcaffeine <- apply(data[,c('coffee', 'tea', 'softdrink',
                                    'energydrink')], 1, mean, na.rm=TRUE)*4

  data$totalnicotine <- apply(data[, c('nicotine', 'nicotineother')],
                             1, mean, na.rm=TRUE)*2

  data$totaltravel <- apply(data[,c('travelcar', 'travelplane', 'traveltrain',
                                  'travelmotorcycle', 'travelship')], 1, mean, na.rm=TRUE)*5

  data$physicalstrains <- apply(data[,c('tiredness', 'neckpain', 'yawning',
                                      'poorconcentration', 'hunger')], 1, mean, na.rm=TRUE)*5

  data$severity <- rep(0, nrow(data))
  data$severity[is.na(data$peakseverity) | is.na(data$headachetype)] <- NA
  data$severity[data$peakseverity=="mild"] <- 1
  data$severity[data$peakseverity=="moderate"] <- 2
  data$severity[data$peakseverity=="severe"] <- 3

  data$severity.migraine <- data$severity
  data$severity.migraine[data$headachetype=='headache' | is.na(data$headachetype)] <- NA

  data$severity.nonmigraine <- data$severity
  data$severity.nonmigraine[data$headachetype=='migraine' | is.na(data$headachetype)] <- NA

  patientInclude<-data$patient[which(data$daystracked==90 & data$day<=120)]
  data<-data[is.element(data$patient,patientInclude),]
  data<-data[data$daystracked<=90,]

  return(data)
}
