# TODO ####
# paralellize loops


# function to read WRAP data
# this should become it's own standalone package or be added to 
# PersonAlytics


library(rjson)
library(summarytools)
library(lubridate)
library(plyr)
library(tidyverse)
library(foreach)
library(doParallel)
library(ranger)
library(kableExtra)


#' Read and extract json WRAP data from multiple participants
#' 
#' @param file a vector json file names from WRAP
#' @param times see \code{\link(add_phase)}.
#' @param labels see \code{\link(add_phase)}.
read_wraps <- function(files, study_start=NULL, times=NULL, labels=NULL)
{
  # read the files
  wraps <- lapply(files, read_wrap, study_start=study_start)
  
  # extract into three data sets
  activity <- do.call(rbind, lapply(wraps, \(x) x$activity))
  stress <- rbind.fill(lapply(wraps, \(x) x$stress))
  sleep <- rbind.fill(lapply(wraps, \(x) x$sleep))
  
  # add phases
  phases <- NULL
  if(!is.null(times))
  {
    activity$phase <- add_phase(activity, times, labels)
    stress$phase <- add_phase(stress, times, labels)
    sleep$phase <- add_phase(sleep, times, labels)
    
    phases <- activity$startTimeInSeconds[which(!duplicated(activity$phase))]
  }
  
  # fix ID variables
  activity$id <- as.numeric(factor(activity$id))
  activity$idlot <- factor(activity$id) 
  
  stress$id <- as.numeric(factor(stress$id)) 
  stress$idlot <- factor(stress$id) 
  
  sleep$id <- as.numeric(factor(sleep$id)) 
  sleep$idlot <- factor(sleep$id) 
  
  # rescale phases
  study_start <- min(min(activity$startTimeInSeconds),
                     min(stress$startTimeInSeconds),
                     min(sleep$startTimeInSeconds))
  p2 <- list()
  p2$phases <- phases
  p2$datetime <- as_datetime(phases)
  p2$time_seconds <- as.numeric(phases-study_start)
  p2$time_minutes <- p2$time_seconds/60
  p2$time_hours <- p2$time_minutes/60
  p2$time_days <- p2$time_hours/24
  
  #
  list(activity=activity, stress=stress, sleep=sleep, phases=p2)
}

#' Read and extract json WRAP data
#' 
#' @param file a character file name for json WRAP data
#' @param study_start a start time for zeroing activity the time variables. 

read_wrap <- function(file, study_start=NULL)
{
  # read the json file
  dat <- fromJSON(file=file)
  
  # extract $data
  temp <- lapply(dat, \(x) x$data)
  resource <- unlist(lapply(dat, \(x) x$resource))
  
  # find the most common data length, ignore the rest for now
  lengths <- freq(unlist(lapply(temp, length)))
  the_length <- as.numeric(attr(lengths, "dimnames")[[1]][1])
  
  # get names for id variables
  nm <- unlist(strsplit(file, "/"))
  nm <- nm[length(nm)]
  nm <- unlist(strsplit(nm, "\\."))[1]
  
  # loop through TODO vectorize or parralelize this
  activity <- sleep <- stress <- list()
  for(i in seq_along(temp))
  {
    if(resource[i]=="epochs") activity[[i]] <- data.frame(temp[[i]]) 
    if(resource[i]=="dailies") 
    {
      temp_steps <- temp[[i]]
      if(length(temp_steps$timeOffsetHeartRateSamples)==0)
      {
        temp_steps$timeOffsetHeartRateSamples[["timeOffsetHeartRateSamples.0"]] <- NA
      }
      stress[[i]] <- data.frame(temp_steps)
      rm(temp_steps)
    }
    if(resource[i]=="sleeps")
    {
      temp_sleep <- lapply(temp[[i]], data.frame)
      for(j in seq_along(temp_sleep))
      {
        if(length(temp_sleep[[j]]==1))
        {
          names(temp_sleep[[j]]) <- names(temp_sleep)[j]
        }
      }
      temp_sleep <- temp_sleep[unlist(lapply(temp_sleep, \(x) length(x) > 0))]
      sleep[[i]] <- do.call(data.frame, temp_sleep) 
    }
  }
  activity <- activity[unlist(lapply(activity, \(x) length(x) > 0))]
  stress   <- stress[unlist(lapply(stress, \(x) length(x) > 0))]
  sleep    <- sleep[unlist(lapply(sleep, \(x) length(x) > 0))]
  
  activity <- rbind.fill(activity)
  stress   <- rbind.fill(stress)
  sleep    <- rbind.fill(sleep)
  
  # start time
  if(is.null(study_start)) study_start <- activity$startTimeInSeconds[1]
  
  # add time objects to activity
  activity$datetime <- as_datetime(activity$startTimeInSeconds)
  activity$period <- seconds_to_period(activity$datetime)
  activity$time_seconds <- as.numeric(activity$startTimeInSeconds-study_start)
  activity$time_minutes <- activity$time_seconds/60
  activity$time_hours <- activity$time_minutes/60
  activity$time_days <- activity$time_hours/24
  
  # add id to activity
  activity$id <- nm
  
  # add time objects to stress
  stress$datetime <- as_datetime(stress$startTimeInSeconds)
  stress$period <- seconds_to_period(stress$datetime)
  stress$time_seconds <- as.numeric(stress$startTimeInSeconds-study_start)
  stress$time_minutes <- stress$time_seconds/60
  stress$time_hours <- stress$time_minutes/60
  stress$time_days <- stress$time_hours/24
  
  # add id to stress
  stress$id <- nm
  
  # add time objects to sleep
  sleep$datetime <- as_datetime(sleep$startTimeInSeconds)
  sleep$period <- seconds_to_period(sleep$datetime)
  sleep$time_seconds <- as.numeric(sleep$startTimeInSeconds-study_start)
  sleep$time_minutes <- sleep$time_seconds/60
  sleep$time_hours <- sleep$time_minutes/60
  sleep$time_days <- sleep$time_hours/24
  
  # add id to sleep
  sleep$id <- nm
  
  # keep last
  activity <- activity[!duplicated(activity$startTimeInSeconds, fromLast = TRUE),]
  sleep <- sleep[!duplicated(sleep$startTimeInSeconds, fromLast = TRUE),]
  stress <- stress[!duplicated(stress$startTimeInSeconds, fromLast = TRUE),]
  
  # return the data
  list(activity=activity, sleep=sleep, stress=stress)
}

read_wrap_par <- function(file, study_start=NULL)
{
  # read the json file
  dat <- fromJSON(file=file)
  
  # extract $data
  temp <- lapply(dat, \(x) x$data)
  resource <- unlist(lapply(dat, \(x) x$resource))
  
  # find the most common data length, ignore the rest for now
  lengths <- freq(unlist(lapply(temp, length)))
  the_length <- as.numeric(attr(lengths, "dimnames")[[1]][1])
  
  # get names for id variables
  nm <- unlist(strsplit(file, "/"))
  nm <- nm[length(nm)]
  nm <- unlist(strsplit(nm, "\\."))[1]
  
  # loop through TODO vectorize or parralelize this
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  out <- foreach(i = seq_along(temp)) %dopar%
  {
    activity <- sleep <- stress <- list()
    if(resource[i]=="epochs") activity[[i]] <- data.frame(temp[[i]]) 
    if(resource[i]=="dailies") 
    {
      temp_steps <- temp[[i]]
      if(length(temp_steps$timeOffsetHeartRateSamples)==0)
      {
        temp_steps$timeOffsetHeartRateSamples[["timeOffsetHeartRateSamples.0"]] <- NA
      }
      stress[[i]] <- data.frame(temp_steps)
      rm(temp_steps)
    }
    if(resource[i]=="sleeps")
    {
      temp_sleep <- lapply(temp[[i]], data.frame)
      for(j in seq_along(temp_sleep))
      {
        if(length(temp_sleep[[j]]==1))
        {
          names(temp_sleep[[j]]) <- names(temp_sleep)[j]
        }
      }
      temp_sleep <- temp_sleep[unlist(lapply(temp_sleep, \(x) length(x) > 0))]
      sleep[[i]] <- do.call(data.frame, temp_sleep) 
    }
    list(activity = activity[[i]], 
         sleep    = sleep   [[i]],
         stress   = stress  [[i]])
  }
  activity <- lapply(out, \(x) x$activity)
  stress   <- lapply(out, \(x) x$stress  )
  sleep    <- lapply(out, \(x) x$sleep   )
  
  activity <- activity[unlist(lapply(activity, \(x) length(x) > 0))]
  stress   <- stress[unlist(lapply(stress, \(x) length(x) > 0))]
  sleep    <- sleep[unlist(lapply(sleep, \(x) length(x) > 0))]
  
  activity <- rbind.fill(activity)
  stress   <- rbind.fill(stress)
  sleep    <- rbind.fill(sleep)
  
  # start time
  if(is.null(study_start)) study_start <- activity$startTimeInSeconds[1]
  
  # add time objects to activity
  activity$datetime <- as_datetime(activity$startTimeInSeconds)
  activity$period <- seconds_to_period(activity$datetime)
  activity$time_seconds <- as.numeric(activity$startTimeInSeconds-study_start)
  activity$time_minutes <- activity$time_seconds/60
  activity$time_hours <- activity$time_minutes/60
  activity$time_days <- activity$time_hours/24
  
  # add id to activity
  activity$id <- nm
  
  # add time objects to stress
  stress$datetime <- as_datetime(stress$startTimeInSeconds)
  stress$period <- seconds_to_period(stress$datetime)
  stress$time_seconds <- as.numeric(stress$startTimeInSeconds-study_start)
  stress$time_minutes <- stress$time_seconds/60
  stress$time_hours <- stress$time_minutes/60
  stress$time_days <- stress$time_hours/24
  
  # add id to stress
  stress$id <- nm
  
  # add time objects to sleep
  sleep$datetime <- as_datetime(sleep$startTimeInSeconds)
  sleep$period <- seconds_to_period(sleep$datetime)
  sleep$time_seconds <- as.numeric(sleep$startTimeInSeconds-study_start)
  sleep$time_minutes <- sleep$time_seconds/60
  sleep$time_hours <- sleep$time_minutes/60
  sleep$time_days <- sleep$time_hours/24
  
  # add id to sleep
  sleep$id <- nm
  
  # keep last
  activity <- activity[!duplicated(activity$startTimeInSeconds, fromLast = TRUE),]
  sleep <- sleep[!duplicated(sleep$startTimeInSeconds, fromLast = TRUE),]
  stress <- stress[!duplicated(stress$startTimeInSeconds, fromLast = TRUE),]
  
  # return the data
  list(activity=activity, sleep=sleep, stress=stress)
}

#' Function to get the earliest time from a set of json WRAP files
#' 
#' @param files a list of json file names. Use full path if files are not in
#' the working directory.
study_starts <- function(files)
{
  times <- list()
  for(i in seq_along(files))
  {
    dat <- fromJSON(file=files[i])
    times[[i]] <- dat[[1]]$data$startTimeInSeconds
  }
  
  list(times = unlist(times), study_start = min(unlist(times)))
}


#' Function to overlay phase information on data
#'
#' @param x a `sleep`, `activity`, or `stress` object produced by `read_wrap`
#' @param times a vector of UTC times that will serve as cutpoints an the `startTimeInSeconds`
#' variable. There should be one fewer `times` than the desired number of phase
#' to be produced.
#' @param labels a vector of the same length as `times` containing labels for the
#' phases. Optional. If `labels` is `NULL` the phase labels will be 1, 2, 3, etc..
add_phase <- function(x, times, labels=NULL)
{
  phase <- cut(x$startTimeInSeconds, breaks = c(-Inf, times, Inf))
  if(is.null(labels)) labels <- c(1:(length(times)+1))
  phase <- factor(phase, labels=labels)
  phase
}


#' Function to add daily data to WRAP data
#' 
#' @param wrap WRAP data produced by `read_wraps` (e.g., sleep, activity, or stress) 
#' @param daily A data frame in wide format with the first column containing the 
#' ID variable and the remaining columns with daily data with names in the format
#' `varName_#` where `#` are integers that overlap with the range of values in
#' the WRAP data variable `time_days`, eg. `range(sleep$time_days)`.
#' @param what Whether the `all` of the WRAP data appended with the daily data
#' should be returned, or only the `new` data should be returned (i.e., the 
#' daily data in long format with variable names required for merging with the
#' WRAP data). The function `append_dailies` loops through a list of daily
#' variables and uses the `new` option to reduce memory use.
append_daily <- function(wrap, daily, what="all")
{
  # check that wrap and daily share an id variable
  id <- names(daily)[1]
  check_id <- id %in% names(wrap)
  if(!check_id) stop("The `daily` ID variable `", id, 
                     "` is not in the WRAP data.\n")
  
  # check wrap day numbers
  wdays <- round(range(unique(wrap$time_days)))
  
  # get day numbers
  dnms <- names(daily)[2:ncol(daily)]
  dnm  <- strsplit(dnms, "_")
  dnum <- as.numeric(unlist(lapply(dnm, \(x) x[2])))
  dnm  <- unlist(lapply(dnm, \(x) x[1]))
  
  # check that the daily variable is not already a WRAP variable
  if(dnm[1] %in% names(wrap)) stop("The daily variable ", dnm[1], 
                                   " is already in the `wrap` data.\n",
                                   "You must first rename the daily variable.")
  
  # checks
  check_dnm <- all(dnm == dnm[1])
  if(!check_dnm) stop("Variable prefix in `daily` are not all the same.\n",
                      "They should all be\n\n", dnm[1], 
                      "\n\nbut the variable prefixes are\n\n",
                      paste(dnm, collapse="\n"))
  check_dnum <- any(dnum %in% wdays[1]:wdays[2])
  if(!check_dnum) stop("The numeric suffixes of the `daily` variable names\n",
                       "do not overlap with the `time_days` in the WRAP data.\n",
                       "The `daily` suffixes are:\n\n", 
                       paste(dnum, collapse="\n"),
                       "\n\nwhile the `time_days` in the WRAP data are:\n\n",
                       paste(wdays[1]:wdays[2], collapse="\n"))
  
  # restructure the daily data
  daily_long <- pivot_longer(daily, 2:ncol(daily))
  daily_long$name <- unlist(lapply(strsplit(daily_long$name, "_"), \(x) x[2]))
  names(daily_long)[which(names(daily_long)=="name")] <- "Day"
  names(daily_long)[which(names(daily_long)=="value")] <- dnm[1]
  
  # create an integer version of the time_days in wrap for merging
  wrap$Day <- round(wrap$time_days)
  
  # merge with wrap 
  dat <- merge(wrap, as.data.frame(daily_long), by=c(id, "Day"), all=TRUE)
  
  # if data are character, convert "" to NA
  dat[[dnm[1]]][dat[[dnm[1]]]==""] <- NA
  
  # return the data
  if(what=="all") return(dat)
  if(what=="new") return(dat[,c(id, "Day", "time_days", dnm[1])])
}

#' Function to add daily data to WRAP data for a list of daily variables
#' 
#' @param data A data frame in wide format containing multiple daily variables.
#' This data must include an ID variable named `id` with integer values starting
#' at 1 to match the default ID variable naming used by WRAP.
#' @param wrap WRAP data produced by `read_wraps` (e.g., sleep, activity, or stress) 
#' @param dailies A list of daily variable names 
#' `e.g., list(c(x_1, x_2, ..., x_t), c(y_1, y_2, ..., y_2))` where names are in 
#' the format `varName_#`, where `#` are integers that overlap with the range of 
#' values in the WRAP data variable `time_days`, e.g. `range(sleep$time_days)`.
append_dailies <- function(data, wrap, dailies)
{
  out <- list()
  for(i in seq_along(dailies))
  {
    vars <- dailies[[i]]
    missing_vars <- vars[!vars %in% names(data)]
    if(length(missing_vars) > 0)
    {
      warning("Missing variables:\n\n", 
              paste(missing_vars, collapse="\n"))
    }
    vars <- vars[vars %in% names(data)]
    daily <- data[,c("id", vars)]
    out[[i]] <- append_daily(wrap, daily, "new")
  }
  out <- Reduce(\(x,y) merge(x, y), out)
  merge(wrap, out)
}


#' Function to add time invariant and demographic data to WRAP data
#' 
#' @param wrap WRAP data produced by `read_wraps` (e.g., sleep, activity, or stress) 
#' @param survey A data frame in wide format with the first column containing the 
#' ID variable and the remaining columns with survey data with names in the format
#' `varName_#` where `#` are integers that overlap with the range of values in
#' the WRAP data variable `time_days`, eg. `range(sleep$time_days)`.
#' @param what Whether the `all` of the WRAP data appended with the survey data
#' should be returned, or only the `new` data should be returned (i.e., the 
#' survey data in long format with variable names required for merging with the
#' WRAP data). The function `append_dailies` loops through a list of survey
#' variables and uses the `new` option to reduce memory use.
append_survey <- function(vars, wrap, survey)
{
  # check that wrap and survey share an id variable
  id <- names(survey)[1]
  check_id <- id %in% names(wrap)
  if(!check_id) stop("The `survey` ID variable `", id, 
                     "` is not in the WRAP data.\n")
  
  check_vars <- vars %in% names(survey)
  if(!all(check_vars)) stop("The following variables in `vars` are not in `survey`:\n\n",
                       paste(vars[! vars %in% names(survey)], collapse="\n"))
  
  merge(wrap, survey[,c(id, vars)])
}

