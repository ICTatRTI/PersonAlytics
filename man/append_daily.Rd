% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_WRAP.R
\name{append_daily}
\alias{append_daily}
\title{Function to add daily data to WRAP data}
\usage{
append_daily(wrap, daily, what = "all")
}
\arguments{
\item{wrap}{WRAP data produced by `read_wraps` (e.g., sleep, activity, or stress)}

\item{daily}{A data frame in wide format with the first column containing the 
ID variable and the remaining columns with daily data with names in the format
`varName_#` where `#` are integers that overlap with the range of values in
the WRAP data variable `time_days`, eg. `range(sleep$time_days)`.}

\item{what}{Whether the `all` of the WRAP data appended with the daily data
should be returned, or only the `new` data should be returned (i.e., the 
daily data in long format with variable names required for merging with the
WRAP data). The function `append_dailies` loops through a list of daily
variables and uses the `new` option to reduce memory use.}
}
\description{
Function to add daily data to WRAP data
}
