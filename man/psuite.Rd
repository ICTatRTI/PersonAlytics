% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psuite.R
\name{psuite}
\alias{psuite}
\title{psuite - apply p.adjust and add output and graphics}
\usage{
psuite(
  DVout,
  ids,
  output,
  method = "BY",
  nbest = NULL,
  alpha = 0.05,
  rawdata = NULL
)
}
\arguments{
\item{DVout}{A data frame created by \code{\link{PersonAlytic}}.}

\item{ids}{The variable in \code{DVout} across which to adjust.}

\item{output}{A character string for labeling the output}

\item{method}{One of \code{\link{p.adjust.methods}}.}

\item{nbest}{The \code{nbest} largest parameters for the target predictor with
\code{p<alpha} to be print to output.}

\item{alpha}{The type I error rate.}

\item{rawdata}{The raw data for plotting (not implemented).}
}
\description{
Apply p.adjust and add output and graphics
}
\author{
Stephen Tueller \email{stueller@rti.org}
}
