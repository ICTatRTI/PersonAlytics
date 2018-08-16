#' PersonAlyticsPro: A package for Ideographic Clinical Trial analyses with high
#' through put.
#'
#' The PersonAlyticsPro package provides the simplified user interface
#' for implementing linear mixed effects models for ideagraphic clinical trial (ICT) data using
#' \code{\link{gamlss}} or \code{\link{lme}} using maximum likelihood (ML) for model comparisons and
#' restricted maximum likelihood (REML) for model results. PersonAlyticsPro automates
#' model comparisons for determining autocorrelation structure (for all patients or for
#' each patient) and for determining the the relationship between time and the outcome,
#' e.g., \code{time^2} vs. \code{time^3} (for all patients or for each patient). Models
#' will be repeated for a list of outcomes, a list of target covariates, or all pairs
#' of both a list of outcomes and a list of target covariates. Options for Type I error
#' correction and false discovery rate adjustments are implemented post-implementation
#' across all analyses and also within individuals if individual level analyses are requested.
#'
#' @details
#' \code{PersonAlyticsPro} builds on \code{\link{PersonAlyticsLite}}
#'
#' The basic model for ICTs is \eqn{dv=time+phase+phase*time} with random intercepts
#' and random slopes for time. The phase variable is optional. Additional covariates
#' can be added.
#'
#' This is done via the \code{\link{PalyticHTP}} object, which
#' offers additional options for advanced users including generalized linear mixed effects
#' models (see \code{\link{gamlss.family}}), user specific correlation structures, and
#' user specified random effects structures.
#'
#'
#' @section PersonAlytics functions
#'
#' @section the PalyticHTP class
#' See \code{\link{PersonAlyticLite::PalyticHTP}}
#'
#' @docType package
#' @name PersonAlyticsPro
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#'
NULL