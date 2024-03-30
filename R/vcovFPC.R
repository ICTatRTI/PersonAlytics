#' FPC function to get the finite population corrected standard errors and
#' update the tTable of an lme object
#'
#' @export
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param object An \code{\link{lme}} object from \code{\link{nlme}}, or a merMod
#' object from \code{\link{lmer}}.
#'
#' @param popsize2 The finite population size at level two, i.e., persons when
#' there are repeated measures at level 1.
#'
FPC <- function(object, popsize2)
{
  if(! any(class(object) %in% c('lme', 'merMod')))
  {
    stop('`object` must be an lme or merMod object.')
  }

  if(inherits(object, "lme"))
  {
    tTable <- data.frame( summary(object)$tTable )
  }

  if(inherits(object, "merMod"))
  {
    tTable <- data.frame( coef(summary(object)) )
    names(tTable) <- c('Value', 'Std.Error', 't-value')
    # the table doesn't have DF, wee need those before we can get tests
    stop('FPC is not fully implemented for lem4/merMod objects.')
  }

  # use the approach of `summary.lme` for t- and p-values
  tTable$Std.Error <- sqrt(diag(vcovFPC(object, popsize2 = popsize2)))
  tTable$t.value   <- tTable$Value/tTable$Std.Error
  tTable$p.value   <- 2 * pt(-abs(tTable$t.value), tTable$DF)

  # return the FPC tTable
  return(tTable)
}

#' fpcCheck
#'
#' @keywords internal
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
fpcCheck <- function(popsize2, n)
{
  if( !exists("popsize2") )
  {
    stop('\nA finite population was requested but no level-2',
         '\nfinite population size `popsize2` was provided.')
  }
  if( ! popsize2 > n )
  {

    stop('\nA finite population was requested but the finite',
         '\npopulation size `popsize2`=', popsize2, ' , which is',
         '\nsmaller than the total sample size n=', n)
  }
}

#' @name vcovFPC
#' @aliases vcovFPC.meMod
#' @aliases vcovFPC.lme
#'
#' @import pbkrtest
#'
#' @title \code{vcovFPC} Obtain finite-population-adjusted standard errors for
#' fixed effects estimates for a fitted multilevel model
#'
#' @export
#'
#' @author Mark H. C. Lai \email{mark.lai@@uc.edu}, Oi-man Kwok, Yu-Yu Hsiao,
#' and Quina Cao. Updated by Stephen Tueller \email{Stueller@@rti.org} for
#' compatibility with \code{nlme}
#' for inclusion in \code{\link{PersonAlytics}}.
#'
#' @references https://psycnet.apa.org/doiLanding?doi=10.1037%2Fmet0000137
#'
#' @param object an R object of class merMod as resulting from lmer() or
#' an object of class lme resulting from lme().
#'
#' @param popsize2 population size at Level-2; if NULL, an infinite Level-2
#' population is assumed.
#'
#' @param popsize1 population size at Level-1; if NULL, an infinite Level-1
#' population is assumed.
#'
#' @param KR Whether Kenward-Roger approximation of standard errors should be used,
#' which is recommended for small number of clusters and average cluster size.
#' Default to FALSE. Not available for lme objects.
#'
#' @return
#' The variance-covariance matrix of the fixed effect estimates, as
#' returned by vcov()
#'
#' @examples
#'
#' library(nlme)
#' library(lme4)
#'
#' # illustrate equivalence in lme and merMod results with a simple ICT model
#' mod.lme <- lme(follicles ~ Time*Phase, data = OvaryICT, random = ~ Time | Mare,
#'                method = 'ML')
#' mod.merMod <- lmer(follicles ~ Time*Phase + (Time | Mare), data = OvaryICT,
#'                    REML = FALSE)
#'
#' tTable.lme <- data.frame( summary(mod.lme)$tTable[,1:2] )
#' tTable.merMod <- data.frame( coef(summary(mod.merMod))[,1:2] )
#'
#' fpc.lme <- vcovFPC(mod.lme, popsize2 = 100)
#' fpc.merMod <- vcovFPC(mod.merMod, popsize2 = 100)
#'
#' tTable.lme$FPC.Std.Error <- sqrt(diag(fpc.lme))
#' tTable.merMod$FPC.Std.Error <- sqrt(diag(fpc.merMod))
#'
#' tTable.lme
#' tTable.merMod
#'
#' # now show how the standard errors change whith an autocorrelation structure
#' mod.lmema2 <- lme(follicles ~ Time*Phase, data = OvaryICT, random = ~ Time | Mare,
#'                method = 'ML', correlation = corARMA(q=2))
#' tTable.lmema2 <- data.frame( summary(mod.lmema2)$tTable[,1:2])
#' fpc.lmema2 <- vcovFPC(mod.lmema2, popsize2 = 100)
#' tTable.lmema2$FPC.Std.Error <- sqrt(diag(fpc.lmema2))
#'
#' # including a good fitting autocorrelation structure reduces the standard
#' # errors relative to an unstructured correlation, and the FPC further
#' # reduces the standard errors (thought the correlation structure has a
#' # much bigger impact than the FPC)
#' tTable.lmema2
#'

vcovFPC <- function(object, ...)
{
  UseMethod("vcovFPC")
}

#' @rdname vcovFPC
#'
#' @export
vcovFPC.merMod <- function(object, popsize2 = NULL,
                    popsize1 = NULL, KR = FALSE)
{

  if (!inherits(object, "merMod")) {
    stop("Wrong input: Not a fitted model from lmer() with class merMod")
  }
  if (length(object@flist) != 1) {
    stop("Wrong input: Only models with two levels are supported")

  }
  if (is.null(popsize1) & is.null(popsize2)) {
    message("No FPC specified; return results from lme4::vcov.merMod()")
    return(vcov(object))
  }
  PR <- object@pp
  N <- unname(object@devcomp$dims["n"])
  nclus <- unname(lme4::ngrps(object))
  if (isTRUE(popsize2 > nclus)) fpc2 <- 1 - nclus / popsize2
  else {
    fpc2 <- 1
    message("No FPC needed at Level-2")
  }
  if (isTRUE(popsize1 > N)) fpc1 <- 1 - N / popsize1
  else {
    fpc1 <- 1
    #message("No FPC needed at Level-1") # never FPC@level1 which is obs w/in person
  }
  if (fpc1 == 1 & fpc2 ==1) {
    message("Return results from lme4::vcov.merMod()")
    return(vcov(object))
  }
  A <- PR$Lambdat %*% PR$Zt # equivalent to: A <- getME(object, "A")
  Astar <- as.matrix( A * sqrt(fpc2) )
  X <- PR$X
  Astar_X <- Astar %*% X
  D <- Matrix::Diagonal(nrow(Astar), fpc1) + tcrossprod(Astar)
  Fisher_I <- (crossprod(X) - crossprod(solve(t(chol(D)), Astar_X))) / fpc1
  Phi <- solve(Fisher_I) * sigma(object)^2
  #Phi <- as(Phi, "dpoMatrix")
  nmsX <- colnames(X)
  dimnames(Phi) <- list(nmsX, nmsX)
  if (!KR) {
    return(Phi)
  } else {
    if (!require("pbkrtest")) {
      stop("Please install the 'pbkrtest' package for the use of Kenward-Roger correction!")
    } else {
      SigmaG <- pbkrtest::get_SigmaG(object, details = 0)
      vcov_kr <- pbkrtest:::vcovAdj16_internal(Phi, SigmaG, X, details = 0)
      vcov_kr <- as(Phi, "dpoMatrix")
      return(vcov_kr)
    }
  }
}

#' @rdname vcovFPC
#'
#' @export
vcovFPC.lme <- function(object, popsize2 = NULL,
                        popsize1 = NULL, KR = FALSE)
{

  if (!inherits(object, "lme")) {
    stop("Wrong input: Not a fitted model from lme() with class lme")
  }
  if (ncol(object$groups) != 1) {
    stop("Wrong input: Only models with two levels are supported")

  }
  if (is.null(popsize1) & is.null(popsize2)) {
    message("No FPC specified; return results from lme4::vcov.merMod()")
    return(vcov(object))
  }
  N <- unname(unlist(object$dims["N"]))
  nclus <- unname(unlist(object$dims["ngrps"])[1])
  if (isTRUE(popsize2 > nclus)) fpc2 <- 1 - nclus / popsize2
  else {
    fpc2 <- 1
    message("No FPC needed at Level-2")
  }
  if (isTRUE(popsize1 > N)) fpc1 <- 1 - N / popsize1
  else {
    fpc1 <- 1
    #message("No FPC needed at Level-1")  # never FPC@level1 which is obs w/in person
  }
  if (fpc1 == 1 & fpc2 ==1) {
    message("Return results from nlme::vcov.lme()")
    return(vcov(object))
  }
  A <- getLambdat(object) %*% getZt(object)
  Astar <- as.matrix( A * sqrt(fpc2) )
  X <- model.matrix(object, data = object$data)
  Astar_X <- Astar %*% X
  D <- Matrix::Diagonal(nrow(Astar), fpc1) + tcrossprod( Astar )
  Fisher_I <- (crossprod(X) - crossprod(solve(t(chol(D)), Astar_X))) / fpc1
  Phi <- solve(Fisher_I) * sigma(object)^2
  #Phi <- as(Phi, "dpoMatrix")
  nmsX <- colnames(X)
  dimnames(Phi) <- list(nmsX, nmsX)
  if (!KR) {
    return(Phi)
  } else {
    stop('kr not implemented for lme objects.')
  }
}

#' getLambdat - get Lambdat from an lme object analogous to
#' getME(object.merMod, "Lambdat"). See \code{\link{getME}}.
#'
#' @export
#'
#' @param object An \code{lme} model
#'
#' @examples
#'
#' # see examples in \code{vcovFPC}
getLambdat <- function(object)
{
  if (!inherits(object, "lme")) {
    stop("Wrong input: Not a fitted model from lme() with class lme")
  }
  lambda  <- chol(getVarCov(object))/object$sigma
  class(lambda) <- "matrix"
  nclus   <- unname(unlist(object$dims["ngrps"])[1])
  Lambdat <- Matrix::bdiag( rep(list(lambda), nclus) )
  return(Lambdat)
}

#' getZt - get Zt from an lme object analogous to
#' getME(object.merMod, "Zt"). See \code{\link{getME}}.
#'
#' @export
#'
#' @param object An \code{lme} model
#'
#' @examples
#'
#' # see examples in \code{vcovFPC}
getZt <- function(object)
{
  id  <- names(object$groups)[1]
  ids <- object$data[[id]]
  # to delete
  #o   <- as.numeric(levels(ids))
  #ido <- as.numeric(levels(ids))[ids]

  # X is just the fixed effects design matrix
  X   <- as.data.frame( model.matrix(object, data = object$data) )
  # discard fixed effects only columns
  X   <- X[,names(ranef(object))]
  # break X up by id
  XL  <- lapply(split(X, ids), as.matrix)
  # transpose the blocks
  XL  <- lapply(XL, t)
  # get the block order from XL
  XLnms <- as.list( names(XL) )
  XLlen <- lapply(XL, ncol)
  Zo <- unlist(
    mapply(function(nms, len){
      rep(nms, len)
    }, nms=XLnms, len=XLlen)
  )
  Zo <- as.numeric(Zo)

  # create Z, which will have the right vertical positions but wrong column block order
  Z   <- Matrix::bdiag(XL)

  # split Z, transposing first b/c split works on rows
  Z   <- split(as.data.frame(t(as.matrix(Z))), Zo) #ids

  Zt  <- lapply(Z, t)
  Zt  <- do.call(cbind, Zt)
  Zt  <- Matrix::bdiag(Zt)
  return(Zt)
}
