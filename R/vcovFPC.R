#' @name vcovFPC
#' @aliases vcovFPC.meMod
#' @aliases vcovFPC.lme
#'
#' @title \code{vcovFPC} Obtain finite-population-adjusted standard errors for
#' fixed effects estimates for a fitted multilevel model
#'
#' @import Matrix
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
#' library(gamlss)
#' library(lme4)
#'
#' # nlme Ovary data
#'
#' mod.lme <- lme(follicles ~ Time, data = Ovary, random = ~ Time | Mare,
#'              method = 'ML')
#'
#' mod.merMod <- lmer(follicles ~ Time + (Time | Mare), data = Ovary,
#'                 REML = FALSE)
#'
#' popsize2 <- seq(100,1000,by=100)
#'
#' fpc.se.merMod <- fpc.se.lme <- list()
#' for(i in seq_along(popsize2))
#' {
#'   fpc.se.merMod[[i]] <- sqrt(diag(vcovFPC(mod.merMod, popsize2 = popsize2[i])))
#'   fpc.se.lme[[i]] <- sqrt(diag(vcovFPC(mod.lme, popsize2 = popsize2[[i]])))
#' }
#'
#' all.equal(
#' do.call(rbind, fpc.se.merMod),
#' do.call(rbind, fpc.se.lme),
#' tolerance = 9e-06
#' )
#'
#' \dontrun{
#'
#' # use getLambdat to create Lambdat and show equivalence to lme4 (within rounding)
#' Lambdat <- getLambdat(mod.lme)
#' all.equal(Lambdat, getME(mod.merMod, "Lambdat"), tolerance = 9e-06)
#'
#' # use getZt to create Zt and show equivalence to lme4
#' Zt <- getZt(mod.lme)
#' all.equal(unname(Zt), unname(getME(mod.merMod, "Zt")))
#'
#' }
#'
#' # nlme rat BodyWeight example
#'
#' mod.lme <- lme(weight ~ Time, data = BodyWeight, random = ~ Time | Rat,
#'              method = 'ML')
#'
#' mod.merMod <- lmer(weight ~ Time + (Time | Rat), data = BodyWeight,
#'                 REML = FALSE)
#'
#' fpc.se.merMod <- fpc.se.lme <- list()
#' for(i in seq_along(popsize2))
#' {
#'   fpc.se.merMod[[i]] <- sqrt(diag(vcovFPC(mod.merMod, popsize2 = popsize2[i])))
#'   fpc.se.lme[[i]] <- sqrt(diag(vcovFPC(mod.lme, popsize2 = popsize2[[i]])))
#' }
#'
#' all.equal(
#' do.call(rbind, fpc.se.merMod),
#' do.call(rbind, fpc.se.lme),
#' tolerance = 9e-06
#' )
#'
#' \dontrun{
#'
#' # use getLambdat to create Lambdat and show equivalence to lme4 (within rounding)
#' Lambdat <- getLambdat(mod.lme)
#' all.equal(Lambdat, getME(mod.merMod, "Lambdat"), tolerance = 9e-06)
#'
#' # use getZt to create Zt and show equivalence to lme4
#' Zt <- getZt(mod.lme)
#' all.equal(unname(Zt), unname(getME(mod.merMod, "Zt")))
#'
#' }
#'
#' \dontrun{
#'
#' # some background tests
#'
#' # random effects covariance matrices
#' getVarCov(mod.lme) # this is the same as VarCorr
#' VarCorr(mod.lme)
#'
#' vcov(mod.merMod) # this is fixed effects, not the same matrix as VarCorr
#' VarCorr(mod.merMod)
#'
#' # not implemented for gamlss
#'
#' mod.gamlss <- gamlss(follicles ~ Time + re(random = ~ Time | Mare,
#'                 method = 'ML'), sigma.fo = ~ 0,
#'                 data = Ovary, )
#' mod.gamlss <- gamlss(weight ~ Time + re(random = ~ Time | Rat,
#'                 method = 'ML'), sigma.fo = ~ 0,
#'                 data = BodyWeight, )
#' }

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
  Astar <- A * sqrt(fpc2)
  X <- PR$X
  Astar_X <- Astar %*% X
  D <- Matrix::Diagonal(nrow(Astar), fpc1) + tcrossprod(Astar)
  Fisher_I <- (crossprod(X) - crossprod(solve(t(chol(D)), Astar_X))) / fpc1
  Phi <- solve(Fisher_I) * sigma(object)^2 # TODO: this used to be 2^, is my correction correct?
    Phi <- as(Phi, "dpoMatrix")
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
  N <- unname(object$dims["N"])
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
  Astar <- A * sqrt(fpc2)
  X <- model.matrix(object, data = object$data)
  Astar_X <- Astar %*% X
  D <- Matrix::Diagonal(nrow(Astar), fpc1) + tcrossprod(Astar)
  Fisher_I <- (crossprod(X) - crossprod(solve(t(chol(D)), Astar_X))) / fpc1
  Phi <- solve(Fisher_I) * sigma(object)^2
    Phi <- as(Phi, "dpoMatrix")
  nmsX <- colnames(X)
  dimnames(Phi) <- list(nmsX, nmsX)
  if (!KR) {
    return(Phi)
  } else {
    stop('kr not implemented for lme objects.')
  }
}

#' getLambdat - get Lambdat from an lme object analogous to
#' getME(object.merMod, "Lambdat"). See \code{\link{getMe}}.
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
#' getME(object.merMod, "Zt"). See \code{\link{getMe}}.
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

# requires equal # of time points per person
if(1==2)
{
  getZt.old <- function(object)
{
  id  <- names(object$groups)[1]
  ids <- object$data[[id]]
  ido <- as.numeric(levels(ids))[ids]
  o   <- order(as.numeric(levels(ids))) # this won't generalize to character id
  X   <- as.data.frame( model.matrix(object, data = object$data) )
  XL  <- lapply(split(X, ids), as.matrix)
  # we break up and reassemble Z b/c the order must match that of the data,
  # but the do.call(cibind, Zt) only works if the nObs per case are equal
  Z   <- Matrix::bdiag(XL)
  Z   <- split(as.data.frame(as.matrix(Z)), ido)
  Z   <- lapply(Z[o], as.matrix)
  Zt  <- lapply(Z, t)
  Zt  <- do.call(cbind, Zt)
  Zt  <- Matrix::bdiag(Zt)
  return(Zt)
}
}



