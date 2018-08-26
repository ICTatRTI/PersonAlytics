#' ICC function from \href{http://davidakenny.net/papers/k&h/R_output.txt}{David Kenny}
#' @author \href{http://davidakenny.net}{David Kenny}
#' @param out An \code{lme} object
#'
#' @export
#' @import multilevel
#' @import gamlss
#' @importFrom gamlss re
#'
#' @examples
#' egaov <- aov(follicles ~ factor(Mare), data = OvaryICT)
#' multilevel::ICC1(egaov)
#' multilevel::ICC2(egaov)
#' eglme <- nlme::lme(follicles ~ 1, data = OvaryICT, random = ~ 1 | Mare,
#'                    method = 'ML')
#' # VarCorr(eglme)
#' ICC(eglme)
#' library(gamlss)
#' eggamlss1 <- gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
#'                            data = OvaryICT,
#'                            sigma.formula = ~ 1)
#' # VarCorr(getSmo(eggamlss1))
#' ICC(eggamlss1)
#' eggamlss2 <- gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
#'                            data = OvaryICT,
#'                            sigma.formula = ~ 0)
#' # VarCorr(getSmo(eggamlss2))
#' ICC(eggamlss2)
ICC <- function(out) # make this a Palytic method
{
  if('lme' %in% class(out)) varests <- as.numeric(nlme::VarCorr(out)[1:2])
  if('gamlss' %in% class(out)) varests <- as.numeric(nlme::VarCorr(gamlss::getSmo(out))[1:2])
  return( varests[1]/sum(varests) )
}


#' eds - testing that an object exists
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @param x Any object
#'
#' @keywords internal
eds <- function(x)
{
  exists( deparse( substitute(x) ) )
}

#' monotone
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @param x Any object
#'
#' @keywords internal
monotone <- function(ids, time, data)
{
  .m <- function(x) all( diff(x) >= 0 )
  monotonic <- by(data[[time[[1]]]], INDICES = data[[ids]], FUN = .m)
  monotonic <- unlist( as.list(monotonic) )
  data.frame(ids=ids, monotonic=monotonic)
}

#' isCorStruct - function to test whether a string resolves into a valid
#' correlation structure
#'
#' @description Note that \code{NULL} is a valiad \code{corStruct}
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @param x Any character string
#'
#' @keywords internal
#'
#' @examples
#' iscorStruct(NULL)
#' iscorStruct("NULL")
#' iscorStruct("corARMA(3,3)")
#' iscorStruct(c(2,2))
iscorStruct <- function(x)
{
  if( ! class(x) %in% c('NULL', 'character') &
      ( length(x)!=2 & !is.numeric(x) )
    )
  {
    stop('`correlation` must be one of\n\n',
         '- A quoted character string of an `nlme` `corStruct`, see `?corStruct`\n',
         '- `"NULL"` or `NULL` to get the default for `lme`\n',
         '- A numeric vector of length 2, `correlation` and `detectAR` in `?PersonAlytic`')
  }
  if( is.numeric(x) )
  {
    if( sum(x) == 0 )
    {
      stop('At least one of `P` or `Q` in `correlation=c(P,Q)` must be > 0. ',
           '\n\nIf you wish to specify a specific `ARMA(p,q)` model, use ',
           '`correlation="corARMA(p,q)" \ninstead of `correlation=c(P,Q)` which ',
           'initializes an automatice search among\n`p=1,...,P` and ',
           '`q=1,...,Q`.')
    }
  }
  if( ! is.null(x) & is.character(x) )
  {
    if(x!="NULL")
    {
      if(! 'corStruct' %in% class( eval( parse( text = x ) ) ) )
      {
        stop('`correlation` is not a valid `corStruct`')
      }
    }
  }
  return(TRUE)
}


#' forms - function to construct and return formulae (fixed, random, formula)
#' whenever a change to any variable that can go into a formula is made
#'
#' Note that only the following combinations work, but these are enforced in
#' \code{\link{Palytitic}} rather than enforced internally
#' 1. Objects except for PalyticObj, fixed, random, and formula (e.g., Palytic$new)
#' 2. PalyticObj and exactly 1 other parameter
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param See \code{\link{Palytic}}
#'
#' @keywords internal
#'
forms <- function(data                  ,
                  PalyticObj   = NULL   ,
                  ids          = NULL   ,
                  dv           = NULL   ,
                  time         = NULL   ,
                  phase        = NULL   ,
                  ivs          = NULL   ,
                  interactions = NULL   ,
                  time_power   = NULL   ,
                  correlation  = NULL   ,
                  family       = NULL   ,
                  fixed        = NULL   ,
                  random       = NULL   ,
                  formula      = NULL   ,
                  method       = "REML" ,
                  dropTime     = FALSE  ,
                  corFromPalyticObj = TRUE)
{
  # since NULL is a valid option for correlation, we must override it using
  #corFromPalyticObj

  if(!is.null(PalyticObj))
  {
    # unpack PalyticObj
    if(is.null(ids         )) ids          <- PalyticObj$ids
    if(is.null(dv          )) dv           <- PalyticObj$dv
    if(is.null(time        )) time         <- PalyticObj$time[[1]]
    if(is.null(phase       )) phase        <- PalyticObj$phase
    if(is.null(ivs         )) ivs          <- PalyticObj$ivs
    if(is.null(interactions)) interactions <- PalyticObj$interactions
    if(is.null(time_power  )) time_power   <- PalyticObj$time_power
    if( corFromPalyticObj)    correlation  <- PalyticObj$correlation
    if(!corFromPalyticObj)    correlation  <- correlation
    if(is.null(family      )) family       <- PalyticObj$family
    if(is.null(method      )) method       <- PalyticObj$method

    # overwrite objects if fixed or random are provided
    if(!is.null(fixed))
    {
      dv  <- all.vars(fixed)[1]
      ivs <- strsplit(as.character(fixed), "~")[[3]]
    }
    if(!is.null(random))
    {
      random.t <- gsub(" ", "", unlist( strsplit(as.character(random), "\\|") ) )
      time     <- random.t[2]
      ids      <- random.t[3]
    }
  }

  if(is.null(PalyticObj))
  {
    if(!is.null(formula))
    {
      theInputs   <- decompFormula(formula)
      ids         <- theInputs$ids
      dv          <- theInputs$dv
      time        <- theInputs$time
      ivs         <- theInputs$ivs
      correlation <- theInputs$correlation
      method      <- theInputs$method
    }
  }

  # update time using time_power
  if(!is.null(time_power))
  {
    if( is.null(time_power) ) time_power <- 1
    if(time_power > 1)
    {
      time <- c(time, paste("I(", time, "^", 2:time_power, ")", sep=''))
    }
    if(time_power == 1) time <- time[[1]]
  }

  theForms <- makeForms(ids          = ids         ,
                        dv           = dv          ,
                        time         = time        ,
                        phase        = phase       ,
                        ivs          = ivs         ,
                        interactions = interactions,
                        correlation  = correlation ,
                        family       = family      ,
                        dropTime     = dropTime    ,
                        method       = method      )
  fixed   <- theForms$fixed
  random  <- theForms$random
  formula <- theForms$formula

  # check that the variables are in the data
  ivs.test <- unlist(lapply(ivs, strsplit, ":"))
  vars <- unique( c(ids, dv, time, phase, ivs.test,
                    all.vars(fixed), all.vars(random), all.vars(formula)) )
  vars <- unique( gsub(" ", "", unlist(lapply(strsplit(vars, '\\*|\\+'), unlist))) )
  wi <- which( substr(vars, 1, 2) == "I(" )
  vars[wi] <- unlist(strsplit(gsub("\\I|\\(|\\)", "", vars[wi]), "\\^"))[1]
  vars <- vars[which(vars!="1")]
  vars <- vars[!is.na(vars)]

  wvars <- which( ! vars %in% names(data) )

  if( length(wvars) > 0 )
  {
    stop( paste('\n`', vars[wvars], '` is not in the data\n'))
  }

  # return
  return( list(  ids          = ids         ,
                 dv           = dv          ,
                 time         = time        ,
                 phase        = phase       ,
                 ivs          = ivs         ,
                 interactions = interactions,
                 time_power   = time_power  ,
                 correlation  = correlation ,
                 family       = family      ,
                 fixed        = fixed       ,
                 random       = random      ,
                 formula      = formula     ,
                 method       = method      ) )

}


#' makeForms - function to construct and return formulae (fixed, random, formula)
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param See \code{\link{Palytic}}
#'
#' @keywords internal
#'
makeForms <- function(ids          = "Mare"                   ,
                      dv           = "follicles"              ,
                      time         = "Time"                   ,
                      phase        = "Phase"                  ,
                      ivs          = "iv1"                    ,
                      interactions = list(c("iv1", "Phase"))  ,
                      correlation  = "NULL"                   ,
                      family       = NO()                     ,
                      dropTime     = FALSE                    ,
                      method       = "REML"                   )
{
  # construct or update the formula objects

  # fixed
  if(  is.null(phase) )
  {
    rhs <- paste(time, collapse = '+')
  }
  if( !is.null(phase) )
  {
    rhs <- paste( paste(time, phase, sep='*'), collapse = '+')
  }
  fixed  <- formula( paste(dv, "~", rhs ) )


  # random
  if(!dropTime)
  {
    random <- formula( paste("~", paste(time, collapse = '+'), "|", ids))
  }
  if( dropTime)
  {
    random <- formula( paste("~", paste(1, collapse = '+'), "|", ids))
  }

  # formula
  formula <- formula( paste(paste(deparse(fixed), collapse=''),
                            "+ re(random = ",
                            #"+ gamlss::re(random = ", # error, use `@importFrom gamlss re`
                            paste(deparse(random), collapse=''),
                            ", method=",
                            deparse(method),
                            ", correlation =",
                            ifelse(is.null(correlation), "NULL", correlation),
                            ")")
                      )

  if(!is.null(ivs) & length(ivs) > 0)
  {
    formula <- update(formula, paste("~ . +", paste(unlist(ivs), collapse='+')))
    fixed   <- update(fixed,   paste("~ . +", paste(unlist(ivs), collapse='+')))
  }

  # interactions
  if(!is.null(interactions))
  {
    ### NEEDS DEBUGGING
    # if phase is present, check whether the phase x time interaction is
    # requested and remove the request, phase x time is currently enforced
    #if(!is.null(phase))
    #{
    #  wpt <- lapply(interactions, function(x) all( c(phase, time) %in% x ))
    #  wpt <- which(unlist(wpt))
    #  interactions[[wpt]] <- NULL
    #}
    for(i in 1:length(interactions))
    {
      formula <- update(formula, paste("~ . +", paste(interactions[[i]], collapse="*")))
      fixed   <- update(fixed,   paste("~ . +", paste(interactions[[i]], collapse="*")))
    }
  }

  return(list(fixed=fixed, random=random, formula=formula))
}

#' decompFormula - function to deconstruct a return formula
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param See \code{\link{Palytic}}
#'
#' @keywords internal
#'
decompFormula <- function(formula=NULL)
{
  if(is.null(formula))
  {
    formula <- formula(y~x+z+a*b+re(random=time + I(time^2)|id,
                       correlation=corARMA(p=1,q=1), method="ML") +
                         k + l)
    if(1==2){
    formula <- follicles ~ Time * Phase + re(random = ~1 | Mare, method = "REML",
                                                correlation = corARMA(p = 1, q = 0))

    formula <- PASAT ~ (re(random = ~1 | ID, method = "REML",
                           correlation = nlme::corARMA(p = 2, q = 2))) +
                          Time2 + Tx + Time2:Tx

    formula <- PASAT ~ (re(random = ~1 | ID, method = "REML", correlation = NULL)) +
      Time2 + Tx + X19.26_888.2232n + Time2:Tx

    formula <- PASAT ~ Time2 + Tx + (re(random = ~1 | ID, method = "REML",
                                        correlation = nlme::corARMA(p = 1, q = 1))) +
      X19.26_888.2232n + Time2:Tx

    formula <- DSST ~ Time2 + Tx + I(Time2^2) + I(Time2^3) +
      (re(random = ~Time2 + I(Time2^2) + I(Time2^3) | ID,
          method = "REML", correlation = nlme::corARMA(p = 3, q = 3))) +
      X3.81_442.2293m.z + Batch + Session2 + Time2:Tx + Tx:I(Time2^2) + Tx:I(Time2^3)

    formula <- 1 ~ TimeSin * Phase + re(random = ~1 | Mare, method = "REML",
                                        correlation = NULL, family = )
    }
  }

  # convert to character
  f.char <- as.character( formula )

  # LHS - get the dependent variable
  dv <- f.char[2]

  # deconstruct the RHS into fixed and re+
  # (where the `+` in re+ may include more fixed effects)
  f.char3 <- paste(unlist(strsplit(f.char[3], "\\(re")), collapse="re")
  f.char3 <- gsub("))", ")", f.char3)
  f.char3 <- unlist( strsplit(f.char3, ' \\+ re') )

  if(length(f.char3)<2)
  {
    f.temp <- unlist( strsplit(f.char3, "\\+") )
    if(! grepl('re', f.temp[1]) )
    {
      stop('Problem parsing random effects in `decompFormula()`')
    }
    f.char3 <- c(paste(f.temp[2:length(f.temp)], collapse='+'),
                f.temp[1])
  }

  # get the ivs
  ivs <- gsub(" ", "", as.list( unlist( strsplit(f.char3[1], '\\+') ) ) )

  # deconstruct the random effects
  re <- unlist( strsplit(f.char3[2], ', ') )

  # find random, correlation, and method, should they exist
  # these will fail if someone passes variable in with partial name
  random.w <- which( grepl('random', re) )
  correl.w <- which( grepl('correlation', re) )
  method.w <- which( grepl('method', re) )

  # random must be present (technically not true, gamlss will run ok w/o)
  if(length(random.w)==0)
  {
    stop('The required `random` component in `formula` is missing.')
  }

  # get LHS and RHS of `|` in random
  random.t <- unlist( strsplit(re[random.w], '=') )[2]
  random.t <- gsub(" ", "", unlist( strsplit(random.t, '\\|') ))

  # get ids
  ids <- random.t[2]

  # get time, should it exist, may just be `1`
  time <- gsub("~", "", random.t[1])

  # get method, set to REML if it doesn't exist
  if(length(method.w)>0)
  {
    method <- ifelse(grepl("REML", re[method.w]), 'REML', 'ML')
  }
  if(length(method.w)==0)
  {
    method <- "REML"
  }

  # Set moreiv to null
  moreiv <- NULL

  # find a full correlation statement
  correl.t <- "NULL"
  if( ! grepl(")", re[correl.w]) )
  {
    if( grepl(")", re[correl.w+1]) )
    {
      correl.t <- paste(re[correl.w], re[correl.w+1], sep=",")
      # check for additional ivs, strip them out
      correl.t <- unlist( strsplit( correl.t, '\\+') )
      if(length(correl.t)>1) moreiv <- gsub(" ", "", correl.t[2:length(correl.t)])
      correl.t <- correl.t[1]
      correl.t <- gsub("))", ")", correl.t)
      correl.t <- gsub(" ", "", correl.t)
    }
    else warning("Cannot extract a correlation structure frome ",
                 re, ", using NULL instead")
  }
  if( grepl(")", re[correl.w]) )
  {
    correl.t <- unlist( strsplit(re[correl.w], ")"))[1]
    correl.t <- gsub(" ", "", unlist( strsplit(correl.t, '='))[2])
  }
  if(correl.t=="NULL") correlation <- correl.t

  # test whether the correlation structure is valid
  if(correl.t!="NULL")
  {
    correl.t <- unlist(strsplit(correl.t, "="))
    if(length(correl.t) > 2)
    {
      correl.t <- paste(correl.t[2:length(correl.t)], collapse="=")
    }
    correl.t <- gsub(" ", "", correl.t)
    correl.test <- try(iscorStruct((correl.t)), silent = TRUE)
    if( "try-error" %in% is(correl.test) )
    {
      correl.test <- iscorStruct(paste("nlme::", correl.t, sep=""))
    }
    if( correl.test) correlation <- correl.t
    if(!correl.test)
    {
      warning(correl.t, " is not a correlation structure, using NULL instead.")
      correlation <- "NULL"
    }
  }

  # add additional ivs
  ivs <- unique( c(ivs, moreiv) )


  return( list(ids=ids, dv=dv, time=time, ivs=ivs, correlation=correlation,
               method=method) )

}

#' frmToChar - convert formula to text
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @param x Any formula
#'
#' @keywords internal
frmToChar <- function(x)
{
  as.character( attr(terms(x), "variables") )[-1L]
}

# this will only work if all side effects apply in the parent environment;
# sourcing code won't work unless you can intelligently find the install folder
#' dohtp -
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @keywords internal
dohtp <- function()
{
  # this is supposedly taboo but I've been unable to work around it b/c we
  # cannot :: the %dopar% operator
  require(foreach)
  # parralelization
  funcs <- c("mean") # c(".eds") -- not importing from PersonAlytic correctly
  cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
  snow::clusterExport(cl, funcs)
  doSNOW::registerDoSNOW(cl)
  pkgs  <- c("gamlss", "nlme")
}

#' alf - attach(as.list(formals(x))), a debugging tool
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @keywords internal
alf <- function(x)
{
  attach( as.list( formals(x) ) )
}


#' to01 - function to convent any variable to the [0,1] range
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @param x A numeric vector.
#' @param na.rm Logical, should missing data be excluded when calculating min/max
#' @references Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods, 11(1), 54.
#' @export
to01 <- function(x, na.rm=TRUE)
{
  (x - min(x, na.rm=na.rm))/(max(x, na.rm=na.rm)-min(x, na.rm=na.rm))
}

#' subdat - function to take a subset of data with column select
#' @author Stephen Tueller \email{Stueller@@rti.org}
#' @param subgroup Logical vector.
#' @param data data.frame.
#' @param formula a formula.
#' @keywords internal
subdat <- function(data, subgroup, formula)
{
  if(is.null(subgroup)) subgroup <- rep(TRUE, nrow(data))
  na.omit( subset(data, subgroup, all.vars(formula)) )
}
