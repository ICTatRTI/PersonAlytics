#' getPalytic, a function to create a Palytic object from PersonAlytics output
#'
#' @param output Character. The path to a csv file produced by the \code{PersonAlytic} function.
#' @param data Data Frame. The same data provided to your \code{PersonAlytic} run.
#' @param rowNum Numeric. The row number you want to create a \code{Palytic} object for. Set to
#'   \code{NULL} if \code{id}, \code{dv}, \code{target_iv}, and \code{ids}
#'   will be provided, noting that all four of these parameters must be specified except
#'   \code{target_iv}, but only if \code{target_iv} is not a column in \code{output}.
#' @param id Character. The value in the 'id' column of \code{output}.
#' @param dv Character. The value in the 'dv' column of \code{output}.
#' @param target_iv Character. The value in the 'target_iv' column of \code{output}. Leave as
#'   \code{NULL} if this column is absent in \code{output}.
#' @param id Character. The value in the 'id' column of \code{output}.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @examples
#' \dontrun{
#'
#' t1 <- PersonAlytic(output          = 'Test1'     ,
#'                    data            = OvaryICT    ,
#'                    ids             = "Mare"      ,
#'                    dvs             = "follicles" ,
#'                    phase           = "Phase"     ,
#'                    time            = "Time"      ,
#'                    package         = "arma"      ,
#'                    individual_mods = TRUE        )
#'  Mare1 <- getPalytic('Test1_PersonAlytic.csv', rowNum=1, data=OvaryICT)
#'
#' }
#'
#' @export
getPalytic <- function(output, data, rowNum,
                       id=NULL,	dv=NULL, target_iv=NULL, ids=NULL)
{
  paout <- read.csv(output)

  if(is.null(rowNum))
  {
    hasTarg <- "target_iv" %in% names(paout)
    if( hasTarg) rowNum <- which(paout$id==id & paout$dv==dv & paout$ids==ids & paout$target_iv==target_iv)
    if(!hasTarg) rowNum <- which(paout$id==id & paout$dv==dv & paout$ids==ids)
  }

  cat( "standardize <- list(", as.character( paout$standardize[rowNum]), ")", file="temp.R")
  source("temp.R"); file.remove("temp.R")

  ids   = as.character( paout$ids[rowNum] )
  dv    = as.character( paout$dv[rowNum]  )
  time  = as.character( paout$time[rowNum] )
  phase = as.character( paout$phase[rowNum] )
  ivs   = as.character( paout$ivs[rowNum])

  fixed  = formula( as.character( paout$fixed[rowNum] ) )
  random = formula( as.character( paout$random[rowNum] ) )

  interactions <- attr(terms(fixed), "term.labels")
  interactions <- interactions[grepl("\\:", interactions)]
  interactions <- strsplit(interactions, ":")

  vars <- unique( c(ids, dv, time, phase, unlist(ivs),
                    all.vars(fixed), all.vars(random),
                    all.vars(formula)) )
  vars <- vars[!is.na(vars)]
  ivs  <- vars[!vars %in% c(ids, dv, time, phase)]

  # for reverse compatibility with versions that don't have `gamlss.family`
  fam <- NO()
  if(!is.null(paout$gamlss.family))
  {
    fam = as.character( paout$gamlss.family[rowNum] )
    # this fails when the scope is inside the function
    #cat( "library(gamlss)\nfam<<-", fam, "()", sep = "", file="temp.R")
    #source("temp.R"); file.remove("temp.R")
    # but this works both ways
    eval(parse(text=paste("fam<-",fam,"()",sep="")))
  }

  # formula -- only implement if fixed and random are insufficient

  correlation  = as.character( paout$correlation[rowNum] )
  corList <- c("corAR1",
               "corARMA",
               "corCAR1",
               "corCompSymm",
               "corExp",
               "corGaus",
               "corLin",
               "corRatio",
               "corSpher",
               "corSymm",
               "NULL")
  if(!correlation %in% corList)
  {
    if(correlation == "See 'call' column'")
    {
      correlation <- as.character( paout$call[rowNum] )
      if(substr(correlation,1,5)=="arima")
      {
        correlation <- NULL
      }
      # place holder for other options if needed
      else correlation <- NULL
    }
  }

  pa <- Palytic$new(
    data         = data                                     ,
    ids          = ids                                      ,
    dv           = dv                                       ,
    time         = time                                     ,
    phase        = phase                                    ,
    ivs          = ivs                                      ,
    interactions = interactions                             ,
    time_power   = paout$time_power[rowNum]                 ,
    correlation  = correlation                              ,
    family       = fam                                      ,
    fixed        = fixed                                    ,
    random       = random                                   ,
    method       = as.character( paout$method[rowNum] )     ,
    standardize  = standardize                              ,
    alignPhase   = as.character( paout$alignPhase[rowNum] )
  )

  return(pa)

}
