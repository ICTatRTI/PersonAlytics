#' simplify exists(deparse(substitute(x)))
#'
#' @param x. The name of any type of object
#'
#' @return Logical. Does the object exist?
#'

eds <- function(x)
{
  exists( deparse( substitute(x) ) )
}
