# you'll need to @export these if you want them viablem by client, rename all
# with a . first

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
frmToChar <- function(x)
{
  as.character( attr(terms(x), "variables") )[-1L]
}
catt <- function(text)
{
  cat(text, file='model.txt', append=TRUE)
}
catu <- function(text)
{
  cat(text, file='uncond.txt', append=TRUE)
}
catd <- function(text)
{
  cat(text, file='uncgrowth.txt', append=TRUE)
}
catt4 <- function(text)
{
  cat(text, file='model4.txt', append=TRUE)
}
catu4 <- function(text)
{
  cat(text, file='uncond4.txt', append=TRUE)
}
catd4 <- function(text)
{
  cat(text, file='uncgrowth4.txt', append=TRUE)
}
tr <- function(expr)
{
  t0 <- Sys.time()
  attempt <- try(expr, silent=TRUE)
  print(Sys.time()-t0)
  return(attempt)
}
