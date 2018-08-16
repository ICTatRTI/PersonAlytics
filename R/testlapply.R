if(1==2)
{
  x = as.list(rep(100,4))
  x = lapply(x, function(x) rnorm(100))
  y = lapply(x, function(x) data.frame(x=x, y=x + rnorm(length(x),0,.3)))

  lapply(y, function(x) lm(x$y~x$x))

  names(y) <- paste('matrix', 1:length(y), sep='')
  y <- lapply(y, function(x) x$name <- names(x))
  #ydf <- do.call()
}


