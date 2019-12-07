PersonAlyticsStartupMessage <- function()
{
  # https://www.askapache.com/online-tools/figlet-ascii/
  msg <- c(paste0(
"
'||'''|,
 ||   ||
 ||...|' .|''|, '||''| ('''' .|''|, `||''|,
 ||      ||..||  ||     `'') ||  ||  ||  ||  ====
.||      `|...  .||.   `...' `|..|' .||  ||.

         /.\      '||`            ||
        // \\      ||             ||     ''
       //...\\     ||  '||  ||` ''||''   ||  .|'', (''''
      //     \\    ||   `|..||    ||     ||  ||     `'')
    .//       \\. .||.      ||    `|..' .||. `|..' `...'
                        ,  |'
                        ''

Version ",
packageVersion("PersonAlytics")),
"\n\nType 'citation(\"PersonAlytics\")' for citing this R package in publications.",
"\n\nType 'vignette(package=\"PersonAlytics\")' for the user's guide.\n")
  return(msg)
}


.onAttach <- function(lib, pkg)
{
  msg <- PersonAlyticsStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'PersonAlytics' version", packageVersion("PersonAlytics"))
  packageStartupMessage(msg)
  invisible()
}
