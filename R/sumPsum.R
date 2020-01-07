#' sumPsum, a function to summarize multiple p-value summaries
#'
#' @export
#'
#' @param wd A path to a directory containing multiple subdirectories, each with
#' one *.txt file produced by the \code{\link{psuite}} function.
#'
#' @return A *.csv file with the names of subfolders of the ‘run’ folders
#' in the first column, the proportion of unadjusted p values < .05 in the
#' second column, and the proportion of adjusted p values < .05 in the third column.
#'
#' #' @author Stephen Tueller \email{stueller@@rti.org}
sumPsum <- function(wd)
{
  dirs <- dir(wd, glob2rx("PersonAlytics output for*"), full.names = TRUE)
  dirn <- dir(wd, glob2rx("PersonAlytics output for*"), full.names = FALSE)
  psum <- list()
  for(i in seq_along(dirs))
  {
    files <- dir(dirs[i], glob2rx("*PAv*.txt"), full.names = TRUE)
    if(length(files) > 1)
    {
      stop("\nThe directory ", dirs[i], " has more than one text file, please remove one of them.\n")
    }
    filec <- scan(files, what = "character", sep="\n", quiet = TRUE)
    psum[[dirn[i]]] <- filec[grepl("Proportion p < alpha =  0.05", filec)]
  }
  psum <- data.frame( do.call(rbind, psum) )
  names(psum) <- c("Before p-value adjustment", "After p-value adjustment")

  write.csv(psum, paste(basename(wd), " p-value summary.csv", sep=""))
}
