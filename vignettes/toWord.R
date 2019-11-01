# render to word
if(1==2)
{
  library(rmarkdown)

  render("vignettes//PersonAlytics_Users_Guide.Rmd", word_document(toc=TRUE))
  render("vignettes//PersonAlytics_Users_Guide.Rmd", html_document(toc=TRUE))

}
