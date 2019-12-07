# http://r-pkgs.had.co.nz/vignettes.html
if(1==2)
{
  # start a vignette
  #devtools::use_vignette("test")

  # render to word for coathur comments, also html for formatting
  library(rmarkdown)
  render("vignettes//PersonAlytics_Users_Guide.Rmd", word_document(toc=TRUE))
  render("vignettes//PersonAlytics_Users_Guide.Rmd", html_document(toc=TRUE))


  devtools::build_vignettes()
  devtools::install()
}
