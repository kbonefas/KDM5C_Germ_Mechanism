#code to generate the pdf from the rmd file

library('tinytex')
library('rmarkdown')
render('submission/manuscript.Rmd')

