library(devtools)
library(roxygen2)
library(available)
devtools::load_all()
#devtools::create('DNEAdev')
available('DNEA')

usethis::use_mit_license('DNEADEV')


devtools::document()
devtools::check()
