library(devtools)
library(roxygen2)
library(available)
library(usethis)
library(BiocCheck)
devtools::load_all("~/Documents/Karnovsky_lab/DNEA/", reset = TRUE, recompile = TRUE)
#devtools::create('DNEAdev')
available('DNEA')

#add packages
usethis::use_mit_license('DNEADEV')
usethis::use_description()
usethis:: use_mit_license()
usethis::use_package("Matrix")
# usethis::use_package("corpcor")
usethis::use_package("dplyr")
# usethis::use_package("furrr")
usethis::use_package("gdata")
usethis::use_package("glasso")
# usethis::use_package("glmnet")
usethis::use_package("igraph")
usethis::use_package("janitor")
usethis::use_package("methods")
usethis::use_package("stats")
# usethis::use_package("pbapply")
usethis::use_package("stringr")
# usethis::use_package("zoo")
usethis::use_package("netgsa")
usethis::use_package("BiocParallel")
usethis::use_package("autoimage")
usethis::use_package("utils")
usethis::use_package("withr", type = "Suggests")
# usethis::use_package("future")
usethis::use_package()

#add data
usethis::use_data(TEDDY, overwrite = TRUE)
usethis::use_data(group_labels, overwrite = TRUE)
usethis::use_data(TEDDYresults, overwrite = TRUE)
usethis::use_data(T1Dmeta)

#add to .Rbuildignore
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/dev_files/parallelTest.R", escape = TRUE)
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/dev_files/TEDDYplasmaIA.csv", escape = TRUE)
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/dev_files/workflow.R", escape = TRUE)
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/dev_files/Dev.R", escape = TRUE)
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/.gitignore", escape = TRUE)
use_build_ignore("~/Documents/Karnovsky_lab/DNEA/dev_files/DNEAdevTesting.R", escape = TRUE)

#add tests

usethis::use_testthat()
mkdir("~/Documents/Karnovsky_lab/DNEA/tests/testthat/testdata/")
saveRDS(T1Dmeta, "~/Documents/Karnovsky_lab/DNEA/tests/testthat/testdata/test-T1Dmeta.rds")
saveRDS(TEDDY, "~/Documents/Karnovsky_lab/DNEA/tests/testthat/testdata/test-TEDDY.rds")
saveRDS(TEDDYresults, "~/Documents/Karnovsky_lab/DNEA/tests/testthat/testdata/test-TEDDYresults.rds")
usethis::use_test("utilities-internals.R")
usethis::use_test("utilities-external.R")
# usethis::use_test("MainFunctions.R")
usethis::use_test("BICtune.R")
usethis::use_test("stabilitySelection.R")
usethis::use_test("getNetworks.R")
usethis::use_test("clusterNet.R")
usethis::use_test("runNetGSA.R")
usethis::use_test("createDNEAobject.R")
usethis::use_test("reduceFeatures.R")

BiocCheck::BiocCheck("~/Documents/Karnovsky_lab/DNEA/")

devtools::test()
devtools::document()
devtools::build_manual()
devtools::check()
devtools::build()







message(paste0(names(unweighted_adjacency_matrices)[[1]]," network specific edges: ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
message(paste0(names(unweighted_adjacency_matrices)[[2]]," network specific edges: ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)
message(paste0("Number of edges shared by both networks: ", sum(edge_list$edge == "Both")))
message(paste0("Total number of edges in dataset: ", nrow(edge_list)))
