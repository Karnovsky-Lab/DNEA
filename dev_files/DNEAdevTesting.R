
BiocManager::install(c("RCy3", "org.Hs.eg.db", "graphite", "graph", "genefilter", "AnnotationDbi", "BiocParallel"))
devtools::install("~/Documents/Karnovsky_lab/DNEAdev/")
library(DNEAdev)
library(BiocParallel)

BP_plan <- MulticoreParam(workers = 4, RNGseed = 101)
set.seed(101)

dat <- read.csv('~/Documents/Karnovsky_lab/DNEAproject//published_files/adjT1DplasmaLastVisitpaired_04252023.csv')
rownames(dat) <- dat$sample
dat<- dat[,-1]

object<-createDNEAobject(project_name = 'testing', expression_data = dat, case = 'DM:case', control = 'DM:control')
object <- BICtune(object = object, BPPARAM = BP_plan)
object <- stabilitySelection(object = object, subSample = FALSE, nreps = 4, BPPARAM = BP_plan)
object <- getNeworks(object = object)
object <- runConsensusCluster(object = object, tau = 0.5)
object <- runNetGSA(object)
