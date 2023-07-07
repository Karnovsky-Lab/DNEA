#dat<- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
#dat <- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
BP_plan <- MulticoreParam(workers = 4, RNGseed = 101)
set.seed(101)
# set.seed(101)
#dat <- read.csv('~/Documents/Karnovsky_lab/published_files/adjT1DfullPlasma10262022.csv')
start <- Sys.time()
dat <- read.csv('~/Documents/Karnovsky_lab/DNEAproject/published_files/adjT1DplasmaLastVisitpaired_04252023.csv')
rownames(dat) <- dat$sample
dat<- dat[,-1]

object<-createDNEAobject(project_name = 'testing', expression_data = dat, case = 'DM:case', control = 'DM:control')
object <- BICtune(object = object, BPPARAM = BP_plan)
object <- stabilitySelection(object = object, subSample = FALSE, nreps = 4, BPPARAM = BP_plan)
finish <- Sys.time()
finish - start

object <- getNeworks(object = object)

object <- filterNetworks(object, pcor = 0.3)
# object <- filterNetworks(object, top_percent_edges = 0.2)
object <- runConsensusCluster(object = object, tau = 0.5)
object <- runNetGSA(object)
plotNetworks(object, type = "group_networks")
plotNetworks(object, type = "subnetworks", subnetwork = 1)
load("~/Documents/Karnovsky_lab/DNEAproject/published_files/test/adjT1DplasmaLastVisitpaired_LOG-SCALED_04252023_BIC_tuning.rda")

load("~/Documents/Karnovsky_lab/DNEAproject/published_files/test/adjT1DplasmaLastVisitpaired_LOG-SCALED_04252023_stable_networks.rda")
sel2 <- vector("list", length(networkGroups(object)))
names(sel2) <- networkGroups(object)

#initiate list for stability selection results converted to probabilities
selp <- vector("list", length(networkGroups(object)))
names(selp) <- networkGroups(object)

for (k in 1:length(sel2)){
  sel2[[k]] <- lapply(stab_guo, function(r) r$mat[[k]])
  sel2[[k]] <- Reduce("+", sel2[[k]])

  if (subSample){
    message(paste0("Calculating selection probabilities WITH subsampling for...", names(sel2)[[k]],"..."), appendLF = TRUE)
    selp[[k]] <- sel2[[k]]/(nreps)
  } else {
    message(paste0("Calculating selection probabilities WITHOUT subsampling for...",names(sel2)[[k]],"..."), appendLF = TRUE)
    selp[[k]] <- sel2[[k]]/(2 * nreps)
  }
}

ed<- read.table("~/Documents/Karnovsky_lab/DNEAproject/published_files/test/adjT1DplasmaLastVisitpaired_LOG-SCALED_04252023_edgelist.txt",
                header = TRUE)
ed <- ed[,c(1:4, 7)]
ed2<-edgeList(object)



dat2<- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/LIPID/TEDDY_POS_DM_LIPID_full_adjusted.csv')
rownames(dat2) <- dat2$sample
dat2<- dat2[,-1]
object<- createDNEAobject(project_name = 'testing', expression_data = dat2, case = 'DM:case', control = 'DM:control')
object<-reduceFeatures(object, method = 'correlation', correlation_threshold = 0.3)
#scale TEDDY
dat<- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
rownames(dat) <- dat$sample
dat<- dat[,-1]
conditions<- c(case = 'IA:case', control = 'IA:control')
# dat_numeric <- data.frame(lapply(dat[,-1], as.numeric))
for( i in conditions){
  dat[which(dat$group == i),-1] <- scale(dat[which(dat$group == i),-1])
}
dat<-cbind.data.frame(rownames(dat),dat)
rownames(dat) <- NULL
colnames(dat)[1] <- 'sample'
write.csv(dat,'~/Documents/Karnovsky_lab/DNEAdev/data/TEDDYplasmaIA.csv', row.names = FALSE)


#stability selection, which requires lastar.guo from the previous step
lastar <- object@BIC$optimizedLambda

#split data by condition
listX = lapply(split_by_condition(dat = scaledExpressionData(object),
                                                    condition_levels = object@Dataset_summary[['condition_levels']],
                                                    condition_by_sample = condition(object)), function(d) t(d))

#initialize static variables to pass to workers
init_param <- stabsel_init(listX = listX, nreps = nreps)







