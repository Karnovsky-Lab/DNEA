#dat<- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
dat <- read.csv('~/Documents/Karnovsky_lab/DNEAdev/data/TEDDYplasmaIA.csv')
rownames(dat) <- dat$sample
dat<- dat[,-1]

object<-createDNEAobject(project_name = 'testing', scaled_expression_data = dat, case = 'IA:case', control = 'IA:control')
object <- BICtune(object = object, nCores = 4)
object<- stabilitySelection(object = object, runParallel = TRUE, subSample = FALSE, nreps = 4, nCores = 4)
object <- getNeworks(object = object)
object <- runConsensusCluster(object = object, tau = 0.5)
object <- runNetGSA(object)





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







