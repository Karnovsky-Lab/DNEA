#test old
dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
ncond = 2
group1 = 'IA:control'
group2 = 'IA:case'

testingX<-readRDS('~/Documents/Karnovsky_lab/Filigree/DNEA-script/trainX.rds')
testingY<-readRDS('~~/Documents/Karnovsky_lab/Filigree/DNEA-script/trainY.rds')
all(testingX == t(do.call(cbind, object@Dataset_summary$scaled_separated_conditions)))
isTRUE(all.equal(testingX, trainX2))
for (i in 1:ncol(trainX2)){
  if(isTRUE(all.equal(trainX2[,i],testingX[,i]))){
    TRUE
  }else{
    print(i)
  }
}
all(testingY == c(rep(1, ncol(object@Dataset_summary$scaled_separated_conditions[[1]])),
                  rep(2, ncol(object@Dataset_summary$scaled_separated_conditions[[2]]))))
testLambda<-readRDS("~/Documents/Karnovsky_lab/Filigree/DNEA-script/testLambda.rds")
all(testLambda == lambda.guo)
bicTest<-readRDS('~/Documents/Karnovsky_lab/Filigree/DNEA-script/bic.guo.rds')
for (i in 1:length(bicTest)){
  if(as.numeric(bicTest[[i]][1]) == as.numeric(object@BIC$bic.guo[[i]][1])){
    print(paste0(i,' BIC same'))
  }else{
    print(paste0(i, 'BIC NOT same'))
  }
  if(as.numeric(bicTest[[i]][2]) == as.numeric(object@BIC$bic.guo[[i]][2])){
    print(paste0(i,' liklihood same'))
  }else{
    print(paste0(i, 'liklihood NOT same'))
  }
}
#test new
dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
#dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/T1D_PLASMA_first_visit_adjusted.csv')
rownames(dat)<-dat$sample
dat<-dat[,-1]

# main.seed <- 101
# set.seed(main.seed)
object <- createDNEAobject(Project.Name = 'TESTING', NormalExpression = dat, control = 'IA:control', case = 'IA:case')
object <- BICtune(object = object)
object<- StabilitySelection(object = object,nreps = 2, main.seed = 101)

# nCores = 2
# cl <- parallel::makeCluster(nCores, type = "PSOCK")
# registerDoParallel(cl)
# clusterMap(cl, x = 1:10, function(x, y) seq_len(x) + y,
#            c(a =  1, b = 2, c = 3), c(A = 10, B = 0, C = -10))
# stopCluster(cl)

load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_BIC_tuning.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_adjacency_matrices.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_dataset_summary.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_stable_networks.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_dataset_summary.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_joint_graph.rda')
load('~/Documents/Karnovsky_lab/adjusted_output/PLASMA/T1D/adjFULL/TEDDY_DM_PLASMA_full_adjusted_netgsa_results.rda')
