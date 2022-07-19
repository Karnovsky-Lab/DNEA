#test old
dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
ncond = 2
group1 = 'IA:control'
group2 = 'IA:case'

testing<-readRDS('~/Documents/Karnovsky_lab/DNEAdev/trainX.rds')
#test new
dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
#dat<-read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/T1D_PLASMA_first_visit_adjusted.csv')
rownames(dat)<-dat$sample
dat<-dat[,-1]

main.seed <- 101
set.seed(main.seed)
object <- createDNEAobject(Project.Name = 'TESTING', NormalExpression = dat, control = 'IA:control', case = 'IA:case')
object <- BICtune(object = object)

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
