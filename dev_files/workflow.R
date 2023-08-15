#dat<- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
#dat <- read.csv('~/Documents/Karnovsky_lab/Datasets/TEDDY/adjusted/PLASMA/IA_PLASMA_first_visit_adjusted_V2.csv')
library(BiocParallel)
BP_plan <- SerialParam(RNGseed = 417)
set.seed(417)

# dat <- read.csv('~/Documents/Karnovsky_lab/DNEAproject/published_files/adjT1DplasmaLastVisitpaired_04252023.csv')
# dat <- read.csv('~/Documents/Karnovsky_lab/DNEAproject/published_files/adjT1DplasmaLastVisitpaired_non-transformed_07122023.csv')

# dat <- read.csv('~/Documents/Karnovsky_lab/DNEAproject/published_files/adjT1DplasmaLastVisitAll-nontransformed-07122023.csv')
#
# dat <- dat[, !grepl("nist", colnames(dat))]
# rownames(dat) <- dat$sample
# group_labels <- dat$group
# names(group_labels) <- dat$sample
# group_labels <- factor(group_labels, levels = c("DM:control", "DM:case"))
# group_labels[1:10]
# dat<- dat[,-c(1,2)]
# dat <- t(dat)
# TEDDY <- dat
# object<-createDNEAobject(project_name = 'testing', expression_data = dat, group_labels = group_labels)

data("TEDDY")
data("T1Dmeta")
if(!(all(rownames(T1Dmeta) == colnames(TEDDY)))) stop("problem!")
group_labels <- T1Dmeta$group
names(group_labels) <- rownames(T1Dmeta)
object <- createDNEAobject(project_name = 'testing', expression_data = TEDDY, group_labels = group_labels)

#test addExpressionData
TEDDY <- t(TEDDY)
dat <- list('DM:control' = TEDDY[T1Dmeta$group == "DM:control",],
            'DM:case' = TEDDY[T1Dmeta$group == "DM:case",])

#log-transform and median center the expression data without scaling
newdat <- NULL
for(cond in dat){
  for(i in 1:ncol(cond)){
    my_median = median(cond[, i], na.rm = TRUE)
    my_range = range(cond[, i], na.rm = TRUE)
    scale_factor = max(abs(my_range-my_median))
    cond[, i] <- (cond[, i] - my_median) / scale_factor
  }
  newdat <- rbind(newdat, cond)
}

#reorder to match TEDDYresults
newdat <- newdat[sampleNames(TEDDYresults), featureNames(TEDDYresults)]
newdat <- t(newdat)
#add data
TEDDYresults <- addExpressionData(object = object, data = newdat)

# #test node collapsing
# TEDDY_groups <- data.frame(features = rownames(expressionData(TEDDYresults, normalized = FALSE)),
#                            groups = rownames(expressionData(TEDDYresults, normalized = FALSE)),
#                            row.names = rownames(expressionData(TEDDYresults, normalized = FALSE)))
#
# TEDDY_groups$groups[TEDDY_groups$groups %in% c("isoleucine", "leucine", "valine")] <- "BCAAs"
# TEDDY_groups$groups[grep("acid", TEDDY_groups$groups)] <- "fatty_acids"
# object <- reduceFeatures(object, method = "knowledge", correlation_threshold = 0.7, feature_groups = TEDDY_groups)
# object <- reduceFeatures(object, method = "correlation", correlation_threshold = 0.9)
# object <- reduceFeatures(object, method = "hybrid", correlation_threshold = 0.7, feature_groups = TEDDY_groups)
object <- BICtune(object = object, BPPARAM = BP_plan)
object <- stabilitySelection(object = object, subSample = FALSE, nreps = 4, BPPARAM = BP_plan)


object <- getNetworks(object = object)

object <- filterNetworks(object, pcor = 0.166)
# object <- filterNetworks(object, top_percent_edges = 0.2)
object <- clusterNet(object = object, tau = 0.5)
object <- runNetGSA(object)
plotNetworks(object, type = "group_networks", subtype = "DM:case", layout_func = layout_nicely,
             label_font = 2, label_size = 0.5, node_size = 7)
# plotNetworks(object, type = "subnetworks", subtype = 1, layout_func = layout_nicely,
#              font = 2, label_size = 0.5, node_size = 20)
plotNetworks(object, type = "subnetworks", subnetwork = 1)

data(TEDDYresults)

#simulate group labels
TEDDY_groups <- data.frame(features = rownames(expressionData(object, normalized = FALSE)),groups = rownames(expressionData(object, normalized = FALSE)),row.names = rownames(expressionData(object, normalized = FALSE)))

TEDDY_groups$groups[TEDDY_groups$groups %in% c("isoleucine", "leucine", "valine")] <- "BCAAs"
TEDDY_groups$groups[grep("acid", TEDDY_groups$groups)] <- "fatty_acids"

object <- reduceFeatures(object, method = "hybrid", correlation_threshold = 0.7, feature_groups = TEDDY_groups)

#check correlations for TEDDY_groups
dat_cor <- cor(t(dat))
dat_cor <- tidyr::pivot_longer(data.frame(metab1 = rownames(dat_cor), dat_cor), cols = colnames(dat_cor), names_to = "metab2",)
dat_cor <- dat_cor[dat_cor$value != 1, ]
dat_cor <- dat_cor[abs(dat_cor$value) >= 0.7,]
dat_cor <- dat_cor[order(dat_cor$metab1, dat_cor$metab2),]

T1Dmeta<-readRDS("~/Documents/Karnovsky_lab/DNEAproject/published_files/T1D_meta.rda")
T1Dmeta <- T1Dmeta[colnames(dat),]
all(rownames(T1Dmeta) == colnames(dat))








