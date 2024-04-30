library(tidyverse)
library(dplyr)
library(zoo)
library(janitor)
library(caret)

#path to save data
filepath <- "[FILEPATH TO SAVE /DNEA/data/scripts/]"
output_filepath <- "[FILEPATH TO OUTPUT FOLDER"
setwd(filepath)
source(paste0(filepath, '/LMadjust.R'))

#######################################################################################################
#*Read in data*#
################

#read in T1D metadata
IAmeta <- readRDS(paste0(filepath, "/IAmeta.rds"))

#read in expression data
plasma_expression_dat <- readRDS(paste0(filepath, "/ST001386_AN002314_Results_and_metabolite_metadata.rds"))

#clean metabolite names and save in variable
TEDDY_metabnames_original <- colnames(plasma_expression_dat)[-c(1,2)]
TEDDY_metabnames <- make_clean_names(TEDDY_metabnames_original)

#replace colnames with cleaned metab names
colnames(plasma_expression_dat) <- c("samples", "factors", TEDDY_metabnames)

#grab the sample names
TEDDY_samplenames <- paste0("sample_",plasma_expression_dat[,1])

#replace rownames with sample names
rownames(plasma_expression_dat) <- TEDDY_samplenames
plasma_expression_dat <- plasma_expression_dat[,-c(1,2)]

#convert expression data to numeric
plasma_expression_dat <-as.data.frame(sapply(plasma_expression_dat, as.numeric))
rownames(plasma_expression_dat) <- TEDDY_samplenames

#remove nist samples
plasma_expression_dat <- plasma_expression_dat[, !grepl("_nist", colnames(plasma_expression_dat))]

#clear space
rm(TEDDY_metabnames, TEDDY_metabnames_original, TEDDY_samplenames)
#######################################################################################################
#*create sampled datasets for analysis*#
########################################

#remove all un-related samples
IAmeta <- IAmeta[!(is.na(IAmeta$group)),]
IAmeta <- IAmeta[IAmeta$group != 'BOTH',]
IAmeta <- IAmeta[complete.cases(IAmeta),]

#grab expression data to match filtered IAmeta
plasma_expression_dat <- plasma_expression_dat[IAmeta$sample,]

#replace NA's with 0
plasma_expression_dat[is.na(plasma_expression_dat)] <- 0

#calculate vector of length ncol(plasma_expression_dat) that contains the # of samples
#with no expression of that metabolite
metab_missingness <- apply(plasma_expression_dat, 2, function(x) sum(x == 0))

#create boolean vector that is TRUE if metabolite is expressed in >= 70% of samples
IAkeep<-unlist(lapply(metab_missingness, function(x) x <= ceiling(nrow(plasma_expression_dat) * 0.3)))

#remove metabolites missing in more than 30% of samples
plasma_expression_dat<-plasma_expression_dat[,IAkeep]

#clear space
rm(metab_missingness, IAkeep)

#impute missing values with median of column
for(i in 1:ncol(plasma_expression_dat)){

  metab_col<-plasma_expression_dat[,i]
  metab_col[metab_col == 0] <- median(metab_col)
  plasma_expression_dat[,i] <- metab_col

}
#clear space
rm(metab_col, i)

#order expression data and metadata to merge
plasma_expression_dat <- plasma_expression_dat[IAmeta$sample,]
if(!all(IAmeta$sample == rownames(plasma_expression_dat))) stop("metadata sample order != expression data sample order!")

#######################################################################################################
#*adjust data for Age and Sex*#
###############################

#create covariates dataframe
IA_covariates <- IAmeta[,c(4,5)]
rownames(IA_covariates) <- IAmeta$sample

adj_plasma_expression_dat <- LMadjust(dataset = plasma_expression_dat,
                                      covariates = IA_covariates,
                                      log_transform = FALSE,
                                      scale_data = FALSE)

#clear space
rm(IA_covariates, plasma_expression_dat)

if(!all(IAmeta$sample == rownames(adj_plasma_expression_dat))) stop("there was an error adjusting data!")
adj_plasma_expression_dat<-cbind.data.frame(IAmeta, adj_plasma_expression_dat)

#save whole dataset
write.csv(adj_plasma_expression_dat[,-c(seq(5))], paste0(output_filepath,'adj_TEDDY_IA_all.csv'), row.names = FALSE)
#######################################################################################################
#*downsample dataset*#
######################

#set seed for downsampling
set.seed(1017)

#set downsample percents
downsample_percent <- c(.9, .75, .50, .25, .1, .05)
names(downsample_percent) <- c("90p", "75p", "50p", "25p", "10p", "5p")

for(down_p in seq(length(downsample_percent))){

  #randomly sample the data
  trainIndex <- createDataPartition(adj_plasma_expression_dat$group, p = downsample_percent[down_p], list = FALSE, times = 1)
  downsampled_data <- adj_plasma_expression_dat[trainIndex,]

  #check to make sure the split remains the same
  sum(adj_plasma_expression_dat$group == 'IA:control') / nrow(adj_plasma_expression_dat)
  sum(downsampled_data$group == 'IA:control') / nrow(downsampled_data)

  write.csv(downsampled_data[,-c(seq(5))], paste0(output_filepath, "/adj_TEDDY_IA_", names(downsample_percent)[down_p],".csv"), row.names = FALSE)
}

#clear space
rm(downsampled_data, trainIndex, down_p, downsample_percent)
#######################################################################################################
#*Take only the first visit*#
#############################

if(!all(IAmeta$sample == rownames(adj_plasma_expression_dat))) stop("meta data and expression data no longer match!")

#subject names in T1D
subjects <- unique(adj_plasma_expression_dat$subject)

#initialize dataframe and loop through data - only keep first visits
IA_first_visit<-data.frame()
for (ia_subject in subjects){

  subject_samples <- adj_plasma_expression_dat[adj_plasma_expression_dat$subject == ia_subject,]
  subject_samples <- subject_samples[order(subject_samples$Age),]
  subject_samples <- subject_samples[1,]

  #add first sample for this subject to the output
  IA_first_visit<-rbind.data.frame(IA_first_visit,subject_samples)
  rm(subject_samples)
}

write.csv(IA_first_visit[, -c(seq(5))], paste0(output_filepath,"/adj_TEDDY_IA_FirstVisit.csv"), row.names = FALSE)

#clear space
rm(IA_first_visit, subjects, ia_subject)
#######################################################################################################
#*Take only the last visit*#
############################

#check order again
if(!all(IAmeta$sample == rownames(adj_plasma_expression_dat))) stop("meta data and expression data no longer match!")

#subject names in T1D
subjects <- unique(adj_plasma_expression_dat$subject)

#initialize dataframe and loop through data - only keep last visits
IA_last_visit<-data.frame()

for (ia_subject in subjects){
  subject_samples <- adj_plasma_expression_dat[adj_plasma_expression_dat$subject == ia_subject,]
  subject_samples <- subject_samples[order(subject_samples$Age),]

  #include control sample so long as they have been in the study for more than a year
  if(all(subject_samples$group == "IA:control")){
    if(subject_samples$Age[nrow(subject_samples)] - subject_samples$Age[1] >= 365) {
      subject_samples <- subject_samples[nrow(subject_samples),]

      #add to output
      IA_last_visit<-rbind.data.frame(IA_last_visit,subject_samples)
    }
  }

  #include case sample as long as they have been in the study more than a year
  if(all(subject_samples$group == "IA:case")){
    if(subject_samples$Age[nrow(subject_samples)] - subject_samples$Age[1] >= 365){
      if(subject_samples$Age[nrow(subject_samples)] - as.numeric(subject_samples$Endpoint1[nrow(subject_samples)]) >= -28){
        subject_samples <- subject_samples[nrow(subject_samples),]
        IA_last_visit<-rbind.data.frame(IA_last_visit,subject_samples)
      }
    }
  }
  rm(subject_samples)
}

#save data
write.csv(IA_last_visit[,-c(seq(5))], paste0(output_filepath,"/adj_TEDDY_IA_LastVisit.csv"), row.names = FALSE)

#clear space
rm(IA_last_visit, subjects, ia_subject)
######################################################################################################







