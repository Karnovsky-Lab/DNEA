#
#
#path to save data
filepath <- "[FILEPATH TO SAVE /DNEA/data/scripts/]"
filepath <- "~/Documents/Karnovsky_lab/DNEAproject/DNEAdata/"
output_filepath <- "[FILEPATH TO OUTPUT FOLDER"
setwd(filepath)



#################################
#**Read in and set up metadata**#
#################################

#meta data file
meta<-readRDS("TEDDY_all_metadat.rds")
meta$Sample.name <- paste0("sample_",meta$Sample.name)

#add straightforward column for T1D subset
t1d<-data.frame(meta)
t1d<-mutate(t1d, group = NA)
for (i in 1:nrow(t1d)){
  if (t1d$t1d_control[i] == '1'){
    t1d$group[i]<-'DM:control'
  }
  if(t1d$t1d_case[i] == '1'){
    t1d$group[i]<-'DM:case'
  }
  if ((t1d$t1d_control[i] == '1') & (t1d$t1d_case[i] == '1')){
    t1d$group[i]<-'BOTH'
  }
}

t1d<-cbind.data.frame(t1d$Subject.name,t1d$t1d_endptage_1, t1d$t1d_endptage_2,t1d$Age.in.Days,t1d$Sex,t1d$Sample.name,t1d$group)
colnames(t1d)<-c('subject','Endpoint1','Endpoint2','Age','Sex','sample','group')
rownames(t1d) <- t1d$sample
saveRDS(t1d, paste0(output_filepath, "/T1Dmeta.rds"))

#add straightforward column for IA subset
IA<-data.frame(meta)
IA<-mutate(IA, group = NA)

for (i in 1:nrow(IA)){
  if (IA$ia_control[i] == '1'){
    IA$group[i]<-'IA:control'
  }
  if(IA$ia_case[i] == '1'){
    if(as.numeric(IA$Age.in.Days[i]) <= as.numeric(IA$ia_endptage_1[i])){
      IA$group[i]<-'IA:case'
    }else{
      NA
    }
  }
  if ((IA$ia_control[i] == '1') & (IA$ia_case[i] == '1')){
    IA$group[i]<-'BOTH'
  }
}

IA <- cbind.data.frame(IA$Subject.name,IA$ia_endptage_1, IA$ia_endptage_2,IA$Age.in.Days,IA$Sex,IA$Sample.name,IA$group)
colnames(IA)<-c('subject','Endpoint1','Endpoint2','Age','Sex','sample','group')
rownames(IA) <- IA$sample
saveRDS(IA, paste0(output_filepath, "/IAmeta.rds"))








