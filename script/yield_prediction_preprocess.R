
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(rrBLUP)

# read data
maizeNE=read.csv("Bugeater2020_merged_v6_preprocess.csv")
maizeMI=read.csv("2020maizeMI_NEname_preprocess.tsv",sep="\t")

# MI flowering time
maizeMI$plantdate="2020-05-26"
maizeMI$Anthesisday=as.numeric(as.Date(as.character(maizeMI$Anthesis),format = "%Y-%m-%d")-as.Date(as.character(maizeMI$plantdate),format="%Y-%m-%d"))
# recode N Y to 0 1
maizeNE$PoorGerm=as.numeric(recode(maizeNE$PoorGerm, N=0,Y=1))
maizeNE$PotentialPollenContamination=as.numeric(recode(maizeNE$PotentialPollenContamination,  N=0,Y=1))

# get genotypes and the numeric columns
maizeNE.num=cbind.data.frame(maizeNE[4],apply(maizeNE[c(5:ncol(maizeNE))],2,as.numeric))
maizeMI.num=cbind.data.frame(maizeMI[6],apply(maizeMI[c(14:(ncol(maizeMI)-3),ncol(maizeMI))],2,as.numeric))

# average by genotype
avgMI=maizeMI.num%>%group_by(GenotypeID=Name)%>%summarise(across(everything(),list(mean=MImean)))
avgNE=maizeNE.num%>%group_by(GenotypeID)%>%summarise(across(everything(),list(mean=NEmean)))
avg_merge=merge(avgMI,avgNE,by="GenotypeID")
write.table(avgNE,"2020NE_avg.csv",sep = ",",quote = F,row.names = F)
write.table(avgMI,"2020MI_avg.csv",sep = ",",quote = F,row.names = F)
write.table(avg_merge,"2020merge_avg.csv",sep = ",",quote = F,row.names = F)

# preprocess NA
avg_merge[colSums(is.na(avg_merge)) > 0.1*nrow(avg_merge)] =NULL
with_na=colnames(avg_merge)[colSums(is.na(avg_merge)) > 0]
for (col in with_na) {
  avg_merge[is.na(avg_merge[,col]),col]=mean(avg_merge[,col],na.rm=T)
}
write.table(avg_merge,"2020merge_avg_narm.csv",sep = ",",quote = F,row.names = F)
