library(tidyverse)

# read data
maizeNE=read.csv("data/Bugeater2020_merged_v6_preprocess.csv")
maizeMI=read.csv("data/2020maizeMI_NEname_preprocess.tsv",sep="\t")

# MI flowering time
maizeMI$plantdate="2020-05-26"
maizeMI$Anthesisday=as.numeric(as.Date(as.character(maizeMI$Anthesis),format = "%Y-%m-%d")-as.Date(as.character(maizeMI$plantdate),format="%Y-%m-%d"))
# recode N Y to 0 1
maizeNE$PoorGerm=as.numeric(recode(maizeNE$PoorGerm, N=0,Y=1))
maizeNE$PotentialPollenContamination=as.numeric(recode(maizeNE$PotentialPollenContamination,  N=0,Y=1))

# get genotypes and the numeric columns
maizeNE.num=cbind.data.frame(maizeNE[4],apply(maizeNE[c(5:ncol(maizeNE))],2,as.numeric))
maizeMI.num=cbind.data.frame(maizeMI[6],apply(maizeMI[c(10,14:(ncol(maizeMI)-3),ncol(maizeMI))],2,as.numeric))

# average MI by plot
maizeMI.num = maizeMI.num %>% transmute(GenotypeID = Name,
                                        EarHeight = rowMeans(maizeMI.num[,c("EarHeightP1", "EarHeightP2")], na.rm = T), 
                                     EarLeaf = rowMeans(maizeMI.num[,c("EarLeafP1", "EarLeafP2")], na.rm =T), 
                                     LeafNumber = rowMeans(maizeMI.num[,c("FlagLeafP1", "FlagLeafP2")], na.rm = T), 
                                     EarPcttoFlag = rowMeans(maizeMI.num[,c("Ear.percent.to.Flag.P1", "Ear.percent.to.Flag.P2")], na.rm =T), 
                                     LargestLeafLength = rowMeans(maizeMI.num[,c("LargestLeafLengthP1", "LargestLeafLengthP2")], na.rm =T), 
                                     LargestLeafWidth = rowMeans(maizeMI.num[,c("LargestLeafWidthP1", "LargestLeafWidthP2")], na.rm = T), 
                                     LargestLeafNumber = rowMeans(maizeMI.num[,c("LargestLeafNumberP1", "LargestLeafNumberP2")], na.rm = T), 
                                     NumberofLeavesLargesttoEar = Number.of.leaves.largest.to.ear, 
                                     TasselHeight = rowMeans(maizeMI.num[,c("TasselHeightP1", "TasselHeightP2")], na.rm = T),
                                     PlantHeight = rowMeans(maizeMI.num[,c("FlagHeightP1", "FlagHeightP2")], na.rm = T), 
                                     WeightofDriedLeafDIsksMG = Weight.of.Dried.Leaf.Disks..mg., 
                                     NumberofDisks = Number.of.disks, 
                                     SpecificLeafWeightMG = Specific.Leaf.Weight..mg.disc., 
                                     Moisture = Moisture...., 
                                     DaysToPollen = Anthesisday, Yield = Total.Weight.of.kernels..g.*2,
                                     SeedsNumber = final...of.seeds.in.packet)
maizeNE.num = maizeNE.num %>% mutate(LeafNumber = rowMeans(maizeNE.num[,c("ExtantLeafNumber1","ExtantLeafNumber2")], na.rm = T)) %>%
  mutate(ExtantLeafNumber1 = NULL, ExtantLeafNumber2 = NULL)

# average by genotype
avgMI=maizeMI.num%>%group_by(GenotypeID)%>%summarise(across(everything(),list(MImean=mean))) # MImean = mean, or mean = MImean
avgNE=maizeNE.num%>%group_by(GenotypeID)%>%summarise(across(everything(),list(NEmean=mean)))
avg_merge=merge(avgMI,avgNE,by="GenotypeID")
write.table(avgNE,"data/2020NE_avg.csv",sep = ",",quote = F,row.names = F)
write.table(avgMI,"data/2020MI_avg.csv",sep = ",",quote = F,row.names = F)
write.table(avg_merge,"data/2020merge_avg.csv",sep = ",",quote = F,row.names = F)

# preprocess NA
avg_merge[colSums(is.na(avg_merge)) > 0.1*nrow(avg_merge)] =NULL # remove columns that contain >10% missing
with_na=colnames(avg_merge)[colSums(is.na(avg_merge)) > 0]
for (col in with_na) {
  avg_merge[is.na(avg_merge[,col]),col]=mean(avg_merge[,col],na.rm=T)
}
write.table(avg_merge,"data/2020merge_avg_narm.csv",sep = ",",quote = F,row.names = F)
