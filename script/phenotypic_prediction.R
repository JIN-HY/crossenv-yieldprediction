library(randomForest)
library(tidyverse)
library(grid)
library(ggpubr)



avgNE=read.csv("data/2020NE_avg.csv")
avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("data/2020merge_avg_narm.csv")


x.idx = apply(matrix(colnames(avg_merge)), 1, grepl,pattern="NE")
x.idx = c(TRUE,x.idx[-1])
xy=avg_merge[,x.idx]
eartraits <- c("EarsPerPlant_NEmean", "CobWeightGrams_NEmean", "HundredKernelMassGrams_NEmean", 
               "EarLengthCM_NEmean", "EarWidthCM_NEmean", "EarFilledLengthCM_NEmean", 
               "KernelRowNumber_NEmean", "KernelsPerRow_NEmean", "PercentFill_NEmean")
xy <- xy[!colnames(xy) %in% eartraits]

xy$Yield_MImean  <-  avg_merge$Yield_MImean
write.table(xy,"data/2020NEMI_RFfulltraining.csv",sep = ",",quote = F,row.names = F)

# RFresult=data.frame(GenotypeID=xy$GenotypeID,Yield_MIpredicted=xy$Yield_MImean)
# find nas
# colnames(xy)[colSums(is.na(xy)) > 0]
# xy[,colnames(xy)[colSums(is.na(xy)) > nrow(xy)*0.1]]=NULL
# 
# for(i in 1:ncol(xy)){
#   xy[is.na(xy[,i]), i] <- mean(xy[,i], na.rm = TRUE)
# }
#############################################################################################

gnumber <- nrow(xy)
random.order = sample(1:gnumber, gnumber)
k = 5
validate = matrix(random.order, ncol = k)
RFcv5.result=RFresult


# for (i in 1:k) {
#   validation.ids = validate[,i] # select the the validation entries
#   training.ids = as.vector(validate[,-i])
#   
#   xy.training = xy[training.ids,-1]
#   #    x.training = xy[training.ids,ncol(xy)]
#   
#   xy.validate = xy[validation.ids,-1]
#   #    x.validate = xy[validate.ids,ncol(xy)]
#   
#   yield.forest <- randomForest(Yield_MImean ~ ., data = xy.training, importance = TRUE)
#   yield.forest
#   
#   RFcv5.result[validation.ids,"Yield_MIpredicted"]=predict(yield.forest,xy.validate)
# }
# 
# RFcv5_MI=merge(xy,RFcv5.result,by="GenotypeID")
# write.table(RFcv5_MI,paste0("data/2020NE_RFcv5_MI.csv",t),quote = F,row.names = F,sep=",")

for (t in 1:20) { # RFCV
  set.seed(t)
  gnumber <- nrow(xy)
  random.order = sample(1:gnumber, gnumber)
  k = 5
  validate = matrix(random.order, ncol = k)
  xyNE <- xy
  RFcv_NE.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_NEpredicted=0)
  xyNE$Yield_MImean <- NULL
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  xyNE$Yield_NEmean <- avg_merge$TotalGrainMassGrams_NEmean
  for (i in 1:k) {
    validation.ids = validate[,i] # select the the validation entries
    training.ids = as.vector(validate[,-i])
    
    xyNE.training = xyNE[training.ids,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    xyNE.validate = xyNE[validation.ids,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_NEmean ~ ., data = xyNE.training, importance = TRUE)
    yield.forest
    
    RFcv_NE.result[validation.ids,"Yield_NEpredicted"]=predict(yield.forest,xyNE.validate)
  }
  write.table(RFcv_NE.result,paste0("data/RF/RF_5cv.csv",t),row.names = F,quote = F,sep=",")
}

for (t in 1:20) { # RFBT
  set.seed(t)
  gnumber <- nrow(xy)
  sp <- sample(1:gnumber, as.integer(0.8*gnumber))
  xyNE <- xy
  RFbs_NE.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_NEpredicted=0)
  RFbs_NE.result <- RFbs_NE.result[sp,]
  xyNE$Yield_MImean <- NULL
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  xyNE$Yield_NEmean <- avg_merge$TotalGrainMassGrams_NEmean
  xyNE <- xyNE[sp,]
  xyNE.training = xyNE[,-1]
  #    x.training = xy[training.ids,ncol(xy)]
  
  xyNE.validate = xyNE[,-1]
  #    x.validate = xy[validate.ids,ncol(xy)]
  
  yield.forest <- randomForest(Yield_NEmean ~ ., data = xyNE.training, importance = TRUE)
  yield.forest
  
  RFbs_NE.result[,"Yield_NEpredicted"]=predict(yield.forest,xyNE.validate)
  write.table(RFbs_NE.result,paste0("data/RF/RF_bs.csv",t),row.names = F,quote = F,sep=",")
}

ggplot(data=RFcv_NE.result,aes(x=Yield_NEpredicted, Yield_NEmean))+
  geom_point()+
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_cor(
    aes(label =  paste( ..r.label.., sep = "~~~~"))) +
  xlab("Predicted Values")+ylab("Observed NE yield")

########################################################################################################################################
# cross-environment

trait_MI <- c("LeafNumber_MImean", "LargestLeafLength_MImean", "LargestLeafWidth_MImean", "PlantHeight_MImean", "DaysToPollen_MImean")
trait_NE <- c("LeafNumber_NEmean", "LeafLengthCM_NEmean", "LeafWidthCM_NEmean", "PlantHeightCM_NEmean", "DaysToPollen_NEmean")

avg_cross <- avg_merge[names(avg_merge) %in% trait_MI | names(avg_merge) %in% trait_NE]
avgNE_cross <- avg_merge[names(avg_merge) %in% trait_NE]
avgMI_cross <- avg_merge[names(avg_merge) %in% trait_MI]

avgNE_cross <- cbind(avg_merge[,c("GenotypeID", "TotalGrainMassGrams_NEmean")], avgNE_cross)
avgMI_cross <- cbind(avg_merge[,c("GenotypeID", "Yield_MImean")], avgMI_cross)


xnames = c("DaysToPollen", "LeafLength", "LeafWidth", "PlantHeight", "LeafNumber")
colnames(avgNE_cross) = c("GenotypeID", "Yield_NEmean", xnames)
colnames(avgMI_cross) = c("GenotypeID", "Yield_MImean", xnames)

avgMI_cross[,c(3:7)] <- sapply(c(3:7), function(col){lm(avgMI_cross[,col]~avgNE_cross[,col])$fitted.values})

for (t in 1:20) {
  set.seed(t)
  
  gnumber <- nrow(avg_cross)
  random.order = sample(1:gnumber, gnumber)
  k = 5
  validate = matrix(random.order, ncol = k)
  
  
  #RFcrossbs.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=NA, Yield_MImean = avgMI_cross$Yield_MImean)
  RFcrosscv.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=avgMI_cross$Yield_MImean, Yield_MImean = avgMI_cross$Yield_MImean)
  RFcv.result = data.frame(GenotypeID=avgNE_cross$GenotypeID,Yield_NEpredicted=avgNE_cross$Yield_NEmean, Yield_NEmean = avgNE_cross$Yield_NEmean)
  
  for (i in 1:k) {
    validate.ids = validate[,i] # select the the validation entries
    training.ids = as.vector(validate[,-i])
    
    avgNE_cross.training = avgNE_cross[training.ids,-1]
    avgMI_cross.training = avgMI_cross[training.ids,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    avgNE_cross.validate = avgNE_cross[validate.ids,-1]
    avgMI_cross.validate = avgMI_cross[validate.ids,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_NEmean ~ ., data = avgNE_cross.training, importance = TRUE)
    yield.forest
    
    #RFcrossbs.result[training.ids, "Yield_MIpredicted"]=predict(yield.forest,avgMI_cross.training)
    RFcrosscv.result[validate.ids,"Yield_MIpredicted"]=predict(yield.forest,avgMI_cross.validate)
    RFcv.result[validate.ids,"Yield_NEpredicted"]=predict(yield.forest,avgNE_cross.validate)
    
    #write.table(RFcrossbs.result, paste0("data/RF/RF_crossenv_bs.csv", k), quote=F, row.names = F, sep=",")
    
    
    #corr_RFcrossbs = cor(RFcrossbs.result$Yield_MIpredicted,RFcrossbs.result$Yield_MImean); print(corr_RFcrossbs)
    corr_RFcrosscv = cor(RFcrosscv.result$Yield_MIpredicted,RFcrosscv.result$Yield_MImean); print(corr_RFcrosscv)
    corr_RFcv = cor(RFcv.result$Yield_NEpredicted, RFcv.result$Yield_NEmean, use = "complete.obs"); print(corr_RFcv)
  }
  
  write.table(RFcrosscv.result, paste0("data/RF/RF_crossenv_cv.csv", t), quote=F, row.names = F, sep=",")
  write.table(RFcv.result, paste0("data/RF/RF_NE_cv.csv", t), quote=F, row.names = F, sep=",")
  
} # for t


