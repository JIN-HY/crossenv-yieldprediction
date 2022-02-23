library(randomForest)
library(tidyverse)
library(rrBLUP)
library(grid)
library(ggpubr)

avgNE=read.csv("data/2020NE_avg.csv")
avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("2020merge_avg_narm.csv")

x.idx = apply(matrix(colnames(avg_merge)), 1, grepl,pattern="NE")
x.idx = c(TRUE,x.idx[-1])
xy=avg_merge[,x.idx]
eartraits <- c("EarsPerPlant_NEmean", "CobWeightGrams_NEmean", "HundredKernelMassGrams_NEmean", 
               "EarLengthCM_NEmean", "EarWidthCM_NEmean", "EarFilledLengthCM_NEmean", 
               "KernelRowNumber_NEmean", "KernelsPerRow_NEmean", "PercentFill_NEmean")
xy <- xy[!colnames(xy) %in% eartraits]

xy$Total.Weight.of.kernels..g._MImean = avg_merge$Total.Weight.of.kernels..g._MImean
write.table(xy,"data/2020NEMI_RFfulltraining.csv",sep = ",",quote = F,row.names = F)

RFresult=data.frame(GenotypeID=xy$GenotypeID,Total.Weight.of.kernels..g._MI=y)
# find nas
# colnames(xy)[colSums(is.na(xy)) > 0]
# xy[,colnames(xy)[colSums(is.na(xy)) > nrow(xy)*0.1]]=NULL
# 
# for(i in 1:ncol(xy)){
#   xy[is.na(xy[,i]), i] <- mean(xy[,i], na.rm = TRUE)
# }
#############################################################################################

y = xy[,"Total.Weight.of.kernels..g._MImean"]
total_pheno = length(y)
set.seed(10) 
random.order = sample(1:total_pheno,total_pheno)
k = 5
validate = matrix(random.order, ncol = k)
cv5.result=RFresult


for (i in 1:k) {
  validation.ids = validate[,i] # select the the validation entries
  training.ids = as.vector(validate[,-i])
  
  xy.training = xy[training.ids,-1]
  #    x.training = xy[training.ids,ncol(xy)]
  
  xy.validate = xy[validation.ids,-1,]
  #    x.validate = xy[validate.ids,ncol(xy)]
  
  yield.forest <- randomForest(Total.Weight.of.kernels..g._MImean ~ ., data = xy.training, importance = TRUE)
  yield.forest
  
  cv5.result[validation.ids,"Total.Weight.of.kernels..g._MI"]=predict(yield.forest,xy.validate)
}

RFcv5_MI=merge(avg_merge,cv5.result,by="GenotypeID")
write.table(RFcv5_MI,"2020NE_RFallcv5_MI.csv",quote = F,row.names = F,sep=",")
ggplot(data=RFcv5_MI,aes(x=Total.Weight.of.kernels..g._MI,Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("Predicted Values")+ylab("Observed MI yield")

########################################################################################################################################

# yield components only

yieldcomp=c("GenotypeID","KernelRowNumber_NEmean","KernelsPerRow_NEmean","HundredKernelMassGrams_NEmean","Total.Weight.of.kernels..g._MImean")
xyyiedcomp=xy[,yieldcomp]

y = xyyieldcomp[,"Total.Weight.of.kernels..g._MImean"]
total_pheno = length(y)
set.seed(10) 
random.order = sample(1:total_pheno,total_pheno)
k = 5
validate = matrix(random.order, ncol = k)
cv5yieldcomp.result=RFresult


for (i in 1:k) {
  validation.ids = validate[,i] # select the the validation entries
  training.ids = as.vector(validate[,-i])
  
  xy.training = xyyieldcomp[training.ids,-1]
  #    x.training = xy[training.ids,ncol(xy)]
  
  xy.validate = xyyieldcom[validation.ids,-1,]
  #    x.validate = xy[validate.ids,ncol(xy)]
  
  yield.forest <- randomForest(Total.Weight.of.kernels..g._MImean ~ ., data = xy.training, importance = TRUE)
  yield.forest
  
  cv5yiedlcomp.result[validation.ids,"Total.Weight.of.kernels..g._MI"]=predict(yield.forest,xy.validate)
}

RFcv5yieldcomp_MI=merge(avg_merge,cv5yieldcomp.result,by="GenotypeID")
write.table(RFcv5yieldcomp_MI,"2020NE_RFallcv5_MI.csv",quote = F,row.names = F,sep=",")
ggplot(data=RFcv5yieldcomp_MI,aes(x=Total.Weight.of.kernels..g._MI,Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("Predicted Values")+ylab("Observed MI yield")

#write.table(predicted_validate_result,"NEyieldcomponent_rrblupcv.csv",quote=F,col.names = T,row.names = F)

varImpPlot(tmp.forest)

################ RF with rrblup dataset #############################

rrblupfull=read.csv("rrblup/2020NEMI_rrblup_full.csv")
rrblup5cv=read.csv("rrblup/2020NEMI_rrblup_5cv.csv")
rrbluponeout=read.csv("rrblup/2020NEMI_rrblup_full.csv")

rrbluponeout_x.idx = apply(matrix(colnames(rrbluponeout)), 1, grepl,pattern="NE")
rrbluponeout_x.idx = c(TRUE,rrbluponeout_x.idx[-1])
rrbluponeout_x = rrbluponeout[,rrbluponeout_x.idx]
rrbluponeout_xy = merge(rrbluponeout_x,RFresult,by="GenotypeID")[,-1]

rrbluponeout_cv5.result = RFresult

for (i in 1:k) {
  validation.ids = validate[,i] # select the the validation entries
  training.ids = as.vector(validate[,-i])
  
  rrbluponeout_xy.training = rrbluponeout_xy[training.ids,-1]
  #    x.training = xy[training.ids,ncol(xy)]
  
  rrbluponeout_xy.validate = rrbluponeout_xy[validation.ids,-1,]
  #    x.validate = xy[validate.ids,ncol(xy)]
  
  yield.forest <- randomForest(Total.Weight.of.kernels..g._MI ~ ., data = rrbluponeout_xy.training, importance = TRUE)
  yield.forest
  
  rrbluponeout_cv5.result[validation.ids,"Total.Weight.of.kernels..g._MI"]=predict(yield.forest,rrbluponeout_xy.validate)
}

RFrrbluponeout_MI=merge(avg_merge,rrbluponeout_cv5.result,by="GenotypeID")

write.table(RFrrbluponeout_MI,"2020NE_RFrrbluponeoutallcv5_MI.csv",quote = F,row.names = F,sep=",")
ggplot(data=RFrrbluponeout_MI,aes(x=Total.Weight.of.kernels..g._MI,Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("Predicted Values")+ylab("Observed MI yield")


