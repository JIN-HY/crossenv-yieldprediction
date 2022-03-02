library(randomForest)
library(tidyverse)
library(rrBLUP)
library(grid)
library(ggpubr)
library(MASS)
library(candisc)
library(data.table)



t <- commandArgs(trailingOnly = T)
t <- as.numeric(t)
set.seed(t) 


setDTthreads(threads = 0, restore_after_fork = NULL, throttle = NULL)

snp_markers  = fread("WiDiv_maf005_725geno.num.txt")
snp_markers = snp_markers[order(snp_markers[,1]),]
snp_markers2 = snp_markers[1:nrow(snp_markers),2:ncol(snp_markers)]


avgNE=read.csv("data/2020NE_avg.csv")
avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("data/2020merge_avg_narm.csv")



#############################################################################################


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



# for (t in 1:20) {
# set.seed(t)
# 
# gnumber <- nrow(avg_cross)
# random.order = sample(1:gnumber, gnumber)
# k = 5
# validate = matrix(random.order, ncol = k)
# 
# 
# #RFcrossbs.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=NA, Yield_MImean = avgMI_cross$Yield_MImean)
# RFcrosscv.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=avgMI_cross$Yield_MImean, Yield_MImean = avgMI_cross$Yield_MImean)
# RFcv.result = data.frame(GenotypeID=avgNE_cross$GenotypeID,Yield_NEpredicted=avgNE_cross$Yield_NEmean, Yield_NEmean = avgNE_cross$Yield_NEmean)
# 
# for (i in 1:k) {
#   validate.ids = validate[,i] # select the the validation entries
#   training.ids = as.vector(validate[,-i])
#   
#   avgNE_cross.training = avgNE_cross[training.ids,-1]
#   avgMI_cross.training = avgMI_cross[training.ids,-1]
#   #    x.training = xy[training.ids,ncol(xy)]
#   
#   avgNE_cross.validate = avgNE_cross[validate.ids,-1]
#   avgMI_cross.validate = avgMI_cross[validate.ids,-1]
#   #    x.validate = xy[validate.ids,ncol(xy)]
#   
#   yield.forest <- randomForest(Yield_NEmean ~ ., data = avgNE_cross.training, importance = TRUE)
#   yield.forest
#   
#   RFcrossbs.result[training.ids, "Yield_MIpredicted"]=predict(yield.forest,avgMI_cross.training)
#   RFcrosscv.result[validate.ids,"Yield_MIpredicted"]=predict(yield.forest,avgMI_cross.validate)
#   RFcv.result[validate.ids,"Yield_NEpredicted"]=predict(yield.forest,avgNE_cross.validate)
#   
#   #write.table(RFcrossbs.result, paste0("data/RF/RF_crossenv_bs.csv", k), quote=F, row.names = F, sep=",")
# 
#   
#   #corr_RFcrossbs = cor(RFcrossbs.result$Yield_MIpredicted,RFcrossbs.result$Yield_MImean); print(corr_RFcrossbs)
#   corr_RFcrosscv = cor(RFcrosscv.result$Yield_MIpredicted,RFcrosscv.result$Yield_MImean); print(corr_RFcrosscv)
#   corr_RFcv = cor(RFcv.result$Yield_NEpredicted, RFcv.result$Yield_NEmean, use = "complete.obs"); print(corr_RFcv)
# }
# 
# write.table(RFcrosscv.result, paste0("data/RF/RF_crossenv_cv.csv", t), quote=F, row.names = F, sep=",")
# write.table(RFcv.result, paste0("data/RF/RF_NE_cv.csv", t), quote=F, row.names = F, sep=",")
# 
# } # for t


# set.seed(t)
# 
# gnumber <- nrow(avg_cross)
# random.order = sample(1:gnumber, gnumber)
# k = 5
# validate = matrix(random.order, ncol = k)
# 
# CANcrosscv.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=avgMI_cross$Yield_MImean, Yield_MImean = avgMI_cross$Yield_MImean)
# CANcv.result = data.frame(GenotypeID=avgNE_cross$GenotypeID,Yield_NEpredicted=avgNE_cross$Yield_NEmean, Yield_NEmean = avgNE_cross$Yield_NEmean)
# 
# for (i in 1:k) {
#   
#   validate.ids = validate[,i] # select the the validation entries
#   training.ids = as.vector(validate[,-i])
#   
#   avgNE_cross.training = avgNE_cross[training.ids,-1]
#   avgMI_cross.training = avgMI_cross[training.ids,-1]
#   #    x.training = xy[training.ids,ncol(xy)]
#   
#   avgNE_cross.validate = avgNE_cross[validate.ids,-1]
#   avgMI_cross.validate = avgMI_cross[validate.ids,-1]
#   
#   # snp_training_data = snp_markers2[training_entries,] # Converting to matrix here causes a vector as result which creates problems in model
#   # snp_validate_data = snp_markers2[validation_entries,]
#   # 
#   # snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
#   # snp_training_data = as.matrix(snp_training_data)
#   # 
#   # snp_validate_data = sapply(snp_validate_data, as.numeric) # Matrix must be numeric
#   # snp_validate_data = as.matrix(snp_validate_data,K=NULL)
#   
#   trained_cancor = cancor(as.matrix(avgNE_cross.training[,-1]),matrix(avgNE_cross.training$Yield_NEmean, ncol=1))
#   can_l = trained_cancor$xcoef[,1]
#   can_u_training = as.matrix(avgNE_cross.training[,-1])%*%matrix(can_l)
#   X_training_data = cbind(can_u_training,1)
#   can_u_validate = as.matrix(avgNE_cross.validate[,-1])%*%matrix(can_l)
#   X_validate_data = cbind(can_u_validate,1)
#   can_u_crossenv_validate = as.matrix(avgMI_cross.validate[,-1])%*%matrix(can_l)
#   X_crossenv_validate = cbind(can_u_crossenv_validate,1)
#   
#   trained_model = mixed.solve(y=matrix(avgNE_cross.training$Yield_NEmean), Z = diag(length(avgNE_cross.training$Yield_NEmean)), X=X_training_data)
# #  trained_model = mixed.solve(y=matrix(avgNE_cross.training$Yield_NEmean),Z = snp_training_data, X=X_training_data)
#   summary(trained_model)
#   
#   marker_effects = as.matrix(trained_model$u) # marker effects center around zero
#   BLUE = as.matrix(trained_model$beta) # BLUE is baseline of which marker effects center
#   
#   predicted_train =  X_training_data%*%BLUE
#   predicted_validate =  X_validate_data%*%BLUE
#   predicted_crossenv_validate = X_crossenv_validate%*%BLUE
#   
#   CANcrosscv.result[validate.ids,"Yield_MIpredicted"]=predicted_crossenv_validate
#   CANcv.result[validate.ids,"Yield_NEpredicted"]=predicted_validate
# 
#   corr_CANcrosscv = cor(CANcrosscv.result$Yield_MIpredicted,CANcrosscv.result$Yield_MImean); print(corr_CANcrosscv)
#   corr_CANcv = cor(CANcv.result$Yield_NEpredicted, CANcv.result$Yield_NEmean, use = "complete.obs"); print(corr_CANcv)
# }


######## cda

avgNE_cda <- cbind(avgNE_cross[,-2], env = "NE")
avgMI_cda <- cbind(avgMI_cross[,-2], env = "MI")

avg_cross_cda <- rbind(avgNE_cda, avgMI_cda)
avg_cross_cda$GenotypeID <- factor(avg_cross_cda$GenotypeID)
avg_cross_cda$env <- factor(avg_cross_cda$env)
avg_cross.mod <- lm(cbind(DaysToPollen, LeafLength, LeafWidth, PlantHeight, LeafNumber) ~ GenotypeID + env, data = avg_cross_cda)
avg_cross.manova <- Anova(avg_cross.mod, type = "III", test = "Wilks")
avg_cross.cda.g <- candisc(avg_cross.mod, term = "GenotypeID")
avg_cross.cda.e <- candisc(avg_cross.mod, term = "env")
plot(avg_cross.cda.g)
plot(avg_cross.cda.e)

l_can <- avg_cross.cda.e$coeffs.raw
can1NE <- avg_cross.cda.e$scores[avg_cross.cda.e$scores$env=="NE",]
can1MI <- avg_cross.cda.e$scores[avg_cross.cda.e$scores$env=="MI",]

set.seed(t)

gnumber <- nrow(avg_cross)
random.order = sample(1:gnumber, gnumber)
k = 5
validate = matrix(random.order, ncol = k)

CDAcrosscv.result = data.frame(GenotypeID=avgMI_cross$GenotypeID,Yield_MIpredicted=avgMI_cross$Yield_MImean, Yield_MImean = avgMI_cross$Yield_MImean)
CDAcv.result = data.frame(GenotypeID=avgNE_cross$GenotypeID,Yield_NEpredicted=avgNE_cross$Yield_NEmean, Yield_NEmean = avgNE_cross$Yield_NEmean)

for (i in 1:k) {
  
  validate.ids = validate[,i] # select the the validation entries
  training.ids = as.vector(validate[,-i])
  
  avgNE_cross.training = avgNE_cross[training.ids,-1]
  avgMI_cross.training = avgMI_cross[training.ids,-1]
  #    x.training = xy[training.ids,ncol(xy)]
  
  avgNE_cross.validate = avgNE_cross[validate.ids,-1]
  avgMI_cross.validate = avgMI_cross[validate.ids,-1]
  
  snp_training_data = snp_markers2[training.ids,] # Converting to matrix here causes a vector as result which creates problems in model
  snp_validate_data = snp_markers2[validate.ids,]

  snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
  snp_training_data = as.matrix(snp_training_data)

  snp_validate_data = sapply(snp_validate_data, as.numeric) # Matrix must be numeric
  snp_validate_data = as.matrix(snp_validate_data,K=NULL)
  
  can1NE_training = can1NE[training.ids,"Can1"]
  can1NE_validate = can1NE[validate.ids,"Can1"]
  can1MI_validate = can1MI[validate.ids,"Can1"]
  
  trained_model = mixed.solve(y=matrix(avgNE_cross.training$Yield_NEmean), Z = snp_training_data, X=cbind(can1NE_training,1))
  #  trained_model = mixed.solve(y=matrix(avgNE_cross.training$Yield_NEmean),Z = snp_training_data, X=X_training_data)
  summary(trained_model)
  
  marker_effects = as.matrix(trained_model$u) # marker effects center around zero
  BLUE = as.matrix(trained_model$beta) # BLUE is baseline of which marker effects center
  
  predicted_validate =  cbind(can1NE_validate,1)%*%BLUE + as.matrix(snp_validate_data) %*% marker_effects
  predicted_crossenv_validate = cbind(-can1MI_validate,1)%*%BLUE + as.matrix(snp_validate_data) %*% marker_effects
  
  CDAcrosscv.result[validate.ids,"Yield_MIpredicted"]=predicted_crossenv_validate
  CDAcv.result[validate.ids,"Yield_NEpredicted"]=predicted_validate
  
}
write.table(CDAbt.result, paste0("data/GP/CDArrblup_crossenv_cv.csv", t), quote = F, row.names = F, sep = ",")
write.table(CDAbt.result, paste0("data/GP/CDArrblup_NE_cv.csv", t), quote = F, row.names = F, sep = ",")

corr_CDAcrosscv = cor(CDAcrosscv.result$Yield_MIpredicted,CDAcrosscv.result$Yield_MImean); print(corr_CDAcrosscv)
corr_CDAcv = cor(CDAcv.result$Yield_NEpredicted, CDAcv.result$Yield_NEmean, use = "complete.obs"); print(corr_CDAcv)


sp <- sample(1:gnumber, gnumber*0.8)
snp_bs <- snp_markers2[sp,]
can1NEbs <- can1NE[sp,"Can1"]
can1MIbs <- can1MI[sp,"Can1"]
avgNE_cross_bs <- avgNE_cross[sp,]
avgMI_cross_bs <- avgMI_cross[sp,]
CDAbt.result <- avgMI_cross_bs

trained_model = mixed.solve(y=matrix(avgNE_cross_bs$Yield_NEmean), Z = as.matrix(snp_bs), X=cbind(can1NEbs,1))

marker_effects = as.matrix(trained_model$u) # marker effects center around zero
BLUE = as.matrix(trained_model$beta)

predicted_crossenv_bs = cbind(-can1MIbs,1)%*%BLUE + as.matrix(snp_bs) %*% marker_effects

CDAbt.result$Yield_MIpredicted <- predicted_crossenv_bs
write.table(CDAbt.result, paste0("data/GP/CDArrblup_NE_bt.csv", t), quote = F, row.names = F, sep = ",")
cor_CDAbt <- cor(CDAbt.result$Yield_MImean, CDAbt.result$Yield_MIpredicted)

