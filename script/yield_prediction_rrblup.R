
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(rrBLUP)

avgNE=read.csv("2020NE_avg.csv")
avgMI=read.csv("2020MI_avg.csv")
avg_merge=read.csv("2020merge_avg_narm.csv")


# rrblup
snp_markers  = read.csv("rrblup/WDPinter.num.txt")

#Remove Header and PIs
snp_markers2 = snp_markers[1:nrow(snp_markers),2:ncol(snp_markers)]

full_train_result=avg_merge

for (i in 2:ncol(avg_merge)) {
  pheno2 = as.matrix(avgyld[,i]) # has to be matrix, even though one column, otherwise throws errors when splitting data
  
  pheno_training_data = as.matrix(pheno2)
  snp_training_data = snp_markers2 # Converting to matrix here causes a vector as result which creates problems in model
  
  
  snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
  snp_training_data = as.matrix(snp_training_data)
  
  trained_model = mixed.solve(y=pheno_training_data,Z = snp_training_data)
  
  marker_effects = as.matrix(trained_model$u) # marker effects center around zero
  BLUE = as.vector(trained_model$beta) # BLUE is baseline of which marker effects center
  
  predicted_train = as.matrix(snp_training_data) %*% marker_effects # rrBLUPS
  full_train_result[i] = as.vector((predicted_train[,1])+BLUE)
}
write.table(full_train_result,"2020NEMI_rrblup_full.csv",quote = F,row.names = F,sep = ",")


##############################################################################
# cross validation

total_pheno = length(pheno2[,1])

set.seed(10) 
random.order = sample(1:total_pheno,total_pheno)
k = 5
validate = matrix(random.order, ncol = k)
predicted_cv_result=avg_merge

for (trait in 2:ncol(avg_merge)) {
  for (i in 1:k) {
    pheno2 = as.matrix(avgyld[,trait])
    
    validation_entries = validate[,i] # select the the validation entries
    training_entries = as.vector(validate[,-i])
    
    pheno_training_data = pheno2[training_entries,]
    snp_training_data = snp_markers2[training_entries,] # Converting to matrix here causes a vector as result which creates problems in model
    
    pheno_validate_data = pheno2[validation_entries,]
    snp_validate_data = snp_markers2[validation_entries,]
    
    snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
    snp_training_data = as.matrix(snp_training_data)
    
    snp_validate_data = sapply(snp_validate_data, as.numeric) # Matrix must be numeric
    snp_validate_data = as.matrix(snp_validate_data,K=NULL)
    
    trained_model = mixed.solve(y=pheno_training_data,Z = snp_training_data)
    summary(trained_model)
    
    marker_effects = as.matrix(trained_model$u) # marker effects center around zero
    BLUE = as.vector(trained_model$beta) # BLUE is baseline of which marker effects center
    
    predicted_train = as.matrix(snp_training_data) %*% marker_effects # rrBLUPS
    predicted_validate = as.matrix(snp_validate_data) %*% marker_effects   # rrBLUPS
    
    predicted_train_result = as.vector((predicted_train[,1])+BLUE)
    predicted_cv_result[validation_entries,trait] = as.vector((predicted_validate[,1])+BLUE)
  }
}
write.table(predicted_cv_result,"2020NEMI_rrblupcv.csv",quote=F,col.names = T,row.names = F)

####################################################################################################

oneout_result=avg_merge

for (trait in 2:ncol(avg_merge)) {
  for (i in 1:nrow(avg_merge)) {
    pheno2 = as.matrix(avg_merge[,trait])
    
    pheno_training_data = as.matrix(pheno2[-i,],ncol=1)
    snp_training_data = snp_markers2[-i,] # Converting to matrix here causes a vector as result which creates problems in model
    
    pheno_validate_data = pheno2[i,]
    snp_validate_data = snp_markers2[-i,]
    
    snp_training_data = sapply(snp_training_data, as.numeric) # Matrix must be numeric
    snp_training_data = as.matrix(snp_training_data)
    
    snp_validate_data = sapply(snp_validate_data, as.numeric) # Matrix must be numeric
    snp_validate_data = as.matrix(snp_validate_data,K=NULL)
    
    trained_model = mixed.solve(y=pheno_training_data,Z = snp_training_data)
    summary(trained_model)
    
    marker_effects = as.matrix(trained_model$u) # marker effects center around zero
    BLUE = as.vector(trained_model$beta) # BLUE is baseline of which marker effects center
    
    predicted_train = as.matrix(snp_training_data) %*% marker_effects # rrBLUPS
    predicted_validate = as.matrix(snp_validate_data) %*% marker_effects   # rrBLUPS
    
    predicted_train_result = as.vector((predicted_train[,1])+BLUE)
    oneout_result[-i,trait] = as.vector((predicted_validate[,1])+BLUE)
  }
}


write.table(oneout_result,"2020NIMI_rrblup_oneout.csv",quote = F,row.names = F,sep = ",")