library(randomForest)
library(tidyverse)
library(grid)
library(ggpubr)



avgNE=read.csv("data/2020NE_avg.csv")
avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("data2/2020merge_avg_narm.csv")


# x.idx = apply(matrix(colnames(avg_merge)), 1, grepl,pattern="NE")
# x.idx = c(TRUE,x.idx[-1])
# xy=avg_merge[,x.idx]
xyids <- c("GenotypeID","TotalGrainMassGrams","Total.Weight.of.kernels..g.",
        "DaysToPollen", "DaysToSilk", "RootLodgingPct","StalkLodgingPct","LeafLengthCM",
        "LeafWidthCM", "PlantHeightCM", "ExtantLeafNumber1", "BranchesPerTassel", 
        "TasselLengthCM", "TasselSpikeLengthCM", "NodesWithBraceRoots", "BranchZoneLengthCM")
xy <- avg_merge[colnames(avg_merge) %in% xyids]
# xycolnames <- c("GenotypeID","Days to pollen", "Days to silk", "Root lodging rate",
#                 "Stalk lodging rate", "Leaf length", "Leaf width", "Plant height", 
#                 "Leaf number", "Number of nodes with brace roots", "Branches per tassel",
#                 "Tassel length", "Tassel branch zone length", "Tassel spike length",
#                 "TotalGrainMassGrams_NEmean","Yield_MImean")
# colnames(xy) <- xycolnames
colnames(xy)[15:16] <- c("TotalGrainMassGrams_NEmean","Yield_MImean")
xy$Yield_MImean <- xy$Yield_MImean*2
write.table(xy,"data2/2020NEMI_RFfulltraining.csv",sep = ",",quote = F,row.names = F)

  
#### correlation among x
library(corrplot)
x = xy[,c(2:14)]
xcorr = cor(x, use="complete.obs")
corrplot(xcorr, method=c( "square"), type = "upper", order="FPC", diag = F, 
         tl.col = "black", col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100),
         tl.srt = 45, tl.cex=0.8)



RFresult=data.frame(GenotypeID=xy$GenotypeID,Yield_MIpredicted=xy$Yield_MImean)
# find nas
colnames(xy)[colSums(is.na(xy)) > 0]
xy[,colnames(xy)[colSums(is.na(xy)) > nrow(xy)*0.1]]=NULL

for(i in 1:ncol(xy)){
  xy[is.na(xy[,i]), i] <- mean(xy[,i], na.rm = TRUE)
}
#############################################################################################

gnumber <- nrow(xy)
random.order = sample(1:gnumber, gnumber)
k = 5
validate = matrix(random.order, ncol = k)
RFcv5.result=RFresult


for (i in 1:k) {
  validation.ids = validate[,i] # select the the validation entries
  training.ids = as.vector(validate[,-i])

  xy.training = xy[training.ids,-1]
  #    x.training = xy[training.ids,ncol(xy)]

  xy.validate = xy[validation.ids,-1]
  #    x.validate = xy[validate.ids,ncol(xy)]

  yield.forest <- randomForest(Yield_MImean ~ ., data = xy.training, importance = TRUE)#, na.action = na.roughfix)
  yield.forest

  RFcv5.result[validation.ids,"Yield_MIpredicted"]=predict(yield.forest,xy.validate)
}

RFcv5_MI=merge(xy,RFcv5.result,by="GenotypeID")
write.table(RFcv5_MI,"data2/2020NE_RFcv5_MI.csv",quote = F,row.names = F,sep=",")

# for (t in 1:30) { # RFCV
#   set.seed(t)
#   gnumber <- nrow(xy)
#   random.order = sample(1:gnumber, gnumber)
#   k = 5
#   validate = matrix(random.order, ncol = k)
# xyNE <- xy
# RFcv_NE.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_NEpredicted=0)
# xyNE$Yield_MImean <- NULL
# xyNE$TotalGrainMassGrams_NEmean <- NULL
# xyNE$Yield_NEmean <- avg_merge$TotalGrainMassGrams_NEmean
# for (i in 1:k) {
#   validation.ids = validate[,i] # select the the validation entries
#   training.ids = as.vector(validate[,-i])
#   
#   xyNE.training = xyNE[training.ids,-1]
#   #    x.training = xy[training.ids,ncol(xy)]
#   
#   xyNE.validate = xyNE[validation.ids,-1]
#   #    x.validate = xy[validate.ids,ncol(xy)]
#   
#   yield.forest <- randomForest(Yield_NEmean ~ ., data = xyNE.training, importance = TRUE)
#   yield.forest
#   
#   RFcv_NE.result[validation.ids,"Yield_NEpredicted"]=predict(yield.forest,xyNE.validate)
# }
# write.table(RFcv_NE.result,paste0("data2/RF/RF_5cv.csv",t),row.names = F,quote = F,sep=",")
# }


ipscores_cv <- data.frame()
for (t in 1:30) { # RFCV
  set.seed(t)
  gnumber <- nrow(xy)
  random.order = sample(1:gnumber, gnumber)
  k = 5
  validate = matrix(random.order, ncol = k)
  xyNE <- xy
  RFcv_MI.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_MIpredicted=0)
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  for (i in 1:k) {
    validation.ids = validate[,i] # select the the validation entries
    training.ids = as.vector(validate[,-i])
    
    xyNE.training = xyNE[training.ids,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    xyNE.validate = xyNE[validation.ids,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_MImean ~ ., data = xyNE.training, importance = TRUE)
    yield.forest
    
    ipscore <- importance(yield.forest, type = 1)
    ipscore <- data.frame(run=t, trait=rownames(ipscore), score = ipscore[,1])
    ipscores_cv <- rbind(ipscores_cv, ipscore)
    
    RFcv_MI.result[validation.ids,"Yield_MIpredicted"]=predict(yield.forest,xyNE.validate)
  }
  write.table(RFcv_MI.result,paste0("data2/RF/RF_5cv_MI.csv",t),row.names = F,quote = F,sep=",")
}

ipscores_cv_summary <- ipscores_cv %>% group_by(trait) %>% summarise(importance_score=mean(score), scoresd=sd(score, na.rm = T))
ipscores_cv_summary=ipscores_cv_summary[order(ipscores_cv_summary$importance_score,decreasing=T),]
ipscores_cv_summary$trait <- c("Branches per tassel", "Days to silk","Days to pollen","Leaf number", "Plant height","Stalk loding rate","Leaf width","Root lodging rate",  
                                "Tassel spike length","Tassel branch zone length", "Leaf length", "Tassel length","Number of nodes with brace roots"
                               )


#ipscores_cv_summary$t=factor(ipscores_cv_summary$trait, levels=paste0("k", 1:14))
ipscores_cv_summary$trait=factor(ipscores_cv_summary$trait, levels=ipscores_cv_summary$trait[order(ipscores_cv_summary$importance_score,decreasing=T)])
ggplot(ipscores_cv_summary, aes(x=trait, y=importance_score)) + geom_bar(stat = "identity", color="black") + 
  geom_errorbar(aes(ymin=importance_score-scoresd, ymax=importance_score+scoresd), width=0.2,) +ylab("Importance score")+
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=13, color="black"), axis.text.y = element_text(size=13),
                     axis.title.y = element_text(size = 15))

# reduce variable
xyNE <- xy
RFcv_MI.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_MIpredicted=0)
for (r in 1:13) {
  red=ipscores_cv_summary$trait[1:r]
  xy <- xy[!colnames(xy) %in% red]

ipscores_cv_reduce=list()
ipscores_cv_reduce[[r]]=data.frame()

for (t in 1:20) { # RFCV
  set.seed(t)
  gnumber <- nrow(xy)
  random.order = sample(1:gnumber, gnumber)
  k = 5
  validate = matrix(random.order, ncol = k)
  xyNE <- xy
  
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  for (i in 1:k) {
    validation.ids = validate[,i] # select the the validation entries
    training.ids = as.vector(validate[,-i])
    
    xyNE.training = xyNE[training.ids,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    xyNE.validate = xyNE[validation.ids,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_MImean ~ ., data = xyNE.training, importance = TRUE)
    yield.forest
    
    ipscore <- importance(yield.forest, type = 1)
    ipscore <- data.frame(run=t, trait=rownames(ipscore), score = ipscore[,1])
    ipscores_cv_reduce[[r]] <- rbind(ipscores_cv_reduce[[r]], ipscore)
    
    RFcv_MI.result[validation.ids,paste0("Yield_MIpredicted_reduce",r)]=predict(yield.forest,xyNE.validate)
  }
  
}
write.table(RFcv_MI.result,paste0("data/RF/RF_5cv_MI_reduce.csv",t),row.names = F,quote = F,sep=",")
}

for (col  in 5:17) {
  acc=cor(RFcv_MI.result$Yield_MImean,RFcv_MI.result[,col])
  print(acc)
}


ez=c("GenotypeID","TotalGrainMassGrams_NEmean","Yield_MImean","LeafLengthCM_NEmean", "RootLodgingPct_NEmean", "LeafWidth_NEmean", "PlantHeightCM_NEmean", "LeafNumber_NEmean", "DaysToPollen_NEmean", 
     "DaysToSilk_NEmean")
xy <- xy[colnames(xy) %in% ez]

for (t in 1:20) { # RFCV
  set.seed(t)
  gnumber <- nrow(xy)
  random.order = sample(1:gnumber, gnumber)
  k = 5
  validate = matrix(random.order, ncol = k)
  xyNE <- xy
  RFcv_MI.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_MIpredicted=0)
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  for (i in 1:k) {
    validation.ids = validate[,i] # select the the validation entries
    training.ids = as.vector(validate[,-i])
    
    xyNE.training = xyNE[training.ids,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    xyNE.validate = xyNE[validation.ids,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_MImean ~ ., data = xyNE.training, importance = TRUE)
    yield.forest
    
    ipscore <- importance(yield.forest, type = 1)
    ipscore <- data.frame(run=t, trait=rownames(ipscore), score = ipscore[,1])
    ipscores_cv <- rbind(ipscores_cv, ipscore)
    
    RFcv_MI.result[validation.ids,"Yield_MIpredicted"]=predict(yield.forest,xyNE.validate)
  }
  write.table(RFcv_MI.result,paste0("data/RF/RF_5cv_MI_ez.csv",t),row.names = F,quote = F,sep=",")
}
acc=cor(RFcv_MI.result$Yield_MImean,RFcv_MI.result$Yield_MIpredicted)
print(acc)


### end of reduced models

ipscores_bs <- data.frame()
for (t in 1:30) { # RFBS
  set.seed(t)
  gnumber <- nrow(xy)
  sp <- sample(1:gnumber, as.integer(0.8*gnumber))
  xyNE <- xy
  RFbs_NE.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_NEpredicted=0)
  RFbs_NE.result <- RFbs_NE.result[sp,]
  xyNE$Yield_MImean <- NULL
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  xyNE$Yield_NEmean <- xy$TotalGrainMassGrams_NEmean
  xyNE <- xyNE[sp,]
    xyNE.training = xyNE[,-1]
    #    x.training = xy[training.ids,ncol(xy)]
    
    xyNE.validate = xyNE[,-1]
    #    x.validate = xy[validate.ids,ncol(xy)]
    
    yield.forest <- randomForest(Yield_NEmean ~ ., data = xyNE.training, importance = TRUE)
    yield.forest
    
    ipscore <- importance(yield.forest, type = 1)
    ipscore <- data.frame(run=t, trait=rownames(ipscore), score = ipscore[,1])
    ipscores_bs <- rbind(ipscores_bs, ipscore)
    
    RFbs_NE.result[,"Yield_NEpredicted"]=predict(yield.forest,xyNE.validate)
  write.table(RFbs_NE.result,paste0("data2/RF/RF_bs.csv",t),row.names = F,quote = F,sep=",")
}


ipscores_bs_summary <- ipscores_bs %>% group_by(trait) %>% summarise(importance_score=mean(score), scoresd=sd(score, na.rm = T))
ipscores_bs_summary$trait <- c("Branches per tassel", "Tassel branch zone length","Days to pollen", "Days to silk", "Leaf length", "Leaf number", "Leaf width",
                               "Number of nodes with brace roots", "Plant height", "Root lodging rate", "Stalk loding rate","Tassel length", "Tassel spike length"
                               )

#ipscores_bs_summary=ipscores_bs_summary[order(ipscores_bs_summary$importance_score,decreasing=T),]
ipscores_bs_summary$trait=factor(ipscores_bs_summary$trait, levels=ipscores_bs_summary$trait[order(ipscores_bs_summary$importance_score,decreasing=T)])
ggplot(ipscores_bs_summary, aes(x=trait, y=importance_score)) + geom_bar(stat = "identity", color="black") + 
  geom_errorbar(aes(ymin=importance_score-scoresd, ymax=importance_score+scoresd), width=0.2,) + ylab("Importance score")+
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=13, color="black"), axis.text.y = element_text(size=13),
                     axis.title.y = element_text(size = 15))


ipscores_bs_summary$prediction <- "NE"
ipscores_cv_summary$prediction <- "MI"
ipscores_summary <- rbind(ipscores_bs_summary, ipscores_cv_summary)




pdf("plot/ipscoresummary.pdf", width = 5, height = 3)
ggplot(ipscores_summary, aes(x=trait, y=importance_score, fill=prediction)) + geom_col(width = 0.6, position=position_dodge(width=0.6)) + 
  geom_errorbar(aes(ymin=importance_score-scoresd, ymax=importance_score+scoresd), width=0.2, position=position_dodge(width=0.6)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = c("NE"="blue","MI"="red"))
dev.off()

pdf("plot/FTvsYLD.svg", width = 3, height = 3)
ggplot(xy) + geom_point(aes(x=DaysToSilk_NEmean, y=TotalGrainMassGrams_NEmean, color="NE", alpha=0.5)) + geom_point(aes(x=DaysToSilk_NEmean, y=Yield_MImean, color="MI", alpha=0.5)) +
  ylab("Average Yield (g)") + stat_smooth(span=1.2,aes(x=DaysToSilk_NEmean, y=TotalGrainMassGrams_NEmean, color="NE")) + stat_smooth(span=1.2,aes(x=DaysToSilk_NEmean, y=Yield_MImean, color="MI")) +
  theme_bw() +scale_color_manual(values = c("NE"="blue","MI"="red")) + theme(legend.position="none")
dev.off()

pdf("plot/NBTvsYLD.pdf", width = 3, height = 3)
ggplot(xy) + geom_point(aes(x=BranchesPerTassel_NEmean, y=TotalGrainMassGrams_NEmean, color="NE",alpha=0.5)) + geom_point(aes(x=BranchesPerTassel_NEmean, y=Yield_MImean, color="MI", alpha=0.5)) +
  ylab("Average Yield (g)") + stat_smooth(span=1.2,aes(x=BranchesPerTassel_NEmean, y=TotalGrainMassGrams_NEmean, color="NE")) + stat_smooth(span=1.2,aes(x=BranchesPerTassel_NEmean, y=Yield_MImean, color="MI")) +
  theme_bw() +scale_color_manual(values = c("NE"="blue","MI"="red"))+ theme(legend.position="none")
dev.off()

pdf("plot/PHvsYLD.pdf", width = 3, height = 3)
ggplot(xy) + geom_point(aes(x=PlantHeightCM_NEmean, y=TotalGrainMassGrams_NEmean, color="NE", alpha=0.5)) + geom_point(aes(x=PlantHeightCM_NEmean, y=Yield_MImean, color="MI", alpha=0.5)) +
  ylab("Average Yield (g)") + stat_smooth(span=1.2,aes(x=PlantHeightCM_NEmean, y=TotalGrainMassGrams_NEmean, color="NE")) + stat_smooth(span=1.2,aes(x=PlantHeightCM_NEmean, y=Yield_MImean, color="MI")) +
  theme_bw()+scale_color_manual(values = c("NE"="blue","MI"="red"))+ theme(legend.position="none")
dev.off()

ggplot(xy) + geom_point(aes(x=TillersPerPlant_NEmean, y=TotalGrainMassGrams_NEmean, color="NE")) + geom_point(aes(x=TillersPerPlant_NEmean, y=Yield_MImean, color="MI")) +
  ylab("Average Yield (g)")+theme_bw()+scale_color_manual(values = c("NE"="blue","MI"="red"))+ theme(legend.position="none")

ggplot(data=RFcv_NE.result,aes(x=Yield_NEpredicted, Yield_NEmean))+
  geom_point()+
  theme_bw()+
  stat_smooth(span=1.2,method = "lm")+
  stat_cor(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("Predicted Values")+ylab("Observed NE yield")

pdf("plot/FTcorrelation.pdf", width = 3, height = 3)
ggplot(data=avg_merge,aes(x=DaysToPollen_NEmean, y=DaysToPollen_MImean))+
  geom_point()+
  theme_bw()+
  stat_smooth(span=1.2,method = "lm")+
  stat_cor(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("Anthesis NE (d)")+ylab("Anthesis MI (d)")
dev.off()

################# rrblup with phenotypic data
library(rrBLUP)


for (t in 1:20) { # RFBT
  set.seed(t)
  gnumber <- nrow(xy)
  sp <- sample(1:gnumber, as.integer(0.8*gnumber))
  xyNE <- xy
  PPRRbs_NE.result=data.frame(GenotypeID=xyNE$GenotypeID,Yield_NEmean=xyNE$TotalGrainMassGrams_NEmean,Yield_MImean=xyNE$Yield_MImean,Yield_NEpredicted=0)
  PPRRbs_NE.result <- PPRRbs_NE.result[sp,]
  xyNE$Yield_MImean <- NULL
  xyNE$TotalGrainMassGrams_NEmean <- NULL
  xyNE$Yield_NEmean <- avg_merge$TotalGrainMassGrams_NEmean
  xyNE <- xyNE[sp,]
  xNE.training = xyNE[,-c(1,ncol(xyNE))]
  #    x.training = xy[training.ids,ncol(xy)]
  yNE.training <- xyNE[,ncol(xyNE)]
  xNE.validate = xyNE[,-c(1,ncol(xyNE))]
  #    x.validate = xy[validate.ids,ncol(xy)]
  
  yield.rr <- mixed.solve(y=yNE.training,Z = xNE.training)
  marker_effects = as.matrix(yield.rr$u) # marker effects center around zero
  BLUE = as.vector(yield.rr$beta) # BLUE is baseline of which marker effects center
  
  #predicted_train = as.matrix(xNE.validate) %*% marker_effects # rrBLUPS
  
  PPRRbs_NE.result[,"Yield_NEpredicted"]=as.matrix(xNE.validate) %*% marker_effects + BLUE
  write.table(PPRRbs_NE.result,paste0("data/RF/PPRR_bs.csv",t),row.names = F,quote = F,sep=",")
}

pdf("plot/PP_rrblup.pdf", height = 3, width = 3)
ggplot(data=PPRRbs_NE.result,aes(x=Yield_NEpredicted, Yield_MImean))+
  geom_point()+
  theme_bw()+
  stat_smooth(span=1.2,method = "lm")+
  stat_cor(
    aes(label =  paste( ..rr.label.., sep = "~~~~")))+
  xlab("Predicted NE yield (g)")+ylab("Observed NE yield (g)")
dev.off()


########################################################################################################################################
# 
library(pls)
X=as.matrix(xy0[,c(2:14)])
test_plsr <- function(n){
  plsr_model <- plsr(xy0$TotalGrainMassGrams_NEmean ~ X, ncomp = n)
  test=predict(plsr_model, X, ncomp=n)
  res=cor(test, xy0$Yield_MImean, use = "complete.obs")
  re=RMSE(xy0$Yield_MImean,test)
  print(res)
  print(re)
  return(res)
}

for (n in 2:13) {
  test_plsr(n)
}

plsr_model <- plsr(xy0$TotalGrainMassGrams_NEmean ~ X, ncomp = 6)
test=predict(plsr_model, X, ncomp=6)
xy0$test=test
ggplot(data=xy0,aes(x=test,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Phenotypic predicted yield (g)")+ylab("Observed MI yield (g)")+
  annotate("text",  x = 0, y = 870,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)


