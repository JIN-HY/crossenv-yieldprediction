library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)


# read data
maizeNE=read.csv("data/Bugeater2020_merged_v6_preprocess.csv")
# maizeMI=read.csv("data/2020maizeMI_NEname_preprocess.tsv",sep="\t")
# avgNE=read.csv("data/2020NE_avg.csv")
# avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("data2/2020merge_avg_narm.csv")
rrblupfull=read.csv("data2/rrblupresult/rrblup_full.csv")
rrblup5cv=read.csv("data2/rrblupresult/NE_rrblupcv.csv1")
RFfull <- read.csv("data2/RF/RF_bs.csv")
#RFcv = read.csv("data/RF/RF_5cv.csv1")
RFcvMI = read.csv("data2/RF/RF_5cv_MI.csv1")

avg_merge$Yield_MImean <- 2*avg_merge$Total.Weight.of.kernels..g.
rmse <- sqrt(mean((avg_merge$Yield_MImean-avg_merge$TotalGrainMassGrams)^2))
# mean vs mean
# svg("plot/NEvsMI.svg")
p2a <- ggplot(data=avg_merge,aes(x=TotalGrainMassGrams,y=Yield_MImean))+geom_point()+
  #xlim(45,95)+
  theme_classic()+xlab("Observed NE yield (g)")+ylab("Observed MI yield (g)")+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20))+
  annotate("text",  x = 0, y = 870,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)
  #stat_regline_equation(
   # aes(label =  paste( ..r.label.., sep = "~~~~")))

# dev.off()

# spread reps
#repMI=maizeMI[,c("Name","Total.Weight.of.kernels..g.","Rep")] %>% group_by(Name, Rep) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% spread(key=Rep,value=Total.Weight.of.kernels..g.)
repNE=maizeNE[,c("GenotypeID","TotalGrainMassGrams","Replicate")] %>% group_by(GenotypeID, Replicate) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% spread(key=Replicate,value=TotalGrainMassGrams)
colnames(repNE)=c("GenotypeID","rep1","rep2")
repNE_avgMI=merge(repNE,avg_merge,by = "GenotypeID" )

svg("plot/NErep1vsMI.svg")
ggplot(data=repNE_avgMI,aes(x=rep1,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+xlim(0,1000)+ylim(0,1000)+
  stat_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=10) + theme(text = element_text(size = 30))+
  xlab("NE rep1 (g)")+ylab("MI yield (g)")
dev.off()

svg("plot/NErep2vsMI.svg")
ggplot(data=repNE_avgMI,aes(x=rep2,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+xlim(0,1000)+ylim(0,1000)+
  stat_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("NE rep2 (g)")+ylab("MI yield (g)")
dev.off()

# rrblup NE full
#rrblupfull_avg= merge(rrblupfull,avg_merge,by = "GenotypeID")
rmse = sqrt(mean((rrblupfull$Yield_MI*2-rrblupfull$predicted_yield)^2))
# svg("plot/NErrblupvsMI.svg")
p2b <- ggplot(data=rrblupfull,aes(x=predicted_yield,y=Yield_MI*2))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Genomic predicted yield (g)")+ylab("Observed MI yield (g)")+
  annotate("text",  x = 0, y = 870,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)
#+
  #ggtitle("Predicting MI with NE rrblup full model")
# dev.off()

# rrblupcv vs MI
svg("plot/NErrblupcvvsMI.svg")
rmse=sqrt(mean((rrblup5cv$Yield_MI*2-rrblup5cv$Yield_NEpredicted)^2))
p3a=ggplot(data=rrblup5cv,aes(x=Yield_NEpredicted,y=Yield_MI*2))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,750)+ylim(0,750)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Genomic predicted yield (g)")+ylab("Observed MI yield (g)")+
  annotate("text",  x = 0, y = 640,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)
#+
  #ggtitle("Predicting MI untested lines with NE rrblup cv")
dev.off()

# rrblupcv vs NE
svg("plot/NErrblupcvvsNE.svg")
rmse=sqrt(mean((rrblup5cv$Yield_NE-rrblup5cv$Yield_NEpredicted)^2))
ggplot(data=rrblup5cv,aes(x=Yield_NEpredicted,y=Yield_NE))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Genomic predicted yield (g)")+ylab("Observed NE yield (g)")+
  annotate("text",  x = 0, y = 870,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)

#+
  #ggtitle("Predicting NE untested lines with NE rrblup cv")
dev.off()

# RFcv vs MI
#nsvg("plot/NERFcvvsMI.svg")
rmse=sqrt(mean((RFcvMI$Yield_MIpredicted-RFcvMI$Yield_MImean)^2))
p3b=ggplot(data=RFcvMI,aes(x=Yield_MIpredicted,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,750)+ylim(0,750)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Phenotypic predicted yield (g)")+ylab("Observed MI yield (g)")+
  annotate("text",  x = 0, y = 640,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)
#+
  #ggtitle("Predicting MI untested lines with NE RF cv")
# dev.off()

# RFcv vs NE ### hasn;t been done!!!!!!!!!!!!!
svg("plot/NERFcvvsNE.svg")
rmse=sqrt(mean((RFcv)))
ggplot(data=RFcv,aes(x=Yield_NEpredicted,y=Yield_NEmean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=10) + theme(text = element_text(size = 30)) +
  xlab("RF predicted yield")+ylab("Observed NE yield")#+
  #ggtitle("Predicting NE untested lines with NE RF cv")
dev.off()

# RFfull vs MI
svg("plot/NERFfullvsMI.svg")
ggplot(data=RFfull,aes(x=Yield_NEpredicted,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Phenomic predicted yield")+ylab("Observed MI yield")#+
  #ggtitle("Predicting MI tested lines with NE RF full")
dev.off()



##### RFbs vs MI
rmse=sqrt(mean((RFfull$Yield_MImean-RFfull$Yield_NEpredicted)^2))
p2c=ggplot(data=RFfull,aes(x=Yield_NEpredicted,y=Yield_MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_classic()+
  stat_smooth(method = "lm")+xlim(0,1000)+ylim(0,1000)+
  stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=5) + theme(text = element_text(size = 20)) +
  xlab("Phenotypic predicted yield (g)")+ylab("Observed MI yield (g)")+
  annotate("text",  x = 0, y = 870,hjust=0,
           label = paste("RMSE =", round(rmse, 1)), size=5)

p2d <- ggplot()
p3c <- ggplot()
panel2 <- grid.arrange(p2a, p2b, p2c, p2d, ncol = 2)
panel3 <- grid.arrange(p3a, p3b, p3c, ncol=3)
