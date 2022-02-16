library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(svglite)

# read data
maizeNE=read.csv("Bugeater2020_merged_v6_preprocess.csv")
maizeMI=read.csv("2020maizeMI_NEname_preprocess.tsv",sep="\t")
avgNE=read.csv("2020NE_avg.csv")
avgMI=read.csv("2020MI_avg.csv")
avg_merge=read.csv("2020merge_avg_narm.csv")
rrblupfull=read.csv("rrblup/2020NEMI_rrblup_full.csv")
rrblup5cv=read.csv("rrblup/2020NEMI_rrblup_5cv.csv")
rrbluponeout=read.csv("rrblup/2020NEMI_rrblup_full.csv")

# mean vs mean
p1=ggplot(data=avg_merge,aes(x=TotalGrainMassGrams_NEmean,y=Total.Weight.of.kernels..g._MImean))+geom_point()+
  #xlim(45,95)+
  theme_bw()+xlab("Average yield NE (g)")+ylab("Average yield MI (g)")+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~")))
ggsave("NEmean_MI.svg",p1)

# spread reps
repMI=maizeMI[,c("Name","Total.Weight.of.kernels..g.","Rep")] %>% group_by(Name, Rep) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% spread(key=Rep,value=Total.Weight.of.kernels..g.)
repNE=maizeNE[,c("GenotypeID","TotalGrainMassGrams","Replicate")] %>% group_by(GenotypeID, Replicate) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% spread(key=Replicate,value=TotalGrainMassGrams)
colnames(repNE)=c("GenotypeID","rep1","rep2")
repNE_avgMI=merge(repNE,avgMI,by = "GenotypeID" )

p2=ggplot(data=repNE_avgMI,aes(x=rep1,y=Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("NE rep1")+ylab("MI yield")
ggsave("NErep1_MI.svg",p2)

p3=ggplot(data=repNE_avgMI,aes(x=rep2,y=Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("NE rep2")+ylab("MI yield")
ggsave("NErep2_MI.svg",p3)

# rrblup NE full
rrblupfull_avg= merge(rrblupfull,avg_merge,by = "GenotypeID")
p4=ggplot(data=rrblupfull_avg,aes(x=TotalGrainMassGrams_NE,y=Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("NE rrblup")+ylab("MI yield")+
  ggtitle("Predicting MI with NE rrblup full model")
ggsave("NErrblupfull_MI.svg",p4)

# rrblup NE oneout
rrbluponeout_avg= merge(rrblupfull,avg_merge,by = "GenotypeID")
p5=ggplot(data=rrbluponeout_avg,aes(x=TotalGrainMassGrams_NE,y=Total.Weight.of.kernels..g._MImean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("NE rrblup")+ylab("MI yield")+
  ggtitle("Predicting MI with NE rrblup one out")
ggsave("NErrbluponeout_MI.svg",p5)
p6=ggplot(data=rrbluponeout_avg,aes(x=TotalGrainMassGrams_NE,y=TotalGrainMassGrams_NEmean))+
  geom_point()+
  #facet_wrap(ncol=3,vars(var),scales = "free") +
  theme_bw()+
  stat_smooth(method = "lm")+
  stat_regline_equation(
    aes(label =  paste( ..rr.label.., sep = "~~~~"))) +
  xlab("NE rrblup")+ylab("MI yield")+
  ggtitle("Predicting 1 unobserved genotype in NE")
ggsave("NErrbluponeout_NE.svg",p6)

