library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(svglite)


rmse <- function(observed, predicted){sqrt(sum((observed-predicted)^2/observed)/length(observed))}
rmdse <- function(observed, predicted){sqrt(median((observed-predicted)^2/observed))}
rmaxse <- function(observed, predicted){sqrt(max((observed-predicted)^2/observed))}

##########################
# rrblup result
ffullnames=paste("data/rrblupresult/2020NEMI_rrblupfull_comparisons.csv",c(1:20),sep = "")
ffulls=lapply(ffullnames, read.csv)

fcvnames=paste("data/rrblupresult/NE_rrblupcv.csv",c(1:20),sep = "")
fcvs=lapply(fcvnames, read.csv)

# tested g untested e
corr_NEobMIbs <- sapply(ffulls, FUN = function(df){corr <- cor(df$MIyield,df$NEyield)})
rmse_NEobMIbs <- sapply(ffulls, FUN = function(df){rmse <- rmse(df$MIyield,df$NEyield)})
corr_NEobMIbs <- data.frame(Accuracy =corr_NEobMIbs, Method="Observed")
corr_NErrblupMIbs <- sapply(ffulls, FUN = function(df){corr <- cor(df$MIyield,df$GPyield)})
rmse_NErrblupMIbs <- sapply(ffulls, FUN = function(df){rmse <- rmse(df$MIyield,df$GPyield)})
corr_NErrblupMIbs <- data.frame(Accuracy =corr_NErrblupMIbs, Method="GP")
#testedg_untestede <- rbind(ob,gp)

# untested e tested g
corr_NErrblupcvNE <- sapply(fcvs, FUN = function(df){corr <- cor(df$TotalGrainMassGrams_NEmean,df$Yield_NEpredicted)})
rmse_NErrblupcvNE <- sapply(fcvs, FUN = function(df){rmse <- rmse(df$TotalGrainMassGrams_NEmean,df$Yield_NEpredicted)})
corr_NErrblupcvNE <- data.frame(Accuracy =corr_NErrblupcvNE, Method="GP")

# untested untested
corr_NErrblupcvMI <- sapply(fcvs, FUN = function(df){corr <- cor(df$Yield_MImean,df$Yield_NEpredicted)})
rmse_NErrblupcvMI <- sapply(fcvs, FUN = function(df){rmse <- rmse(df$Yield_MImean,df$Yield_NEpredicted)})
corr_NErrblupcvMI <- data.frame(Accuracy =corr_NErrblupcvMI, Method="GP")


# RF 2 datasets
RFNEcvs <- paste("data/RF/RF_NE_cv.csv",c(1:20),sep="")
RFNEcvs <- lapply(RFNEcvs, read.csv)
corr_RFNEcvs <- sapply(RFNEcvs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
rmse_RFNEcvs <- sapply(RFNEcvs, function(df){corr <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})
corr_RFNEcvs <- data.frame(Accuracy =corr_RFNEcvs, Method="PP")


RFcrosscvs <- paste("data/RF/RF_crossenv_cv.csv",c(1:20),sep="")
RFcrosscvs <- lapply(RFcrosscvs, read.csv)
corr_RFcrosscvs <- sapply(RFcrosscvs, function(df){corr <- cor(df$Yield_MIpredicted,df$Yield_MImean)})
rmse_RFcrosscvs <- sapply(RFcrosscvs, function(df){rmse <- rmse(df$Yield_MImean,df$Yield_MIpredicted)})
corr_RFcrosscvs <- data.frame(Accuracy =corr_RFNEcrosscvs, Method="PP")

# RF 1 dataset
RFNE5cvs <- paste("data/RF/RF_5cv.csv",c(1:20),sep="")
RFNE5cvs <- lapply(RFNE5cvs, read.csv)
corr_RFNE5cvs <- sapply(RFNE5cvs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
rmse_RFNE5cvs <- sapply(RFNE5cvs, function(df){corr <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})
corr_RFNE5cvs <- data.frame(Accuracy =corr_RFNE5cvs, Method="PP")

corr_RFcross5cvs <- sapply(RFNE5cvs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_MImean)})
rmse_RFcross5cvs <- sapply(RFNE5cvs, function(df){rmse <- rmse(df$Yield_MImean,df$Yield_NEpredicted)})
corr_RFcross5cvs <- data.frame(Accuracy =corr_RFcross5cvs, Method="PP")

RFNEbs <- paste("data/RF/RF_bs.csv",c(1:20),sep="")
RFNEbs <- lapply(RFNEbs, read.csv)
corr_RFNEbs <- sapply(RFNEbs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
rmse_RFNEbs <- sapply(RFNEbs, function(df){corr <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})
corr_RFNEbs <- data.frame(Accuracy =corr_RFNEbs, Method="PP")

corr_RFcrossbs <- sapply(RFNEbs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_MImean)})
rmse_RFcrossbs <- sapply(RFNEbs, function(df){rmse <- rmse(df$Yield_MImean,df$Yield_NEpredicted)})
corr_RFcrossbs <- data.frame(Accuracy =corr_RFcrossbs, Method="PP")

# rrblup CDA
CDArrblup_crossenv_cv=paste("data/GP/CDArrblup_crossenv_cv.csv",c(1:20),sep = "")
CDArrblup_crossenv_cv=lapply(CDArrblup_crossenv_cv, read.csv)
corr_CDAcrosscvs <- sapply(CDArrblup_crossenv_cv, function(df){corr <- cor(df$Yield_MIpredicted,df$Yield_MImean)})
rmse_CDAcrosscvs <- sapply(CDArrblup_crossenv_cv, function(df){rmse <- rmse(df$Yield_MImean,df$Yield_MIpredicted)})

CDArrblup_NE_cv=paste("data/GP/CDArrblup_NE_cv.csv",c(1:20),sep = "")
CDArrblup_NE_cv=lapply(CDArrblup_NE_cv, read.csv)
corr_CDANEcvs <- sapply(CDArrblup_NE_cv, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
rmse_CDANEcvs <- sapply(CDArrblup_NE_cv, function(df){rmse <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})


df2 <- rbind(corr_NEobMIbs, corr_NErrblupMIbs, corr_RFcrossbs)
comparisons2 <- list(c("Observed", "GP"), c("Observed", "PP"), c("GP", "PP"))
df3 <- rbind(corr_NErrblupcvNE, corr_RFNE5cvs)
df4 <- rbind(corr_NErrblupcvMI, corr_RFcross5cvs)

svg("compare2.svg")
ggboxplot(df2, x="Method", y="Accuracy",
          palette = "jco",
          add = "jitter")+
  stat_compare_means(comparisons = comparisons2, label.y = c(0.73, 0.765, 0.74))
dev.off()

svg("compare3.svg")
ggboxplot(df3, x="Method", y="Accuracy",
          palette = "jco",
          add = "jitter")+
  stat_compare_means(method = "t.test")
dev.off()

svg("compare4.svg")
ggboxplot(df4, x="Method", y="Accuracy",
          palette = "jco",
          add = "jitter")+
  stat_compare_means(method = "t.test")
dev.off()
