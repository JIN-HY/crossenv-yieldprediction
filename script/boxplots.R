library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(cowplot)


RMSE <- function(observed, predicted){sqrt(sum((observed-predicted)^2)/length(observed))}
#rmdse <- function(observed, predicted){sqrt(median((observed-predicted)^2/observed))}
#rmaxse <- function(observed, predicted){sqrt(max((observed-predicted)^2/observed))}

##########################
# rrblup result
ffullnames=paste("data2/rrblupresult/NE_rrblupfullbs.csv",c(1:30),sep = "")
ffulls=lapply(ffullnames, read.csv)

fcvnames=paste("data2/rrblupresult/NE_rrblupcv.csv",c(1:30),sep = "")
fcvs=lapply(fcvnames, read.csv)

# tested g untested e
corr_NEobMIbs <- sapply(ffulls, FUN = function(df){corr <- cor(df$MIyield,df$NEyield)})
rmse_NEobMIbs <- sapply(ffulls, FUN = function(df){rmse <- RMSE(df$MIyield*2,df$NEyield)})
corr_NEobMIbs <- data.frame(`R.squared` =corr_NEobMIbs^2, Method="Observed NE")
RMSE_NEobMIbs <- data.frame(RMSE =rmse_NEobMIbs, Method="Observed NE")
corr_NErrblupMIbs <- sapply(ffulls, FUN = function(df){corr <- cor(df$MIyield,df$GPyield)})
rmse_NErrblupMIbs <- sapply(ffulls, FUN = function(df){rmse <- RMSE(df$MIyield*2,df$GPyield)})
corr_NErrblupMIbs <- data.frame(R.squared =corr_NErrblupMIbs^2, Method="Genomic prediction")
RMSE_NErrblupMIbs <- data.frame(RMSE =rmse_NErrblupMIbs, Method="Genomic prediction")
#testedg_untestede <- rbind(ob,gp)

# # untested e tested g
# corr_NErrblupcvNE <- sapply(fcvs, FUN = function(df){corr <- cor(df$TotalGrainMassGrams_NEmean,df$Yield_NEpredicted)})
# rmse_NErrblupcvNE <- sapply(fcvs, FUN = function(df){rmse <- rmse(df$TotalGrainMassGrams_NEmean,df$Yield_NEpredicted)})
# corr_NErrblupcvNE <- data.frame(Accuracy =corr_NErrblupcvNE, Method="Genomic prediction")
# RMSE_NErrblupcvNE <- data.frame(RMSE =RMSE_NErrblupcvNE, Method="Genomic prediction")

# untested untested
corr_NErrblupcvMI <- sapply(fcvs, FUN = function(df){corr <- cor(df$Yield_MI,df$Yield_NEpredicted)})
rmse_NErrblupcvMI <- sapply(fcvs, FUN = function(df){rmse <- RMSE(df$Yield_MI*2,df$Yield_NEpredicted)})
corr_NErrblupcvMI <- data.frame(R.squared =corr_NErrblupcvMI^2, Method="Genomic prediction")
RMSE_NErrblupcvMI <- data.frame(RMSE =rmse_NErrblupcvMI, Method="Genomic prediction")


# RF 1 dataset
RFNE5cvs <- paste("data2/RF/RF_5cv_MI.csv",c(1:30),sep="")
RFNE5cvs <- lapply(RFNE5cvs, read.csv)
# corr_RFNE5cvs <- sapply(RFNE5cvs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
# rmse_RFNE5cvs <- sapply(RFNE5cvs, function(df){corr <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})
# corr_RFNE5cvs <- data.frame(Accuracy =corr_RFNE5cvs, Method="Phenotypic prediction")

corr_RFcross5cvs <- sapply(RFNE5cvs, function(df){corr <- cor(df$Yield_MIpredicted,df$Yield_MImean)})
rmse_RFcross5cvs <- sapply(RFNE5cvs, function(df){rmse <- RMSE(df$Yield_MImean,df$Yield_MIpredicted)})
corr_RFcross5cvs <- data.frame(R.squared =corr_RFcross5cvs^2, Method="Phenotypic prediction")
RMSE_RFcross5cvs <- data.frame(RMSE =rmse_RFcross5cvs, Method="Phenotypic prediction")


RFNEbs <- paste("data2/RF/RF_bs.csv",c(1:30),sep="")
RFNEbs <- lapply(RFNEbs, read.csv)
# corr_RFNEbs <- sapply(RFNEbs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_NEmean)})
# rmse_RFNEbs <- sapply(RFNEbs, function(df){corr <- rmse(df$Yield_NEmean,df$Yield_NEpredicted)})
# corr_RFNEbs <- data.frame(Accuracy =corr_RFNEbs, Method="Phenotypic prediction")


corr_RFcrossbs <- sapply(RFNEbs, function(df){corr <- cor(df$Yield_NEpredicted,df$Yield_MImean)})
rmse_RFcrossbs <- sapply(RFNEbs, function(df){rmse <- RMSE(df$Yield_MImean,df$Yield_NEpredicted)})
corr_RFcrossbs <- data.frame(R.squared =corr_RFcrossbs^2, Method="Phenotypic prediction")
RMSE_RFcrossbs <- data.frame(RMSE =rmse_RFcrossbs, Method="Phenotypic prediction")



# Example ANOVA model without AR(1) covariate
anova_model <- lm(`R.squared` ~ Method, data = df2)

# Plot the residuals and fitted values
plot(anova_model$residuals, type = "o", col = "blue", pch = 16, main = "AR(1) Model for Residuals")
#lines(fitted(ar_model), col = "red", type = "o", pch = 16)

# Diagnostic plots
par(mfrow = c(2, 2))
#plot(ar_model)






df2 <- rbind(corr_NEobMIbs, corr_NErrblupMIbs, corr_RFcrossbs)
comparisons2 <- list(c("Observed NE", "Genomic prediction"), c("Observed NE", "Phenotypic prediction"), c("Genomic prediction", "Phenotypic prediction"))
#df3 <- rbind(corr_NErrblupcvNE, corr_RFNE5cvs)
df4 <- rbind(corr_NErrblupcvMI, corr_RFcross5cvs)
comparison4 <- list(c("Genomic prediction", "Phenotypic prediction"))

library(car)

anova_model <- lm(`R.squared` ~ Method, data = df2)
dw <- durbinWatsonTest(anova_model, reps = 30, method="resample")
ar_model <- arima(anova_model$residuals, order = c(1, 0, 0))
summary(ar_model)
anova_ar_model <- aov(df2$`R.squared` ~ df2$Method + ar_model$residuals)
summary(anova_ar_model)
TukeyHSD(anova_ar_model)

anova_model <- lm(`R.squared` ~ Method, data = df4)
dw <- durbinWatsonTest(anova_model, reps = 30, method="resample")
ar_model <- arima(anova_model$residuals, order = c(1, 0, 0))
summary(ar_model)
anova_ar_model <- aov(df4$`R.squared` ~ df4$Method + ar_model$residuals)
summary(anova_ar_model)



p2d <- ggboxplot(df2, x="Method", y="R.squared",
          palette = "jco",
          add = "jitter")+ylab('R-squared')+ylim(0.43,0.58)+
  stat_compare_means(comparisons = comparisons2, label.y = c(0.52, 0.555, 0.54),label = "p.signif",size=5)+
  xlab(' ')+theme(text = element_text(size = 20),axis.text.y = element_text(size=14),
                  axis.text.x = element_text(angle = 25, hjust=0.5, vjust=0.5))


p3c <- ggboxplot(df4, x="Method", y="`R.squared`",
          palette = "jco",
          add = "jitter")+ylab('R-squared')+ylim(0.1,0.45)+
  #stat_compare_means(method = "t.test")+xlab(" ")
  stat_compare_means(comparisons=comparison4, label = "p.signif",size=5)+xlab(" ")+
  theme(text = element_text(size = 20),axis.text.y = element_text(size=14),
        axis.text.x = element_text(angle = 25, hjust=0.5, vjust=0.5))









df2E <- rbind(RMSE_NEobMIbs, RMSE_NErrblupMIbs, RMSE_RFcrossbs)
#comparisons2 <- list(c("Observed", "GP"), c("Observed", "PP"), c("GP", "PP"))
df4E <- rbind(RMSE_NErrblupcvMI, RMSE_RFcross5cvs)


anova_model <- lm(RMSE ~ Method, data = df2E)
dw <- durbinWatsonTest(anova_model, reps = 30, method="resample")
ar_model <- arima(anova_model$residuals, order = c(1, 0, 0))
summary(ar_model)
anova_ar_model <- aov(df2E$RMSE ~ df2E$Method + ar_model$residuals)
summary(anova_ar_model)
TukeyHSD(anova_ar_model)

anova_model <- lm(RMSE ~ Method, data = df4E)
dw <- durbinWatsonTest(anova_model, reps = 30, method="resample")
ar_model <- arima(anova_model$residuals, order = c(1, 0, 0))
summary(ar_model)
anova_ar_model <- aov(df4E$RMSE ~ df4E$Method + ar_model$residuals)
summary(anova_ar_model)
TukeyHSD(anova_ar_model)

p2e <- ggboxplot(df2E, x="Method", y="RMSE",
          palette = "jco",
          add = "jitter")+ylim(110,138)+
  stat_compare_means(comparisons = comparisons2, label.y = c(132, 135, 123), label = "p.signif", size=5)+
  xlab(' ')+theme(text = element_text(size = 20),axis.text.y = element_text(size=14),
                  axis.text.x = element_text(angle = 25, hjust=0.5, vjust=0.5))





p3d <- ggboxplot(df4E, x="Method", y="RMSE",
          palette = "jco",
          add = "jitter")+ylim(120,145)+
  stat_compare_means(comparisons = comparison4, label = "p.signif", size=5)+
  xlab(" ")+theme(text = element_text(size = 20),axis.text.y = element_text(size=14),
                  axis.text.x = element_text(angle = 25, hjust=0.5, vjust=0.5))


panel2.4 <- plot_grid(p2d, p2e, ncol = 2)
p2.1 <- plot_grid(p2a, p2b, p2c, ncol=3)
panel2 <- plot_grid(p2.1, panel2.4, ncol = 1)
panel3 <- plot_grid(p3a, p3b, p3c, p3d, ncol=2)
