library(rpart)
library(tidyverse)
library(lares)
source("script/x2yfunctions.R")

avgNE=read.csv("data/2020NE_avg.csv")
avgMI=read.csv("data/2020MI_avg.csv")
avg_merge=read.csv("data/2020merge_avg_narm.csv")

# svg("plot/NEvsMI.svg")
# ggplot(data=avg_merge2,aes(x=value,y=Yield_MImean))+geom_point()+
#   #xlim(45,95)+
#   theme_bw()+xlab("Average NE yield (g)")+ylab("Average MI yield (g)")+
#   stat_smooth(method = "lm")+
#   stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),size=10) + theme(text = element_text(size = 30))
# #stat_regline_equation(
# # aes(label =  paste( ..r.label.., sep = "~~~~")))
# dev.off()

x2yNE <- x2y(avgNE, target_y = T, target = "TotalGrainMassGrams_NEmean",confidence=T, top=NULL, max_cat=30)
x2yMI <- x2y(avgMI, target_y = T, target = "Yield_MImean",confidence=T)

write.table(x2yNE, "x2yNE2.csv", quote = F, row.names = F, sep = ",")
write.table(x2yMI, "x2yMI.csv", quote = F, row.names = F, sep = ",")

plot.x2ypreds <- function(x, corr = FALSE, ...) {
  if (!inherits(x, "x2y_preds")) stop("Object must be class x2y_preds")
  p <- ggplot(x, aes(x = .data$x)) +
    geom_point(aes(y = .data$y), size = 0.5) +
    geom_line(aes(y = .data$p),
              colour = names(lares_pal()[[2]])[2],
              size = 0.8, alpha = 0.7
    ) +
    scale_color_brewer(name = NULL) +
    labs(
      title = "x's predictive power over y",
      subtitle = sprintf("x2y: %s", x2y_metric(x$x, x$y)$x2y)
    ) +
    theme_bw()
  if (corr & is.numeric(x$x) & is.numeric(x$y)) {
    p <- p + labs(caption = paste("Correlation:", signif(cor(x$x, x$y), 1))) +
      geom_smooth(aes(y = .data$y), method = "lm", formula = "y ~ x", size = 0.5)
  }
  return(p)
}


predictors0 = c("DaysToPollen_NEmean",        "DaysToSilk_NEmean"         , "RootLodgingPct_NEmean",     
               "StalkLodgingPct_NEmean",     "LeafLengthCM_NEmean"      ,  "LeafWidthCM_NEmean"    ,    
               "PlantHeightCM_NEmean" ,      "NodesWithBraceRoots_NEmean", "TillersPerPlant_NEmean" ,   
                "BranchesPerTassel_NEmean",   "TasselLengthCM_NEmean"    ,  "BranchZoneLengthCM_NEmean", 
               "TasselSpikeLengthCM_NEmean", "LeafNumber_NEmean" )

predictors = c("Days to pollen",        "Days to silk"         , "Root lodging rate",     
               "Stalk lodging rate",     "Leaf length"      ,  "Leaf width"    ,    
               "Plant height" ,      "Number of nodes with brace roots", "Tillers per plant" ,   
               "Branches per tassel",   "Tassel length"    ,  "Tassel branch zone length", 
               "Tassel spike length", "Leaf number" )

colnames(avg_merge)[colnames(avg_merge) %in% predictors0] <- predictors

for(i in predictors){
  px2yNE <- plot.x2ypreds(x2y_preds(avg_merge[,i], avg_merge$TotalGrainMassGrams_NEmean), corr=T)
  px2yNE <- px2yNE+labs(x=i, y="Average NE yield (g)")+
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  ggsave(paste0(i,"x2yNE.svg"), plot=px2yNE, width=3.5, height=3.5, units="in")
  px2yMI <- plot.x2ypreds(x2y_preds(avg_merge[,i], avg_merge$Yield_MImean), corr=T)
  px2yMI <- px2yMI+labs(x=i, y="Average MI yield (g)")+
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  ggsave(paste0(i,"x2yMI.svg"), plot=px2yMI, width=3.5, height=3.5, units="in")
}

plot.scattercurve <- function(x, corr = FALSE, ...) {
  if (!inherits(x, "x2y_preds")) stop("Object must be class x2y_preds")
  p <- ggplot(x, aes(x = .data$x)) +
    geom_point(aes(y = .data$y), size = 0.5) +
    geom_smooth(aes(y = .data$y),
              colour = names(lares_pal()[[2]])[2],
              size = 0.8, alpha = 0.7
    ) +
    scale_color_brewer(name = NULL) +
    labs(
      title = "x's predictive power over y",
      subtitle = sprintf("x2y: %s", x2y_metric(x$x, x$y)$x2y)
    ) +
    theme_bw()
  if (corr & is.numeric(x$x) & is.numeric(x$y)) {
    p <- p + labs(caption = paste("Correlation:", signif(cor(x$x, x$y), 1))) +
      geom_smooth(aes(y = .data$y), method = "lm", formula = "y ~ x", size = 0.5)
  }
  return(p)
}

for(i in predictors){
  px2yNE <- plot.x2ypreds(x2y_preds(avg_merge[,i], avg_merge$TotalGrainMassGrams_NEmean), corr=T)
  px2yNE <- px2yNE+labs(x=i, y="Average NE yield (g)")+
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  ggsave(paste0(i,"x2yNEscatter.svg"), plot=px2yNE, width=3.5, height=3.5, units="in")
  px2yMI <- plot.x2ypreds(x2y_preds(avg_merge[,i], avg_merge$Yield_MImean), corr=T)
  px2yMI <- px2yMI+labs(x=i, y="Average MI yield (g)")+
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  ggsave(paste0(i,"x2yMIscatter.svg"), plot=px2yMI, width=3.5, height=3.5, units="in")
}
