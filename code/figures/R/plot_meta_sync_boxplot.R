library(ggplot2)
library(ggsignif)
library(cowplot)
meta <- read.csv("~/OneDrive/projects/flow/analysis/figures/python/meta_01608.csv")
sync <- read.csv("~/OneDrive/projects/flow/analysis/figures/python/sync_01608.csv")

sub <- meta[,1]
meta <- meta[,-1]
sync <- sync[,-1]
metastability <- unlist(meta)
synchrony <- unlist(sync)
brain_region <- rep(c('Global', 'FPCN', 'Reward', 'FPCN&Reward'), each=3*35)
Condition <- rep(rep(c('Boredom', 'Flow', 'Frustration'), each=35), 4)
subject <- rep(sub, 12)

meta_data_df = data.frame(subject, metastability, brain_region, Condition)
meta_data_df$brain_region <- factor(meta_data_df$brain_region , levels=c("Global", "FPCN", "Reward", "FPCN&Reward"))
sync_data_df = data.frame(subject, metastability, brain_region, Condition)
sync_data_df$brain_region <- factor(sync_data_df$brain_region , levels=c("Global", "FPCN", "Reward", "FPCN&Reward"))


meta_plot <-ggplot(meta_data_df, aes(x=brain_region, y=metastability, fill=Condition)) + 
  geom_boxplot() + 
  ggtitle("Boxplot of Metastability")+xlab("Brain Region")+ylab("Metastability")+
  geom_signif(
    y_position = c(0.26, 0.26,
                   0.26, 0.26,
                   0.26, 0.26,
                   0.26, 0.26), 
          xmin = c(0.7, 1.05,
                   1.7, 2.05,
                   2.7, 3.05,
                   3.7, 4.05),
          xmax = c(0.95, 1.3,
                   1.95, 2.3,
                   2.95, 3.3,
                   3.95, 4.3),
    
    annotation = c("NS", "NS",
                   "NS", "NS",
                   "NS", "NS",
                   "NS", "*"), 
    textsize = 5,
    
    tip_length = 0.02
  ) + theme(text = element_text(size = 20))  + theme(legend.position = c(0.8, 0.15))


sync_plot <- ggplot(sync_data_df, aes(x=brain_region, y=synchrony, fill=Condition)) + 
  geom_boxplot() + 
  ggtitle("Boxplot of Synchrony")+xlab("Brain Region")+ylab("Synchrony")+
  theme(text = element_text(size = 20))+
  geom_signif(
    y_position = c(0.6, 0.6,
                   0.6, 0.6,
                   0.6, 0.6,
                   0.6, 0.6), 
    xmin = c(0.7, 1.05,
             1.7, 2.05,
             2.7, 3.05,
             3.7, 4.05),
    xmax = c(0.95, 1.3,
             1.95, 2.3,
             2.95, 3.3,
             3.95, 4.3),
    
    annotation = c("*", "NS",
                   "NS", "*",
                   "NS", "NS",
                   "NS", "*"), 
    textsize = 5,
    
    tip_length = 0.02
  ) + theme(legend.position = "none")



plot_grid(sync_plot, meta_plot, labels = "AUTO")


