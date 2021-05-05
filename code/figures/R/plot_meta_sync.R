library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(R.matlab)

# Load the csv files for each meta sync values for each filters
data_06125 <- read.csv('06125.csv', header=FALSE)
data_0407 <- read.csv('0407.csv', header=FALSE)
data_01608 <- read.csv('01608.csv', header=FALSE)

# Load the subject ID
sub_number <- c("sub-005",
                "sub-006",
                "sub-007",
                "sub-008",
                "sub-009",
                "sub-011",
                "sub-012",
                "sub-013",
                "sub-014",
                "sub-015",
                "sub-016",
                "sub-017",
                "sub-018",
                "sub-019",
                "sub-020",
                "sub-021",
                "sub-022",
                "sub-023",
                "sub-024",
                "sub-025",
                "sub-026",
                "sub-028",
                "sub-029",
                "sub-030",
                "sub-031",
                "sub-032",
                "sub-033",
                "sub-034",
                "sub-035",
                "sub-036",
                "sub-037",
                "sub-039",
                "sub-040",
                "sub-041",
                "sub-042")

# Construct the dataframe for synchrony 
sync_06125 <- unlist(data_06125[,1:12])
meta_06125 <- unlist(data_06125[,13:24])

sync_0407 <- unlist(data_0407[,1:12])
meta_0407 <- unlist(data_0407[,13:24])

sync_01608 <- unlist(data_01608[,1:12])
meta_01608 <- unlist(data_01608[,13:24])

data_sync <- data.frame("Sync" = c(sync_06125, sync_0407, sync_01608), "Participant" = as.factor(rep(sub_number, 12*3)),
                        "Condition" = as.factor(rep(rep(1:3, each=35),4*3)),
                        "Brain" = as.factor(rep(rep(1:4, each = 35*3),3)),
                        "Bandpass_Filter" = as.factor(rep(1:3, each=35*3*4)))

levels(data_sync$Condition) <- c("Boredom", "Flow", "Frustration")
levels(data_sync$Brain) <- c("Global", "FPCN", "Reward", "FPCN&Reward")
levels(data_sync$Bandpass_Filter) <- c("06-125 Hz", "04-07 Hz", "016-08 Hz")
data_sync_groupby <- data_sync %>% group_by(Bandpass_Filter, Brain)
splited_sync <- group_split(data_sync_groupby)

# Construct the dataframe for metastability 
data_meta <- data.frame("Metastability" = c(meta_06125, meta_0407, meta_01608), "Participant" = as.factor(rep(sub_number, 12*3)),
                        "Condition" = as.factor(rep(rep(1:3, each=35),4*3)),
                        "Brain" = as.factor(rep(rep(1:4, each = 35*3),3)),
                        "Bandpass_Filter" = as.factor(rep(1:3, each=35*3*4)))

levels(data_meta$Condition) <- c("Boredom", "Flow", "Frustration")
levels(data_meta$Brain) <- c("Global", "FPCN", "Reward", "FPCN&Reward")
levels(data_meta$Bandpass_Filter) <- c("06-125 Hz", "04-07 Hz", "016-08 Hz")
data_meta_groupby <- data_meta %>% group_by(Bandpass_Filter, Brain)
splited_meta <- group_split(data_meta_groupby)

# fit the anova models for synchrony
sync_models = c()
for (sync_df in splited_sync) {
  model <- anova_test(data = sync_df,
                      dv = sync, wid = participant, within = condition)
  sync_models = c(sync_models, get_anova_table(model)$p)
}
sync_models

# fit the anova models for metastability
meta_models = c()
for (meta_df in splited_meta) {
  model <- anova_test(data = meta_df,
                      dv = meta, wid = participant, within = condition)
  meta_models = c(meta_models, get_anova_table(model)$p)
}
meta_models

# Plot the synchrony boxplot
p_sync <- data_sync %>% 
  ggplot(aes(x= Brain, y=dv,  fill=Condition)) +
  geom_boxplot() +
  geom_jitter(width=0.3,alpha=0.2) +
  xlab("")+ 
  ylab("Synchrony")+ 
  facet_wrap(~Bandpass_Filter) +
  theme(axis.text.x=element_blank())  + 
  ggtitle("Synchrony & Metastability across condition, brain region, and bandpass filter")

# Plot the metastability boxplot
p_meta <- data_meta %>% 
  ggplot(aes(x= Brain, y=dv,  fill=Condition), show.legend = FALSE) +
  geom_boxplot() +
  geom_jitter(width=0.3,alpha=0.2) +
  xlab("Bandpass Filter")+ 
  ylab("Metastability")+ 
  facet_wrap(~Bandpass_Filter) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Arrange the two plots and show them
grid.arrange(p_sync,p_meta,nrow = 2)