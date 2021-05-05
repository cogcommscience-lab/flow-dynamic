# this file is used for create the plot for the dense weighted graph analysis

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(R.matlab)

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

# Load the dataset
indi_mod_flex_df <- readMat("dense_individual_mod_flex.mat")
indi_mod_flex_df <- indi_mod_flex_df$individual.mod.fex

flex_df_ind <- indi_mod_flex_df[,-c(1:3)]
mod_df_ind <- indi_mod_flex_df[,1:3]

# Create the flexibility dataframe
data_mod_flex <- data.frame("flex" = c(flex_df_ind), "Participant" = as.factor(rep(sub_number, 12)),
                            "Condition" = as.factor(rep(rep(1:3, each=35),4)),
                            "Brain" = as.factor(rep(rep(1:4, each = 35*3))))
levels(data_mod_flex$Condition) <- c('bore', 'flow', 'frus')
levels(data_mod_flex$Brain) <- c('GB', 'FPCN', 'RN', 'FPCN&RN')

# plot the flexibility boxplot
p_flex_ind <- data_mod_flex %>% 
  ggplot(aes(x= Brain, y=flex,  fill=Condition), show.legend = FALSE) +
  geom_boxplot() +
  geom_jitter(width=0.3,alpha=0.2) +
  xlab("Brain Region")+ 
  ylab("Flexibility") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p_flex_ind

# Create the modularity dataframe
data_mod_ind <- data.frame("mod" = c(mod_df_ind), "Participant" = as.factor(rep(sub_number, 3)),
                           "Condition" = as.factor(rep(rep(1:3, each=35))))
levels(data_mod_ind$Condition) <- c('bore', 'flow', 'frus')

# plot the modularity boxplot
p_mod_ind <- data_mod_ind %>% 
  ggplot(aes(x = Condition, y=mod)) +
  geom_boxplot() +
  geom_jitter(width=0.3,alpha=0.2) +
  xlab("experimental condition")+ 
  ylab("Modularity")
p_mod_ind

