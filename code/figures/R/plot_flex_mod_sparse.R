# This file is used for creating the plot for the flexibility and modularity at different thretholds
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(R.matlab)

# Load the dataset
mod_flex_df <- readMat("flex_df_threthold.mat")
mod_df <- mod_flex_df[['modularity.mean.threth.mean']]
flex_df <- mod_flex_df[['flex.df']]

# Construct the modularity dataframe
mod <- data.frame('mod' = c(mod_df),
                  'condition' = rep(c('Low-Difficulty', 'Balanced-Difficulty', 'High-Difficulty') ,5),
                  'threthold' = rep(c('10%', '20%', '30%', '40%', '50%'), each=3))
mod <- mod %>% filter(threthold %in% c('10%', '20%', '30%', '40%', '50%'))
mod$condition <- factor(mod$condition, levels = c('Low-Difficulty', 'Balanced-Difficulty', 'High-Difficulty'))

# Construct the flexibility dataframe
flex_01 <- unlist(flex_df[[1]][[1]])
flex_02 <- unlist(flex_df[[2]][[1]])
flex_03 <- unlist(flex_df[[3]][[1]])
flex_04 <- unlist(flex_df[[4]][[1]])
flex_05 <- unlist(flex_df[[5]][[1]])

flex <- data.frame('flex' = c(flex_01, flex_02, flex_03, flex_04, flex_05),
                   'condition' = rep(rep(c('Low-Difficulty', 'Balanced-Difficulty', 'High-Difficulty'), each=4),5),
                   'brain' = rep(c('GB', 'FPCN', 'RN', 'FPRN'), 3*5),
                   'threthold' = rep(c('10%', '20%', '30%', '40%', '50%'), each=3*4))

flex$condition <- as.factor(flex$condition)
flex$brain <- as.factor(flex$brain)
flex$threthold <- as.factor(flex$threthold)

# plot modularity
p_mod <- mod %>% mutate(condition = fct_relevel(condition,
                                                'Low-Difficulty',
                                                'Balanced-Difficulty',
                                                'High-Difficulty')) %>%
  ggplot(aes(x=condition, y=mod, group=threthold)) +
  geom_line(aes(color=threthold))+
  ylab('Modularity') +
  xlab('Condition') + 
  geom_point()+ 
  ggtitle('Figure 2A') +
  theme(legend.position = "none")

# plot flexibility
p_flex <- flex %>% mutate(condition = fct_relevel(condition,
                                                  'Low-Difficulty',
                                                  'Balanced-Difficulty',
                                                  'High-Difficulty')) %>%
  ggplot(aes(x=condition, y=flex, group=threthold)) +
  geom_line(aes(color=threthold))+
  ylab("Flexibility") +
  geom_point() + 
  xlab('Condition') + 
  facet_wrap(~brain) + 
  ggtitle('Figure 2B') +
  theme(axis.text.x = element_text(angle = 13, hjust = 1),
        legend.title = element_blank()) 

# plot the two fitures and save figure as png file
grid.arrange(p_mod,p_flex,ncol = 2, widths = c(1.2,2))
g <- arrangeGrob(p_mod,p_flex,ncol = 2, widths = c(1.2,2))
ggsave("Figure_flex_mod.png", g, dpi=300, dev='png', height=4, width=8, units="in")
