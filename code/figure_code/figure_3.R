## Figure 3 (Plotting the min-max growth rate)
# Written by: Cyriel Huijer
# Last updated: 30 March 2025

# Load in packages:
library(ggbreak) # ggbreak_0.1.2
library(tidyverse) # tidyverse_2.0.0

# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)

# Panel A: 0% error

meas_0 <- read.csv("data/growth_rate_simulations/output/measured_fluxes_0_bounds.csv")

meas_0$Condition <- fct_rev(factor(meas_0$Condition, levels = unique(meas_0$Condition)))

ggplot(meas_0, aes(x = Condition)) +
  geom_errorbar(aes(ymin = MinSolution, ymax = MaxSolution), 
                width = 0.2, size = 1, alpha = 1) +
  geom_point(aes(y = GrowthRate), color = "red", size = 3) +  # Overlapping growth rate point
  coord_flip() +
  scale_y_continuous(name = "Growth Rate", labels = scales::comma_format()) +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank()) + 
  ylim(0,0.8) +
  ylab("Growth Rate")+
  scale_y_break(c(0.25,0.7), ticklabels=c(0.7,0.75,0.8))

ggsave("figures/Fig3a.pdf", width = 6, height = 6, units = "in")



# Panel B:  20% error

meas_20 <- read.csv("data/growth_rate_simulations/output/measured_fluxes_20_bounds.csv")

meas_20$Condition <- fct_rev(factor(meas_20$Condition, levels = unique(meas_20$Condition)))

ggplot(meas_20, aes(x = Condition)) +
  geom_errorbar(aes(ymin = MinSolution, ymax = MaxSolution), 
                width = 0.2, size = 1, alpha = 1) +
  geom_point(aes(y = GrowthRate), color = "red", size = 3) +  # Overlapping growth rate point
  coord_flip() +
  scale_y_continuous(name = "Growth Rate", labels = scales::comma_format()) +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank()) + 
  ylim(0,0.8) +
  ylab("Growth Rate")+
  scale_y_break(c(0.25,0.7), ticklabels=c(0.7,0.75,0.8))

ggsave("figures/Fig3b.pdf", width = 6, height = 6, units = "in")



# Panel C: 700% bounds: 

meas_700 <- read.csv("data/growth_rate_simulations/output/measured_fluxes_700_bounds.csv")

meas_700$Condition <- fct_rev(factor(meas_20$Condition, levels = unique(meas_700$Condition)))

ggplot(meas_700, aes(x = Condition)) +
  geom_errorbar(aes(ymin = MinSolution, ymax = MaxSolution), 
                width = 0.2, size = 1, alpha = 1) +
  geom_point(aes(y = GrowthRate), color = "red", size = 3) +  # Overlapping growth rate point
  coord_flip() +
  scale_y_continuous(name = "Growth Rate", labels = scales::comma_format()) +
  theme_bw(base_size = 16) +
  theme(axis.title.y = element_blank()) +
  ylim(0,0.8) +
  ylab("Growth Rate") +
  scale_y_break(c(0.25,0.7), ticklabels=c(0.7,0.75,0.8))

ggsave("figures/Fig3c.pdf", width = 6, height = 6, units = "in")




