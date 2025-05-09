## Figure 1c
# Written by: Rosemary Yu
# Last updated: 30/03/2025

# Load in packages:

library(tidyverse) # tidyverse_2.0.0
library(readr) # readr_2.1.5
library(ggplot2) # ggplot2_3.5.1

# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)

# Load in data

cell_counts <- read_csv("data/cell_count/IncuCyte_cell_counts_MCF10A.csv")
cell_counts$replicate <- as.factor(cell_counts$replicate)
metab <- read_csv("data/MCF10A_metabolomics_data_raw.csv")
metab$replicate <- as.factor(metab$replicate)


cell_count_t0 <- 3e5 #number of seeded cells
gDW_t0 <- cell_count_t0 * 397e-12 / 1.5e-3 # number of seeded cells * cell dry weight / culture volume
gly_t0 <- 0.25 #mM in DMEM/F12
glu_t0 <- 0.05 #mM in DMEM/F12


#panel A: growth curve
cell_counts <- rbind(c(0,3e5,1),
                     c(0,3e5,2),
                     cell_counts)
count_means <- cell_counts %>% 
  group_by(timepoint) %>%
  summarize(mean_count = mean(cell_count)) %>%
  mutate(ln_mean_count = log(mean_count)) %>%
  mutate(phase = c(rep('lag',15),rep('exp',20),rep('stn',15)))

ggplot(count_means) + 
  geom_line(aes(x = timepoint, y=ln_mean_count),color='grey', size=3) + 
  geom_point(aes(x=timepoint, y=ln_mean_count, shape=phase), size=3, stroke=0.5, fill='grey') + 
  scale_shape_manual(values=c(21,21,21)) +
  geom_smooth(data=subset(count_means,phase=="exp"), 
              aes(x=timepoint, y=ln_mean_count, color=phase), 
              method=lm, se=F, fullrange=T, size=2) +
  annotate('rect', xmin=0, xmax=12.5, ymin=12.25, ymax=14.25, alpha=.1, fill='blue') +
  annotate('rect', xmin=17.5, xmax=40, ymin=12.25, ymax=14.25, alpha=.1, fill='red') +
  xlim(0,55) + ylim(12.25,14.25) + 
  theme_bw() + theme(legend.position='none') +
  labs(x='time', y='ln(mean cell count)')

ggsave("figures/Fig2_growth.pdf", width=3.7, height=3.7)



#panel B: glutamine
glu <- metab[metab$metabolite=="Glu",2:4]
glu <- rbind(c(0,3,glu_t0/gDW_t0), glu)

#set up for the dashed lines
b <- summary(lm(unlist(glu[2:6,3]) ~ unlist(glu[2:6,1])))$coefficients[1,1]
m <- summary(lm(unlist(glu[2:6,3]) ~ unlist(glu[2:6,1])))$coefficients[2,1] #first value b, second value m
glu_t15 <- m*15+b
glu_t35 <- m*35+b
remove(m,b)

#plot
ggplot(glu) + 
  geom_point(aes(x=timepoint, y=mmol_gDW, shape=replicate), size=3, stroke=1, fill='grey') + 
  scale_shape_manual(values=21) +
  xlim(0,55) + ylim(0,5.5) +
  annotate('rect', xmin=0, xmax=12.5, ymin=0, ymax=5.5, alpha=.1, fill='blue') +
  annotate('rect', xmin=17.5, xmax=40, ymin=0, ymax=5.5, alpha=.1, fill='red') +
  annotate('segment', x=0, xend=15, y=unlist(glu[1,3]), yend=glu_t15, 
           alpha=.5, color='blue',size=2,linetype='dashed') +
  annotate('segment', x=15, xend=35, y=glu_t15, yend=glu_t35, 
           alpha=.5, color='red',size=2) +
  theme_bw() + theme(legend.position='none') +
  labs(x='time', y='mmol glutamate in culture media, \nnormalized to cell count and cell dry weight)')

ggsave("figures/Fig2_glu.pdf", width=3.7, height=3.7)

#panel C: glycine
gly <- metab[metab$metabolite=="Gly",2:4]
gly <- rbind(c(0,2,gly_t0/gDW_t0),
             c(0,3,gly_t0/gDW_t0), gly)


#set up for the dashed lines
b <- summary(lm(unlist(gly[3:12,3]) ~ unlist(gly[3:12,1])))$coefficients[1,1]
m <- summary(lm(unlist(gly[3:12,3]) ~ unlist(gly[3:12,1])))$coefficients[2,1] #first value b, second value m
gly_t15 <- m*15+b
gly_t35 <- m*35+b
remove(m,b)

#plot
ggplot(gly) + 
  geom_point(aes(x=timepoint, y=mmol_gDW, shape=replicate), size=3, stroke=1, fill='grey') + 
  scale_shape_manual(values=c(21,24)) +
  xlim(0,55) + ylim(0,3.5) +
  annotate('rect', xmin=0, xmax=12.5, ymin=0, ymax=3.5, alpha=.1, fill='blue') +
  annotate('rect', xmin=17.5, xmax=40, ymin=0, ymax=3.5, alpha=.1, fill='red') +
  annotate('segment', x=0, xend=15, y=unlist(gly[1,3]), yend=gly_t15, 
           alpha=.5, color='blue',size=2,linetype='dashed') +
  annotate('segment', x=15, xend=35, y=gly_t15, yend=gly_t35, 
           alpha=.5, color='red',size=2) +
  theme_bw() + theme(legend.position='none') +
  labs(x='time', y='mmol glycine in culture media, \nnormalized to cell count and cell dry weight)')

ggsave("figures/Fig2_gly.pdf", width=3.7, height=3.7)
