## Figure 4 (FVA heatmap):
# Written by: Cyriel Huijer
# Last updated: 30 March 2025

# Load in libraries:
library(tidyverse) # tidyverse_2.0.0
library(ComplexHeatmap) # ComplexHeatmap_2.15.4
library(circlize) # circlize_0.4.16

# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)


# Load in data:
z_scores <- read.csv('data/fva_heatmap/z_scores_heatmap_20percent_GP.csv')

# Convert to numeric matrix and transpose
z_scores_matrix <- as.matrix(z_scores)
z_scores_matrix <- t(z_scores_matrix)  

z_scores_matrix_capped <- pmax(pmin(z_scores_matrix, 3), -3)
pdf("figures/Fig4.pdf", width = 15, height = 15)
Heatmap(z_scores_matrix_capped, 
        cluster_columns = FALSE, 
        col = colorRamp2(c(-3, 0, 3), c("#228be6", "white", "#f03e3e")))  
dev.off()


