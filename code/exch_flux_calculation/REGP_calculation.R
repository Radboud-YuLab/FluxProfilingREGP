## Script calculation of REGP exchange fluxes:
# Written by: Cyriel Huijer
# Last updated: 31 March 2025

# Load in packages:
library(tidyverse)
# Set working directory:

root <- "" # This should be the base directory of the github repo
setwd(root)

# Functions used:
calculate_slope <- function(data) {
  fit <- lm(mmol_cell ~ timepoint, data = data)
  slope <- coef(fit)["timepoint"]
  R2 <- summary(fit)$r.squared
  return(list(slope = slope, R2 = R2))
}

# Load in data:
metab_concentration <- read.csv("data/metabolomics/calculated_metabolite_concentrations.csv")
manual_inspection <- read.csv("data/metabolomics/manual_inspection_curves.csv") # Manual inspection was done to determine low quality samples
cell_count <- read.csv("data/cell_count/CyQUANT_cellCount.csv") # CyQUANT derived cell counts

metab_concentration <- merge(metab_concentration, manual_inspection, by = c("batch", "amino_acid"))
metab_concentration <- subset(metab_concentration, keep_manual_inspection != "n")

# cells were grown in 1.5 mL, calculate amount of metabolite in total 
metab_concentration$final_concentration_mM <- metab_concentration$final_concentration_uM/1000
metab_concentration$volume_corr <- metab_concentration$final_concentration_mM *0.0015
# cells were grown in 6-well plate, correct for this:
cell_count$cell_number_6_well <- cell_count$Cell_Number * (9.6/0.32)
merged_df <- left_join(metab_concentration, cell_count, by = c("timepoint"))
merged_df$mmol_cell <- merged_df$volume_corr / merged_df$cell_number_6_well
# Calculate mmol/gDW, cell dry weight was determined to be 397 pg
merged_df$mmol_gDW <- merged_df$mmol_cell * (1 / (397 * 10^-12))

# make condition column
merged_df <- merged_df %>%
  filter(cell_line %in% c("E","M","T","TM"),dox_conc %in% c(0, 100,1000))
merged_df$condition <- paste0(merged_df$cell_line,"_",as.character(merged_df$dox_conc))

slopes_df <- merged_df %>%
  group_by(condition, amino_acid, replicate) %>%
  summarize(
    slope = calculate_slope(cur_data())$slope,
    R2 = calculate_slope(cur_data())$R2,
    .groups = 'drop'
  ) %>%
  filter(R2 > 0.7) # Slope has to be over to be included in further analysis
slopes_df$cell_line <- gsub("_.*", "", slopes_df$condition)
slopes_df$dox_conc <- sub(".*_", "", slopes_df$condition)

# Used for further CORE/REGP comparison, in mmol/cell/hr, converted later ######
write_csv(slopes_df,"data/metabolomics/output/REGP_exchange_fluxes.csv")



