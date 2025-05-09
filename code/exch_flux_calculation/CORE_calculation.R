## Script calculation of CORE exchange fluxes:
# Written by: Cyriel Huijer
# Last updated: 31 March 2025

# Load packages:
library(tidyverse)
# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)

# Load data:
auc_data <- read.csv("data/cell_count/auc_cell_count_curves_plate_correction.csv")
metabolomics_data <- read.csv("data/metabolomics/calculated_metabolite_concentrations_M0_included.csv")
manual_inspection <- read.csv("data/metabolomics/manual_inspection_curves.csv") # Manual inspection was done to determine low quality samples

# Convert g/L to uM for glucose & lactate (to compare against measured media concentrations raw measurement is used)
metabolomics_data <- metabolomics_data %>%
  mutate(
    concentration_sample_uM = case_when(
      amino_acid == "Glc" ~ (concentration_sample_uM * 1e6) / 180.16,  # Glucose conversion
      amino_acid == "Lac" ~ (concentration_sample_uM * 1e6) / 90.08,   # Lactate conversion
      TRUE ~ concentration_sample_uM  # Keep original if not Glc or Lac
    )
  )

# Calculate concentration in mM
metabolomics_data$concentration_sample_mM <- metabolomics_data$concentration_sample_uM / 1000

# Subset dataframe with only sample measurements
sample_df <- metabolomics_data %>%
  filter(sample!= 'M0')
# Filter out 
manual_selection_df <- merge(sample_df, manual_inspection, by = c("batch", "amino_acid"))
manual_selection_filtered_df <- subset(manual_selection_df, keep_manual_inspection != "n")


df_m0_only <- metabolomics_data %>%
  filter(sample == 'M0') %>%
  dplyr::select(amino_acid, batch, concentration_sample_mM)
mean_concentration_m0 <- df_m0_only %>%
  group_by(amino_acid, batch) %>%
  summarize(mean_concentration_mM = mean(concentration_sample_mM, na.rm = TRUE))

merged_df <- manual_selection_filtered_df %>%
  left_join(mean_concentration_m0, by = c("amino_acid", "batch"))
merged_df <- merged_df %>%
  mutate(substraction = 0.0015 * (concentration_sample_mM - mean_concentration_mM))
# Merge merged_df with auc_data based on metab_tp and timepoint
final_df <- merged_df %>%
  left_join(auc_data, by = c("timepoint" = "metab_tp"))

final_df <- final_df %>%
  mutate(CORE_value = substraction / mean_auc)


# Save file for three way comparison with Jain method (Robinson data, our data)
write.csv(final_df, file = "data/metabolomics/output/CORE_exchange_fluxes.csv",row.names=F)




