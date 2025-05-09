## Script calculation of  metabolite concentrations
# Written by: Cyriel Huijer
# Last updated: 31 March 2025

# Load in packages:
library(tidyverse)
library(ggplot2)

# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)


# Load in data:
data <- read.csv(file = "data/metabolomics/feed_into_R.csv")
media_composition <- read.csv(file="data/metabolomics/DMEM_F12_media_composition_mM.csv")
media_composition$AA_concentration_uM <- media_composition$AA_concentration_mM * 1000

# Select only media measurements (denoted by m0)

m0_means <- data %>%
  filter(sample == "M0") %>%
  group_by(batch, amino_acid) %>%
  summarise(mean_m0 = mean(concentration_sample_uM, na.rm = TRUE), .groups = "drop")

# Add media blank media concentration to normalize (measurements are normalized to measured media concentration)

data_with_m0_means <- data %>%
  left_join(m0_means, by = c("batch", "amino_acid"))

data_with_ratios <- data_with_m0_means %>%
  mutate(ratio_to_m0 = if_else(sample != "M0",
                               concentration_sample_uM / mean_m0,
                               NA_real_))  # Leave NA for "M0" rows themselves

# The lactate values are also normalized, however, since lactate is not present in blank medium I normalize for the glucose concentration

glucose_m0_means <- data %>%
  filter(sample == "M0", amino_acid == "Glc") %>%
  group_by(batch) %>%
  summarise(mean_glc_m0 = mean(concentration_sample_uM, na.rm = TRUE), .groups = "drop")
reference_glc_m0 <- 3.151
glucose_m0_means <- glucose_m0_means %>%
  mutate(correction_factor = reference_glc_m0 / mean_glc_m0)

data_with_ratios <- data_with_ratios %>%
  left_join(media_composition, by = "amino_acid")

final_data <- data_with_ratios %>%
  left_join(glucose_m0_means, by = "batch") %>%
  mutate(final_concentration = if_else(
    amino_acid == "Lac",  # Lactate not in media, so use glucose correction factor
    concentration_sample_uM * correction_factor,  # Apply correction
    ratio_to_m0 * AA_concentration_uM  # Normal calculation for all other amino acids
  ))


# Correct glucose & lactate concentration:

final_data <- final_data %>%
  mutate(
    final_concentration_uM = case_when(
      amino_acid == "Glc" ~ (final_concentration * 1000) / 180.16,  # Glucose conversion (mg/L to uM)
      amino_acid == "Lac" ~ (final_concentration * 1e6) / 90.08,   # Lactate conversion (g/L to uM)
      TRUE ~ final_concentration  # Keep original if not Glc or Lac
    )
  )

# Filter out rows for media sample (media sample is not necessary anymore for REGP calculation)
final_data_filtered <- final_data %>%
  filter(sample != "M0")

# For CORE, the media sample needs to be included:
write.csv(final_data, "data/metabolomics/calculated_metabolite_concentrations_M0_included.csv", row.names = FALSE)

# For REGP
write.csv(final_data_filtered, "data/metabolomics/calculated_metabolite_concentrations.csv", row.names = FALSE)


