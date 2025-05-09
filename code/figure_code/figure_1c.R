## Figure 1c
# Written by: Cyriel Huijer
# Last updated: 30/03/2025

# Load in packages:
library(tidyverse) # tidyverse_2.0.0

# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)

# Load in data:
mcf10a_core <- read.csv(file = "data/exchange_fluxes/MCF10A_CORE.csv")
mcf10a_regp <- read.csv(file = "data/exchange_fluxes/MCF10A_REGP.csv")
nci60 <- read.csv(file = "data/exchange_fluxes/exchangeFluxes_NCI60_formatted.txt", sep = "\t") 

# Format & combine MCF10A data:

# Combine datafiles: 
mcf10a_combined <- rbind(mcf10a_core, mcf10a_regp)
colnames(mcf10a_combined)[2] <- "metabolite"
# Convert data to mmol/gDW/hr:
mcf10a_combined$mmol_gDW_hr <- mcf10a_combined$mmol_cell_hr * (1 / (397 * 10^-12))
mcf10a_combined <- mcf10a_combined %>%
  dplyr::select(metabolite,condition,mmol_gDW_hr,method)

# Format nci60 data:
nci_60_data <- nci60 %>%
  pivot_longer(cols = HT29:SR,         
               names_to = "condition", 
               values_to = "value")
nci_60_data <- nci_60_data %>%
  dplyr::select(metabolite,condition,value)
colnames(nci_60_data) <- c("metabolite","condition","mmol_gDW_hr")
nci_60_data$method <- "NCI-60"

# Select overlapping metabolites:

mets_in_nci60 <- unique(nci_60_data$metabolite)
mets_in_mcf10a_exp_gp <- unique(mcf10a_regp$amino_acid)
mets_in_mcf10a_core <- unique(mcf10a_core$amino_acid)

overlap_mcf10a <- intersect(mets_in_mcf10a_exp_gp, mets_in_mcf10a_core)
overlap_datasets <- intersect(overlap_mcf10a,mets_in_nci60)

nci_60_data <- nci_60_data %>%
  filter(metabolite %in% overlap_datasets)
mcf10a_combined <- mcf10a_combined %>%
  filter(metabolite %in% overlap_datasets)

# Combine dataframes:
all_data <- rbind(nci_60_data,mcf10a_combined)
# Factorize for plot:
all_data <- all_data %>%
  mutate(method = factor(method, levels = c("NCI-60", "CORE", "REGP")))

# Plot:
ggplot(all_data, aes(x = method, y = mmol_gDW_hr, fill = method)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(), alpha = 0.7) +  # Mean bar
  geom_jitter(aes(color = method), width = 0.2, size = 2, alpha = 0.8) +  # Jittered points
  facet_wrap(~metabolite, scales = "free_y") +  # One plot per metabolite
  theme_minimal() +
  labs(x = "Method", y = "mmol/gDW/hr", title = "Exchange Flux per Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/Fig1c_exch_flux_comparison.pdf", width=10, height=10)



