## Script calculation of CORE exchange fluxes:
# Written by: Cyriel Huijer
# Last updated: 1  April2025

# Load packages:
library(tidyverse)
# Set working directory:
root <- "" # This should be the base directory of the github repo
setwd(root)

# MCF10A EV0 data:
mcf10a_core <- read.csv(file = "data/metabolomics/output/CORE_exchange_fluxes.csv")
mcf10a_regp <- read.csv(file = "data/metabolomics/output/REGP_exchange_fluxes.csv")

# Load in NCI-60 cell line data:
nci60_core <- read.csv(file = "data/exchange_fluxes/exchangeFluxes_NCI60_formatted.txt", sep = "\t")

# Format the dataframes:

mcf10a_core$condition <- paste0(mcf10a_core$cell_line, "_", mcf10a_core$dox_conc)
mcf10a_core <- mcf10a_core %>%
  filter(timepoint == 24)
mcf10a_core <- mcf10a_core %>% 
  dplyr::select(condition,amino_acid,replicate,CORE_value)
mcf10a_regp <- mcf10a_regp %>%
  dplyr::select(condition,amino_acid,replicate,slope)
colnames(mcf10a_core)[4] <- "mmol_cell_hr"
colnames(mcf10a_regp)[4] <- "mmol_cell_hr"

# Select only emtpy vector 0 condition
mcf10a_core <- mcf10a_core %>%
  filter(condition == "E_0")
mcf10a_regp <- mcf10a_regp %>%
  filter(condition == "E_0")

# Add column for method:
mcf10a_core$method <- "CORE"
mcf10a_regp$method <- "REGP"

# Output files for figure 1C:
#write.csv(mcf10a_regp, file ="data/exchange_fluxes/MCF10A_REGP.csv",row.names = F)
#write.csv(mcf10a_core, file ="data/exchange_fluxes/MCF10A_CORE.csv",row.names = F)


# Combine datafiles: 
mcf10a_combined <- rbind(mcf10a_core, mcf10a_regp)
colnames(mcf10a_combined)[2] <- "metabolite"
# Convert data to mmol/gDW/hr:
mcf10a_combined$mmol_gDW_hr <- mcf10a_combined$mmol_cell_hr * (1 / (397 * 10^-12)) # Dry weight of 397 pg per cell
mcf10a_combined <- mcf10a_combined %>%
  dplyr::select(metabolite,condition,mmol_gDW_hr,method)


# Format Robinson data:
nci60_core <- nci60_core %>%
  pivot_longer(cols = HT29:SR,         
               names_to = "condition", 
               values_to = "value")
nci60_core <- nci60_core %>%
  dplyr::select(metabolite,condition,value)
colnames(nci60_core) <- c("metabolite","condition","mmol_gDW_hr")
nci60_core$method <- "NCI-60"

# Select overlapping metabolites:

mets_in_nci60 <- unique(nci60_core$metabolite)
mets_in_mcf10a_exp_gp <- unique(mcf10a_regp$amino_acid)
mets_in_mcf10a_core <- unique(mcf10a_core$amino_acid)

overlap_mcf10a <- intersect(mets_in_mcf10a_exp_gp, mets_in_mcf10a_core)
overlap_datasets <- intersect(overlap_mcf10a,mets_in_nci60)

nci60_core <- nci60_core %>%
  filter(metabolite %in% overlap_datasets)
mcf10a_combined <- mcf10a_combined %>%
  filter(metabolite %in% overlap_datasets)

# Combine dataframes:
all_data <- rbind(nci60_core,mcf10a_combined)
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

# Calculate mean exchange flux per metabolite, this is used for metabolic modeling

all_data_means <-all_data %>%
  group_by(metabolite,method) %>%
  summarise(
    mean = mean(mmol_gDW_hr, na.rm = TRUE),
    sd = sd(mmol_gDW_hr, na.rm = TRUE)
  )

mean_bounds <- all_data_means %>%
  ungroup() %>%
  mutate(column_name = method) %>%
  dplyr::select(metabolite, column_name,mean) %>%
  pivot_wider(names_from = column_name, values_from = mean)

# Format so from robinson all cell lines are included (before the mean of all cell lines calculated, but we want the exch. flux for every cell line):
nci60_core <- read.csv(file = "data/exchange_fluxes/exchangeFluxes_NCI60_formatted.txt", sep = "\t")
mean_bounds <- merge(mean_bounds, nci60_core, by.x = "metabolite", by.y = "metabolite")
mean_bounds <- mean_bounds %>%
  dplyr::select(-c('NCI-60',Rxn.ID))

# change names for metabolite (this matches HumanGEM metabolite identifiers):
mean_bounds$metabolite <- c('arginine','asparagine','glucose','glutamine','glutamate','glycine','L-lactate','methionine','phenylalanine','proline','serine','threonine','tryptophan','tyrosine','valine')

# Save results, this data is used for running the GEMs. 
write_csv(mean_bounds,'data/exchange_fluxes/output/exch_fluxes_MCF10A_NCI60.csv')






