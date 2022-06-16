#!/usr/bin/env Rscript --vanilla

# PURPOSE: To perform corrections to selected data. This is to avoid correcting
# the initial data and do clustering or other compute intensive processes all over again. 
# Corrections done:
# 1. Add continuous age.
# 2. Correct gender labels.

library(tidyverse)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Selected data. 
selected_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))

# Read pediatric subject data. 
# NOTE: this data was not in the original args file hence was not included
# in FlowSom clustering in the Robinson's pipeline.
subj_dat <- read_tsv(file.path("meta", "pediatric_subject_data.txt"))

# Add continuous age at sample to the FlowSom clustering data.
resp_map <- selected_data$response_mapping_main %>%
  dplyr::select(!gender)
selected_data$response_mapping_main <- left_join(resp_map, dplyr::select(subj_dat, patient_id, gender, age_at_sample))

# Write updated data back to selected_dall_all.rds
#saveRDS(selected_data, file.path(results_folder, "selected_data_all.rds"))
