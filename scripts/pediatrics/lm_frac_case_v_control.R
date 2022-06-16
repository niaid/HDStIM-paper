#!/usr/bin/env Rscript --vanilla

# PURPOSE: Build linear models for case v control.

library(tidyverse)

results_folder <- file.path("results", "pediatrics")
figures_folder <- file.path("figures", "pediatrics")

if(!dir.exists(results_folder)){
  loginfo("Creating %s folder", results_folder)
  dir.create(results_folder)
} else {
  loginfo("%s folder already exists. Output will be over written.", results_folder)
}

if(!dir.exists(figures_folder)){
  loginfo("Creating %s folder", figures_folder)
  dir.create(figures_folder)
} else {
  loginfo("%s folder already exists. Output will be over written.", figures_folder)
}

# Load per subject fraction data.
df_frac <- read_tsv(file.path(results_folder, "responding_unstim_per_subj_frac_m_f.tsv")) %>% 
  dplyr::rename("frac" = freq)
sub_dat <- read_tsv(file.path("meta", "mini", "MINI_subject_data.txt"))

# Add condition/group to to fraction data.
df_lm <- left_join(df_frac, dplyr::select(sub_dat, patient_id, group), by = "patient_id")
df_lm_case <- dplyr::filter(df_lm, group == "Case")
df_lm_control <- dplyr::filter(df_lm, group == "Control")
