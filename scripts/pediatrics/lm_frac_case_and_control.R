#!/usr/bin/env Rscript --vanilla

# PURPOSE: Build linear models for responding unstimulated fractions separately for case and control
# samples.

library(tidyverse)

results_folder <- file.path("results", "pediatrics")
figures_folder <- file.path("figures", "pediatrics")

# CALCULATE PER SUBJECT FRACTIONS
# Selected data. 
selected_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))

# Calculate total number of cells per subject for the stim samples.
total_count <- selected_data$response_mapping_main %>%
  group_by(comb_no, gender, cell_population, stim_type, patient_id) %>%
  droplevels() %>%
  count(name = "total_count") 

## Per combination select stim type label.
comb_stim <- ungroup(total_count) %>% 
  dplyr::select(comb_no, stim_type) %>% 
  unique %>% filter(stim_type != "U") %>%
  droplevels() %>%
  dplyr::rename("stim_type_u" = stim_type)

# Calculate number of cells per subject for the resp. unstim. samples.
resp_count <- selected_data$response_mapping_main %>% 
  group_by(comb_no, gender, cell_population, stim_type, patient_id) %>%
  dplyr::filter(response_status == "Resp. Unstim.") %>%
  droplevels() %>%
  count(name = "resp_count")

# Calculate fraction of resp cells.
df_resp_stim <- merge(resp_count, total_count) %>%
  dplyr::mutate("freq" = resp_count / total_count) %>%
  as_tibble()

## Join stim_type and stim_type_u columns.
plot_dat <- left_join(df_resp_stim, comb_stim) %>%
  tidyr::unite("stim_type", stim_type,stim_type_u, sep = "_")

# Add age column and save the plot data for correlation analysis. 
age <- selected_data$response_mapping_main %>%
  dplyr::select(patient_id, age_at_sample) %>%
  unique()

#write_tsv(left_join(plot_dat, age), file.path(results_folder, "responding_unstim_per_subj_frac_m_f.tsv"))

# LINEAR MODELING
# Load per subject fraction data.
df_frac <- read_tsv(file.path(results_folder, "responding_unstim_per_subj_frac_m_f.tsv")) %>% 
  dplyr::rename("frac" = freq)
sub_dat <- read_tsv(file.path("meta", "pediatrics_subject_data.txt"))

# Add condition/group to to fraction data.
df_lm <- left_join(df_frac, dplyr::select(sub_dat, patient_id, group), by = "patient_id")
df_lm_case <- dplyr::filter(df_lm, group == "Case")
df_lm_control <- dplyr::filter(df_lm, group == "Control")

# Build linear models.
## Case
df_groupped <- dplyr::group_by(df_lm_case, cell_population, stim_type)
df_split <- dplyr::group_split(df_groupped)
df_lm_case_out <- data.frame()
for(i in 1:length(df_split)){
  lm_dat <- df_split[[i]] %>%
    dplyr::select(!comb_no) %>%
    dplyr::mutate("frac_rank" = rank(frac))
  st <- as.character(unique(lm_dat$stim_type))
  cp <- as.character(unique(lm_dat$cell_population))
  form <- as.formula(paste0("frac_rank ~ age_at_sample * gender"))
  lm_res <- lm(form, data = lm_dat)
  lm_sum <- summary(lm_res) %>%
    broom::tidy()
  lm_glance <- broom::glance(lm_res)
  colnames(lm_glance) <- paste("mod",colnames(lm_glance), sep="_")
  df_temp <- data.frame("cell_population" = cp, "stim_type" = st, lm_sum, lm_glance)
  df_lm_case_out <- rbind(df_lm_case_out, df_temp)
}
#write_tsv(df_lm_case_out, file.path(results_folder, "resp_unstim_frac_lm_rank_case.tsv"))

## Control
df_groupped <- dplyr::group_by(df_lm_control, cell_population, stim_type)
df_split <- dplyr::group_split(df_groupped)
df_lm_control_out <- data.frame()
for(i in 1:length(df_split)){
  lm_dat <- df_split[[i]] %>%
    dplyr::select(!comb_no) %>%
    dplyr::mutate("frac_rank" = rank(frac))
  st <- as.character(unique(lm_dat$stim_type))
  cp <- as.character(unique(lm_dat$cell_population))
  form <- as.formula(paste0("frac_rank ~ age_at_sample * gender"))
  lm_res <- lm(form, data = lm_dat)
  lm_sum <- summary(lm_res) %>%
    broom::tidy()
  lm_glance <- broom::glance(lm_res)
  colnames(lm_glance) <- paste("mod",colnames(lm_glance), sep="_")
  df_temp <- data.frame("cell_population" = cp, "stim_type" = st, lm_sum, lm_glance)
  df_lm_control_out <- rbind(df_lm_control_out, df_temp)
}
#write_tsv(df_lm_control_out, file.path(results_folder, "resp_unstim_frac_lm_rank_control.tsv"))
