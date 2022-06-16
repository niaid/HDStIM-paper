#!/usr/bin/env Rscript --vanilla

# PURPOSE: Build linear models for case v control.

library(tidyverse)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load per subject fraction data.
df_frac <- read_tsv(file.path(results_folder, "responding_unstim_per_subj_frac_m_f.tsv")) %>% 
  dplyr::rename("frac" = freq)
sub_dat <- read_tsv(file.path("meta", "pediatric_subject_data.txt"))

# Add condition/group to to fraction data.
df_lm <- left_join(df_frac, dplyr::select(sub_dat, patient_id, group), by = "patient_id")
# df_lm_case <- dplyr::filter(df_lm, group == "Case")
# df_lm_control <- dplyr::filter(df_lm, group == "Control")

df_groupped <- dplyr::group_by(df_lm, cell_population, stim_type)
df_split <- dplyr::group_split(df_groupped)
df_lm_cc_out <- data.frame()
for(i in 1:length(df_split)){
  lm_dat <- df_split[[i]] %>%
    dplyr::select(!comb_no) %>%
    dplyr::mutate("frac_rank" = rank(frac))
  st <- as.character(unique(lm_dat$stim_type))
  cp <- as.character(unique(lm_dat$cell_population))
  form <- as.formula(paste0("frac_rank ~ group"))
  lm_res <- lm(form, data = lm_dat)
  lm_sum <- summary(lm_res) %>%
    broom::tidy()
  lm_glance <- broom::glance(lm_res)
  colnames(lm_glance) <- paste("mod",colnames(lm_glance), sep="_")
  df_temp <- data.frame("cell_population" = cp, "stim_type" = st, lm_sum, lm_glance)
  df_lm_cc_out <- rbind(df_lm_cc_out, df_temp)
}
#write_tsv(df_lm_cc_out, file.path(results_folder, "resp_unstim_frac_lm_rank_case_v_control.tsv"))
