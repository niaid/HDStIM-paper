#!/usr/bin/env Rscript --vanilla

# PURPOSE: Build linear models for responding unstimulated fractions.

library(tidyverse)
library(ggfortify)
library(sjPlot)
library(ggpubr)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load per subject fraction data.
df_frac <- read_tsv(file.path(results_folder, "responding_unstim_per_subj_frac_m_f.tsv")) %>% 
  dplyr::rename("frac" = freq)
sub_dat <- read_tsv(file.path("meta", "pediatric_subject_data.txt"))

# Add condition/group to to fraction data.
df_lm <- left_join(df_frac, dplyr::select(sub_dat, patient_id, group), by = "patient_id")

# Figure path.
fig_path <- file.path(figures_folder, "resp_unstim_lm_rank")
# dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

fig_path_marg <- file.path(figures_folder, "resp_unstim_lm_marginal_effects_rank")
# dir.create(fig_path_marg, recursive = TRUE, showWarnings = TRUE)


# Build linear model.
df_groupped <- dplyr::group_by(df_lm, cell_population, stim_type)
df_split <- dplyr::group_split(df_groupped)
df_lm_out <- data.frame()
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
  df_lm_out <- rbind(df_lm_out, df_temp)
  
  # Generate diagnostic plots.
  file_name <- paste0("resp_unstim_",cp,"_",st,".png")
  #   png(filename = file.path(fig_path, file_name), units = "in", width = 10, height = 8, res = 300)
  #   print(autoplot(lm_res))
  #   dev.off()
  
  # Generate marginal effects plots.
  p1 <- sjPlot::plot_model(lm_res, type="pred", term = c("age_at_sample", "gender"), show.data = TRUE)
  p2 <- sjPlot::plot_model(lm_res, type="pred", term = c("gender", "age_at_sample"))
  p <- ggpubr::ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)
  #   file_name <- paste0("resp_unstim_",cp,"_",st,"_marg_eff.png")
  #   ggsave(file_name, plot = p, path = fig_path_marg, 
  #          width = 14, height = 5, device = "png", dpi = 300)
}

# write_tsv(df_lm_out, file.path(results_folder, "resp_unstim_frac_lm_rank.tsv"))


