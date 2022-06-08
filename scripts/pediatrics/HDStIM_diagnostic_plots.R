#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate HDStIM diagnostic plots.

library(HDStIM)
library(tidyverse)

results_folder <- file.path("results", "pediatrics")
figures_folder <- file.path("figures", "pediatrics")

# Selected data. 
selected_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))

# Write Fisher's and k-means stats.
# df_f_all <- selected_data$all_fisher_p_val
# write_tsv(df_f_all, file.path(results_folder, "F_stats_all.tsv")) 

# df_k_all <- selected_data$all_k_means_dat
# write_tsv(df_k_all, file.path(results_folder, "K_stats_all.tsv"))

# Stacked bar plots.
pk <- plot_K_Fisher(selected_data, path = file.path(figures_folder, "k_plots_all"), verbose = TRUE)

# UMAPs
pu <- plot_umap(selected_data, path = file.path(figures_folder, "umap_plots_all"), verbose = TRUE)

# KDE plots.
pe <- plot_exprs(selected_data, path = file.path(figures_folder, "exprs_plots_all"), verbose = TRUE)



