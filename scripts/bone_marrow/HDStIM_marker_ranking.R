#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run marker ranking function from HDStIM.

library(HDStIM)
library(tidyverse)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

# Load results from all marker HDStIM.
selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))

br_res_path <- file.path(figures_folder, "boruta_plots")
attribute_stats <- marker_ranking_boruta(selected_data, path = br_res_path,
                                         max_runs = 1000, seed_val = 123,
                                         verbose = 2)

# write_tsv(attribute_stats, file.path(results_folder, "K_clust_boruta_stats.tsv"))

