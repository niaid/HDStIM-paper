#!/usr/bin/env Rscript --vanilla

# PURPOSE: To plot HDStIM diagnostic plots.
# Figure 1 C & D
# Note: These figures were generated using the development version of the 
# package. Figures folder only contain figures that were used in the paper.  

library(tidyverse)
library(HDStIM)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

# Mapped data. 
#selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))

# Stacked bar plots.
# pkf <- plot_K_Fisher(selected_data, path = file.path(figures_folder, "sbp_plots_all"), verbose = TRUE)

# UMAPs
# pu <- plot_umap(selected_data, path = file.path(figures_folder, "umap_plots_all"), verbose = TRUE)

# KDE plots.
#pe <- plot_exprs(selected_data, path = file.path(figures_folder, "exprs_plots_all"), verbose = TRUE)


