#!/usr/bin/R

# Purpose: To run HDStIM on CyTOF data.

library(tidyverse)
library(HDStIM)
library(arrow)

results_folder <- file.path("results", "pregnancy_baseline_fj")
figures_folder <- file.path("figures", "pregnancy_baseline_fj")
dir.create(figures_folder, recursive = TRUE)


# Load data.
cytof_dat <- arrow::read_feather(file.path(results_folder,"cytof-dat-asinh-transform.feather"))

# Read pannel info
panel <- read_tsv(file.path("meta", "pregnancy_phospho_panel.txt"), show_col_types = FALSE)

stims <- setdiff(unique(cytof_dat$stim_type), "Unstim")
unstim <- "Unstim"
state_marker <- panel %>% dplyr::filter(marker_class == "state") %>%
        dplyr::pull(antigen)
state_marker <- str_replace_all(state_marker, "-", "_")


cytof_dat_filt <- dplyr::filter(cytof_dat, !cell_population %in% c("Mononuclear cells", "Leukocytes"))

cytof_mapped_data <-  HDStIM(cytof_dat_filt, state_marker,
                             "cell_population", stims, unstim, seed_val = 123, 
                             umap = TRUE, umap_cells = 5000, 
                             verbose = TRUE)

saveRDS(cytof_mapped_data, file.path(results_folder, "cytof-mapped-data.rds"))

# Generate diagnostic plots
cytof_k_plots <- plot_K_Fisher(cytof_mapped_data, path = file.path(figures_folder, "cytof-k-Fisher-plots") , verbose = TRUE)

cytof_u_plots <- plot_umap(cytof_mapped_data, path = file.path(figures_folder, "cytof-u-plots"), verbose = TRUE)

cytof_e_plots <- plot_exprs(cytof_mapped_data, path = file.path(figures_folder, "cytof-exprs-plots"), verbose = TRUE)

stop()
# Run marker ranking
ranking <- marker_ranking_boruta(cytof_mapped_data,
  path = figures_folder,
  n_cells = 50000,
  max_runs = 100,
  seed_val = 123,
  verbose = 1
)
saveRDS(ranking, file.path(results_folder, "marker-ranking.rds"))

