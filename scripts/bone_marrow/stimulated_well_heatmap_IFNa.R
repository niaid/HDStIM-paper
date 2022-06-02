#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a heamap with and without the z-score calculated on
# median marker intensities in the simulated wells for IFNa.
# Supplementary figure 1 D

library(tidyverse)
library(ComplexHeatmap)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))
dat <- selected_data$response_mapping_main
state_markers <- selected_data$state_markers

plot_dat <- dat %>% dplyr::filter(stim_type == "IFNa") %>%
  dplyr::select(stim_type, cell_population, all_of(state_markers)) %>%
  droplevels() %>%
  group_by(stim_type, cell_population) %>%
  summarise_at(all_of(state_markers), mean)

mat <- matrix(nrow = length(state_markers), ncol = nrow(plot_dat), dimnames = list(state_markers, plot_dat$cell_population))
for(i in 1:nrow(plot_dat)) {
  cell_pop <- as.character(plot_dat[[i, "cell_population"]])
  for(state in state_markers){
    mat[state, cell_pop] <- plot_dat[[i, state]]
  }
}

ht <- Heatmap(mat, column_dend_side = "bottom", 
              column_dend_height = unit(3, "cm"),
              column_dend_reorder = TRUE,
              heatmap_legend_param = list(title = "Mean Exp."))
# png(file.path(figures_folder, "pre_scs_stim_mean_IFNa.png"), width = 8, height = 10, units = "in", res = 300)
# draw(ht)
# dev.off()

# Z score plot.
df_z_out <- matrix(nrow = length(state_markers), ncol =0)
for(coln in colnames(mat)){
  z_scores <- (mat[,coln ] - mean(mat[,coln], na.rm = TRUE)) / sd(mat[,coln], na.rm = TRUE)
  z_scores <- data.frame(z_scores)
  colnames(z_scores) <- coln
  df_z_out <- cbind(df_z_out, z_scores)
}

ht <- Heatmap(as.matrix(df_z_out), column_dend_side = "bottom", 
              column_dend_height = unit(3, "cm"),
              column_dend_reorder = TRUE,
              column_title = "IFNa",
              heatmap_legend_param = list(title = "Z-Score"))
png(file.path(figures_folder, "pre_scs_stim_z_IFNa.png"), width = 5, height = 8, units = "in", res = 300)
draw(ht)
dev.off()
