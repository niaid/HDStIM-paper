#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a heamap with the fold chain between stimulated and unstimulated
# wells for IFNa.
# Supplementary figure 1 D.

library(tidyverse)
library(ComplexHeatmap)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))
dat <- selected_data$response_mapping_main
state_markers <- selected_data$state_markers

mean_dat_ifna <- dat %>% dplyr::filter(stim_type == "IFNa") %>%
  dplyr::select(stim_type, cell_population, all_of(state_markers)) %>%
  droplevels() %>%
  group_by(stim_type, cell_population) %>%
  summarise_at(all_of(state_markers), mean)

mean_dat_basal1 <- dat %>% dplyr::filter(stim_type == "Basal1") %>%
  dplyr::select(stim_type, cell_population, all_of(state_markers)) %>%
  droplevels() %>%
  group_by(stim_type, cell_population) %>%
  summarise_at(all_of(state_markers), mean)

merge_dat <- merge(mean_dat_ifna, mean_dat_basal1, by = "cell_population")

cell_pop <- as.character(unique(merge_dat$cell_population))

mat_fc_out <- matrix(nrow = length(state_markers), ncol = length(cell_pop), dimnames = list(state_markers, cell_pop))
for(cp in cell_pop){
  for(state in state_markers){
    st.x <- paste0(state, ".x")
    st.y <- paste0(state, ".y")
    fc <- merge_dat[merge_dat$cell_population == cp, st.x] / merge_dat[merge_dat$cell_population == cp, st.y]
    mat_fc_out[state, cp] <- fc
  }
}

ht <- Heatmap(mat_fc_out, column_dend_side = "bottom", 
              column_dend_height = unit(3, "cm"),
              column_dend_reorder = TRUE,
              column_title = "IFNa",
              heatmap_legend_param = list(title = "Fold Change"))
png(file.path(figures_folder, "pre_scs_stim_unstim_fc_IFNa.png"), width = 5, height = 8, units = "in", res = 300)
draw(ht)
dev.off()
