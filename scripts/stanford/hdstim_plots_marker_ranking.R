#!/usr/bin/R

# Purpose: To run HDStIM on CyTOF data.

library(tidyverse)
library(HDStIM)
library(ComplexHeatmap)

results_folder <- file.path("results", "stanford")
figures_folder <- file.path("figures", "stanford")
dir.create(figures_folder, recursive = TRUE)


# Load data.
cytof_mapped_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))

# Generate diagnostic plots
#cytof_k_plots <- plot_K_Fisher(cytof_mapped_data, path = file.path(figures_folder, "cytof-k-Fisher-plots") , verbose = TRUE)

#cytof_u_plots <- plot_umap(cytof_mapped_data, path = file.path(figures_folder, "cytof-u-plots"), verbose = TRUE)

#cytof_e_plots <- plot_exprs(cytof_mapped_data, path = file.path(figures_folder, "cytof-exprs-plots"), verbose = TRUE)

stop()
# Run marker ranking
ranking <- marker_ranking_boruta(cytof_mapped_data,
  path = file.path(figures_folder, "marker-ranking"),
  n_cells = 50000,
  max_runs = 100,
  seed_val = 123,
  verbose = 1
)
saveRDS(ranking, file.path(results_folder, "marker-ranking.rds"))

# Plot marker ranking heatmap
cytof_m_ranks <- readRDS(file.path(results_folder, "marker-ranking.rds"))

cytof_mat <- cytof_m_ranks$attribute_stats %>%
  dplyr::select(stim_type, cell_population, state_marker, meanImp) %>%
  dplyr::group_by(stim_type, cell_population) %>%
  mutate(min_max = (meanImp - min(meanImp)) /(max(meanImp)-min(meanImp))) %>%
  dplyr::select(stim_type, cell_population, state_marker, min_max) %>%
  dplyr::mutate("stim_pop" = paste0(stim_type, "_", cell_population)) %>%
  dplyr::ungroup()

mat <- matrix(nrow = length(unique(cytof_mat$stim_pop)), ncol = length(unique(cytof_mat$state_marker)))
colnames(mat) <- unique(cytof_mat$state_marker)
rownames(mat) <- unique(cytof_mat$stim_pop)
for(i in 1:nrow(cytof_mat)){
        rowi <- cytof_mat$stim_pop[i]
        coli <- as.character(cytof_mat$state_marker[i])
        val <- cytof_mat$min_max[i]
        mat[rowi, coli] <- val 
}

library(circlize)
col_fun = colorRamp2(c(0, 1), c("black", "yellow"))
col_fun(seq(0, 1))

r_breaks <-  as.character(data.frame(str_split(rownames(mat), "_"))[1,])

hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = FALSE, cluster_columns = TRUE, 
                                col = col_fun, row_split = r_breaks,
                                row_title = NULL,
                                heatmap_legend_param = list(title = "Mean\nImportance"),
                                column_names_gp = gpar(fontsize = 10),
                                row_names_gp = gpar(fontsize = 10))
png(filename = file.path(figures_folder, "marker_ranking_heatmap.png") , width = 6, 
    height = 14, units = "in", pointsize = 12, bg = "white", res = 600)
draw(hmap)
dev.off()

