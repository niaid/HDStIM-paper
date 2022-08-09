#!/usr/bin/env Rscript --vanilla

# PURPOSE: To rank phospho markers on how well they separate k-means clusters in HDStIM.

library(HDStIM)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(DescTools)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Figures path.
# path <- file.path(figures_folder, "resp_noresp_marker_ranking")

# Selected data. 
# selected_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))

# Run Boruta (it's time taking). This will also generate importance score plots.
# att_stats <- marker_ranking_boruta(selected_data, path = path, n_cells = 5000, max_runs = 100, seed_val = 123, verbose = 1)
# write_tsv(att_stats, file.path(results_folder, "boruta_att_stats.tsv"))

# Read Boruta run results and plot heatmaps with importance scores across cell-population and state markers.
# For each stimulation.
att_stats <- read_tsv( file.path(results_folder, "boruta_att_stats.tsv"))

# Typo correction
att_stats$cell_population <- stringr::str_replace(att_stats$cell_population, "Niave cytotoxic T cells", "Naive cytotoxic T cells")
att_stats$cell_population <- stringr::str_replace(att_stats$cell_population, "Niave T helper cells", "Naive T helper cells")

stim_a <- att_stats %>% dplyr::filter(stim_type == "A")
stim_t <- att_stats %>% dplyr::filter(stim_type == "T")
stim_l <- att_stats %>% dplyr::filter(stim_type == "L")
stim_g <- att_stats %>% dplyr::filter(stim_type == "G")

stim_a <- stim_a[order(stim_a$medianImp),] %>% group_by(cell_population) %>% mutate(imp_idx = dense_rank(desc(medianImp)))
stim_t <- stim_t[order(stim_t$medianImp),] %>% group_by(cell_population) %>% mutate(imp_idx = dense_rank(desc(medianImp)))
stim_l <- stim_l[order(stim_l$medianImp),] %>% group_by(cell_population) %>% mutate(imp_idx = dense_rank(desc(medianImp)))
stim_g <- stim_g[order(stim_g$medianImp),] %>% group_by(cell_population) %>% mutate(imp_idx = dense_rank(desc(medianImp)))

stim_a_in <- xtabs(imp_idx ~ state_marker + cell_population, stim_a)
stim_t_in <- xtabs(imp_idx ~ state_marker + cell_population, stim_t)
stim_l_in <- xtabs(imp_idx ~ state_marker + cell_population, stim_l)
stim_g_in <- xtabs(imp_idx ~ state_marker + cell_population, stim_g)

col_fun = colorRamp2(c(0, 10), c("yellow", "black"))

stim_a_ht <- Heatmap(as.matrix.xtabs(stim_a_in),
     column_title = "A",
        heatmap_legend_param = list(title = "Marker\nRank"),
          col = col_fun,
         show_heatmap_legend = FALSE,
          show_column_dend = FALSE,
     column_names_gp = gpar(fontsize = 13))

stim_t_ht <- Heatmap(as.matrix.xtabs(stim_t_in),
     column_title = "T",
        heatmap_legend_param = list(title = "Marker\nRank"),
          col = col_fun, show_heatmap_legend = FALSE,
     show_column_dend = FALSE,
     column_names_gp = gpar(fontsize = 13))

stim_l_ht <- Heatmap(as.matrix.xtabs(stim_l_in),
     column_title = "L",
        heatmap_legend_param = list(title = "Marker\nRank"),
          col = col_fun, show_heatmap_legend = FALSE,
     show_column_dend = FALSE,
     column_names_gp = gpar(fontsize = 13))

stim_g_ht <- Heatmap(as.matrix.xtabs(stim_g_in),
     column_title = "G",
        heatmap_legend_param = list(title = "Marker\nRank"),
          col = col_fun,
     show_column_dend = FALSE,
     column_names_gp = gpar(fontsize = 13)
     )

ht <- stim_a_ht + stim_t_ht + stim_g_ht + stim_l_ht
png(filename = file.path(figures_folder, "resp_noresp_boruta_all_v2.png"), width = 7, height = 5, units = "in", res = 600)
draw(ht, padding = unit(c(10, 5, 5, 5), "mm"), show_row_dend = FALSE)
dev.off()

