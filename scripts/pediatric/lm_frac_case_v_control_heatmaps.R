#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate heatmaps for case vs control linear regressions.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# CASE vs CONTROL
# STIMULATED
# Load fraction lm data based on ranks. 
frac_lm <- read_tsv(file.path(results_folder, "resp_stim_frac_lm_rank_case_v_control.tsv"))
frac_age <- frac_lm %>% dplyr::filter(term == "groupControl") %>%
  dplyr::mutate("p.value.log" = -log10(p.value))

c_names <- as.character(unique(frac_age$cell_population))
r_names <- as.character(unique(frac_age$stim_type))
mat <- matrix(nrow = length(r_names), ncol = length(c_names) )
rownames(mat) <- r_names
colnames(mat) <- c_names
for(i in 1:nrow(frac_age)){
  st <- frac_age[[i, "stim_type"]]
  cp <- frac_age[[i, "cell_population"]]
  mat[st,cp] <- frac_age[[i, "p.value.log"]]
}

file_n <- file.path(figures_folder, "resp_stim_frac_lm_case_vs_control_heatmap.png")
mat <- t(mat)
col_fun = colorRamp2(c(0, 1, 3), c("blue", "white", "red"))
col_fun(seq(0, 3))
png(filename = file_n, width = 5, height = 6, units = "in", res = 300)
hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
                                heatmap_legend_param = list(title = "-log10(p.value)"),
                                cell_fun = function(j, i, x, y, width, height, fill){
                                  if(!is.na(mat[i, j])){
                                    if(mat[i, j] > -log10(0.05)) {
                                      grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                                    }
                                  } 
                                },
                                column_names_gp = gpar(fontsize = 12),
                                row_names_gp = gpar(fontsize = 12),
                                col = col_fun
)
draw(hmap)
dev.off()

# BASELINE
# Load fraction lm data based on ranks. 
frac_lm <- read_tsv(file.path(results_folder, "resp_unstim_frac_lm_rank_case_v_control.tsv"))
frac_age <- frac_lm %>% dplyr::filter(term == "groupControl") %>%
  dplyr::mutate("p.value.log" = -log10(p.value))

c_names <- as.character(unique(frac_age$cell_population))
r_names <- as.character(unique(frac_age$stim_type))
mat <- matrix(nrow = length(r_names), ncol = length(c_names) )
rownames(mat) <- r_names
colnames(mat) <- c_names
for(i in 1:nrow(frac_age)){
  st <- frac_age[[i, "stim_type"]]
  cp <- frac_age[[i, "cell_population"]]
  mat[st,cp] <- frac_age[[i, "p.value.log"]]
}


file_n <- file.path(figures_folder, "resp_unstim_frac_lm_case_vs_control_heatmap.png")
mat <- t(mat)
col_fun = colorRamp2(c(0, 1, 3), c("blue", "white", "red"))
col_fun(seq(0, 3))
png(filename = file_n, width = 5, height = 6, units = "in", res = 300)
hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
                                heatmap_legend_param = list(title = "-log10(p.value)"),
                                cell_fun = function(j, i, x, y, width, height, fill){
                                  if(!is.na(mat[i, j])){
                                    if(mat[i, j] > -log10(0.05)) {
                                      grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                                    }
                                  } 
                                },
                                column_names_gp = gpar(fontsize = 12),
                                row_names_gp = gpar(fontsize = 12),
                                col = col_fun
)
draw(hmap)
dev.off()

