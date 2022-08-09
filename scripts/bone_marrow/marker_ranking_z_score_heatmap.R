#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a heamap with the z-scores calculated on
# median importance score calculated by Boruta. 
# Stim - Cluster on the columns and state markers on the rows.
# Figure 1 F

library(tidyverse)
library(ComplexHeatmap)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

# Read panel.
panel <- read_tsv(file.path("meta", "paper_stim_panel.txt"))
state_markers <- panel[panel$marker_class == "state", ][["antigen"]]

# Read args file.
args_file <- read_tsv(file.path("meta", "paper_fcs_info.txt"))
all_stims <- unique(as.character(args_file$stim_type))
all_stims <- all_stims[6:18]
unstim <- "Basal1"

clusters <- unique(args_file$cluster)

# Generate a matrix with median importance score from Boruta.
df_out_all <- matrix( nrow = 13 * 24 ,ncol = 18)
stim_cluster <- as.character(sapply(all_stims, FUN = function(x) paste0(x, "_", clusters)))
colnames(df_out_all) <- state_markers
rownames(df_out_all) <- stim_cluster

## Load statistics from Boruta run.
b_res_all <- read_tsv(file.path(results_folder, "K_clust_boruta_stats.tsv"))
for(i in 1:nrow(b_res_all)){
  stim <- b_res_all$stim_type[i]
  state <- b_res_all$state_marker[i]
  clust <- b_res_all$cell_population[i]
  med_imp <- b_res_all$medianImp[i]
  
  stm_clust <- paste0(stim, "_", clust)
  df_out_all[stm_clust, state] <- med_imp
}

df_out_all_write <- data.frame(df_out_all)
df_out_all_write <- rownames_to_column(df_out_all_write, "stim_cluster")
# write_tsv(df_out_all_write, file.path(results_folder, "boruta_stats_martix_stim_cluster.tsv"))


# Calculate z-scores on median importance values calculated by Boruta. 
# Load the above generated dataframe from the tsv file.
# b_stat_st_cl <- read_tsv(file.path(results_folder, "boruta_stats_martix_stim_cluster.tsv"))
b_stat_st_cl <- df_out_all_write
cols_ <- colnames(b_stat_st_cl)
b_dat <- b_stat_st_cl[, cols_[2:19]]
b_dat <- as.matrix(b_dat)
rownames(b_dat) <- b_stat_st_cl$stim_cluster 
b_mat <- b_dat[rowSums(is.na(b_dat)) != ncol(b_dat), ] # remove rows with all NAs.

df_z_out <- matrix(nrow = 18, ncol =0)
t_b_mat <- t(b_mat)
for(coln in colnames(t_b_mat)){
  z_scores <- (t_b_mat[,coln ] - mean(t_b_mat[,coln], na.rm = TRUE)) / sd(t_b_mat[,coln], na.rm = TRUE)
  z_scores <- data.frame(z_scores)
  colnames(z_scores) <- coln
  df_z_out <- cbind(df_z_out, z_scores)
}

# Generate the heatmap on the z-score matrix calculated above.
in_mat <- as.matrix(df_z_out)

# Generate the heatmap on the z-score matrix calculated above with per stim partitioning.
in_mat <- as.matrix(df_z_out)
c_breaks <- colnames(in_mat)
c_breaks <- as.character(data.frame(strsplit(c_breaks, "_"))[1,])
colnames(in_mat) <- as.character(data.frame(strsplit(colnames(in_mat), "_"))[2,])
ht <- Heatmap(in_mat, show_column_dend = FALSE, 
              show_row_dend = FALSE,
              column_split = c_breaks,
              heatmap_legend_param = list(title = "Z-Score",fontsize = 20,
                                          legend_height = unit(4, "cm"),
                                          legend_width = unit(3, "cm"),
                                          title_position = "leftcenter-rot"),
              row_names_gp = gpar(fontsize = 16),
              column_names_gp = gpar(fontsize = 15),
              column_title_gp = gpar(fontsize = 20))
#png(file.path(figures_folder, "boruta_z_heatmap_horizontal_stim_breaks_no_dendo.png"), width = 30, height = 8, units = "in", res = 300)
#draw(ht)
#dev.off()

# Plot a heatmap for only INFa.
inf_dat <- df_z_out %>% dplyr::select(starts_with("IFNa"))
inf_mat <- as.matrix(inf_dat)
c_breaks <- colnames(inf_mat)
c_breaks <- as.character(data.frame(strsplit(c_breaks, "_"))[1,])
colnames(inf_mat) <- as.character(data.frame(strsplit(colnames(inf_mat), "_"))[2,])
infht <- Heatmap(inf_mat, show_column_dend = FALSE, 
              show_row_dend = FALSE,
              heatmap_legend_param = list(title = "Z-Score",fontsize = 14,
                                          title_position = "leftcenter-rot"),
              row_names_gp = gpar(fontsize = 14),
              column_names_gp = gpar(fontsize = 14),
              column_title = "IFNa",
              column_title_gp = gpar(fontsize = 16))
png(file.path(figures_folder, "boruta_z_heatmap_IFNa.png"), width = 4.5, height = 5.5, units = "in", res = 600)
draw(infht)
dev.off()

