#!/usr/bin/env Rscript --vanilla

# PURPOSE: Run Boruta on marker magnitude to predict age in responding unstimulated samples.

library(tidyverse)
library(Boruta)
library(ComplexHeatmap)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

stop()
# Run Boruta (time consuming).
# Load per subject marker magnitude data.
f_mag <- read_tsv(file.path(results_folder, "resp_unstim_per_subj_phospho_mag_m_f.tsv")) 
state_markers <- colnames(df_mag)[6:15]

# Figure path.
fig_path <- file.path(figures_folder, "resp_unstim_phospo_mag_boruta")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

df_groupped <- dplyr::group_by(df_mag, cell_population, stim_type)
df_split <- dplyr::group_split(df_groupped)

max_runs <- 100
verbose <- 2
seed_val <- 123
n_cells = 5000
path <- fig_path
df_stats_out <- data.frame()
for(i in 1:length(df_split)){
  dat_boruta <- df_split[[i]] %>%
    dplyr::select(!comb_no)
  stim <- as.character(unique(dat_boruta$stim_type))
  clust <- as.character(unique(dat_boruta$cell_population))
  
  form_boruta <- as.formula(paste0("age_at_sample ~ gender + ",paste0(state_markers, collapse = " + ")))
  
  # Run Boruta.
  set.seed(seed_val)
  res_boruta <- Boruta(form_boruta, data = dat_boruta, doTrace = verbose, maxRuns = max_runs)
  att_stats <- attStats(res_boruta)
    
  # Estimate variable importance and generate output data frame.
  df_imp <- tibble::rownames_to_column(att_stats, "state_marker")
  df_imp <- df_imp[order(df_imp$medianImp),]
  df_imp$state_marker <- factor(df_imp$state_marker, levels = df_imp$state_marker)
  df_imp$decision <- as.character(df_imp$decision)
  df_imp <- cbind("stim_type" = stim, "cell_population" = clust, df_imp)
  df_stats_out <- rbind(df_stats_out, df_imp)
    
  # Plot
  plot_title <- paste0("Stimulation: ", stim,"\nCell Population: ", clust)
  my_cols <- c("Confirmed" = "darkgreen", "Tentative" = "darkblue", "Rejected" = "darkred")
  
  att_plot <- ggplot(df_imp, aes(x = state_marker, y = medianImp, col = decision)) +
    geom_point() +
    geom_errorbar(aes(ymin=minImp, ymax=maxImp), width = 0.2) +
    labs(x = "Predictor", y = "Importance", title = plot_title, col = "Decision") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_colour_manual(values = my_cols)
  
    att_file <- paste0("imp_", stim, "_", clust, ".png")
    ggsave(att_file, plot = att_plot, path = path,
           device = "png", dpi = 300, width = 7, height = 5, units = "in")
}

# write_tsv(df_stats_out, file.path(results_folder, "resp_unstim_mag_boruta.tsv"))

stop()
# Generate heatmap.
df_stats_out <- read_tsv(file.path(results_folder, "resp_unstim_mag_boruta.tsv"))
df_stats_out <- df_stats_out[order(df_stats_out$stim_type),]

# Typo correction
df_stats_out$cell_population <- stringr::str_replace(df_stats_out$cell_population, "Niave cytotoxic T cells", "Naive cytotoxic T cells")
df_stats_out$cell_population <- stringr::str_replace(df_stats_out$cell_population, "Niave T helper cells", "Naive T helper cells")
df_stats_out$state_marker <- stringr::str_replace(df_stats_out$state_marker, "pERK1_2", "pERK1/2")


c_names <- as.character(unique(paste0(df_stats_out$cell_population," ~ ", df_stats_out$stim_type)))
r_names <- as.character(unique(df_stats_out$state_marker))
mat <- matrix(nrow = length(r_names), ncol = length(c_names) )
rownames(mat) <- r_names
colnames(mat) <- c_names
df_stats_out <- df_stats_out %>% mutate(cell_pop_st = paste0(df_stats_out$cell_population," ~ ", df_stats_out$stim_type))
for(i in 1:nrow(df_stats_out)){
  st <- df_stats_out[[i, "state_marker"]]
  cp <- df_stats_out[[i, "cell_pop_st"]]
  mat[st,cp] <- df_stats_out[[i, "meanImp"]]
}

mat_decision <- matrix(nrow = length(r_names), ncol = length(c_names) )
rownames(mat_decision) <- r_names
colnames(mat_decision) <- c_names
for(i in 1:nrow(df_stats_out)){
  st <- df_stats_out[[i, "state_marker"]]
  cp <- df_stats_out[[i, "cell_pop_st"]]
  mat_decision[st,cp] <- df_stats_out[[i, "decision"]]
}

file_n <- file.path(figures_folder, "resp_unstim_mag_ranking_heatmap_v3.png")

c_breaks <- c(rep("U_IFNa", 14), rep("U_IFNg", 4), rep("U_LPS", 2), rep("U_TCR", 8))
png(filename = file_n, width = 9, height = 4.5, units = "in", res = 600)
colnames(mat) <- str_split_fixed(colnames(mat), " ~ ",2)[,1]
hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = TRUE, cluster_columns = TRUE,
                                heatmap_legend_param = list(title = "Mean\nImportance"),
                                column_split = c_breaks,
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if(mat_decision[i, j] == "Confirmed"){
                                    grid.text("C", x, y, gp = gpar(fontsize = 9))
                                  }else if(mat_decision[i, j] == "Tentative"){
                                    grid.text("T", x, y, gp = gpar(fontsize = 9))
                                  }
                                  
                                },
                                column_names_gp = gpar(fontsize = 12),
                                row_names_gp = gpar(fontsize = 11),
                                show_column_dend = FALSE,
                                show_row_dend = FALSE)
draw(hmap)
dev.off()

