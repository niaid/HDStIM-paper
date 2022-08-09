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
selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))

# Stacked bar plots.
# pkf <- plot_K_Fisher(selected_data, path = file.path(figures_folder, "sbp_plots_all"), verbose = TRUE)

# UMAPs
# pu <- plot_umap(selected_data, path = file.path(figures_folder, "umap_plots_all"), verbose = TRUE)

# KDE plots.
#pe <- plot_exprs(selected_data, path = file.path(figures_folder, "exprs_plots_all"), verbose = TRUE)
pe <- plot_exprs(selected_data, path = NULL, verbose = TRUE)

ke_out <- tibble()
for(i in 1:length(pe)){
  cp <- as.character(unique(ggplot_build(pe[[i]])$plot$data$cell_population))
  st <- as.character(unique(ggplot_build(pe[[i]])$plot$data$stim_type))[1]
  temp <- tibble("entry" = i, "cp" = cp, "st" = st)
  ke_out <- bind_rows(ke_out, temp)
}

kplt <- function(ggobj){
  ggplt <- ggobj +
    geom_rect(data = subset(ggplot_build(ggobj$plot$data), state_marker %in% c("pP38", "pSTAT3")), 
              fill = NA, colour = "red", xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf) +
    theme(plot.title = element_text(size = 25)) +
    theme(axis.text = element_text(size = 20)) +
    theme(axis.title.x = element_text(size = 22)) +
    theme(axis.title.y = element_text(size = 22)) +
    theme(legend.text = element_text(size = 22))+
    theme(legend.title = element_text(size = 22))
  
    return(ggplt)
}

pop1 <- dplyr::filter(ke_out, cp == "Immature B" & st == "PVO4") %>%
  dplyr::pull(entry)

pop1_fig <- kplt(pe[[pop1]])
ggsave("exprs_Immature_B_PVO4_v2.png", path = figures_folder, plot = pop1_fig,
       width = 12, height = 10, units = "in", dpi = 600)

