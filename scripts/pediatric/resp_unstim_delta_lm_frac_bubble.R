#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a bubble plot for linear models on fractions of 
# responding cells at the baseline and delta. 

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load fraction lm data based on ranks. 
frac_lm_unstim <- read_tsv(file.path(results_folder, "resp_unstim_frac_lm_rank.tsv")) %>% 
  dplyr::filter(term == "age_at_sample")
frac_lm_delta <- read_tsv(file.path(results_folder, "delta_frac_lm_rank.tsv")) %>%
  dplyr::filter(term == "age_at_sample")

comb_dat <- bind_rows(tibble("dataset" = "unstim", frac_lm_unstim),
                      tibble("dataset" = "delta", frac_lm_delta))

# Typo correction
comb_dat$cell_population <- stringr::str_replace(comb_dat$cell_population, "Niave cytotoxic T cells", "Naive cytotoxic T cells")
comb_dat$cell_population <- stringr::str_replace(comb_dat$cell_population, "Niave T helper cells", "Naive T helper cells")
comb_dat$stim_type <- stringr::str_replace_all(comb_dat$stim_type, c("A" = "IFNa", "T" = "TCR", 
                                                                     "G" = "IFNg", "L" = "LPS"))


# Bubble plot with estimate is denoted by the size of the bubble
# and p.value by the color.
comb_dat$stim_type <- factor(comb_dat$stim_type, levels = c("U_IFNa", "U_TCR", 
                                                            "U_IFNg", "U_LPS",
                                                            "IFNa", "TCR", 
                                                            "IFNg", "LPS"))
bplt <- ggplot(comb_dat, aes(x=stim_type, y=cell_population, size=estimate, color=p.value)) +
  geom_point() +
  labs(x = "Experimental Condition", y = "Cell Population", size = "Estimate", color = "P-Value") +
  scale_colour_gradient(low = "red", high = "white", limits=c(0, 0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 12))
  #ggsave("resp_unstim_delta_lm_rank_bubble_plot_v2.png", plot = bplt, path = figures_folder,
   #       width = 6, height = 4, dpi = 300)


# Convert data frame to matrix and plot a heatmap instead of a bubble plot.
  mate <- dplyr::select(comb_dat, cell_population, stim_type, estimate) %>%
    spread(key = stim_type, value = estimate) %>%
    column_to_rownames(var = "cell_population") %>%
    as.matrix()
  mate <- mate[c(5,1, setdiff(1:nrow(mate), c(1,5))),]
  matp <- dplyr::select(comb_dat, cell_population, stim_type, p.value) %>%
    spread(key = stim_type, value = p.value) %>%
    column_to_rownames(var = "cell_population") %>%
    as.matrix()
  matp <- matp[c(5,1, setdiff(1:nrow(matp), c(1,5))),]
  col_fun = colorRamp2(c(0, 1, 3), c("blue", "white", "red"))
  col_fun(seq(0, 3))
  cbreaks <- c(rep("Baseline/Tonic\nActivation", 4), rep("In vitro\nResponse", 4))
  hmap <- ComplexHeatmap::Heatmap(mate, cluster_rows = FALSE, cluster_columns = FALSE,
                                  heatmap_legend_param = list(title = "Effect\nSize"),
                                  cell_fun = function(j, i, x, y, width, height, fill){
                                    if(!is.na(matp[i, j])){
                                      if(matp[i, j] < 0.05) {
                                        grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 14 ))
                                      }
                                    }
                                  },
                                  column_names_gp = gpar(fontsize = 13),
                                  column_title_gp = gpar(fontsize = 16),
                                  row_names_gp = gpar(fontsize = 13),
                                  column_split = cbreaks)
  file_n <- file.path(figures_folder, "resp_unstim_delta_lm_rank_heatmap_v2.png")
  png(filename = file_n, width = 7, height = 5, units = "in", res = 600)
  draw(hmap)
  dev.off()
  