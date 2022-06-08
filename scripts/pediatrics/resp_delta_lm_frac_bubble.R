#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a bubble plot for linear models on fractions of 
# responding cells at the baseline and delta. 

library(tidyverse)

results_folder <- file.path("results", "pediatrics")
figures_folder <- file.path("figures", "pediatrics")

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

# Bubble plot with estimate is denoted by the size of the bubble
# and p.value by the color.
comb_dat$stim_type <- factor(comb_dat$stim_type, levels = c("U_A", "U_T", "U_G", "U_L",
                                                            "A", "T", "G", "L"))
bplt <- ggplot(comb_dat, aes(x=stim_type, y=cell_population, size=estimate, color=p.value)) +
  geom_point() +
  labs(x = "Experimental Condition", y = "Cell Population", size = "Estimate", color = "P-Value") +
  scale_colour_gradient(low = "red", high = "white", limits=c(0, 0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 12))
ggsave("resp_unstim_delta_lm_rank_bubble_plot.png", plot = bplt, path = figures_folder,
       width = 6, height = 4, dpi = 300)


