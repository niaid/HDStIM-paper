
library(tidyverse)

results_folder <- file.path("results", "stanford")
figures_folder <- file.path("figures", "stanford")

# Load HDStIM data
hdstim <- readRDS(file.path(results_folder,"selected_data_all.rds"))
hdstim <- hdstim$all_fisher_p_val

t_test <- read_tsv(file.path(results_folder, "marker-by-marker-t-tests.tsv"))

t_res <- t_test %>% dplyr::group_by(cell_population, stim_type) %>%
          summarize(t_count = sum(p.value < 0.05, na.rm = TRUE)) %>%
          ungroup()

hdstim_res <- hdstim %>% dplyr::select(cell_population, stim_type, p.value) %>%
          mutate("f_significance" = ifelse(p.value < 0.05, 1, 0)) %>%
          mutate("f_log10_pval" = -log10(p.value)) %>%
          as_tibble()

plot_dat <- inner_join(t_res, hdstim_res)

max_m <- length(unique(t_test$marker))

mfp <- ggplot(plot_dat, aes(x = t_count, y = f_log10_pval)) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  labs(x = "Individual Marker Analysis\nMarkers Significant Per SPC", 
       y = "HD Analysis\nFisher's -log10(p-value)") +
  scale_x_continuous(breaks = seq(from = 0, to = max_m, by = 1 )) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 14))
ggsave("t_test_vs_hdstim.png", path = figures_folder, plot = mfp,
       width = 7, height = 4, units = "in", dpi = 600)


one_to_3 <- dplyr::filter(t_res, t_count >= 1 & t_count <=3) %>%
        inner_join(hdstim_res) %>%
        dplyr::filter(f_significance == 1)

zero <- dplyr::filter(t_res, t_count ==0) %>%
        inner_join(hdstim_res) %>%
        dplyr::filter(f_significance == 1)
