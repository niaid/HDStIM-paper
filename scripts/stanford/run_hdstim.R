# Purpose: read the concatenated feather file and fix stimtype or cell population inconsistencies

library(arrow)
library(tidyverse)
library(HDStIM)

results_folder <- file.path("results", "stanford")


dat <- arrow::read_feather(file.path(results_folder, "arcsinh-data-stim-pop-fixed-for-hdstim.feather")) %>%
        dplyr::filter(!stim_type == "IL")


# Fetch state markers.
state_markers <- dat %>% dplyr::select(where(is.numeric)) %>%
        colnames()

# Define column name with annotated cell populations.
cluster_col <- "cell_population"

# Define stimuli used in the assay.
stim_lab <- setdiff(unique(dat$stim_type), "US")

# Define sample lable for unstimulated samples.
unstim_lab <- "US"


selected_data <-   HDStIM(dat, state_markers, cluster_col, stim_lab, unstim_lab,
       seed_val = 123, umap = TRUE, umap_cells = 10000,
       verbose = TRUE)

saveRDS(selected_data, file.path(results_folder, "selected_data_all.rds"))


