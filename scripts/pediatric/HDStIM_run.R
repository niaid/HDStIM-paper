#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run HDStIM using all the state markers on the pediatric dataset
# NOTE: Annotation for cell populations 5 and 7 were changed at this step.

library(HDStIM)
library(tidyverse)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load exported data from CyTOF workflow after manual merging of clusters. 
dat <- readRDS(file.path(results_folder, "dat_all_meta.rds"))
dat <- as_tibble(dat)

# Load panel file.
panel <- read.table(file.path("meta", "pediatric_stim_panel.txt"), sep = "\t", header = TRUE)

# Load args file.
args_file <- read.table(file.path("meta", "pediatric_fcs_info.txt"), sep = "\t", header = TRUE)

# Fetch state markers.
state_markers <- as.character(panel[panel$marker_class == "state", ][["antigen"]])

# Define column name with annotated cell populations.
cluster_col <- "merging1"

# Define stimuli used in the assay.
stim_lab <- c("A", "G", "L", "T")

# Define sample lable for unstimulated samples.
unstim_lab <- "U"

# Correct population 5 and 7. 
dat$merging1 <- stringr::str_replace(dat$merging1, "Naive total T cells", "45RA+ double negative T cells")
dat$merging1 <- stringr::str_replace(dat$merging1, "Total T cells", "45RA- double negative T cells")

selected_data <-   HDStIM(dat, state_markers, cluster_col, stim_lab, unstim_lab,
       seed_val = 123, umap = TRUE, umap_cells = 10000,
       verbose = TRUE)

# saveRDS(selected_data, file.path(results_folder, "selected_data_all.rds"))


