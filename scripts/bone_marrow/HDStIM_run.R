#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run HDStIM using all the state markers.

library(HDStIM)

results_folder <- file.path("results", "bone_marrow")
figures_folder <- file.path("figures", "bone_marrow")

# Load exported data after manual merging of clusters. . 
dat <- readRDS(file.path(results_folder, "bendall_24_clust.rds"))

# Load panel file.
panel <- read.table(file.path("meta", "paper_stim_panel.txt"), sep = "\t", header = TRUE)

# Load args file.
args_file <- read.table(file.path("meta", "paper_fcs_info.txt"), sep = "\t", header = TRUE)

state_markers <- as.character(panel[panel$marker_class == "state", ][["antigen"]])
cluster_col <- "cluster"
all_unstims <- c("Basal1", "Basal2", "Basal3", "Basal4", "Basal5")
all_stims <- unique(as.character(args_file$stim_type))
stim_lab <- setdiff(all_stims, all_unstims)[1:13]
unstim_lab <- "Basal1"

# Select only donor1.
donor_1 <- dat[dat$patient_id == "marrow1",]
donor_1$patient_id <- droplevels(donor_1$patient_id)

selected_data <-   HDStIM(donor_1, state_markers, cluster_col, stim_lab, unstim_lab,
       seed_val = 123, umap = TRUE, umap_cells = 10000, verbose = TRUE)

# saveRDS(selected_data, file.path(results_folder, "selected_data_all_k_clust.rds"))


