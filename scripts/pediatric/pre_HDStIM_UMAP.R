#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate UMAPs on pre-HDStIM data and color the cells 
# as case or control.

library(tidyverse)
library(uwot)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load exported data after manual merging of clusters. 
dat <- readRDS(file.path(results_folder, "dat_all_meta.rds"))
dat <- as_tibble(dat)

# Load panel file.
panel <- read.table(file.path("meta", "pediatric_stim_panel.txt"), sep = "\t", header = TRUE)

# Load args file.
args_file <- read.table(file.path("meta", "pediatric_fcs_info.txt"), sep = "\t", header = TRUE)

# Fetch state markers.
state_markers <- as.character(panel[panel$marker_class == "state", ][["antigen"]]) %>%
  stringr::str_replace( "-", "_")

# Define column name with annotated cell populations.
cluster_col <- "merging1"

# Define stimuli used in the assay.
stim_lab <- c("A", "G", "L", "T")

# Define sample lable for unstimulated samples.
unstim_lab <- "U"

# Correct population 5 and 7. 
dat$merging1 <- stringr::str_replace(dat$merging1, "Naive total T cells", "45RA+ double negative T cells")
dat$merging1 <- stringr::str_replace(dat$merging1, "Total T cells", "45RA- double negative T cells")

set.seed(123)
a <- dplyr::filter(dat, stim_type == "A") %>%
  droplevels() %>%
  sample_n(10000)
t <- dplyr::filter(dat, stim_type == "T") %>%
  droplevels() %>%
  sample_n(10000)
l <- dplyr::filter(dat, stim_type == "L") %>%
  droplevels() %>%
  sample_n(10000)
g <- dplyr::filter(dat, stim_type == "G") %>%
  droplevels() %>%
  sample_n(10000)
uu <- dplyr::filter(dat, stim_type == "U") %>%
  droplevels() %>%
  sample_n(10000)

u_a <- uwot::umap(dplyr::select(a, all_of(state_markers)))
u_a <- data.frame(u_a) %>%
  cbind("Groups" = a$condition)
u_t <- uwot::umap(dplyr::select(t, all_of(state_markers)))
u_t <- data.frame(u_t) %>%
  cbind("Groups" = t$condition)
u_l <- uwot::umap(dplyr::select(l, all_of(state_markers)))
u_l <- data.frame(u_l) %>%
  cbind("Groups" = l$condition)
u_g <- uwot::umap(dplyr::select(g, all_of(state_markers)))
u_g <- data.frame(u_g) %>%
  cbind("Groups" = g$condition)
u_uu <- uwot::umap(dplyr::select(uu, all_of(state_markers)))
u_uu <- data.frame(u_uu) %>%
  cbind("Groups" = g$condition)

# Plots
p_a <- ggplot(u_a, aes(x = X1, y = X2, color = Groups)) +
  geom_point() +
  labs(x= "UMAP1", y = "UMAP2", title = "Stimulation: A")
ggsave(file = "a_umap.png", plot = p_a, path = figures_folder, width = 7, height = 5, dpi = 300)

p_t <- ggplot(u_t, aes(x = X1, y = X2, color = Groups)) +
  geom_point() +
  labs(x= "UMAP1", y = "UMAP2", title = "Stimulation: T")
ggsave(file = "t_umap.png", plot = p_t, path = figures_folder, width = 7, height = 5, dpi = 300)

p_l <- ggplot(u_l, aes(x = X1, y = X2, color = Groups)) +
  geom_point() +
  labs(x= "UMAP1", y = "UMAP2", title = "Stimulation: L")
ggsave(file = "l_umap.png", plot = p_l, path = figures_folder, width = 7, height = 5, dpi = 300)

p_g <- ggplot(u_g, aes(x = X1, y = X2, color = Groups)) +
  geom_point() +
  labs(x= "UMAP1", y = "UMAP2", title = "Stimulation: G")
ggsave(file = "g_umap.png", plot = p_g, path = figures_folder, width = 7, height = 5, dpi = 300)

p_uu <- ggplot(u_uu, aes(x = X1, y = X2, color = Groups)) +
  geom_point() +
  labs(x= "UMAP1", y = "UMAP2", title = "Unstimulated")
ggsave(file = "unstimulated_umap.png", plot = p_uu, path = figures_folder, width = 7, height = 5, dpi = 300)

