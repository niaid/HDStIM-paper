#!/usr/bin/env Rscript --vanilla

# PURPOSE: Generate heatmap from the spline regression data. 

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Load fraction lm data based on ranks. 
frac_lm <- read_tsv(file.path(results_folder, "resp_unstim_frac_spline_rank_no_gender.tsv"))

# Typo correction
frac_lm$cell_population <- stringr::str_replace(frac_lm$cell_population, "Niave cytotoxic T cells", "Naive cytotoxic T cells")
frac_lm$cell_population <- stringr::str_replace(frac_lm$cell_population, "Niave T helper cells", "Naive T helper cells")

frac_age <- frac_lm %>% dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate("p.value.log" = -log10(p.value)) %>% 
  dplyr::select(cell_population, stim_type, term, p.value.log)
group <- rep(c(1,2,3), nrow(frac_age)/3)
frac_age <- frac_age %>% dplyr::mutate(group = group)
frac_age <- dplyr::filter(frac_age, p.value.log > -log10(0.05))
frac_age$cell_population <- paste0(frac_age$cell_population, " [", frac_age$group, "]")

c_names <- as.character(unique(frac_age$cell_population))
r_names <- c("U_A", "U_T", "U_G")
mat <- matrix(nrow = length(r_names), ncol = length(c_names) )
rownames(mat) <- r_names
colnames(mat) <- c_names
for(i in 1:nrow(frac_age)){
  st <- frac_age[[i, "stim_type"]]
  cp <- frac_age[[i, "cell_population"]]
  mat[st,cp] <- frac_age[[i, "p.value.log"]]
}

# Stimulation specific heatmaps. 
# Load lm model data.
frac_all <- frac_lm %>% dplyr::filter(term != "(Intercept)")
frac_all <- frac_all %>% dplyr::mutate("groups" = case_when(
  term == "bs(age_at_sample, knots = c(8, 14), degree = 1)1" ~ "<8",
  term == "bs(age_at_sample, knots = c(8, 14), degree = 1)2" ~ "8-14",
  term == "bs(age_at_sample, knots = c(8, 14), degree = 1)3" ~ ">14"))

frac_a <- frac_all %>% dplyr::filter(stim_type == "U_A") %>%
  dplyr::filter(cell_population %in% c("Activated monocytes", "Myeloid DC", "pDC")) %>% droplevels()
frac_t <- frac_all %>% dplyr::filter(stim_type == "U_T") %>% 
  dplyr::filter(cell_population %in% c("45RA+ double negative T cells",
                                       "Activated T helper cells",
                                       "Cytotoxic NKT cells",
                                       "Naive cytotoxic T cells",
                                       "Naive T helper cells",
                                       "T helper cells")) %>%
                  droplevels()

cols <- c("<8", "8-14", ">14")
rows_a <- as.character(unique(frac_a$cell_population))
mat_a <- matrix(nrow = length(rows_a), ncol = length(cols), dimnames = list(rows_a, cols))
rows_t <- as.character(unique(frac_t$cell_population))
mat_t <- matrix(nrow = length(rows_t), ncol = length(cols), dimnames = list(rows_t, cols))

for(i in 1:nrow(frac_a)){
  gr <- frac_a[[i, "groups"]]
  cp <- frac_a[[i, "cell_population"]]
  mat_a[cp, gr] <- frac_a[[i, "estimate"]]
}

for(i in 1:nrow(frac_t)){
  gr <- frac_t[[i, "groups"]]
  cp <- frac_t[[i, "cell_population"]]
  mat_t[cp, gr] <- frac_t[[i, "estimate"]]
}


# Plot for U_A and U_T
col_fun = colorRamp2(c(-20, 0, 25), c("blue", "white", "red"))
col_fun(seq(-20, 25))

mat_all <- rbind(mat_a, mat_t)
r_breaks <- c(rep("U_IFNa", 3), rep("U_TCR", 6))
ht <- Heatmap(mat_all, 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_split = r_breaks,
              column_title = "Age Group",
              column_title_side = "bottom",
              column_names_side = "top",
              column_names_rot =0,
              column_names_centered = TRUE,
              width = ncol(mat_all)*unit(20, "mm"),
              heatmap_legend_param = list(title = "Effect\nSize",fontsize = 14,
                                          legend_height = unit(2, "cm"),
                                          legend_width = unit(1, "cm"),
                                          legend_position = "bottom"),
              row_names_gp = gpar(fontsize = 14),
              column_names_gp = gpar(fontsize = 15),
              column_title_gp = gpar(fontsize = 16), 
              show_heatmap_legend = FALSE)
png(file.path(figures_folder, "spline_heatmap_v2.png"), width = 7, 
    height = 5, units = "in", res = 600)
draw(ht)
dev.off()



# OLD CODE
# hmap_a <- ComplexHeatmap::Heatmap(mat_a, cluster_rows = FALSE, cluster_columns = FALSE,
#                                 row_title = "U_A",
#                                 row_title_side = "left",
#                                 column_title = "Age Group",
#                                 column_title_side = "bottom",
#                                 column_names_side = "top",
#                                 column_names_rot = 0,
#                                 column_names_centered = TRUE,
#                                 show_heatmap_legend = FALSE)
# 
# #Plot for U_T
# hmap_t <- ComplexHeatmap::Heatmap(mat_t, cluster_rows = FALSE, cluster_columns = FALSE,
#                                   row_title = "U_T",
#                                   row_title_side = "left",
#                                   column_title = "Age Group",
#                                   column_title_side = "bottom",
#                                   column_names_side = "top",
#                                   column_names_rot = 0,
#                                   column_names_centered = TRUE,
#                                   heatmap_legend_param = list(title = "Effect\nSize"))
# 
# 
# file_a <- file.path(figures_folder, "spline_heatmap_all.png")
# png(filename = file_a, width = 7, height = 5, units = "in", res = 600)
# hmap_all <- hmap_a %v% hmap_t
# draw(hmap_all)
# dev.off()

#  heatmap_legend_param = list(title = "Effect\nSize", title_gp = gpar(fontsize = 4), labels_gp = gpar(fontsize = 4),  at = c(-2, 0, 2)),

