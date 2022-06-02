#!/usr/bin/env Rscript --vanilla

# PURPOSE: Generate plots to compare median based analysis from Bendall etal 
# and re-analysis by HDStIM.
# Figure 1 B & E.
# Supplementary Figure 1 A & B.

library(tidyverse)
library(circlize)
library(ggrepel)
library(ComplexHeatmap)
library(diptest)
library(reshape2)
library(scales)

results_folder <- file.path("results", "bone_marrow" )
figures_folder <- file.path("figures", "bone_marrow" )

# Read all Fisher's exact data. 
all_fisher <- read_tsv(file.path(results_folder, "all_fisher_p_val.tsv"))
tabs3 <- read_tsv(file.path("meta", "paper_donor1_table_s3.tsv")) %>%
  rename("cell_population" = cluster)

tabs3_count <- tabs3 %>% group_by(stim_type, cell_population) %>%
  count(name = "count_marker")

p_m_count <- full_join(all_fisher, tabs3_count)
p_m_count[is.na(p_m_count$count_marker),"count_marker"] <- 0

max_m <- max(p_m_count$count_marker) + 1

three <- p_m_count %>% dplyr::filter(count_marker == 3 & p.value > 0.05)
three_b <- p_m_count %>% dplyr::filter(count_marker == 3 & p.value < 0.05)
two <- p_m_count %>% dplyr::filter(count_marker == 2 & p.value > 0.05)
two_b <- p_m_count %>% dplyr::filter(count_marker == 2 & p.value < 0.05)
one <- p_m_count %>% dplyr::filter(count_marker == 1 & p.value > 0.05)
one_b <- p_m_count %>% dplyr::filter(count_marker == 1 & p.value < 0.05)
zero <- p_m_count %>% dplyr::filter(count_marker == 0 & p.value > 0.05)
zero_b <- p_m_count %>% dplyr::filter(count_marker == 0 & p.value < 0.05)

# Generate heatmap with boruta rankings for all the significant in column 1
# and put an astrix on the marker that median based method found.
bstats <-read_tsv(file.path(results_folder, "K_clust_boruta_stats.tsv"))
mergedat <- merge(one_b, bstats, by = c("stim_type", "cell_population")) %>%
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))

rown <- unique(mergedat$state_marker)
coln <- unique(mergedat$cellpop_stim)
bmat <- matrix(nrow = length(rown), ncol = length(coln), dimnames = list(rown, coln))
for(i in 1:nrow(mergedat)){
  cp <- mergedat[[i, "cellpop_stim"]]
  m <- mergedat[[i, "state_marker"]]
  bmat[m, cp] <- mergedat[[i, "meanImp"]]
}

tabs3_one <- merge(one_b, tabs3, by = c("stim_type", "cell_population")) %>% 
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))
tmat <- matrix(0, nrow = length(rown), ncol = length(coln), dimnames = list(rown, coln))
for(i in 1:nrow(tabs3_one)){
  cp <- tabs3_one[[i, "cellpop_stim"]]
  m <- tabs3_one[[i, "state_marker"]]
  tmat[m, cp] <- 1
}

bmat <- rescale(bmat, to = c(0,1))

col_fun = colorRamp2(c(min(bmat), mean(bmat), max(bmat)), c("#56B4E9", "#009E73", "#E69F00"))
col_fun(seq(min(bmat), mean(bmat), max(bmat)))
bht <- Heatmap(bmat, cluster_columns = FALSE, cluster_rows = FALSE,
               heatmap_legend_param = list(title = "Importance"),
               col = col_fun,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(tmat[i, j] == 1){
                   grid.text("*", x, y, gp = gpar(fontsize = 14))
                 }
               })
# png(file.path(figures_folder, "boruta_column1_sig.png"), width = 10, height = 8, units = "in", res = 300)
# draw(bht)
# dev.off()

## In the heatmap generated above find out how many cellpop-stims have markers 
## ranked higher than the one found significant by the median based method. 
cpst <- colnames(bmat)
df_mgt <- data.frame()
for(cp in cpst){
  df_temp <- data.frame(bmat[,cp]) %>%
    rownames_to_column(var = "state_marker") 
  mi <- data.frame(tmat[,cp] == 1) %>% dplyr::filter(tmat...cp.....1 == TRUE) %>% rownames()
  mi_sc <- df_temp[df_temp$state_marker == mi, "bmat...cp."]
  mgt <- dplyr::filter(df_temp, bmat...cp. > mi_sc) %>%
    nrow
  df_mgt <- rbind(df_mgt, data.frame("cell_population" = cp, "mgt" = mgt))
}

24 - nrow(dplyr::filter(df_mgt, mgt == 0))

# Heatmap with cellpopulation on the columns and stims on the rows
# cells containig propotion of markers that contributed to activated cells
# with multimodal distribution in the parent stimulated samples.
selected_data <- readRDS(file.path(results_folder, "selected_data_all_k_clust.rds"))
stim_types <- selected_data$stim_lab
state_markers <- selected_data$state_markers
selected_data <- selected_data$response_mapping_main %>%
  dplyr::filter(stim_type %in% stim_types) %>%
  droplevels()
cell_pop <- as.character(unique(selected_data$cell_population))

dipmat <- matrix(nrow = length(stim_types), ncol = length(cell_pop), dimnames = list(stim_types, cell_pop))

stim_state <- selected_data %>% group_by(stim_type, cell_population) %>%
  dplyr::select(stim_type, cell_population, all_of(state_markers)) %>%
  group_split()

df_dip_out <- data.frame()
for(z in 1:length(stim_state)){
  dat <- stim_state[[z]]
  st <- as.character(unique(dat$stim_type))
  cp <- as.character(unique(dat$cell_population))
  for(state in state_markers) {
    dipres <- dip.test(dat[,state][[1]])  
    df_temp <- data.frame("stim_type" = st, "cell_population" = cp, "state_marker" = state, "dip_p" = dipres$p.value)
    df_dip_out <- rbind(df_dip_out, df_temp)
  }
}

dip_bstats <- merge(df_dip_out, bstats) %>%
  dplyr::filter(decision == "Confirmed")

df_dip_count <- dplyr::filter(dip_bstats, dip_p < 0.05) %>%
  dplyr::group_by(stim_type, cell_population) %>% count(name = "sig_mark")

df_bstat_count <- dplyr::filter(dip_bstats, decision == "Confirmed") %>%
  dplyr::group_by(stim_type, cell_population) %>% count(name = "decision")

dip_bstats_frac <- merge(df_dip_count, df_bstat_count) %>%
  dplyr::mutate("frac" = sig_mark/decision)

for(i in 1:nrow(dip_bstats_frac)){
  st <- dip_bstats_frac[[i, "stim_type"]]
  cp <- dip_bstats_frac[[i, "cell_population"]]
  dipmat[st, cp] <- dip_bstats_frac[[i, "frac"]]
}

# dipht <- Heatmap(dipmat, cluster_columns = FALSE, cluster_rows = FALSE,
#                heatmap_legend_param = list(title = "Multimodal Fraction"))
# png(file.path(figures_folder, "dip_test_heatmap.png"), width = 8, height = 5, units = "in", res = 300)
# draw(dipht)
# dev.off()

# Generate first scatter plot but now coloured with multimodal fractions.
mfp_dat <- left_join(p_m_count, dip_bstats_frac)

mfp <- ggplot(mfp_dat, aes(x = count_marker, y = -log10(p.value), color = frac)) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  labs(x = "Individual Marker Analysis\nMarkers Significant Per Cell Population-Stimulation", 
       y = "HD Analysis\nFisher's -log10(p-value)",
       color = "Multimodal\nFraction") +
  scale_x_continuous(breaks = seq(from = 0, to = max_m, by = 1 )) +
  scale_color_gradient(low="blue", high="red", na.value="black", breaks = seq(0,1,0.1))
ggsave("no.mark_vs_f_p.value_frac.png", path = figures_folder, plot = mfp,
       width = 7, height = 5, units = "in", dpi = 300)

## Genereate a graph with multiple cutoffs for the diptest. 
dip_bstats <- merge(df_dip_out, bstats) %>%
  dplyr::filter(decision == "Confirmed")

df_dip_mult_out <- data.frame()
for(p in c(0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001)){
    nr <- dplyr::filter(dip_bstats, dip_p < p) %>%
    dplyr::select(stim_type, cell_population) %>%
    unique() %>%
    nrow()
    df_dip_mult_out <- rbind(df_dip_mult_out, data.frame("dip_p" = p, "no_cp_stim" = nr))    
}

df_dip_mult_out$dip_p <- as.character(df_dip_mult_out$dip_p)
dippt <- ggplot(df_dip_mult_out, aes(x = dip_p, y = no_cp_stim)) +
  geom_point() +
  geom_label_repel(aes(label = no_cp_stim)) +
  scale_x_discrete(breaks = c("0.1", "0.05", "0.04", "0.03", "0.02", "0.01", "0.001")) +
  labs(x = "Dip Test P-value Cutoff", y = "No. of Cell Pop. - Stim.\nwith Multimodal Marker Distribution")
ggsave("dip_mult_cutoff.png", path = figures_folder, plot = dippt,
       width = 5, height = 3, units = "in", dpi = 300)


# Generate heatmap with boruta rankings for all the significant in column 1, 2, and 3
# and put an astrix on the marker that median based method found.
bstats <-read_tsv(file.path(results_folder, "K_clust_boruta_stats.tsv"))
onedat <- merge(one_b, bstats, by = c("stim_type", "cell_population")) %>%
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))
twodat <- merge(two_b, bstats, by = c("stim_type", "cell_population")) %>%
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))
threedat <- merge(three_b, bstats, by = c("stim_type", "cell_population")) %>%
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))

tabs3_one <- merge(one_b, tabs3, by = c("stim_type", "cell_population")) %>% 
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))
tabs3_two <- merge(two_b, tabs3, by = c("stim_type", "cell_population")) %>% 
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))
tabs3_three <- merge(three_b, tabs3, by = c("stim_type", "cell_population")) %>% 
  dplyr::mutate("cellpop_stim" = paste0(cell_population, " ~ ", stim_type))

df_1_out <- data.frame()
tmp <- as.character(unique(tabs3_one$cellpop_stim))
for(cpstu in tmp){
  mrk <- dplyr::filter(tabs3_one, cellpop_stim == cpstu)[["state_marker"]]
  imp <- dplyr::filter(onedat, cellpop_stim == cpstu, state_marker == mrk)[["meanImp"]]
  if(nrow(dplyr::filter(onedat, cellpop_stim == cpstu, meanImp > imp)) > 0){
    df_temp <- dplyr::filter(onedat, cellpop_stim == cpstu)
    df_1_out <- rbind(df_1_out, data.frame("group" = 1, "med_mark" = mrk, df_temp))
  }
}

df_2_out <- data.frame()
tmp <- as.character(unique(tabs3_two$cellpop_stim))
for(cpstu in tmp){
  mrk <- dplyr::filter(tabs3_two, cellpop_stim == cpstu)[["state_marker"]]
  imp <- dplyr::filter(twodat, cellpop_stim == cpstu, state_marker %in% mrk)[["meanImp"]]
  if(nrow(dplyr::filter(twodat, cellpop_stim == cpstu, meanImp > min(imp))) > 1){
    df_temp <- dplyr::filter(twodat, cellpop_stim == cpstu)
    df_2_out <- rbind(df_2_out, data.frame("group" = 2, "med_mark" = paste0(mrk, collapse = ", "), df_temp))
  } else{
    print("false")
  }

}

df_3_out <- data.frame()
tmp <- as.character(unique(tabs3_three$cellpop_stim))
for(cpstu in tmp){
  mrk <- dplyr::filter(tabs3_three, cellpop_stim == cpstu)[["state_marker"]]
  imp <- dplyr::filter(threedat, cellpop_stim == cpstu, state_marker %in% mrk)[["meanImp"]]
  if(nrow(dplyr::filter(threedat, cellpop_stim == cpstu, meanImp > min(imp))) > 2){
    df_temp <- dplyr::filter(threedat, cellpop_stim == cpstu)
    df_3_out <- rbind(df_3_out, data.frame("group" = 3, "med_mark" = paste0(mrk, collapse = ", "), df_temp))
  }else{
    print("false")
  }
}

# Generate heatmap.
df123 <- rbind(df_1_out, df_2_out, df_3_out) %>%
  dplyr::mutate("cellpop_stim_grp" = paste0(cellpop_stim, "#", group))
coln <- as.character(unique(df123$cellpop_stim_grp))
mat123 <- reshape2::acast(df123, state_marker ~ cellpop_stim_grp, value.var="meanImp", )
col_breaks <- str_split(colnames(mat123), "#")
col_breaks <- as.numeric(data.frame(col_breaks)[2,])

# Genereate binary matrix with 1 indicating significant in median based method.
tmat123 <- matrix(0, nrow = nrow(mat123), ncol = ncol(mat123), 
                  dimnames = list(rownames(mat123), colnames(mat123)))
tabs3_123 <- rbind(tabs3_one, tabs3_two, tabs3_three) %>%
  dplyr::mutate("cellpop_stim_grp" = paste0(cellpop_stim, "#", count_marker)) %>%
  dplyr::filter(cellpop_stim_grp %in% colnames(tmat123))
for(i in 1:nrow(tabs3_123)){
  cpstg <- tabs3_123[i,"cellpop_stim_grp"]
  st <- tabs3_123[i,"state_marker"]
  tmat123[st, cpstg] <- 1
}

mat123 <- rescale(mat123, to = c(0,1))

col_fun = colorRamp2(c(min(mat123), mean(mat123), max(mat123)), c("#56B4E9", "#009E73", "#E69F00"))
col_fun(seq(min(mat123), mean(mat123), max(mat123)))

mat123ht <- Heatmap(mat123, cluster_columns = FALSE, cluster_rows = FALSE,
               heatmap_legend_param = list(title = "Importance"),
               col = col_fun,
               column_split = col_breaks,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(tmat123[i, j] == 1){
                   grid.text("*", x, y, gp = gpar(fontsize = 14))
                 }
               }
               )
png(file.path(figures_folder, "boruta_column123_sig.png"), width = 14, height = 8, units = "in", res = 300)
draw(mat123ht)
dev.off()


# Generate a scatter plot with Bendall's p-value on x axis
# and Boruta importance score on the y-axis.
b_stats <- read_tsv(file.path(results_folder, "K_clust_boruta_stats.tsv"))
tabs3_p <- read_tsv(file.path("meta", "tableS3WithPValues.tsv"))

bt_dat <- merge(b_stats, tabs3_p)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
bt_dat$meanImp <- range01(bt_dat$meanImp)
btplt <- ggplot(bt_dat, aes(x = adj_p_val, y = meanImp)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  labs(x = "Adjusted P-Value",
       y = "Importance Score") +
  theme(text = element_text(size=20))
ggsave("p_vs_impScore.png", plot = btplt, path = figures_folder, width = 7, height = 5, dpi = 300)

cor_res <- cor.test(bt_dat$meanImp, bt_dat$adj_p_val)
cor_res$p.value
cor_res$estimate


