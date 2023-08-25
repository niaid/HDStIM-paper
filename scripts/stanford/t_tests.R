# Purpose: read the concatenated feather file and fix stimtype or cell population inconsistencies

library(arrow)
library(tidyverse)

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


test_dat <- dat %>% dplyr::group_by(subject, stim_type, cell_population) %>%
        summarize(across(where(is.numeric), median)) %>%
        ungroup() %>%
        dplyr::filter(subject != "ParentName")

cell_pops <- unique(test_dat$cell_population)

df_t_out <- tibble()
for(i in 1:length(cell_pops)){
        pop <- cell_pops[i]
        print(pop)
        for(j in 1:length(stim_lab)){
                stim <- stim_lab[j]
                print(stim)
                for(k in 1:length(state_markers)){
                        marker <- state_markers[k]
                        x <- test_dat %>% dplyr::filter(cell_population == pop & stim_type == stim) %>%
                                dplyr::pull(marker)
                        y <- test_dat %>% dplyr::filter(cell_population == pop & stim_type == unstim_lab) %>%
                                dplyr::pull(marker)
                        t_res <- t.test(x, y) %>%
                                broom::tidy()
                        df_temp <- tibble("cell_population" = pop, "stim_type" = stim, "marker" = marker, t_res)
                        df_t_out <- bind_rows(df_t_out, df_temp)
                }
        }
}

write_tsv(df_t_out, file.path(results_folder, "marker-by-marker-t-tests.tsv"))












