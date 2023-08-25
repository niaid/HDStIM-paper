# Purpose: concatenate all the feather files into a single dataframe,
# do arcsinh transformation for HDStIM and change column names to antigen symbols.

library(arrow)
library(tidyverse)
library(doMC)


feather_folder <- file.path("results", "stanford", "wsp-to-cell-pop-expr")
results_folder <- file.path("results", "stanford")

samp_info <- read_tsv(file.path(results_folder, "wsp-to-fcs-filtered-with-stim-subj.tsv")) %>%
        dplyr::select(name, stim.type, subject)

fea_files <- list.files(feather_folder)

df_all_fea <- data.frame()
registerDoMC(12)
df_all_fea <- foreach(i = 1:length(fea_files), .combine = bind_rows) %dopar% {
        print(fea_files[i])
        dat_in <- arrow::read_feather(file.path(feather_folder, fea_files[i]))
        print(ncol(dat_in))
        dat_in
}

state_markers <-  list("Nd144Di" = "pPLCG2",
                       "Nd150Di" = "pSTAT5",
                       "Eu153Di" = "pSTAT1",
                       "Gd158Di" = "pSTAT3", 
                       "Gd156Di" = "pP38",
                       "Dy164Di" = "IkBa",
                       "Ho165Di" = "pCREB",
                       "Yb171Di" = "pERK1_2",
                       "Yb172Di" = "Ki67")

# Change column names based on the key-value pairs of the state_markers list.
new_colnames <- names(df_all_fea)
for (old_name in names(df_all_fea)) {
  if (old_name %in% names(state_markers)) {
    new_colnames[new_colnames == old_name] <- state_markers[old_name]
  }
}
colnames(df_all_fea) <- new_colnames

# Filter data for NAs and other parameters.
df_state <- df_all_fea %>% dplyr::select(wsp, fcs_file, gate, all_of(unname(unlist(state_markers)))) %>%
        na.omit() %>%
        dplyr::filter(!grepl("_22-", fcs_file)) %>%
        as_tibble()

rm(df_all_fea)
gc()

# Function to apply arcsinh transformation to a vector
arcsinh_transform <- function(x) {
  asinh(x/5)
}

## Add stim type column
fcs_files <- unique(df_state$fcs_file)
registerDoMC(6)
stims_subj_out <- tibble()
stims_subj_out <- foreach(j = 1:length(fcs_files), .combine = bind_rows) %dopar% {
        file_n <- fcs_files[j]
        print(file_n)
        df_temp <- dplyr::filter(df_state, fcs_file == file_n) %>%
                mutate(cell_population = paste0(sub(".*/([^/]+)/[^/]+$", "\\1", gate), "/", sub(".*/([^/]+)$", "\\1", gate))) %>%
                mutate(across(where(is.numeric), arcsinh_transform))

        
        # Fetch stim type
        underscore_position <- max(gregexpr("_", df_temp[[1, "fcs_file"]])[[1]])
        stim <- substr(df_temp[[1, "fcs_file"]], underscore_position + 1, nchar(df_temp[[1,"fcs_file"]]) - 4)

        # Fetch subject id
        underscore_positions <- gregexpr("_", file_n)[[1]][1:2]
        subj <- substr(file_n, underscore_positions[1] + 1, underscore_positions[2] - 1)

        stim_ <- bind_cols(tibble("stim_type" = stim), df_temp)
        bind_cols(tibble("subject" = subj), stim_)
}

arrow::write_feather(stims_subj_out, file.path(results_folder, "arcsinh-data-for-hdstim.feather"))
