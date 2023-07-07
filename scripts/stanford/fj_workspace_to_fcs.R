#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a table of all the FlowJo workspaces and their associated FCS files.
# Find the FCS files that are missing etc.

library(tidyverse)
library(flowCore)
library(CytoML)

data_folder <- file.path("data", "stanford", "wsp-fcs")
results_folder <- file.path("results", "stanford")
dir.create(results_folder, recursive = TRUE)

wsp_files <- list.files(data_folder, pattern = ".wsp")

all_out <- tibble()
for(f in wsp_files){
# Load the FlowJo workspace
        workspace <- CytoML::open_flowjo_xml(file.path(data_folder, f))

# Get the names of the FCS files in the workspace
        fcs_files <- fj_ws_get_samples(workspace) %>%
                as_tibble()
        fcs_files <- bind_cols("wsp" = f, fcs_files)
        bool_out <- tibble()
        for(i in 1:nrow(fcs_files)){
                bool <- file.exists(file.path(data_folder, "FCS", fcs_files[[i, "name"]]))
                bool_out <- bind_rows(bool_out, tibble("file.exists" = bool))
        }

        fcs_files <- bind_cols(fcs_files, bool_out)
        all_out <- bind_rows(all_out, fcs_files)
}

# Write unfiltered data to a tsv file.
write_tsv(all_out, file.path(results_folder, "wsp-to-fcs-unfiltered-data.tsv"))

# Add stim type column
# First filter all_out to get rid of all the missing FCS files. That should get rid of all the technical controls
# and unassigned events.

filt_out <- all_out %>% dplyr::filter(file.exists == TRUE)

stims_out <- tibble()
for(j in 1:nrow(filt_out)){
        underscore_position <- max(gregexpr("_", filt_out[[j, "name"]])[[1]])
        stim <- substr(filt_out[[j, "name"]], underscore_position + 1, nchar(filt_out[[j,"name"]]) - 4)
        stims_out <- bind_rows(stims_out, tibble("stim.type" = stim))
}
filt_out <- bind_cols(filt_out, stims_out)

# Add subject column
subj_out <- tibble()
for(k in 1:nrow(filt_out)){
        underscore_positions <- gregexpr("_", filt_out[[k, "name"]])[[1]][1:2]
        subj <- substr(filt_out[[k, "name"]], underscore_positions[1] + 1, underscore_positions[2] - 1)
        subj_out <- bind_rows(subj_out, tibble("subject" = subj))
}
filt_out <- bind_cols(filt_out, subj_out)

write_tsv(filt_out, file.path(results_folder, "wsp-to-fcs-filtered-with-stim-subj.tsv"))
