#!/usr/bin/env Rscript --vanilla

# PURPOSE: To generate a table of all the FlowJo workspaces and their associated FCS files.
# Find the FCS files that are missing etc.

library(tidyverse)
library(flowCore)
library(CytoML)
library(flowWorkspace)
library(doMC)

data_folder <- file.path("data", "stanford", "wsp-fcs")
results_folder <- file.path("results", "stanford")

wsp_files <- read_tsv(file.path(results_folder, "wsp-to-fcs-filtered-with-stim-subj.tsv" )) %>%
        dplyr::pull("wsp") %>%
        unique()

# Check if all the worspace have the same gating scheme and population names
all_gates <- tibble()
registerDoMC(12)
all_gates <- foreach(f = 1:length(wsp_files), .combine = rbind) %dopar% {
# Load the FlowJo workspace
        print(wsp_files[f])
        ws <- CytoML::open_flowjo_xml(file.path(data_folder, wsp_files[f]))
        
        #fj_ws_get_sample_groups(ws)
        gs <- flowjo_to_gatingset(ws, name = "All Samples", path = file.path(data_folder, "FCS"), execute = FALSE)

        gs_stats <- gs_pop_get_stats(gs, xml = TRUE) %>%
                pull(pop) %>%
                unique()

        temp_df <- tibble("wsp" = wsp_files[f], "gates" = gs_stats)
        temp_df
}

# Count the number of files with the same gating labeling. 
gate_counts <- all_gates %>% dplyr::group_by(gates) %>% summarize("count"= n())
write_tsv(gate_counts, file.path(results_folder, "filtered-wsp-gate-counts.tsv"))

# There were 5 wsp that had a gate "Live" instaed of "Live Singles" make them different from the
# remaining 96.
wsp_with_live <- all_gates %>% dplyr::filter(gates == "/DNA/Intact cells/Live")

# Using the filtered-ws-gate-counts.tsv generate (manually by Excel) two files with the 28 terminal population
# one for 96 wsp and the other for the 5 wsp that were different in gate names.

### DONE SO far ###








