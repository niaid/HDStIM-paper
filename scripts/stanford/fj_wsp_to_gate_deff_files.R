#!/usr/bin/env Rscript --vanilla

# PURPOSE: To export cell populations from workspaces.

library(tidyverse)
library(flowCore)
library(CytoML)
library(flowWorkspace)
library(doMC)
library(arrow)

data_folder <- file.path("data", "stanford", "wsp-fcs")
results_folder <- file.path("results", "stanford")

wsp_files <- read_tsv(file.path(results_folder, "wsp-to-fcs-filtered-with-stim-subj.tsv" )) %>%
        dplyr::pull("wsp") %>%
        unique()

gates_96_samp <- read_tsv(file.path(results_folder, "28-gates-to-use-96-wsp.tsv"))
gates_5_samp <- read_tsv(file.path(results_folder, "28-gates-to-use-5-wsp.tsv"))

wsp_5_live <- read_tsv(file.path(results_folder, "5-wsp-with-live.tsv")) %>%
        dplyr::pull("wsp")



fetch_expr <- function(ws, wsp_file, gates_x_samp){
        print("FlowJo to gating set")
        gs <- flowjo_to_gatingset(ws, name = "All Samples", path = file.path(data_folder, "FCS"), execute = TRUE)
        fcs_files <- sampleNames(gs)

        gate_fcs_expr_out <- data.frame()
        for(g in 1:nrow(gates_x_samp)){
                gate <- gates_x_samp[[g, "gates"]]
                print(gate)
                fcs_expr_out <- data.frame()
                for(s in 1:length(fcs_files)){
                        fcs_file <- fcs_files[s]
                        print(fcs_file)
                        fcs_position <- max(gregexpr(".fcs_", fcs_file)[[1]])
                        fcs_file <- substr(fcs_file, 1, fcs_position + 3)
                        expr_dat <- flowCore::exprs(gs_pop_get_data(gs, gate, inverse.transform = T)[[s]]) %>%
                                data.frame()
                        if(nrow(expr_dat) == 0){
                                df_na <- data.frame(matrix(NA, ncol = ncol(expr_dat), nrow = 1))
                                colnames(df_na) <- colnames(expr_dat)
                                expr_dat <- df_na
                        }

                        fcs_expr_out <- rbind(fcs_expr_out, data.frame("wsp" = wsp_file, "fcs_file" = fcs_file, "gate" = gate, expr_dat))
                }
                gate_fcs_expr_out <- rbind(gate_fcs_expr_out, fcs_expr_out)
        }
        return(gate_fcs_expr_out)
}

registerDoMC(8)
output_folder <- file.path(results_folder, "wsp-to-cell-pop-expr")
dir.create(output_folder, recursive = TRUE)

foreach(f = 1:length(wsp_files)) %dopar% {
# Load the FlowJo workspace
        print(wsp_files[f])
        wsp_file <- wsp_files[f]
        if(!(wsp_files[f] %in% wsp_5_live)){
                ws <- CytoML::open_flowjo_xml(file.path(data_folder, wsp_file))
                wsp_expr_dat <- fetch_expr(ws, wsp_file, gates_96_samp)
                arrow::write_feather(wsp_expr_dat, file.path(output_folder, paste0(wsp_file,".feather")))
        }else{
                ws <- CytoML::open_flowjo_xml(file.path(data_folder, wsp_file))
                wsp_expr_dat <- fetch_expr(ws, wsp_file, gates_5_samp)
                arrow::write_feather(wsp_expr_dat, file.path(output_folder, paste0(wsp_file,".feather")))
        }
}

stop()
# Temp
#temp_files <- as.character(data.frame(str_split(list.files(output_folder),".feather"))[1,])
#wsp_files <- setdiff(wsp_files, temp_files)
