# Purpose: to compare the list of FCS files from the disk, the SDM data sent by Holden and the files exported from the workspaces.

library(tidyverse)
library(readxl)
library(ggVennDiagram)

disk_list <- list.files(file.path("data/stanford/wsp-fcs/FCS"), "*.fcs")

sdm_list <- read_excel(file.path("results", "stanford", "fcs_list_2014pCyTOFalldata.xlsx")) %>%
        dplyr::pull(fcs_file)

wsp_list <- read_excel(file.path("results", "stanford","wsp-to-fcs-filtered-with-stim-subj.xlsx")) %>%
        dplyr::pull(name)

x <- list(disk_list = disk_list,
        sdm_list = sdm_list,
        wsp_list = wsp_list)

plt_ven <- ggVennDiagram(x)
ggsave("fcs-diff-sources.png", plot = plt_ven, path = file.path("figures", "stanford"),
width = 7, height = 7, dpi = 300)
