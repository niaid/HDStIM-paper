# Purpose: read the concatenated feather file and fix stimtype or cell population inconsistencies

library(arrow)
library(tidyverse)

feather_folder <- file.path("results", "stanford", "wsp-to-cell-pop-expr")
results_folder <- file.path("results", "stanford")


dat <- arrow::read_feather(file.path(results_folder, "arcsinh-data-for-hdstim.feather"))

# Fix stim type
dat$stim_type <- str_to_upper(dat$stim_type)
dat$stim_type <- str_replace_all(dat$stim_type, "LPA", "LPS")

# Fix cell populations
dat$cell_population <- str_replace_all(dat$cell_population, "Live Singles/Basophils", "Basophils")
dat$cell_population <- str_replace_all(dat$cell_population, "Live/Basophils", "Basophils")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4 T cels/Central Memory", "CD4+ T Central Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4 T cels/Effector", "CD4+ T Effector")
dat$cell_population <- str_replace_all(dat$cell_population, 'CD4 T cels/HLA\\-DR\\+CD38\\+', "CD4+ T HLA-DR+CD38+")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4 T cels/Naive", "CD4+ T Naive")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4 T cels/T regs", "CD4+ T Regs")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8 T cells/Central Memory", "CD8+ T Central Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8 T cells/Effector", "CD8+ T Effector")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8 T cells/HLA\\-DR\\+CD38\\+", "CD8+ T HLA-DR+CD38+")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8 T cells/Naive", "CD8+ T Naive")
dat$cell_population <- str_replace_all(dat$cell_population, "CD3\\+/non CD4 non CD8", "CD4-CD8- T Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/IgA\\+", "B Cell IgA+")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/IgD\\+ memory", "B Cell IgD+ Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/IgD\\-CD27\\-", "B Cell IgD-CD27-")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/Naive", "B Cell Naive")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/Switched memory", "B Cell Switched Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "B cells/Transitional B cells", "B Cell Transitional")
dat$cell_population <- str_replace_all(dat$cell_population, "mDC/type 2 mDC", "mDC Type 2")
dat$cell_population <- str_replace_all(dat$cell_population, "DC/pDC", "pDC")
dat$cell_population <- str_replace_all(dat$cell_population, "NK cells/CD16\\+ NK cells", "CD16+ NK Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "NK cells/CD16\\- CD56hi NK cells", "CD16-CD56hi NK Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "Lymphocytes/NK T cells", "NKT Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "Monocytes/CD16hi monocytes", "CD16hi Monocyte")
dat$cell_population <- str_replace_all(dat$cell_population, "Monocytes/CD16low monocytes", "CD16low Monocyte")
dat$cell_population <- str_replace_all(dat$cell_population, "Plasmablasts/Q1: IgA\\- , Ki67\\+", "Plasmablasts IgA-Ki67+")
dat$cell_population <- str_replace_all(dat$cell_population, "Plasmablasts/Q2: IgA\\+ , Ki67\\+", "Plasmablasts IgA+Ki67+")
dat$cell_population <- str_replace_all(dat$cell_population, "Plasmablasts/Q3: IgA\\+ , Ki67\\-", "Plasmablasts IgA+Ki67-")
dat$cell_population <- str_replace_all(dat$cell_population, "Plasmablasts/Q4: IgA\\- , Ki67\\-", "Plasmablasts IgA-Ki67-")
dat$cell_population <- str_replace_all(dat$cell_population, "170Er_CD3\\+/CD4\\-CD8\\-", "CD4-CD8- T Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4/Central Memory", "CD4+ T Central Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4/Effector", "CD4+ T Effector")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4/HLA\\-DR\\+CD38\\+", "CD4+ T HLA-DR+CD38+")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4/Naive", "CD4+ T Naive")
dat$cell_population <- str_replace_all(dat$cell_population, "CD4/T regs", "CD4+ T Regs")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8/Central Memory", "CD8+ T Central Memory")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8/Effector", "CD8+ T Effector")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8/HLA\\-DR\\+CD38\\+", "CD8+ T HLA-DR+CD38+")
dat$cell_population <- str_replace_all(dat$cell_population, "CD8/Naive", "CD8+ T Naive")
dat$cell_population <- str_replace_all(dat$cell_population, "NK cells/CD16- CD56 high NK cells", "CD16-CD56hi NK Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "NK cells/CD16- NK cells", "CD16- NK Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "NK cells/CD6\\+ NK cells", "CD16+ NK Cell")
dat$cell_population <- str_replace_all(dat$cell_population, "Monocytes/CD16 high", "CD16hi Monocyte")
dat$cell_population <- str_replace_all(dat$cell_population, "Monocytes/CD16 low", "CD16low Monocyte")

arrow::write_feather(dat, file.path(results_folder, "arcsinh-data-stim-pop-fixed-for-hdstim.feather"))
