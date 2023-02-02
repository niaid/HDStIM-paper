library(CytoML)
library(flowWorkspace)


dataDir <- system.file("extdata",package="flowWorkspaceData")
gs_dir <- list.files(dataDir, pattern = "gs_manual",full = TRUE)
gs <- load_gs(gs_dir)
head(flowCore::exprs(gs_pop_get_data(gs, "CD4")[[1]])[,5:7])

wsFile <- file.path(file.path("results", "pregnancy-baseline.wps.wsp"))
ws <- CytoML::open_flowjo_xml(wsFile)
gs <- flowjo_to_gatingset(ws, name = 1, path = file.path("data", "pregnancy_baseline"), execute = TRUE)

