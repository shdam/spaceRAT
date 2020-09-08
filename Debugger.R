source("R/createEset.R")
load("data/exprs_dmap.rda")
load("data/pData_dmap.rda")
load("data/exprs_ilaria.rda")
load("data/pData_ilaria.rda")
eset_dmap <- createEset(exprs_dmap,pData_dmap)


source("R/findDEGenes.R")
DEgenes <- findDEGenes(eset_dmap,"cell_types",0.05,2)

source("R/buildScaffold.R")
g <- buildScaffold(exprs_dmap,pData_dmap,"cell_types")
plot(g)

source("R/projectSample.R")
projectSample(g,exprs_ilaria,pData_ilaria,"cancer_type")


