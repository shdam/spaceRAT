source("func_create_eset_dmap.R")
eset_dmap <- create_eset_dmap("data_DMap.gct")

source("R/find_de_genes.R")
DE_genes <- find_de_genes(eset_dmap, "cell_types",0.05, 2)

source("R/buildScaffold.R")
g <- buildScaffold(eset_dmap,"cell_types")
print(g)

source("func_create_eset_ilaria.R")
eset_ilaria <- create_eset_ilaria("data_Ilaria_counts.rds","data_Ilaria_pd.rds")
idx_AML <- grep("AML",pData(eset_ilaria)$cancer_type)
eset_ilaria <- eset_ilaria[,idx_AML]

source("R/projectSample.R")
projectSample(g,eset_ilaria,"cancer_type")
