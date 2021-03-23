# if the package is properly built, then there is no need to source() or load()
load("data/exprs_dmap.rda")
load("data/pData_dmap.rda")
load("data/exprs_ilaria.rda")
load("data/pData_ilaria.rda")
load("data/gene_name_mapper_hs.rda")
source("R/createEset.R")
source("R/class_scaffoldSpace.R")
source("R/convertGeneName.R")
source("R/createContrast.R")
source("R/convertGeneName.R")
source("R/mapGene.R")
source("R/findDEGenes.R")
source("R/buildScaffold.R")
source("R/loadingPlot.R")
source("R/projectSample.R")

# test createEset()
eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types")
exprs_dmap[1,] <- NA
exprs_dmap[5,10:15] <- NA
pData_dmap[19:24,1] <-NA
eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types")


# test convertGeneName()
dmap <- read.table("data_DMap.gct",skip=2,sep="\t",header=T,quote = "",colClasses = c("character","character",rep("numeric", 211)))
library(dplyr)
dmap <- dmap %>% group_by(NAME) %>% filter(row_number()==1) %>% tibble::column_to_rownames(var="NAME")
counts <- convertGeneName(dmap,to="hgnc_symbol")
counts <- convertGeneName(dmap,to="entrezgene_id")
counts <- convertGeneName(dmap,to="refseq_mrna")
#test convert transcript
transcripts <- unique(gene_name_mapper_hs$ensembl_transcript_id)[1:5000]
dat <- matrix(1,nrow=5000,ncol=3)
rownames(dat) <- transcripts
convertGeneName(dat)
# test removing version number
count_blue <- read.table("data_blue_counts.txt",header=T)
rownames(count_blue) <- count_blue$gene
count_blue <- convertGeneName(count_blue)


# test mapGene()
mapGene(dmap$Description[1:5],to="entrezgene_id")
mapGene(dmap$Description[1:5],to="hgnc_symbol")
mapGene(dmap$Description[1:5],to="ensembl_gene_id")
mapGene(dmap$Description[1:5],to="refseq_mrna")


# test findDEGenes()
DEgenes <- findDEGenes(eset_dmap,0.05,2)


# test buildScaffold()
g <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",plot_mode ="tiny_label")
DMAP_scaffold <- buildScaffold(exprs_dmap,pData_dmap,"cell_types")
# save(DMAP_scaffold,file="DMAP_scaffold.rda")
# test use of prebuilt
g <- buildScaffold("prebuilt_DMAP")


# test loadingPlot
loadingPlot(g)
loadingPlot(g,num_genes = 3)


# test projectSample
projectSample(g,exprs_ilaria,pData_ilaria,"cancer_type")
g2 <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",pcs=c(3,4))
projectSample(g2,exprs_ilaria,pData_ilaria,"cancer_type")


