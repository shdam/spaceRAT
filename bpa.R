# load human data
dmap <- read.table("data_DMap.gct",skip=2,sep="\t",header=T,quote = "",colClasses = c("character","character",rep("numeric", 211)))
dmap <- dmap[,-grep("NK|BCELL|TCELL|DEND",colnames(dmap))]
rownames(dmap) <- dmap[,1]
dmap <- dmap[,-c(1,2)]
dmap <- as.matrix(dmap)

# do human GO
library(GSEABase)
gene_sets_human <- getGmt("data/c5.go.bp.v7.3.entrez.gmt")
library(singscore)
bpa_res <- multiScore(rankData = rankGenes(dmap), upSetCol= gene_sets_human, knownDirection = TRUE)
library(dplyr)
human_scores <- bpa_res$Scores %>% as.data.frame() %>% drop_na()
library(RAT)
source("R/buildScaffold.R")
buildScaffold(log(human_scores),pData_dmap,"cell_types")


load("data/pData_dmap.rda")
pca <- prcomp(t(human_scores))
PC1 <- pca$x[,1]
PC2 <- pca$x[,2]
scaffold_group <- pData_dmap$cell_types
df <- data.frame(PC1,PC2,scaffold_group)

ggplot2::ggplot()+
        ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=scaffold_group))+
        ggplot2::ggtitle("PCA from BPA score")


# load mouse data
library(MouseGastrulationData)
sce <- EmbryoAtlasData(samples = 21)

# convert MOUSE gene set from rds to gmt
gene_sets_mm <- readRDS("data/Mm.c5.all.v7.1.entrez.rds")
library(MEGENA)
output.geneSet.file(gene_sets_mm,"Mm.c5.all.v7.1.entrez.gmt")

# read gmt file
gene_sets_mm <- getGmt("data/Mm.c5.all.v7.1.entrez.gmt")
