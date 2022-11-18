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
human_scores <- bpa_res$Scores
human_scores <- human_scores[complete.cases(human_scores),]
library(Biobase)
eset <- ExpressionSet(human_scores, phenoData=AnnotatedDataFrame(pData_dmap))
cell_types <- Biobase::pData(eset)[,1]
design <- model.matrix(~0+ cell_types)
colnames(design) <- gsub("cell_types","",colnames(design))
cm <- createContrast(eset)
fit <- limma::lmFit(eset,design)
fit <- limma::contrasts.fit(fit, contrast=cm)
fit <- limma::eBayes(fit)
print(fit)
de_genes <- sapply(colnames(cm), function(name){
  rownames(limma::topTable(fit, coef=name, number = Inf, p.value =0.05, adjust.method="fdr"))
})
debp <- unique(unlist(de_genes))
human_eset <- human_eset[debp,]

# Prepare a dataframe for ggplot2
pca <- prcomp(t(exprs(human_eset)),scale=T)
PC1 <- pca$x[,1]
PC2 <- pca$x[,2]
scaffold_group <- Biobase::pData(human_eset)[["cell_types"]]
df <- data.frame(PC1,PC2,scaffold_group)

# calculate centroids
centroids_df <- suppressMessages(df %>% dplyr::group_by(scaffold_group) %>% dplyr::summarise(mean_PC1=mean(PC1),mean_PC2=mean(PC2)))

# define color scheme
total_types <- length(unique(scaffold_group))
my_col <- RColorBrewer::brewer.pal(total_types,"Paired")

# ggplot2 for dot mode
g <- ggplot2::ggplot()+
    ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=scaffold_group))+
    ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_PC1,mean_PC2,label=scaffold_group,color=scaffold_group),fontface="bold",show.legend = FALSE)
g

# load mouse data
library(MouseGastrulationData)
sce <- EmbryoAtlasData(samples = 21)

# convert MOUSE gene set from rds to gmt
gene_sets_mm <- readRDS("data/Mm.c5.all.v7.1.entrez.rds")
library(MEGENA)
output.geneSet.file(gene_sets_mm,"Mm.c5.all.v7.1.entrez.gmt")

# read gmt file
gene_sets_mm <- getGmt("data/Mm.c5.all.v7.1.entrez.gmt")
