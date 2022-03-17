# load all first

#---------------------------------------
# Update gene name converter regularly
#---------------------------------------
load("data/gene_id_converter_hs.rda")
dim(gene_id_converter_hs)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_id_converter_hs <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","refseq_mrna", "entrezgene_id","hgnc_symbol"),
                           mart = ensembl)

colnames(gene_id_converter_hs) <- c("ensembl_gene", "ensembl_transcript", "refseq_mrna", "entrez", "hgnc_symbol")

# clean up
any(is.na(gene_id_converter_hs$ensembl_gene))
"" %in% gene_id_converter_hs$ensembl_gene

any(is.na(gene_id_converter_hs$ensembl_transcript))
"" %in% gene_id_converter_hs$ensembl_transcript

any(is.na(gene_id_converter_hs$refseq_mrna))
"" %in% gene_id_converter_hs$refseq_mrna
idx <- which(gene_id_converter_hs$refseq_mrna=="")
gene_id_converter_hs$refseq_mrna[idx] <- NA

any(is.na(gene_id_converter_hs$entrez))
"" %in% gene_id_converter_hs$entrez

any(is.na(gene_id_converter_hs$hgnc_symbol))
"" %in% gene_id_converter_hs$hgnc_symbol
idx <- which(gene_id_converter_hs$hgnc_symbol=="")
gene_id_converter_hs$hgnc_symbol[idx] <- NA

usethis::use_data(gene_id_converter_hs)

#------------------------------------
# internal
# test createEset()
#------------------------------------
eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types",classes=c("HSC","MONO"))
exprs_dmap[1,] <- NA
exprs_dmap[5,10:15] <- NA
pData_dmap[19:24,1] <-NA
eset_dmap <- createEset(exprs_dmap,pData_dmap,"cell_types")


#------------------------------------
# internal
# test convertGeneName()
#------------------------------------
dmap <- read.table("data_DMap.gct",skip=2,sep="\t",header=T,quote = "",colClasses = c("character","character",rep("numeric", 211)))
library(dplyr)
dmap <- dmap %>% group_by(NAME) %>% filter(row_number()==1) %>% tibble::column_to_rownames(var="NAME")
counts <- convertGeneName(dmap[,2:212],to="hgnc_symbol")
counts <- convertGeneName(counts,to="entrez")
counts <- convertGeneName(counts,to="refseq_mrna")
counts <- convertGeneName(counts,to="ensembl_gene")

#test convert transcript
transcripts <- unique(gene_id_converter_hs$ensembl_transcript)[1:5000]
dat <- matrix(1,nrow=5000,ncol=3)
rownames(dat) <- transcripts
convertGeneName(dat)

counts <- convertGeneName(dmap[,-1],to="refseq_mrna")
counts <- convertGeneName(counts,to="ensembl_gene_id")

# test removing version number
count_blue <- read.table("data_blue_counts.txt",header=T)
rownames(count_blue) <- count_blue$gene
count_blue <- convertGeneName(count_blue)

#---------------------------------------------
# test mapGene()
#---------------------------------------------------
df <- mapGene(dmap$Description,to="entrez")
df <- mapGene(dmap$Description,to="hgnc_symbol")
df <- mapGene(dmap$Description,to="ensembl_gene")
df <- mapGene(dmap$Description,to="refseq_mrna")

#-------------------------------------------
# test findDEGenes()
#-------------------------------------------
DEgenes <- findDEGenes(eset_dmap,0.05,2)

#-------------------------------------------
# user interface
# test buildScaffold()
#-------------------------------------------
library(leukemiasEset)
data("leukemiasEset")
g1 <- buildScaffold(exprs(leukemiasEset),pData(leukemiasEset),"LeukemiaType",plot_mode ="dot",pcs=c(1,2))
g2 <- buildScaffold(exprs_dmap,pData_dmap,"cell_types",plot_mode="tiny_label",pcs=c(1,2))
g3 <- buildScaffold(exprs_ilaria,pData_ilaria,"cancer_type",pcs=c(1,2))
DMAP_scaffold <- buildScaffold(exprs_dmap,pData_dmap,"cell_types")
# usethis::use_data(DMAP_scaffold)
# test use of prebuilt

g <- buildScaffold("prebuilt_DMAP")
g <- buildScaffold("prebuilt_DMAP",plot_mode="tiny_label",pcs=c(3,4))

#---------------------------------------------------------
# user interface
# test loadingPlot
#--------------------------------------------------------
loadingPlot(g)
loadingPlot(g2)
loadingPlot(g3)
loadingPlot(g3,gene_name ="ensembl_gene_id")
loadingPlot(g,num_genes = 1)

#----------------------------------------------------------
# user interface
# test projectSample()
#---------------------------------------------------------
projectSample(g,exprs_ilaria,pData_ilaria,"cancer_type") # really nice alignment with AML
projectSample(g2,exprs_ilaria,pData_ilaria,"cancer_type")
projectSample(g3,exprs_dmap,pData_dmap,"cell_types")

g <- buildScaffold("prebuilt_DMAP")
projectSample(g,exprs(leukemiasEset),pData(leukemiasEset),"LeukemiaType")

# test projectSample without phenotype table
projectSample(g,exprs(leukemiasEset))


#

# GTEX
plot_samples(map_samples())



