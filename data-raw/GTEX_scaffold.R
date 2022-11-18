## code to prepare `GTEX_scaffold` dataset goes here


counts_scaffold <- readr::read_csv("~/Downloads/Archive/gtex_representativeSetLog2cpm.csv")
pheno_scaffold <- readr::read_csv("~/Downloads/Archive/gtex_representativeSet_metadata.csv")

colname <- "gtex.smts"

counts_scaffold <- tibble::column_to_rownames(counts_scaffold, colnames(counts_scaffold)[1])
pheno_scaffold <- tibble::column_to_rownames(pheno_scaffold, colnames(pheno_scaffold)[1])

# remove genes with total count<10
idx <- which(rowSums(2**(counts_scaffold))<10)
counts_scaffold <- counts_scaffold[-idx,]

eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname, to="ensembl_gene")

pval_cutoff = 0.05
lfc_cutoff = 2

# find DE genes
DEgenes <- findDEGenes(eset_scaffold,pval_cutoff,lfc_cutoff)


usethis::use_data(GTEX_scaffold, overwrite = TRUE)


# TCGA

tcga_counts_scaffold <- readr::read_csv("~/Downloads/Archive/tcga_representativeSetLog2cpm.csv")
tcga_pheno_scaffold <- readr::read_csv("~/Downloads/Archive/tcga_representativeSet_metadata.csv")

colname <- "study"#"tcga.gdc_cases.project.primary_site" # tcga.gdc_cases.project.name

tcga_counts_scaffold <- tibble::column_to_rownames(tcga_counts_scaffold, colnames(tcga_counts_scaffold)[1])
tcga_pheno_scaffold <- tibble::column_to_rownames(tcga_pheno_scaffold, colnames(tcga_pheno_scaffold)[1])

# remove genes with total count<10
idx <- which(rowSums(2**(tcga_counts_scaffold))<10)
tcga_counts_scaffold <- tcga_counts_scaffold[-idx,]

tcga_eset_scaffold <- createEset(tcga_counts_scaffold,tcga_pheno_scaffold,colname, to="ensembl_gene")

pval_cutoff = 0.05
lfc_cutoff = 2

# find DE genes
tcga_DEgenes <- findDEGenes(tcga_eset_scaffold,pval_cutoff,lfc_cutoff)

tcga_eset_scaffold <- tcga_eset_scaffold[tcga_DEgenes,]

# rank
tcga_scaffold_rank<- apply(Biobase::exprs(tcga_eset_scaffold),2,rank)

# PCA
tcga_pca <- stats::prcomp(t(tcga_scaffold_rank),scale=pca_scale)

# record standard data in scaffoldSpace class
tcga_space <- methods::new("scaffoldSpace",
                      DEgene=tcga_DEgenes,
                      label=as.character(Biobase::pData(tcga_eset_scaffold)[,1]),
                      pca=tcga_pca,
                      pcs=pcs,
                      plot_mode=plot_mode)

plotScaffold(tcga_space,title)
