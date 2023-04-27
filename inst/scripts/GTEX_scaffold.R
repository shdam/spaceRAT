## code to prepare `GTEX_scaffold` dataset goes here


counts_scaffold <- readr::read_csv("../spaceRAT_data/gtex_representativeSetLog2cpm.csv")
pheno_scaffold <- readr::read_csv("../spaceRAT_data/gtex_representativeSet_metadata.csv")

colname <- "gtex.smts"

# counts_scaffold <- tibble::column_to_rownames(counts_scaffold, colnames(counts_scaffold)[1])
# pheno_scaffold <- tibble::column_to_rownames(pheno_scaffold, colnames(pheno_scaffold)[1])
# sce <- SingleCellExperiment::SingleCellExperiment(assays=list("counts"=counts_scaffold), colData = S4Vectors::DataFrame(pheno_scaffold))

GTEX_scaffold <- buildScaffold(
    object = counts_scaffold,
    pheno_scaffold = pheno_scaffold,
    colname = colname,
    data = "logged",
    dims = c(1,2),
    dim_reduction = "PCA",
    plot_mode = "dot",
    classes = NULL,
    pval_cutoff = 0.05,
    lfc_cutoff = 2,
    title = "GTEX PCA scaffold",
    pca_scale = FALSE,
    auto_plot = FALSE,
    annotation = "ensembl_gene"
    )

plotScaffold(GTEX_scaffold,"GTEX PCA scaffold")

usethis::use_data(GTEX_scaffold, overwrite = TRUE)



# TCGA

tcga_counts_scaffold <- readr::read_csv("~/Downloads/Archive/tcga_representativeSetLog2cpm.csv")
tcga_pheno_scaffold <- readr::read_csv("~/Downloads/Archive/tcga_representativeSet_metadata.csv")

colname <- "study"#"tcga.gdc_cases.project.primary_site" # tcga.gdc_cases.project.name

TCGA_scaffold <- buildScaffold(
    object = tcga_counts_scaffold,
    pheno_scaffold = tcga_pheno_scaffold,
    colname = colname,
    data = "logged",
    dims = c(1,2),
    plot_mode = "dot",
    classes = NULL,
    pval_cutoff = 0.05,
    lfc_cutoff = 2,
    title = "TCGA PCA scaffold",
    pca_scale = FALSE,
    auto_plot = FALSE,
    annotation = "ensembl_gene"
    )

usethis::use_data(TCGA_scaffold, overwrite = TRUE)


plotScaffold(TCGA_scaffold,"TCGA PCA scaffold")
