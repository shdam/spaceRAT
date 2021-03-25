#' Find Differentially Expressed Genes
#'
#' This function performs DE analysis using limma to find differentially expressed genes between cell X and non-X for all types of cells that has at least 2 samples in the count matrix.
#'
#' @importFrom Biobase pData
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @param eset an ExpressionSet
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @export
#' @return A vector of names of differentially expressed genes
#' @examples
#' find_de_genes(eset_dmap,"cell_types")

findDEGenes <- function(eset,pval_cutoff=0.05,lfc_cutoff=2){
        # remove cell types with less than 2 cells
        count_table <- data.frame(table(Biobase::pData(eset)[,1]))
        cell_types <- count_table[count_table$Freq>1,]$Var1
        eset <- eset[,pData(eset)[,1] %in% cell_types]

        # create contrast matrix
        n <- length(cell_types)
        cm <- matrix(-1/(n-1),n,n)
        diag(cm) <- 1
        rownames(cm) <- cell_types
        colnames(cm) <- paste(cell_types,"others",sep="_")

        # limma
        cell_types <- Biobase::pData(eset)[,1]
        design <- model.matrix(~0+ cell_types)
        colnames(design) <- gsub("cell_types","",colnames(design))
        fit <- limma::lmFit(eset,design)
        fit <- limma::contrasts.fit(fit, contrast=cm)
        fit <- limma::eBayes(fit)
        de_genes <- sapply(colnames(cm), function(name){
                rownames(limma::topTable(fit, coef=name, number = Inf, lfc =lfc_cutoff, p.value = pval_cutoff, adjust.method="fdr"))
        })
        unique(unlist(de_genes))

}
