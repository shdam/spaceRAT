#' Returns either a prebuilt scaffoldSpace object,
#' or a new one calculated from input expression matrix
#'
#' This function returns a scaffold space as a list
#' object that contains all parameters needed for plotting
#' the PCA/UMAP transformed sample space.
#' To use prebuilt scaffold space, simply call:
#' \code{buildScaffold("scaffold_name")}. Get a list of scaffolds with:
#' \code{\link[spaceRATScaffolds:listScaffolds]{listScaffolds()}}
#'
#' To build your own scaffold space, pass as arguments a count matrix,
#' a phenotype table, an a column name of the phenotype table to the function.
#' Missing values are allowed and automatically processed by the function.
#' But please make sure that the all missing values are represented by NA,
#' not "", " "(blank space), "-" or any other character.
#'
#' When building a user-defined scaffold space,
#' this function first preprocesses
#' count matrix by removing genes with total counts less than 10.
#' Then it performs differential expression analysis to select
#' differentially expressed genes (DE genes),
#' subsets the scaffold dataset to contain only the DE genes,
#' ranks the genes within each sample, finally performs
#' principle component analysis (PCA) using ranks.
#'
#'
#'
#' @param object An count or expression matrix of class
#' \code{\link[SummarizedExperiment:SummarizedExperiment]{SummarizedExperiment}}
#' ,
#' \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
#' , \code{matrix}, or a \code{data.frame}.
#' Column names are sample names.
#' Row names are gene names, which can be Ensembl Gene ID, HGNC Symbol,
#' Entrez Gene ID Ensembl, Transcript ID or Refseq mRNA.
#' All gene ID will be automatically converted to Ensembl Gene ID.
#' If you want to retain your gene identifier that is not currently supported,
#' please specify \code{annotation=NA}. See parameter
#' \code{annotation} for more information.
#' Counts of several transcript ID corresponding to same gene will
#' be added and recorded as the counts of the gene.
#'
#' @inheritParams limma::topTable
#' @param pheno Phenotype table corresponding to
#' the expression matrix.
#' Row names are sample names, identical to column names of \code{object}.
#'
#' @param colname A column name of \code{pheno}.
#' Cells will be grouped by this column of values.
#' Thus, differential expression analysis will be performed using this
#' column of phenotype as independent variables.
#'
#' @param data A character indicating whether the matrix is
#' log-transformed or raw counts. By default \code{data="exprs"}.
#' If the count matrix is raw counts, please specify \code{data="counts"}
#' @param subset_deg (Default: TRUE) Inform whether the scaffold should be
#' built using all genes or only the differentially expressed ones in the input.
#' @param assay (Default: NULL) The assay slot to use with your
#' Bioconductor object.
#' @param threshold (Default: 10) Prefiltering threshold for count row sums.
#' @param n_genes (Default: Inf) The number of significant genes to extract
#' from limma output per group in the scaffold
#' @param pval_cutoff A cutoff value for p value when selecting
#' differentially expressed genes. By default \code{pval_cutoff=0.05}.
#' @param lfc_cutoff A cutoff value for logFC when selecting
#' differentially expressed genes. By default \code{lfc_cutoff=2}.
#' @param pca_scale A logical variable determining whether to
#' normalize rows when plotting PCA
#' @param rank_scale A logical determining if ranks should be scaled on min
#' @param add_umap Add a UMAP space to the scaffold
#' @param annotation Type of gene identifier to use for scaffold.
#' Currently "ensembl_gene", "ensembl_transcript", "entrez", "hgnc_symbol",
#' and "refseq_mrna" are supported.
#' @param classes Cell types to use in plot
#' For example, set \code{annotation="hgnc_symbol"} will convert the row names
#' (gene identifiers) of \code{counts_scaffold} to hgnc symbol,
#' so will the expression data and the resulting PCA scaffold.
#' If this attempted translation fails, or your desired gene identifier is
#' not supported (especially when you are analyzing non-human data),
#' please set \code{annotation="hgnc_symbol"} to avoid translation.
#' In this case, please manually make sure that the row names
#' (gene identifiers) of \code{counts_scaffold} and
#' \code{counts_sample} are the same.
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom spaceRATScaffolds listScaffolds
#' @importFrom utils packageVersion
#' @export
#' @return A scaffold space
#' @examples
#' utils::data("DMAP_exprs", "DMAP_pData", package = "spaceRATScaffolds")
#' scaffold <- buildScaffold(DMAP_exprs,DMAP_pData,"cell_types", data = "exprs",
#' pval_cutoff=0.01,pca_scale=TRUE)
buildScaffold <- function(
        object,
        pheno = NULL,
        colname = NULL,
        assay = NULL,
        data = NULL,
        subset_deg = TRUE,
        threshold = 10,
        add_umap = FALSE,
        classes = NULL,
        pval_cutoff = 0.05,
        lfc_cutoff = 2,
        pca_scale = TRUE,
        rank_scale = FALSE,
        n_genes = Inf, 
        sort.by = "B",
        annotation = "ensembl_gene"){

    # Check prebuilt
    if(is(object, "character")) {
        object <- checkPrebuilt(object, classes)
        if(!is(object$pca, "NULL")) return(object)
        else{ # Extract modified prebuilt scaffold
            pheno <- object$pheno
            colname <- object$colname
            object <- object$object
        }}

    stopifnot("Please ensure unique column names in data." = all(
        !duplicated(colnames(object))))
    if (!is(pheno, "NULL") & is(colname, "NULL")) stop(
        "Please specify colname for pheno data")

    if (add_umap && !requireNamespace("uwot")) {
      stop("To add UMAP, please install uwot:\n",
           "install.packages(\"uwot\")")
      }
    
    
    # Preprocessing
    object <- preprocess(
        object,
        colname = colname,
        pheno = pheno,
        assay = assay,
        data = data,
        annotation = annotation,
        classes = classes
        )
    pheno <- object$pheno
    mat <- object$mat
    rm(object)

    # Define scaffold space
    scaffold <- list(
      "label" = pheno[, colname],
      "annotation" = annotation)

    # Find DE genes
    if(subset_deg){
        scaffold$DEgenes <- findDEGenes(
            mat, scaffold$label,
            pval_cutoff = pval_cutoff, lfc_cutoff = lfc_cutoff,
            n_genes = n_genes, sort.by = sort.by)
        mat <- mat[unique(unlist(scaffold$DEgenes)), ]
    }

    # rank
    mat <- ranking(mat, rank_scale = rank_scale)
    scaffold$rank <- mat
    scaffold$rank_scale <- rank_scale

    # dimension reduction
    message("Reducing dimensions.")
    scaffold$pca <- stats::prcomp(t(mat), scale. = pca_scale)
    scaffold$pca_scale <- pca_scale
    if (add_umap) scaffold$umap <- uwot::umap(t(mat), ret_model = TRUE)

    scaffold$spaceRATVersion <- packageVersion("spaceRAT")
    message("Scaffold is built.")
    return(scaffold)
}


