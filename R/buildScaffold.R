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
#' @param pval_cutoff A cutoff value for p value when selecting
#' differentially expressed genes. By default \code{pval_cutoff=0.05}.
#' @param lfc_cutoff A cutoff value for logFC when selecting
#' differentially expressed genes. By default \code{lfc_cutoff=2}.
#' @param pca_scale A logical variable determining whether to
#' normalize rows when plotting PCA
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
#' @usage
#' buildScaffold(
#'     object,
#'     pheno = NULL,
#'     colname = NULL,
#'     assay = NULL,
#'     data = NULL,
#'     subset_deg = TRUE,
#'     threshold = 10,
#'     add_umap = FALSE,
#'     classes = NULL,
#'     pval_cutoff = 0.05,
#'     lfc_cutoff = 2,
#'     pca_scale = FALSE,
#'     annotation = "ensembl_gene"
#'     )
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom uwot umap
#' @importFrom spaceRATScaffolds listScaffolds
#' @export
#' @return A scaffold space
#' @examples
#' utils::data("DMAP_exprs", "DMAP_pData", package = "spaceRATScaffolds")
#' buildScaffold(DMAP_exprs,DMAP_pData,"cell_types", data = "exprs",
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
        ranking = TRUE,
        pval_cutoff = 0.05,
        lfc_cutoff = 2,
        pca_scale = TRUE,
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
    if(!is(pheno, "NULL") & is(colname, "NULL")) stop(
        "Please specify colname for pheno data")

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
    pheno <- object$pheno; mat <- object$mat; rm(object)

    # Define scaffold space
    scaffold <- list("label" = pheno[, colname])

    # Find DE genes
    if(subset_deg){
        message("Finding differentially expressed genes")
        # Remove cell types with less than 2 cells
        count_table <- data.frame(table(scaffold$label))
        keep_labels <- count_table[count_table$Freq > 1, 1]
        stopifnot(
            "All cells have unique phenotype information.
            Please group cells by phenotype." = length(keep_labels) > 0)
        if(!all(scaffold$label %in% keep_labels)){
            warning(
                "The samples:",
                scaffold$label[!(scaffold$label %in% keep_labels)],
                "were removed because there are too few.")
        }
        keep <- scaffold$label %in% keep_labels
        mat <- mat[, keep]
        scaffold$label <- scaffold$label[keep]

        scaffold$DEgenes <- lapply(unique(scaffold$label), function(group){
            findDEGenes2(
                mat = mat, group = group, labels = scaffold$label,
                pval_cutoff = pval_cutoff, lfc_cutoff = lfc_cutoff)
        })
        names(scaffold$DEgenes) <- unique(scaffold$label)

        mat <- mat[unique(unlist(scaffold$DEgenes)), ]
    }

    # Iterate over groups, calculate PCA, and store eigenvalues
    mat <- lapply(names(scaffold$DEgenes), function(group) {
        group_genes <- scaffold$DEgenes[[group]]
        if(length(group_genes) > 0) {
            pca <- prcomp(t(mat[group_genes, ]), scale. = TRUE)
            eigenvalues <- pca$sdev^2
            scaled_mat <- as.matrix(mat[group_genes, ] / eigenvalues[1])
            rownames(scaled_mat) <- paste(group_genes, group, sep = "_")
            return(list(scaled_mat = scaled_mat, eigenvalues = eigenvalues))
        } else {
            scaled_mat <- as.matrix(mat)
            rownames(scaled_mat) <- paste(rownames(scaled_mat), group, sep = "_")
            return(list(scaled_mat = scaled_mat, eigenvalues = 1))
        }
    })

    # Extract eigenvalues and matrices from the results
    scaffold$eigenvalues <- sapply(mat, `[[`, "eigenvalues")
    names(scaffold$eigenvalues) <- names(scaffold$DEgenes)
    mat <- do.call(rbind, lapply(mat, `[[`, "scaled_mat"))


    # rank
    if(ranking) mat <- apply(mat, 2, rank)
    scaffold$rank <- mat
    # scaffold$DEgenes <- unique(unlist(scaffold$DEgenes))
    # dimension reduction
    message("Reducing dimensions.")
    scaffold$pca <- stats::prcomp(t(mat), scale. = TRUE)
    if (add_umap) scaffold$umap <- uwot::umap(t(mat), ret_model = TRUE)

    message("Scaffold is built.")
    return(scaffold)
}


