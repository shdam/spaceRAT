#' Returns either a prebuilt scaffoldSpace object,
#' or a new one calculated from input expression matrix
#'
#' This function returns a \code{\link{scaffoldSpace-class}}
#' object that contains all parameters needed for plotting
#' the PCA transformed sample space.
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
#' By default this function plots the resulting
#' \code{\link{scaffoldSpace-class}} object before returning it,
#' but this autoplot mode can be turned off by specifying
#' \code{auto_plot=FALSE}.
#'
#'
#' @param object An expression matrix of class
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
#' @param pheno_scaffold Phenotype table corresponding to
#' the expression matrix.
#' Row names are sample names, identical to column names of \code{object}.
#'
#' @param colname A column name of \code{pheno_scaffold}.
#' Cells will be grouped by this column of values.
#' Thus, differential expression analysis will be performed using this
#' column of phenotype as independent variables.
#'
#' @param data A character indicating whether the count matrix is
#' log-transformed or raw. By default \code{data="logged"}.
#' If the count matrix is raw counts, please specify \code{data="raw"}
#'
#' @param assay (Default: "counts") The assay slot to use with your
#' Bioconductor object.
#' @param dims A numeric vector containing 2 numbers, indicating
#' which two principle components to plot.
#' @param plot_mode A character indicating whether to add tiny
#' labels to each data point.
#' By default \code{plot_mode="dot"} and tiny labels will not be attached.
#' If more than 12 cell types are to be displayed, setting
#' \code{plot_mode="tiny_label"} may yield better visualization.
#' Shorter names for phenotypes (e.g. cell types) is strongly recommended
#' in "tiny_label" mode.
#' @param dim_reduction A character indicating the method for
#' dimensionality reduction. Currently "PCA" and "UMAP" are supported.
#' @param threshold (Default: 10) Prefiltering threshold for count row sums.
#' @param pval_cutoff A cutoff value for p value when selecting
#' differentially expressed genes. By default \code{pval_cutoff=0.05}.
#' @param lfc_cutoff A cutoff value for logFC when selecting
#' differentially expressed genes. By default \code{lfc_cutoff=2}.
#' @param title Title for the plot
#' @param pca_scale A logical variable determining whether to
#' normalize rows when plotting PCA
#' @param auto_plot A logical variable determining whether
#' to plot the resulting scaffold space when calling the function.
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
#'     pheno_scaffold = NULL,
#'     colname = NULL,
#'     assay = "counts",
#'     data = "logged",
#'     threshold = 10,
#'     dim_reduction = c("PCA", "UMAP"),
#'     dims = c(1, 2),
#'     plot_mode = "dot",
#'     classes = NULL,
#'     pval_cutoff = 0.05,
#'     lfc_cutoff = 2,
#'     title = "Scaffold Plot",
#'     pca_scale = FALSE,
#'     auto_plot = TRUE,
#'     annotation = "ensembl_gene"
#'     )
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom uwot umap
#' @importFrom spaceRATScaffolds listScaffolds
#' @export
#' @return A scaffoldSpace object
#' @examples
#' utils::data("exprs_dmap", "pData_dmap", package = "spaceRATScaffolds")
#' buildScaffold(exprs_dmap,pData_dmap,"cell_types",
#' pval_cutoff=0.01,pca_scale=TRUE, auto_plot = FALSE)
buildScaffold <- function(
        object,
        pheno_scaffold = NULL,
        colname = NULL,
        assay = "counts",
        data = "logged",
        threshold = 10,
        dim_reduction = c("PCA", "UMAP"),
        dims = c(1,2),
        plot_mode = "dot",
        classes = NULL,
        pval_cutoff = 0.05,
        lfc_cutoff = 2,
        title = "Scaffold Plot",
        pca_scale = FALSE,
        auto_plot = TRUE,
        annotation = "ensembl_gene"){

    dim_reduction <- match.arg(dim_reduction)
    # prebuilt_DMAP no samples removed
    if(
        is(object, "character") &&
        object %in% listScaffolds() &&
        is(classes, "NULL")
        ){
        space <- loadData(object)
        space@dims <- dims
        space@plot_mode <- plot_mode
        return(space)
        # prebuilt DMAP samples removed
    } else if(
        is(object, "character") &&
        object == "DMAP_scaffold" &&
        !is(classes, "NULL")
        ){
        object <- loadData("exprs_dmap")
        pheno_scaffold <- loadData("pData_dmap")
        colname <- "cell_types"
    } else if(is(object, "character")){
        stop(
        "Incorrectly named prebuilt scaffold. The available are: ",
        listScaffolds()
        )
    }

    # Preprocessing ----

    if(is(pheno_scaffold, "NULL")) {
        warning("No annotation data provided.
                Expression data colnames are used instead.")
        pheno_scaffold <- data.frame(colnames(object))
        colname <- NULL
    }

    object <- preprocess(
        object,
        colname = colname,
        pheno = pheno_scaffold,
        assay = assay,
        data = data,
        annotation = annotation
        )

    counts_scaffold <- object[[1]]
    cell_types <- object[[2]]
    rm(object)

    # Find DE genes
    DEgenes <- findDEGenes(
        counts_scaffold, cell_types, pval_cutoff, lfc_cutoff)

    # subset
    counts_scaffold <- counts_scaffold[DEgenes, ]

    # rank
    scaffold_rank <- apply(counts_scaffold, 2, rank)

    # dimension reduction
    message("Reducing dimensions.")
    if (dim_reduction == "PCA"){
        reduced_dims <- stats::prcomp(t(scaffold_rank), scale = pca_scale)
    } else if (dim_reduction == "UMAP"){
        reduced_dims <- uwot::umap(t(scaffold_rank))
    }


    # record standard data in scaffoldSpace class
    space <- methods::new(
        "scaffoldSpace",
        DEgene = DEgenes,
        label = as.character(cell_types),
        model = reduced_dims,
        dims = dims,
        plot_mode = plot_mode
        )
    message("Done.")

    return(space)
}


