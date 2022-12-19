#' Project new sample(s) onto existing scaffold PCA plot
#'
#' This function takes in a \code{\link{scaffoldSpace}} objects, subsets
#' the new dataset to contain only the genes that defines the scaffold space,
#' ranks the subset within new sample(s),
#' finally projects the new sample(s) onto the existing scaffold PCA plot.
#'
#' @inheritParams buildScaffold
#' @param space A scaffoldSpace object, returned by function \code{\link{buildScaffold}}
#' @param pheno Phenotype data corresponding to \code{object}. If not specified, the output plot will not show legends for new samples.
#' @param colname A column name of \code{pheno}. This column of values will be used to annotate projected samples.
#' This argument should be set together with \code{pheno}.
#' If \code{pheno} is not specified, this argument will be ignored, and the output plot will not show legends for new samples.
#' @param annotation Annotation type to use for scaffold. counts_scaffold rownames using alternative identifiers will be translated.
#' Currently ensembl_gene, entrez, hgnc_symbol, and refseq_mrna are supported. set to "NA", to avoid translation (both scaffold and projected samples must be the same)
#' @param title Title of the plot.
#' @param verbose A logical vector indicating whether to report the number of genes imputed to make \code{counts_sample} compatible with \code{counts_scaffold}
#' @export
#' @return A ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' utils::data("DMAP_scaffold", "exprs_ilaria", "pData_ilaria", package = "spaceRAT")
#' space <- DMAP_scaffold # or create your own scaffoldSpace
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")


projectSample <- function(space,
                          object,
                          pheno = NULL,
                          colname = "cancer_type",
                          assay = "counts",
                          title = "Samples projected onto scaffold PCA",
                          verbose = TRUE,
                          annotation = "ensembl_gene"){


    # Preprocess data
    object <- preprocess(
        object,
        colname = colname,
        pheno = pheno,
        assay = assay,
        annotation = annotation,
        threshold = NULL
        )

    counts <- object[[1]]
    cell_types <- object[[2]]
    rm(object)


    # add absent genes then subset eset_sample so counts_scaffold and counts_project contain same genes
    absent_genes <- space@DEgene[!(space@DEgene %in% rownames(counts))]
    absent_exprs <- matrix(0,length(absent_genes), ncol(counts), dimnames = list(absent_genes, colnames(counts)))
    rownames(absent_exprs) <- absent_genes
    counts_sample <- rbind(counts,absent_exprs)

    if (verbose){
        message(paste0(length(absent_genes)," genes are added to count matrix with imputed expression level 0."))
    }

    if (length(absent_genes)/length(space@DEgene)>1/4){
        warning("More than 1/4 genes are added to sample with imputed expression level 0!")
    }

    #subset
    counts_sample <- counts_sample[space@DEgene, ]


    # rank and transform exprs_project and multiply with the percent of missing values, to retain comparable numeric range
    ranked_sample <- apply(counts_sample, 2, rank)*(1+(length(absent_genes)/length(space@DEgene)))

    # PCA transform the sample data

    transformed_sample <- stats::predict(space@model, newdata = t(ranked_sample))

    # Prepare dataframe for ggplot
    PC1_sample <- transformed_sample[, space@dims[1]]
    PC2_sample <- transformed_sample[, space@dims[2]]

    graph <- plotScaffold(space, title = title)
    if (is.null(cell_types)){
        # calculate correct color scale
        p <- ggplot2::ggplot_build(graph)$data[[2]]
        cols <- p[["colour"]]
        label <- as.character(p[["label"]])
        idx <- rank(c(label, "New_samples"))[length(label)+1]
        cols <- append(cols,"#000000", after = idx-1)
        df_sample <- data.frame(PC1_sample, PC2_sample)

        g <- graph +
            ggplot2::geom_point(data = df_sample,
                                mapping = ggplot2::aes(PC1_sample, PC2_sample, color = "New_samples")) +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::ggtitle(title) +
            ggplot2::coord_fixed()
    return (g)
        }

    New_samples <- cell_types
    df_sample <- data.frame(PC1_sample, PC2_sample, New_samples)

    # project points
    g <- graph +
        ggplot2::geom_point(data = df_sample,
                            mapping = ggplot2::aes(PC1_sample, PC2_sample, shape = New_samples),
                            color = "black") +
        ggplot2::scale_shape_manual(values = seq_along(unique(New_samples))) +
        ggplot2::ggtitle(title) +
        ggplot2::coord_fixed()
    return(g)
}
