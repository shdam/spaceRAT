#' Project new sample(s) onto existing scaffold PCA plot
#'
#' This function takes in a \code{\link{scaffoldSpace-class}} objects, subsets
#' the new dataset to contain only the genes that defines the scaffold space,
#' ranks the subset within new sample(s),
#' finally projects the new sample(s) onto the existing scaffold PCA plot.
#'
#' @inheritParams buildScaffold
#' @param space A scaffoldSpace object, returned by function
#' \code{\link{buildScaffold}}
#' @param sample sample data to project. A `matrix`, `data.frame`, or
#' `SummarizedExperiment`
#' @param pheno Phenotype data corresponding to \code{sample}.
#' If not specified, the output plot will not show legends for new samples.
#' @param colname A column name of \code{pheno}.
#' This column of values will be used to annotate projected samples.
#' This argument should be set together with \code{pheno}.
#' If \code{pheno} is not specified, this argument will be ignored,
#' and the output plot will not show legends for new samples.
#' @param annotation Annotation type to use for scaffold.
#' counts_scaffold rownames using alternative identifiers will be translated.
#' Currently ensembl_gene, entrez, hgnc_symbol, and refseq_mrna are supported.
#' set to "NA", to avoid translation
#' (both scaffold and projected samples must be the same)
#' @param title Title of the plot.
#' @param verbose A logical vector indicating whether to report the number of
#' genes imputed to make \code{sample} compatible with
#' \code{scaffold}
#' @usage
#' projectSample(
#'     space,
#'     sample,
#'     pheno = NULL,
#'     colname = "cancer_type",
#'     assay = "counts",
#'     title = "Samples projected onto scaffold PCA",
#'     verbose = TRUE,
#'     annotation = "ensembl_gene"
#'     )
#' @export
#' @importFrom stats predict
#' @import ggplot2
#' @return ggplot object with new samples projected to existing scaffold plot
#' @examples
#' utils::data("DMAP_scaffold", "exprs_ilaria", "pData_ilaria",
#' package = "spaceRATScaffolds")
#' space <- DMAP_scaffold # or create your own scaffoldSpace
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")
projectSample <- function(
        space, sample,
        pheno = NULL,
        colname = NULL,
        assay = "counts",
        dims = c(1, 2),
        title = "Samples projected onto scaffold PCA",
        verbose = TRUE,
        annotation = "ensembl_gene"){
    # Preprocess data
    sample <- preprocess(
        sample,
        colname = colname,
        pheno = pheno,
        assay = assay,
        annotation = annotation,
        threshold = NULL
        )

    # add absent genes then subset eset_sample so counts_scaffold
    # and counts_project contain same genes
    absent_genes <- space$DEgene[!(space$DEgene %in% rownames(sample))]
    absent_counts <- matrix(
        0,length(absent_genes), ncol(sample),
        dimnames = list(absent_genes, colnames(sample)))
    rownames(absent_counts) <- absent_genes
    counts_sample <- rbind(assay(sample, assay), absent_counts)

    if (verbose){
        message(
        length(absent_genes)," genes are added to count matrix ",
        "with imputed expression level 0.")
    }

    if (length(absent_genes)/length(space$DEgene)>1/4){
        warning("More than 1/4 genes are added to sample with imputed ",
                "expression level 0!")
    }

    # Subset sample to DEgenes
    counts_sample <- counts_sample[space$DEgene, ]


    # rank and transform exprs_project and multiply with the percent of
    # missing values, to retain comparable numeric range
    ranked_sample <- apply(counts_sample, 2, rank) *
        (1+(length(absent_genes)/length(space$DEgene)))

    # PCA transform the sample data

    transformed_sample <- stats::predict(
        space$pca, newdata = t(ranked_sample))

    # Prepare dataframe for ggplot
    PC1_sample <- transformed_sample[, dims[1]]
    PC2_sample <- transformed_sample[, dims[2]]

    graph <- plotScaffold(space, title = title)

    if (ncol(colData(sample)) == 0){
        # calculate correct color scale
        p <- ggplot2::ggplot_build(graph)$data[[2]]
        cols <- p[["colour"]]
        label <- as.character(p[["label"]])
        idx <- rank(c(label, "New_samples"))[length(label)+1]
        cols <- append(cols,"#000000", after = idx-1)
        df_sample <- data.frame(PC1_sample, PC2_sample)

        suppressMessages(g <- graph +
            ggplot2::geom_point(data = df_sample,
                                mapping = ggplot2::aes(
                                    .data$PC1_sample,
                                    .data$PC2_sample,
                                    color = "New_samples")) +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::ggtitle(title) +
            ggplot2::coord_fixed())
        return (g)
        }

    New_samples <- colData(sample)[, colname]
    df_sample <- data.frame(PC1_sample, PC2_sample, New_samples)

    # project points
    suppressMessages(g <- graph +
        ggplot2::geom_point(data = df_sample,
                            mapping = ggplot2::aes(
                                .data$PC1_sample,
                                .data$PC2_sample,
                                shape = .data$New_samples),
                            color = "black") +
        ggplot2::scale_shape_manual(values = seq_along(unique(New_samples))) +
        ggplot2::ggtitle(title) +
        ggplot2::coord_fixed())
    return(g)
}
