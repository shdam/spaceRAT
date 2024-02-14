#' Project new sample(s) onto existing scaffold PCA plot
#'
#' This function takes in a scaffold space objects, subsets
#' the new dataset to contain only the genes that defines the scaffold space,
#' ranks the subset within new sample(s),
#' finally projects the new sample(s) onto the existing scaffold PCA plot.
#'
#' @inheritParams buildScaffold
#' @inheritParams plotScaffold
#' @param scaffold A scaffold space object, returned by function
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
#'     scaffold, sample,
#'     pheno = NULL,
#'     colname = NULL,
#'     assay = NULL,
#'     dimred = "PCA",
#'     dims = c(1, 2),
#'     plot_mode = "dot",
#'     title = "Samples projected onto scaffold PCA",
#'     verbose = TRUE,
#'     annotation = "ensembl_gene"
#'     )
#' @export
#' @importFrom stats predict
#' @importFrom uwot umap_transform
#' @import ggplot2
#' @return ggplot object with new samples projected to existing scaffold plot
#' @examples
#' utils::data("ilaria_counts", "ilaria_pData", package = "spaceRATScaffolds")
#' scaffold <- buildScaffold("DMAP")
#' # or create your own scaffoldSpace
#' projectSample(scaffold,ilaria_counts,ilaria_pData,"cancer_type")
projectSample <- function(
        scaffold, sample,
        pheno = NULL,
        colname = NULL,
        assay = NULL,
        dimred = "PCA",
        pca_scale = TRUE,
        subset_intersection = FALSE,
        dims = c(1, 2),
        plot_mode = "dot",
        title = "Samples projected onto scaffold PCA",
        verbose = TRUE,
        annotation = "ensembl_gene"){
    if(is(scaffold, "NULL") | is(sample, "NULL")){
        warning("No scaffold or sample provided.")
        return(NULL)}
    # Preprocess data
    sample <- preprocess(
        sample,
        colname = colname,
        pheno = pheno,
        assay = assay,
        annotation = annotation,
        threshold = NULL
        )
    pheno <- sample$pheno
    sample <- sample$mat

    # add absent genes then subset eset_sample so counts_scaffold
    # and counts_project contain same genes
    scaffold_genes <- unique(unlist(scaffold$DEgenes))
    sample_genes <- rownames(sample)
    overlap_genes <- intersect(scaffold_genes, sample_genes)#all_genes[(all_genes %in% rownames(sample))]
    absent_genes <- scaffold_genes[!(scaffold_genes %in% sample_genes)]
    absent_counts <- matrix(
        0,length(absent_genes), ncol(sample),
        dimnames = list(absent_genes, colnames(sample)))
    rownames(absent_counts) <- absent_genes
    sample <- sample[overlap_genes, ]
    if (!subset_intersection) sample <- rbind(sample, absent_counts)

    if (verbose & !subset_intersection){
        message(
        length(absent_genes)," genes are added to count matrix ",
        "with imputed expression level 0.")
    }

    if (subset_intersection & (length(absent_genes)/length(scaffold$DEgenes)>1/4) ){
        warning("More than 1/4 genes are added to sample with imputed ",
                "expression level 0!")
    }

    # Subset sample to DEgenes
    #sample <- sample[scaffold$DEgenes, ]
    
    ranked_sample <- apply(sample, 2, rank)
    
    # Rebuild scaffold
    if (subset_intersection){
        scaffold$rank <- apply(scaffold$rank[overlap_genes, ], 2, rank)
        scaffold$pca <- stats::prcomp(t(scaffold$rank), scale. = pca_scale)
    } else{
        # rank and transform sample and multiply with the percent of
        # missing values, to retain comparable numeric range
        ranked_sample <- ranked_sample *
            (1 + (length(absent_genes) / length(scaffold$DEgenes)))
    }

    
    # ranked_sample <- apply(sample, 2, rank) #*
        #(1+(length(absent_genes)/length(scaffold$DEgenes))); rm(sample)

    if (toupper(dimred) == "PCA"){
        # PCA transform the sample data
        transformed_sample <- stats::predict(
            scaffold$pca, newdata = t(ranked_sample))
    } else if (toupper(dimred) == "UMAP"){
        transformed_sample <- uwot::umap_transform(
            t(ranked_sample), scaffold$umap)
    }
    rm(ranked_sample)
    # Prepare dataframe for ggplot
    df_sample <- data.frame(
        "Dim1" = transformed_sample[, dims[1]],
        "Dim2" = transformed_sample[, dims[2]],
        "shape" = "19", "Scaffold_group" = "New_samples")
    rm(transformed_sample)

    graph <- plotScaffold(
        scaffold = scaffold, title = title,
        dimred = dimred, plot_mode = plot_mode,
        dims = dims)


    if (is(pheno, "NULL")){
        shapes <- 19
        # Hide shape from legend
        graph <- graph + ggplot2::guides(shape = "none")
    } else{
        df_sample$shape <- pheno[, colname]
        shapes <- seq_along(unique(df_sample$shape))
    }
    # calculate correct color scale
    p <- ggplot2::ggplot_build(graph)$data[[2]]
    cols <- unique(p[["colour"]])
    label <- as.character(p[["label"]])
    cols <- append(cols, "#000000")

    # Add projection to scaffold plot
    suppressMessages(
        g <- graph +
            ggplot2::geom_point(
                data = df_sample,
                mapping = ggplot2::aes(
                    x = .data$Dim1,
                    y = .data$Dim2,
                    color = .data$Scaffold_group,
                    shape = .data$shape)) +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::scale_shape_manual(values = shapes) +
            ggplot2::coord_fixed() +
            ggplot2::labs(
                shape = colname,
                title = title)
        )

    return(g)
}
