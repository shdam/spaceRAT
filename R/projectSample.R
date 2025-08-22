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
#' @param annotation (DEPRECATED) Please set the annotation when building the
#' scaffold
#' @param subset_intersection (Default: FALSE) Setting this to TRUE will
#' rebuild the scaffold on the overlapping genes between scaffold and sample
#' @param sample_name Name of new sample label in legend. Default "New_samples"
#' @param title Title of the plot.
#' @param verbose A logical vector indicating whether to report the number of
#' genes imputed to make \code{sample} compatible with
#' \code{scaffold}
#' @export
#' @importFrom stats predict
#' @import ggplot2
#' @return ggplot object with new samples projected to existing scaffold plot
#' @examples
#' utils::data("ilaria_counts", "ilaria_pData", package = "spaceRATScaffolds")
#' scaffold <- buildScaffold("DMAP")
#' # or create your own scaffoldSpace
#' p <- projectSample(scaffold,ilaria_counts,ilaria_pData,"cancer_type")
projectSample <- function(
        scaffold, sample,
        pheno = NULL,
        colname = "Classes",
        assay = NULL,
        dimred = "PCA",
        classes = NULL,
        subset_intersection = FALSE,
        dims = c(1, 2),
        sample_name = "New samples",
        plot_mode = "dot",
        title = "Samples projected onto scaffold PCA",
        verbose = TRUE,
        annotation = "ensembl_gene"){

    if (is(scaffold, "NULL") | is(sample, "NULL")){
        warning("No scaffold or sample provided.")
        return(NULL)
        }
    if (is(scaffold$annotation, "NULL")) {
        scaffold$annotation <- annotation
    }
    if (toupper(dimred) == "UMAP" && !requireNamespace("uwot")) {
        stop("To use UMAP space, please install uwot:\n",
             "install.packages(\"uwot\")")
    }

    # Preprocess data
    sample <- spaceRAT:::preprocess(
        sample,
        colname = colname,
        pheno = pheno,
        assay = assay,
        annotation = scaffold$annotation,
        threshold = NULL
        )
    pheno <- sample$pheno
    sample <- sample$mat

    if (is(scaffold$rank_scale, "NULL")) scaffold$rank_scale <- FALSE
    if (is(scaffold$group_scale, "NULL")) scaffold$group_scale <- FALSE

    # Subset to classes
    if(!is(classes, "NULL")){
        scaffold$DEgenes <- scaffold$DEgenes[classes]
        keep <- scaffold$label %in% classes
        scaffold$label <- scaffold$label[keep]
        scaffold$rank <- scaffold$rank[, keep]
    }

    # add absent genes then subset eset_sample so counts_scaffold
    # and counts_project contain same genes
    scaffold_genes <- unique(unlist(scaffold$DEgenes))
    sample_genes <- rownames(sample)
    overlap_genes <- intersect(scaffold_genes, sample_genes)

    sample <- sample[overlap_genes, ]
  
    
    # if (!subset_intersection) {
    #     absent_genes <- scaffold_genes[!(scaffold_genes %in% sample_genes)]
    #     absent_counts <- matrix(
    #         0,length(absent_genes), ncol(sample),
    #         dimnames = list(absent_genes, colnames(sample)))
    #     rownames(absent_counts) <- absent_genes
    #     sample <- rbind(sample, absent_counts)
    #     # sample <- sample[scaffold_genes, ] # get genes in the same order
    # 
    #     if (verbose) message(
    #         length(absent_genes)," genes are added to count matrix ",
    #         "with imputed expression level 0.")
    #     if ((length(absent_genes)/length(scaffold_genes))>1/4) warning(
    #         "More than 1/4 genes are added to sample with imputed ",
    #         "expression level 0!")
    # }


    if (!is.null(names(scaffold$DEgenes))) {
        sample <- lapply(names(scaffold$DEgenes), function(group) {
          
            group_genes <- scaffold$DEgenes[[group]]
            group_genes <- group_genes[group_genes %in% overlap_genes]

            if(length(group_genes)>0){
                if (scaffold$group_scale) sample <- as.matrix(sample[group_genes, ] / scaffold$eigenvalues[[group]][1])
                else sample <- as.matrix(sample[group_genes, ])
                rownames(sample) <- paste(group_genes, group, sep = "_")
            } else{
                # print(group)
                sample <- as.matrix(sample)
                rownames(sample) <- paste(rownames(sample), group, sep = "_")

            }
            sample[rownames(sample) %in% rownames(scaffold$pca$rotation), ]
        })
        sample <- do.call(rbind, sample)
    }
    
    # Moved from before the group scaling part 
    
    if (!subset_intersection) {
        pca_genes <- rownames(scaffold$pca$rotation)
        sample_genes <- rownames(sample)
        absent_genes <- pca_genes[!(pca_genes %in% sample_genes)]
        absent_counts <- matrix(
          0,length(absent_genes), ncol(sample),
          dimnames = list(absent_genes, colnames(sample)))
        rownames(absent_counts) <- absent_genes
        sample <- rbind(sample, absent_counts)
        sample <- sample[pca_genes, ] # get genes in the same order

        if (verbose) message(
            length(absent_genes)," genes are added to count matrix ",
            "with imputed expression level 0.")
        if ((length(absent_genes)/length(scaffold_genes))>1/4) warning(
            "More than 1/4 genes are added to sample with imputed ",
            "expression level 0!")
    }
    
    # Move end
   
    # Subset sample to DEgenes
    #sample <- sample[scaffold$DEgenes, ]
    ranked_sample <- ranking(sample, rank_scale = scaffold$rank_scale)

    # Rebuild scaffold
    if (subset_intersection & !is.null(scaffold$rank)){
        scaffold$rank <- ranking(scaffold$rank[rownames(sample), ], rank_scale = scaffold$rank_scale)
        # scaffold$pca <- stats::prcomp(t(scaffold$rank), scale. = scaffold$pca$scale) # OLD
        scaffold$pca <- stats::prcomp(t(scaffold$rank), scale. = scaffold$pca$scale[rownames(sample)]) # HW
    } else if (!scaffold$rank_scale){ # or scale ranks
        ranked_sample <- ranked_sample *
            (1 + (length(absent_genes) / length(scaffold$DEgenes)))
    }
    if (subset_intersection & is.null(scaffold$rank)) {
        # Subset PCA if cannot recompute
        scaffold$pca$rotation <- scaffold$pca$rotation[overlap_genes,]
        scaffold$pca$center <- scaffold$pca$center[overlap_genes]
    }

    # Transform the sample data
    if (toupper(dimred) == "PCA"){
        transformed_sample <- stats::predict(
            scaffold$pca, newdata = t(ranked_sample))
    } else if (toupper(dimred) == "UMAP"){
        scaffold$umap <- uwot::umap(t(scaffold$rank), ret_model = TRUE)
        transformed_sample <- uwot::umap_transform(
            t(ranked_sample), scaffold$umap)
    }

    # Prepare dataframe for ggplot
    transformed_sample <- as.data.frame(transformed_sample)
    transformed_sample$shape <- "19"
    transformed_sample$scaffoldGroup <- sample_name
    dim1 <- colnames(transformed_sample)[dims[1]]
    dim2 <- colnames(transformed_sample)[dims[2]]

    # Create scaffold plot
    projection <- plotScaffold(
        scaffold = scaffold, title = title,
        dimred = dimred, plot_mode = plot_mode,
        class_name = colname,
        dims = dims)

    if (is(pheno, "NULL")){
        shapes <- 19
        # Hide shape from legend
        projection <- projection + ggplot2::guides(shape = "none")
    } else{
      transformed_sample$shape <- pheno[, colname]
        shapes <- seq_along(unique(transformed_sample$shape))
    }

    # calculate correct color scale
    p <- ggplot2::ggplot_build(projection)$data[[2]]
    cols <- unique(p[["colour"]])
    label <- as.character(p[["label"]])
    cols <- append(cols, "#000000")

    # Add projection to scaffold plot
    suppressMessages(
        projection <- projection +
            ggplot2::geom_point(
                data = transformed_sample,
                mapping = ggplot2::aes(
                    x = .data[[dim1]],
                    y = .data[[dim2]],
                    color = .data$scaffoldGroup,
                    shape = .data$shape)) +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::scale_shape_manual(values = shapes) +
            ggplot2::coord_fixed() +
            ggplot2::labs(
                shape = colname,
                color = colname,
                title = title)
        )

    # Add sample data to projection output
    projection$sample <- transformed_sample
    projection$scaffold <- scaffold
    projection$rankedSample <- ranked_sample
    projection$pheno <- pheno
    projection$colnames <- colname


    return(projection)
}
