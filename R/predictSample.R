#' Measures the distances to the centroids of existing scaffold labels to new sample(s)
#'
#' This function takes in a scaffold space objects, subsets
#' the new dataset to contain only the genes that defines the scaffold space,
#' ranks the subset within new sample(s),
#' finally projects the new sample(s) onto the existing scaffold PCA plot.
#'
#' @inheritParams buildScaffold
#' @inheritParams plotScaffold
#' @inheritParams projectSample
#' @param projection A scaffold space object, returned by function
#' \code{\link{projectSample}}. Recommended to have run with \code{subset_interaction = TRUE}
#' @param centroid_calculation Calculate centroids by \code{'mean'} or \code{'median'} of class expressions.
#' @param temperature Parameter for softmax function used transform distance to probability. 
#' Larger temperature → smoother probabilities.
#' Smaller temperature → peakier, closer to winner-takes-all.
#' @param title Title of the plot.
#' @param display_numbers Display numbers in tiles of pheatmap (see \code{\link{pheatmap::pheatmap}}).  
#' @param exclude_small_prob Logical. If TRUE, labels with low overall 
#'   probabilities will be excluded from the results. Requires 
#'   \code{min_prob_sum} to be specified.
#' @param min_prob_sum Numeric in the range [0, 1]. Threshold for excluding 
#'   labels when \code{exclude_small_prob = TRUE}. Labels with a 
#'   total summed probability across all samples less than this value will 
#'   be removed. Defaults to \code{NULL}, and must be provided if 
#'   \code{exclude_small_prob = TRUE}.
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
predictSample <- function(
        projection, 
        # sample,
        # pheno = NULL,
        # colname = NULL,
        # assay = NULL,
        # dimred = "PCA",
        # classes = NULL,
        # subset_intersection = FALSE,
        centroid_calculation = 'mean', 
        dims = c(1, 2),
        temperature = 1e3, 
        display_numbers = FALSE, 
        exclude_small_prob = FALSE, 
        min_prob_sum = NULL,
        # plot_mode = "dot",
        title_centroid_projection = "Centroids projected onto scaffold PCA",
        title_distance_pheatmap = "Distance of Sample to Scaffold Centroids",
        title_probability_pheatmap = "Softmax of Rank Distances of Sample to Scaffold Centroids",
        verbose = TRUE,
        annotation = "ensembl_gene"
        ){
    # if (is(scaffold, "NULL") | is(sample, "NULL")){
    #     warning("No scaffold or sample provided.")
    #     return(NULL)}
  # if (is(scaffold$annotation, "NULL")) {
  #   scaffold$annotation <- annotation
  # }
  # if (toupper(dimred) == "UMAP" && !requireNamespace("uwot")) {
  #   stop("To use UMAP space, please install uwot:\n",
  #        "install.packages(\"uwot\")")
  # }
  
    # LOADING FOR TEST - DELETE THIS
    # gtex_scaffold <- buildScaffold("GTEx.v1")
    # data("ilaria_counts", package = "spaceRATScaffolds")
    # data("ilaria_pData", package = "spaceRATScaffolds")
    # 
    # source('spaceRAT/R/projectSample.R')
    # projection <- projectSample(tcga_scaffold, sample = ilaria_counts, pheno = ilaria_pData, colname = 'cancer_type', subset_intersection = TRUE)
    # 
    # rm(gtex_scaffold, ilaria_counts, ilaria_pData)
    # 
    # dims <- c(1,2)
    # temperature <- 1e5
    # display_numbers <- FALSE
    # 
    # # KEEP FROM HERE
    
    # Initiate output 
    prediction <- list()
    
    # Calculate centroids
    class_sample <- projection$scaffold$label 
    centroids <- list()
    for (i in unique(class_sample)) {
      
      if (centroid_calculation == 'mean'){
        centroids[[i]] <- rowMeans(projection$scaffold$rank[,class_sample==i])
      } else if (centroid_calculation == 'median'){
        centroids[[i]] <- apply(projection$scaffold$rank[,class_sample==i], 1, median)
      }
      
    }
    
    # rm(class_sample)
    
    # source('spaceRAT/R/ranking.R')
    # Rank centroids 
    centroids_df <- as.data.frame(centroids) %>% ranking(rank_scale = projection$scaffold$rank_scale)
    
    
    # centroids_pheno <- data.frame(tissue = names(centroids))
    # rownames(centroids_pheno) <- colnames(centroids_df)
    
    # projectSample(projection$scaffold, sample = centroids_df, pheno = centroids_pheno, colname = 'tissue') 
    # projection$scaffold$pca$rotation <- projection$scaffold$pca$rotation[overlap_genes,]
    # projection$scaffold$pca$center <- projection$scaffold$pca$center[overlap_genes]
    # Transform the sample data
    # if (toupper(dimred) == "PCA"){
    transformed_centroids_df <- stats::predict(
      projection$scaffold$pca, newdata = t(centroids_df))
    # } else if (toupper(dimred) == "UMAP"){
    #   transformed_sample <- uwot::umap_transform(
    #     t(ranked_sample), scaffold$umap)
    # }
    
    # Prepare dataframe for ggplot
    transformed_centroids_df <- as.data.frame(transformed_centroids_df)
    transformed_centroids_df$shape <- '21'
    transformed_centroids_df$scaffoldGroup <- rownames(transformed_centroids_df)
    dim1 <- colnames(transformed_centroids_df)[dims[1]]
    dim2 <- colnames(transformed_centroids_df)[dims[2]]
    

    # Scaffold with centroids projected 
    
    p <- ggplot2::ggplot_build(projection)$data[[2]]
    cols <- unique(p[["colour"]])
    label <- as.character(p[["label"]])
    cols <- append(cols, "#000000")
    
    prediction$centroid_on_scaffold <- (plotScaffold(
      scaffold = projection$scaffold, title = title_centroid_projection) + 
      ggplot2::geom_point(
        data = transformed_centroids_df,
        mapping = ggplot2::aes(
          x = !!sym(dim1),
          y = !!sym(dim2),
          fill = scaffoldGroup), 
          shape = 21,
          color = 'black') + 
      ggplot2::scale_fill_manual(values = cols)
    )
    
    ## Calculates distance to centroid for each sample
    sample_distance_df <- c()
    nc_pred <- c()
    for(i in 1:ncol(projection$rankedSample)) {
      centroid_distances <- c()
      # for(ii in 1:length(centroids)) {
      for(ii in 1:length(colnames(centroids_df))) {
        # if (rank == FALSE) {
        #   centroid_distances <- c(centroid_distances, dist(rbind(sample[,i], centroids[[ii]])))
        # } else if (rank == TRUE){
        # centroid_distances <- c(centroid_distances, rankdist::DistancePair(projection$rankedSample[,i], centroids[[ii]]))
        # centroid_distances <- c(centroid_distances, rankdist::DistancePair(projection$rankedSample[,i], centroids_df[,ii]))
        centroid_distances <- c(centroid_distances, dist(rbind(projection$rankedSample[,i], centroids_df[,ii])))
        # }
      }
      # names(centroid_distances) <- names(centroids)
      names(centroid_distances) <- colnames(centroids_df)
      sample_distance_df <- rbind(sample_distance_df, c(colnames(projection$rankedSample)[i], centroid_distances))
      nc_pred <- c(nc_pred, names(which.min(centroid_distances)))
    }
    
    sample_distance_df <- as.data.frame(sample_distance_df)
    sample_distance_df <- tibble::column_to_rownames(sample_distance_df, 'V1')
    
    # Finetune output
    names(nc_pred) <- rownames(sample_distance_df)
    prediction$prediction_label <- nc_pred
    
    # Transform distances to numeric
    options(digits = 10)
    distance_mat <- sapply(sample_distance_df, function(x) as.numeric(x))
    rownames(distance_mat) <- rownames(sample_distance_df)
    
    # Prep pheno for pheatmap 
    if ('pheno' %in% names(projection)){
      pheno <- projection$pheno %>% dplyr::arrange(!!sym(projection$colnames)) %>% dplyr::select(!!sym(projection$colnames))
      distance_mat <- distance_mat[rownames(pheno), ]
    } else {
      pheno <- NA
    }
    
    prediction$distance_pheatmap <- ggplotify::as.ggplot(
      function() { pheatmap::pheatmap(distance_mat,
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    annotation_row = pheno,
                                    display_numbers = display_numbers,
                                    main = title_distance_pheatmap,
                                    silent = FALSE) }
      )
    
    # prediction$distance_pheatmap <- distance_pheatmap
    
    # Distance to probability using softmax
    # Larger temperature → smoother probabilities.
    # Smaller temperature → peakier, closer to winner-takes-all.
    exp_neg <- exp(-(distance_mat - apply(distance_mat, 1, min)) / temperature)
    prob_mat <- exp_neg / rowSums(exp_neg)
    
    if (exclude_small_prob) {
      if (is.null(min_prob_sum)) {
        stop("You must supply a 'min_prob_sum' value when exclude_small_prob = TRUE")
      }
      if (min_prob_sum < 0 || min_prob_sum > 1) {
        stop("'min_prob_sum' must be between 0 and 1")
      }
      prob_mat <- prob_mat[, colSums(prob_mat) > min_prob_sum, drop = FALSE] 
    }
    
    prediction$probability_pheatmap <- ggplotify::as.ggplot(
      function() { pheatmap::pheatmap(prob_mat, 
                                    cluster_rows = FALSE, 
                                    cluster_cols = FALSE, 
                                    annotation_row = pheno, 
                                    display_numbers = display_numbers, 
                                    main = title_probability_pheatmap, 
                                    breaks = seq(0, 1, length.out = 100),
                                    silent = FALSE) }
      )
    
    distance_mat <- tibble::rownames_to_column(.data = as.data.frame(distance_mat), var = 'Sample')
    prediction$distance_mat <- distance_mat
    
    prob_mat <- tibble::rownames_to_column(.data = base::as.data.frame(prob_mat), var = 'Sample')
    prediction$probability_mat <- prob_mat
    
    # prediction$probability_pheatmap <- probability_pheatmap
    
    return(prediction)
}

