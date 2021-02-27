#' Project new sample(s) onto the existing scaffold PCA plot
#'
#' This function takes in a scaffoldSpace objects, subsets the new dataset, ranks the subset, finally projects the new sample(s) onto the existing scaffold PCA plot.
#' @importFrom Biobase exprs pData
#' @param space a scafoldSpace object, returned by function \code{buildScaffold()}
#' @param exprs_sample expression matrix of new sample
#' @param pData_sample phenotype data corresponding to exprs_sample.
#' @param group_sample a column name of pData_sample
#' @param title Title of the plot
#' @param verbose a logical vector indicating whether to report the number of genes added to eset_sample to make it compatible with eset_scaffold
#' @export
#' @return a ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")


projectSample <- function(space,
                          exprs_sample,
                          pData_sample,
                          group_sample,
                          title = "Samples projected on scaffold PCA",
                          y = "PC2",
                          x = "PC1",
                          verbose = T){
        # check gene name concordance
        if (all(! space@DEgenes %in% rownames(exprs_sample))){
                stop("Row names (i.e. gene names) of scaffold expression matrix and new sample expression matrix do not match! ")
        }

        # create eset
        eset_sample <- createEset(exprs_sample,pData_sample)

        # add absent genes then subset eset_sample so eset_scaffold and eset_project contain same genes
        exprs_sample <- Biobase::exprs(eset_sample)
        absent_genes <- space@DEgenes[! space@DEgenes %in% rownames(exprs_sample)]
        absent_exprs <- matrix(0,length(absent_genes),ncol(exprs_sample))
        rownames(absent_exprs) <- absent_genes
        exprs_sample <- rbind(exprs_sample,absent_exprs)

        if (verbose){
                message(paste0(length(absent_genes)," genes are added to sample, with imputed expression level 0."))
        }

        if (length(absent_genes)/length(space@DEgenes)>1/3){
                warning("Over 1/3 genes are added to sample with imputed expression level 0!")
        }

        #subset
        exprs_sample <- exprs_sample[space@DEgenes,]

        # rank and transform exprs_project
        ranked_sample<- apply(exprs_sample,2,rank)

        # PCA transform the sample data
        transformed_sample <- predict(space@pca,newdata=t(ranked_sample))

        # Prepare dataframe for ggplot
        PC1_sample <- transformed_sample[,1]
        PC2_sample <- transformed_sample[,2]
        sample_group <- Biobase::pData(eset_sample)[[group_sample]]
        df_sample <- data.frame(PC1_sample,PC2_sample,sample_group)


        # project points
        g <- space@graph +
                ggplot2::geom_point(data = df_sample,
                                    mapping = ggplot2::aes(PC1_sample,
                                                           PC2_sample,
                                                           shape = sample_group,
                                                           color = "New samples")) +
                ggplot2::scale_shape_manual(values = 1:length(unique((df_sample$sample_group)))) +
                ggplot2::labs(title = title,
                              y = y,
                              x = x) +
                ggplot2::coord_fixed() +
                ggplot2::theme_bw()
        return(g)
}
