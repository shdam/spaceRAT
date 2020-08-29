#' @title Project new sample(s) onto the existing scaffold PCA plot
#' @description This function takes in objects returned from function \code{buildScaffold()}, subsets the new dataset, ranks the subset, finally projects the new sample(s) onto the existing scaffold PCA plot.
#' @import Biobase
#' @param list_from_buildScaffold the list returned by function \code{buildScaffold()}
#' @param eset_sample the sample(s) to project, of class ExpressionSet
#' @param group_sample a column name of pData(eset_sample)
#' @param title Title of the plot
#' @param verbose a logical vector indicating whether to report the number of genes added to eset_sample to make it compatible with eset_scaffold
#' @export
#' @return a ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' projectSample(l,eset_sample,"cancer_type")


projectSample <- function(list_from_buildScaffold, eset_sample, group_sample,title="Samples projected on scaffold PCA",verbose=T){
        # unpack list
        g <- list_from_buildScaffold[[1]]
        pca <- list_from_buildScaffold[[2]]
        DE_genes <- rownames(pca$rotation)

        # add absent genes then subset eset_sample so eset_scaffold and eset_project contain same genes
        exprs_sample <- exprs(eset_sample)
        absent_genes <- DE_genes[! DE_genes %in% rownames(exprs(eset_sample))]
        absent_exprs <- matrix(0,length(absent_genes),ncol(exprs_sample))
        rownames(absent_exprs) <- absent_genes
        exprs_sample <- rbind(exprs_sample,absent_exprs)

        if (verbose){
                print(paste0(length(absent_genes)," genes are added to sample, with imputed expression level 0."))
        }

        #subset
        exprs_sample <- exprs_sample[DE_genes,]

        # rank and transform exprs_project
        ranked_sample<- apply(exprs_sample,2,rank)

        # PCA transform the sample data
        transformed_sample <- predict(pca,newdata=t(ranked_sample))

        # Prepare dataframe for ggplot
        PC1_sample <- transformed_sample[,1]
        PC2_sample <- transformed_sample[,2]
        sample_group <- pData(eset_sample)[[group_sample]]
        df_sample <- data.frame(PC1_sample,PC2_sample,sample_group)

        # project points
        g <- g+ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,shape=sample_group,color="New_samples"))+ggplot2::ggtitle(title)
        print(g)
        return(g)

}
