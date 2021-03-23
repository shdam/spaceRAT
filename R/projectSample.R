#' Project new sample(s) onto the existing scaffold PCA plot
#'
#' This function takes in a scaffoldSpace objects, subsets the new dataset, ranks the subset, finally projects the new sample(s) onto the existing scaffold PCA plot.
#' @importFrom Biobase exprs pData
#' @param space a scafoldSpace object, returned by function \code{buildScaffold()}
#' @param counts_sample expression matrix of new sample
#' @param pheno_sample phenotype data corresponding to exprs_sample.
#' @param colname a column name of pData_sample
#' @param title Title of the plot
#' @param verbose a logical vector indicating whether to report the number of genes added to eset_sample to make it compatible with eset_scaffold
#' @export
#' @return a ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")


projectSample <- function(space, counts_sample, pheno_sample=NULL,colname,title="New samples projected on scaffold PCA",verbose=T){

        # create eset
        eset_sample <- createEset(counts_sample,pheno_sample,colname)

        # add absent genes then subset eset_sample so eset_scaffold and eset_project contain same genes
        counts_sample <- Biobase::exprs(eset_sample)
        absent_genes <- space@DEgenes[! space@DEgenes %in% rownames(counts_sample)]
        absent_exprs <- matrix(0,length(absent_genes),ncol(counts_sample))
        rownames(absent_exprs) <- absent_genes
        counts_sample <- rbind(counts_sample,absent_exprs)

        if (verbose){
                message(paste0(length(absent_genes)," genes are added to count matrix, with imputed expression level 0."))
        }

        if (length(absent_genes)/length(space@DEgenes)>1/3){
                warning("Over 1/3 genes are added to sample with imputed expression level 0!")
        }

        #subset
        counts_sample <- counts_sample[space@DEgenes,]

        # rank and transform exprs_project
        ranked_sample<- apply(counts_sample,2,rank)

        # PCA transform the sample data
        transformed_sample <- predict(space@pca,newdata=t(ranked_sample))

        # Prepare dataframe for ggplot
        PC1_sample <- transformed_sample[,space@pcs[1]]
        PC2_sample <- transformed_sample[,space@pcs[2]]
        new_samples <- Biobase::pData(eset_sample)[[colname]]
        df_sample <- data.frame(PC1_sample,PC2_sample,new_samples)


        # project points
        g <- space@graph+
                ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,shape=new_samples),color="black")+
                ggplot2::scale_shape_manual(values=1:length(unique((new_samples))))+
                ggplot2::ggtitle(title)+
                ggplot2::coord_fixed()
        return(g)
}
