#' Project new sample(s) onto the existing scaffold PCA plot
#'
#' This function takes in a scaffoldSpace objects, subsets the new dataset, ranks the subset, finally projects the new sample(s) onto the existing scaffold PCA plot.
#' @importFrom Biobase exprs pData
#' @param space a scafoldSpace object, returned by function \code{buildScaffold()}
#' @param counts_sample expression matrix of new sample
#' @param pheno_sample phenotype data corresponding to counts_sample. If not specified, the output plot will not show legends for new samples.
#' @param colname a column name of pheno_sample. This argument should be set together with pheno_sample. If pheno_sample is not specified, this argument will be ignored, thus the output plot will not show legends for new samples.
#' @param title Title of the plot
#' @param verbose a logical vector indicating whether to report the number of genes added to eset_sample to make it compatible with eset_scaffold
#' @export
#' @return a ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")


projectSample <- function(space, counts_sample, pheno_sample=NULL,colname=NULL, title="Samples projected on scaffold PCA",verbose=T){

        # create eset
        if (!is.null(pheno_sample)){
                eset_sample <- createEset(counts_sample,pheno_sample,colname)
                counts_sample <- Biobase::exprs(eset_sample)
        }

        # add absent genes then subset eset_sample so eset_scaffold and eset_project contain same genes
        absent_genes <- space@DEgene[! space@DEgene %in% rownames(counts_sample)]
        absent_exprs <- matrix(0,length(absent_genes),ncol(counts_sample))
        rownames(absent_exprs) <- absent_genes
        counts_sample <- rbind(counts_sample,absent_exprs)

        if (verbose){
                message(paste0(length(absent_genes)," genes are added to count matrix, with imputed expression level 0."))
        }

        if (length(absent_genes)/length(space@DEgene)>1/4){
                warning("Over 1/4 genes are added to sample with imputed expression level 0!")
        }

        #subset
        counts_sample <- counts_sample[space@DEgene,]

        # rank and transform exprs_project
        ranked_sample<- apply(counts_sample,2,rank)

        # PCA transform the sample data
        transformed_sample <- predict(space@pca,newdata=t(ranked_sample))

        # Prepare dataframe for ggplot
        PC1_sample <- transformed_sample[,space@pcs[1]]
        PC2_sample <- transformed_sample[,space@pcs[2]]

        graph <- plotScaffold(space,title=title)
        if (is.null(pheno_sample)){
                # calculate correct color scale
                p <- ggplot2::ggplot_build(graph)$data[[2]]
                cols <- p[["colour"]]
                label <- as.character(p[["label"]])
                idx <- rank(c(label,"New_samples"))[length(label)+1]
                cols <- append(cols,"#000000",after = idx-1)
                df_sample <- data.frame(PC1_sample,PC2_sample)
                g <- graph+
                        ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,color="New_samples"))+
                        ggplot2::scale_color_manual(values = cols)
                        ggplot2::ggtitle(title)+
                        ggplot2::coord_fixed()
                return (g)
        }

        new_samples <- Biobase::pData(eset_sample)[[colname]]
        df_sample <- data.frame(PC1_sample,PC2_sample,new_samples)

        # project points
        g <- graph+
                ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,shape=new_samples),color="black")+
                ggplot2::scale_shape_manual(values=1:length(unique((new_samples))))+
                ggplot2::ggtitle(title)+
                ggplot2::coord_fixed()
        return(g)
}
