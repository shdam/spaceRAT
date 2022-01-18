#' Project new sample(s) onto existing scaffold PCA plot
#'
#' This function takes in a \code{\link{scaffoldSpace}} objects, subsets the new dataset to contain only the genes that defines the scaffold space,
#' ranks the subset within new sample(s),
#' finally projects the new sample(s) onto the existing scaffold PCA plot.
#'
#' @importFrom Biobase exprs pData
#' @param space A scaffoldSpace object, returned by function \code{\link{buildScaffold}}
#' @param counts_sample Expression matrix of new sample, can be either logged or raw data without specification.
#' @param pheno_sample Phenotype data corresponding to \code{counts_sample}. If not specified, the output plot will not show legends for new samples.
#' @param colname A column name of \code{pheno_sample}. This column of values will be used to annotate projected samples.
#' This argument should be set together with \code{pheno_sample}.
#' If \code{pheno_sample} is not specified, this argument will be ignored, and the output plot will not show legends for new samples.
#' @param annotation Annotation type to use for scaffold. counts_scaffold rownames using alternative identifyers will be attemped translated. 
#' Currently ensembl_gene_id,entrezgene_id,hgnc_symbol, and refseq_mrna are supported. set to "NA", to avoid translation (both scaffold and projected sampels must be the same)

#' @param title Title of the plot.
#' @param verbose A logical vector indicating whether to report the number of genes imputed to make \code{counts_sample} compatible with \code{counts_scaffold}
#' @export
#' @return A ggplot object with new samples projected to existing scaffold PCA plot
#' @examples
#' projectSample(space,exprs_ilaria,pData_ilaria,"cancer_type")


projectSample <- function(space,
                          counts_sample,
                          pheno_sample=NULL,
                          # classes = NULL,
                          colname="cancer_type",
                          title="Samples projected onto scaffold PCA",
                          verbose=TRUE,
                          annotation="ensembl_gene_id"){

        # create eset
        if (!is.null(pheno_sample)){
                eset_sample <- createEset(counts_sample,pheno_sample,colname, annotation)
                counts_sample <- Biobase::exprs(eset_sample)
        }

        # add absent genes then subset eset_sample so counts_scaffold and counts_project contain same genes
        absent_genes <- space@DEgene[! space@DEgene %in% rownames(counts_sample)]
        absent_exprs <- matrix(0,length(absent_genes),ncol(counts_sample), dimnames = list(absent_genes, colnames(counts_sample)))
        # rownames(absent_exprs) <- absent_genes
        counts_sample <- rbind(counts_sample,absent_exprs)

        if (verbose){
                message(paste0(length(absent_genes)," genes are added to count matrix, with imputed expression level 0."))
        }

        if (length(absent_genes)/length(space@DEgene)>1/4){
                warning("Over 1/4 genes are added to sample with imputed expression level 0!")
        }

        #subset
        counts_sample <- counts_sample[space@DEgene,]

        # rank and transform exprs_project and multiply with the percent of missing values, to retain comparable numeric range
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
                suppressMessages(g <- graph+
                        ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,color="New_samples"))+
                        ggplot2::scale_color_manual(values = cols)+
                        ggplot2::ggtitle(title)+
                        ggplot2::coord_fixed()
                )
                return (g)
        }

        New_samples <- Biobase::pData(eset_sample)[[colname]]
        df_sample <- data.frame(PC1_sample,PC2_sample,New_samples)

        # project points
        suppressMessages(
        g <- graph+
                ggplot2::geom_point(data=df_sample, mapping=ggplot2::aes(PC1_sample,PC2_sample,shape=New_samples),color="black")+
                ggplot2::scale_shape_manual(values=1:length(unique((New_samples))))+
                ggplot2::ggtitle(title)+
                ggplot2::coord_fixed()
        )
        return(g)
}
