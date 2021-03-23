#' Returns a ggplot and PCA parameters for ExpressionSet object
#'
#' This function performs differential expression analysis to select DE gene, then subset the scaffold dataset, rank the subset, finally produce a PCA plot using ranks.
#' @import dplyr
#' @import  ggplot2
#' @importFrom Biobase pData exprs
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @param counts_scaffold an expression matrix of class matrix, or data frame that can be converted to matrix. To use our prebuilt scaffold, you can set count_scaffold="prebuilt_DMAP", or "prebuilt_GTEX".
#' @param pheno_scaffold a phenotype table corresponding to the expression matrix.
#' @param data a character indicating whether the count matrix is log-tansformed or raw. To perform DE analysis, dataset need to be log-transformed, so setting "data=raw" will automatically perform log transformation before DE analysis.
#' @param colname a column name of pData_scaffold. Differential analysis will be performed using this column of phenotype as independent variables.
#' @param pval_cutoff a cutoff value for p value when selecting differentially expressed genes. By default 0.05
#' @param lfc_cutoff a cutoff value for logFC when selecting differentially expressed genes. By default 2
#' @param pca_scale a logical variable determining whether to normalize rows when plotting PCA
#' @export
#' @return A scaffoldSpace object
#' @examples
#' buildScaffold(exprs_dmap,pData_dmap"cell_types")
#' buildScaffold(exprs_dmap,pData_dmap,"cell_types", pval_cutoff=0.01,pca_scale=T)

buildScaffold <- function(counts_scaffold,
                          pheno_scaffold,
                          colname,
                          data="logged",
                          pcs=c(1,2),
                          plot_mode="dot",
                          pval_cutoff=0.05,
                          lfc_cutoff=2,
                          title="Scaffold PCA Plot",
                          pca_scale=FALSE){

        if(is.character(counts_scaffold) && counts_scaffold=="prebuilt_DMAP"){         # not implemented
                load("data/DMAP_scaffold.rda")
                scaffold <- DMAP_scaffold
                scaffold@pcs <- pcs
                return(scaffold)
        }

        # create eset
        if (data=="logged"){
                # remove genes with total count<10
                idx <- which(rowSums(exp(counts_scaffold))<10)
                if (length(idx)>0) counts_scaffold <- counts_scaffold[-idx,]
                eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname)
        }
        if (data=="raw"){
                # remove genes with total count<10
                idx <- which(rowSums(counts_scaffold)<10)
                counts_scaffold <- counts_scaffold[-idx,]
                counts_scaffold <- log(counts_scafold+1)
                eset_scaffold <- createEset(counts_scaffold,pheno_scaffold,colname)
        }

        # find DE genes
        DEgenes <- findDEGenes(eset_scaffold,pval_cutoff,lfc_cutoff)

        # subset
        eset_scaffold <- eset_scaffold[DEgenes,]

        # rank
        scaffold_rank<- apply(Biobase::exprs(eset_scaffold),2,rank)

        # PCA
        pca <- prcomp(t(scaffold_rank),scale=pca_scale)

        # calculate variance explained and add percentage to axis labels
        var_sum <- sum(pca$sdev^2)
        var1 <- round(pca$sdev[pcs[1]]^2/var_sum*100,2)
        var2 <- round(pca$sdev[pcs[2]]^2/var_sum*100,2)
        xlabel <- paste0( "PC",as.character(pcs[1]), " (",as.character(var1), "%", ")" )
        ylabel <- paste0( "PC",as.character(pcs[2]), " (",as.character(var2), "%", ")" )


        # Prepare a dataframe for ggplot2
        PC1 <- pca$x[,pcs[1]]
        PC2 <- pca$x[,pcs[2]]
        scaffold_group <- Biobase::pData(eset_scaffold)[[colname]]
        df <- data.frame(PC1,PC2,scaffold_group)

        # calculate centroids
        library(dplyr)
        centroids_df <- suppressMessages(df %>% dplyr::group_by(scaffold_group) %>% dplyr::summarise(mean_PC1=mean(PC1),mean_PC2=mean(PC2)))

        # define color scheme
        total_types <- length(unique(scaffold_group))
        my_col <- RColorBrewer::brewer.pal(total_types,"Paired")
        if (total_types>12){
                pal <- grDevices::colorRampPalette(my_col)
                my_col <- pal(total_types)
                if (plot_mode=="dot"){
                        warning("More than 12 cell types are to be displayed. Setting 'plot_mode='tiny_label' may yield better visualization.")
                }

        }
        # ggplot2 for dot mode
        if (plot_mode=="dot"){
                g <- ggplot2::ggplot()+
                        ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=scaffold_group))+
                        ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_PC1,mean_PC2,label=scaffold_group,color=scaffold_group),fontface="bold",show.legend = FALSE)+
                        ggplot2::ggtitle(title)+
                        ggplot2::scale_color_manual(colname,values=my_col)+
                        ggplot2::xlab(xlabel)+
                        ggplot2::ylab(ylabel)+
                        ggplot2::theme_bw()+
                        ggplot2::coord_fixed()
        }

        # ggplot2 for label mode
        if (plot_mode=="tiny_label"){
                message("plot_mode='tiny_label', shorter names for cell types in phenotype table yields better visualization.")
                g <- ggplot2::ggplot()+
                        ggplot2::geom_point(data=df,mapping=ggplot2::aes(PC1,PC2,color=scaffold_group),size = 0.8)+
                        ggplot2::geom_text(data=df,ggplot2::aes(PC1,PC2,label=scaffold_group,color=scaffold_group),alpha=0.5,size=3,show.legend = FALSE)+
                        ggplot2::geom_label(data=centroids_df,ggplot2::aes(mean_PC1,mean_PC2,label=scaffold_group,color=scaffold_group,fontface="bold"),show.legend = FALSE)+
                        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2)))+
                        ggplot2::ggtitle(title)+
                        ggplot2::xlab(xlabel)+
                        ggplot2::ylab(ylabel)+
                        ggplot2::scale_color_manual(colname,values=my_col)+
                        ggplot2::theme_bw()+
                        ggplot2::coord_fixed()
        }

        # record standard data in scaffoldSpace class
        print(g)
        scaffold <- new("scaffoldSpace",
                        graph=g,
                        pca=pca,
                        DEgenes=DEgenes,
                        pcs=pcs)

        return(scaffold)
}
