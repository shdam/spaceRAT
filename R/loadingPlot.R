#' Generate a loading plot of scaffold PCA
#'
#' This function takes a scaffoldSpace object and generate a loading plot to visualize the genes that contribute most to selected principle components.
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @param space a scaffoldSpace object created by function buildScaffold()
#' @param num_genes an optional parameter indicating number of genes to be shown
#' @export
#' @return a ggplot object of loading plot.
#' @examples
#' space <- buildScaffold(exprs_dmap,pData_dmap"cell_types")
#' loadingPlot(space)

loadingPlot <- function(space,num_genes=3,gene_name="hgnc_symbol"){
        pca <- space@pca
        pcs <- space@pcs

        # calculate variance explained and add percentage to axis labels
        var_sum <- sum(pca$sdev^2)
        var1 <- round(pca$sdev[pcs[1]]^2/var_sum*100,2)
        var2 <- round(pca$sdev[pcs[2]]^2/var_sum*100,2)
        xlabel <- paste0( "PC",as.character(pcs[1]), " (",as.character(var1), "%", ")" )
        ylabel <- paste0( "PC",as.character(pcs[2]), " (",as.character(var2), "%", ")" )

        # determine most important PCs
        datapc <- as.data.frame(pca$rotation[,pcs])
        datapc <- tibble::rownames_to_column(datapc, "ensembl_gene_id")[,c(2,3,1)]
        pc1_ordered <- datapc[order(datapc[,1], decreasing = T),]
        pc2_ordered <- datapc[order(datapc[,2], decreasing = T),]
        pos1 <- data.frame(head(pc1_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[1])),num_genes))
        neg1 <- data.frame(tail(pc1_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[1])),num_genes))
        pos2 <- data.frame(head(pc2_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[2])),num_genes))
        neg2 <- data.frame(tail(pc2_ordered,num_genes),class=rep(paste0("Top ",num_genes," genes for PC",as.character(pcs[2])),num_genes))
        df <- rbind(pos1,neg1,pos2,neg2)

        gene_mapper <- mapGene(df$ensembl_gene_id,to=gene_name)
        from <- colnames(gene_mapper)[1]
        df <- merge(gene_mapper,df)

        xmin <- min(df[3])
        xmax <- max(df[3])
        ymin <- min(df[4])
        ymax <- max(df[4])

        # ggplot
        g <-ggplot2::ggplot(data=df)+
            ggplot2::geom_text(ggplot2::aes_string(colnames(df)[3], colnames(df)[4],label=colnames(df)[2],color=colnames(df)[5]),size = 3, angle=30,show.legend = FALSE)+
            ggplot2::geom_segment(ggplot2::aes_string(x=0, y=0, xend=colnames(df)[3], yend=colnames(df)[4],color=colnames(df)[5]), arrow=ggplot2::arrow(length=ggplot2::unit(0.1,"cm")), alpha=0.75)+
            ggplot2::theme_bw()+
            ggplot2::theme(legend.text=ggplot2::element_text(hjust=1))+
            ggplot2::xlim(xmin*1.5,xmax*1.5)+
            ggplot2::ylim(ymin*1.5,ymax*1.5)+
            ggplot2::xlab(xlabel)+ggplot2::ylab(ylabel)+
            ggplot2::coord_fixed()+
            ggplot2::ggtitle("Loading plot")

        print(g)

        return(df[order(df[,3],decreasing=T),])
}

