#' Convert between Ensembl gene, Ensembl transcript, Entrez, RefSeq and Gene symbol for vector
#'
#' This function takes in a charactor vector of gene names, determine which gene identifier it is, then converts the gene names to the specified type, finally returns list containing 2 objects.
#' @import dplyr
#' @param vec charactor vector of gene names
#' @param to a character specifying the gene id to convert to. Options are "ensembl_gene_id", "ensembl_transcript_id", "entrezgene_id", "hgnc_symbol" and "refseq_mrna".
#' @export
#' @return a list, containing a character indicating current type of gene name, and a data frame with orginal gene names as first column, converted gene names as second column.

mapGene <- function(vec,to){

        from <- NULL
        # determine which gene id is used
        load("data/gene_name_mapper_hs.rda") # this line shouldn't be needed if the package is built correctly
        if(all(startsWith(vec,"ENSG"))){
                from <- "ensembl_gene_id"
        }else if(all(startsWith(vec,"ENST"))){
                from <- "ensembl_transcript_id"
        }else if(all(startsWith(vec,"NM_"))){
                from <- "refseq_mrna"
        }else if(all(grepl("^[0-9]+$",vec))){
                from <- "entrezgene_id"
        }else{
                idx <- which(gene_name_mapper_hs$hgnc_symbol %in% vec)
                if (length(idx)>0) from <- "hgnc_symbol"
        }

        if (from==to) {
                return(list(from,data.frame(vec,vec)))
        }else if (is.null(from)) {
                stop("Gene names are not recognized. Please manually convert the gene names to ",to, ".")
        }else if(from!=to){
                message("Convert gene names of count matrix from ",from," to ",to,".")
        }

        # subset gene_name_mapper_hs
        library(dplyr)
        if(from!= "hgnc_symbol") idx <- which(gene_name_mapper_hs[[from]] %in% vec)
        df <- gene_name_mapper_hs[idx,c(from,to),drop=F] %>%
                group_by_at(from) %>%
                filter(row_number()==1)%>%    # select the first match, ensure 1-1 matching
                as.data.frame()

        return(list(from,df))
}
