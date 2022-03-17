#' Convert between Ensembl gene, Ensembl transcript, Entrez, RefSeq and Gene symbol for count ?matrix
#'
#' This function takes in a count matrix, and converts its row names to Ensembl (by default), or any specified gene identifier.
#' @import magrittr
#' @importFrom dplyr select group_by_at across everything summarise filter one_of
#' @importFrom tibble column_to_rownames
#' @param counts An expression matrix of class matrix
#' @param to A character specifying the gene id to convert to. Options are "ensembl_gene" (default), "ensembl_transcript", "entrez", "hgnc_symbol" and "refseq_mrna".
#' @return A count matrix with transformed gene names as row names.
#' @noRd
#' @examples
#' convertGeneName(exprs_dmap[1:20,1:10],to="hgnc_symbol")

convertGeneName <- function(counts,to="ensembl_gene"){
        # remove version number
        cur_genes <- gsub("\\.[0-9]+$","",rownames(counts))

        # get gene mapper data frame
        gene_mapper <- mapGene(cur_genes,to=to)

        # case1: mapGene fails to infer gene id
        if (is.null(gene_mapper)){
                stop("Could not infer gene identifiers from row names of expression matrix.
                     Please set annotation = NA to stop gene id conversion and make sure manually that new samples have same gene ids as scaffold dataset")
        }

        from <- colnames(gene_mapper)[1]

        # case 2: from==to. Only need to remove version number before return
        if (from==to) {
                rownames(counts) <- cur_genes
                return(as.matrix(counts))
        }

        # case 3: conversion between gene id and gene id
        if (from %in% c("ensembl_gene","entrez","hgnc_symbol")){
                rownames(gene_mapper) <- gene_mapper[,1]
                counts <- counts[rownames(gene_mapper),]
                rownames(counts) <- gene_mapper[,2]
                return(as.matrix(counts))
        }

        # case 4: conversion from transcripts to genes
        if (from %in% c("ensembl_transcript","refseq_mrna")){ # add all matches
                counts <- merge(gene_mapper,counts,by.x=from,by.y=0,all=F) %>%
                        dplyr::select(-dplyr::one_of(from))%>%
                        dplyr::group_by_at(to) %>%
                        dplyr::summarise(dplyr::across(dplyr::everything(),sum)) %>%
                        tibble::column_to_rownames(var=to)
                return (as.matrix(counts))
        }

        # all unexpected cases
        stop("Could not infer gene identifiers from row names of expression matrix.
                    Please set annotation = NA to stop gene id conversion and make sure manually that new samples have same gene ids as scaffold dataset")
}
